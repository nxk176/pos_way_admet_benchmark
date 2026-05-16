from __future__ import annotations

import argparse
import json
import math
import sqlite3
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

try:
    from rdkit import Chem
except Exception as exc:  # noqa: BLE001 - this extractor standardizes structures with RDKit.
    raise SystemExit(
        "ERROR: extract_chembl_experimental_measurements.py requires RDKit. "
        "Run it with the project environment, for example .\\myenv311\\Scripts\\python.exe"
    ) from exc


MeasurementNormalizer = Callable[[sqlite3.Row], tuple[float, str] | None]


def canon_smiles(smiles: str | None) -> str | None:
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def lower_text(value: Any) -> str:
    return str(value or "").strip().lower()


def normalize_molar_log10(row: sqlite3.Row) -> tuple[float, str] | None:
    value = float(row["standard_value"])
    if value <= 0:
        return None
    standard_type = lower_text(row["standard_type"]).replace(" ", "")
    units = lower_text(row["standard_units"])
    if standard_type in {"logs", "logs0", "logsolubility", "log(kinetic_solubility)"} or "log s" in lower_text(row["standard_type"]):
        return value, "log10_mol_L"
    if units in {"nm", "nanomolar"}:
        return math.log10(value * 1e-9), "log10_mol_L"
    if units in {"um", "µm", "μm", "micromolar"}:
        return math.log10(value * 1e-6), "log10_mol_L"
    if units in {"mm", "millimolar"}:
        return math.log10(value * 1e-3), "log10_mol_L"
    if units in {"m", "molar"}:
        return math.log10(value), "log10_mol_L"
    return None


def normalize_papp_log10(row: sqlite3.Row) -> tuple[float, str] | None:
    value = float(row["standard_value"])
    if value <= 0:
        return None
    standard_type = lower_text(row["standard_type"]).replace(" ", "")
    units = lower_text(row["standard_units"])
    if standard_type in {"logpapp", "logpe"}:
        return value, "log10_cm_s"
    if units in {
        "10^-6 cm/s",
        "10'-6 cm/s",
        "1e-6 cm/s",
        "10e-6 cm/s",
        "ucm/s",
        "cm/s * 10e6",
        "10'6cm/s",
    }:
        return math.log10(value * 1e-6), "log10_cm_s"
    if units in {"10'-7 cm/s", "10^-7 cm/s"}:
        return math.log10(value * 1e-7), "log10_cm_s"
    if units in {"10'-5 cm/s", "10^-5cm/s", "10'-5cm/s"}:
        return math.log10(value * 1e-5), "log10_cm_s"
    if units in {"nm/s", "nm s-1"}:
        return math.log10(value * 1e-7), "log10_cm_s"
    if units in {"cm/s", "cm s-1"}:
        return math.log10(value), "log10_cm_s"
    return None


def normalize_pic50(row: sqlite3.Row) -> tuple[float, str] | None:
    value = float(row["standard_value"])
    units = lower_text(row["standard_units"])
    if value <= 0:
        return None
    if units == "nm":
        return -math.log10(value * 1e-9), "pIC50"
    if units in {"um", "µm", "μm"}:
        return -math.log10(value * 1e-6), "pIC50"
    if units == "mm":
        return -math.log10(value * 1e-3), "pIC50"
    return None


def normalize_clearance(row: sqlite3.Row) -> tuple[float, str] | None:
    value = float(row["standard_value"])
    units = lower_text(row["standard_units"])
    if value < 0:
        return None
    if units in {"ml.min-1.kg-1", "ml/min/kg", "ml min-1 kg-1"}:
        return value, "mL_min_kg"
    return None


def normalize_half_life(row: sqlite3.Row) -> tuple[float, str] | None:
    value = float(row["standard_value"])
    units = lower_text(row["standard_units"])
    if value < 0:
        return None
    if units in {"hr", "h", "hour", "hours"}:
        return value * 60.0, "min"
    if units in {"min", "minute", "minutes"}:
        return value, "min"
    return None


ENDPOINT_RULES: dict[str, dict[str, Any]] = {
    "logS_mol_L": {
        "description_patterns": ["%solub%"],
        "standard_types": ["solubility", "log s", "logs0", "log(kinetic_solubility)"],
        "standard_type_patterns": ["%solub%", "%log s%", "%logs%"],
        "standard_flags": [1],
        "normalizer": normalize_molar_log10,
    },
    "Caco2_logPapp_cm_s": {
        "description_patterns": ["%caco%", "%papp%"],
        "standard_types": ["papp", "permeability", "logpapp"],
        "standard_type_patterns": ["%papp%", "%permeab%", "%logpapp%", "%log pe%"],
        "standard_flags": [0, 1],
        "normalizer": normalize_papp_log10,
    },
    "hERG_pIC50": {
        "description_patterns": ["%herg%"],
        "standard_types": ["ic50", "ki", "potency"],
        "standard_type_patterns": ["%ic50%", "%ki%", "%potency%"],
        "standard_flags": [1],
        "normalizer": normalize_pic50,
    },
    "microsomal_clearance_mL_min_kg": {
        "description_patterns": ["%microsom%", "%clearance%"],
        "standard_types": ["cl", "clint", "clh", "cl/f"],
        "standard_type_patterns": ["%cl%", "%clearance%"],
        "standard_flags": [1],
        "normalizer": normalize_clearance,
    },
    "half_life_min": {
        "description_patterns": ["%half-life%", "%half life%", "%t1/2%"],
        "standard_types": ["t1/2", "t50", "plasma half-life", "plasma half life"],
        "standard_type_patterns": ["%t1/2%", "%t50%", "%half%"],
        "standard_flags": [1],
        "normalizer": normalize_half_life,
    },
}


BASE_QUERY = """
select
  a.activity_id,
  a.standard_type,
  a.standard_relation,
  a.standard_value,
  a.standard_units,
  a.standard_flag,
  a.pchembl_value,
  a.data_validity_comment,
  a.potential_duplicate,
  md.chembl_id as molecule_chembl_id,
  cs.canonical_smiles,
  cs.standard_inchi_key,
  assays.assay_id,
  assays.chembl_id as assay_chembl_id,
  assays.description as assay_description,
  assays.assay_type,
  assays.assay_organism,
  assays.assay_cell_type,
  assays.assay_category,
  docs.year as document_year,
  docs.doi,
  docs.pubmed_id,
  docs.title as document_title
from activities a
join assays on a.assay_id = assays.assay_id
join molecule_dictionary md on a.molregno = md.molregno
join compound_structures cs on a.molregno = cs.molregno
left join docs on a.doc_id = docs.doc_id
where a.standard_value is not null
  and a.standard_flag in ({standard_flag_placeholders})
  and a.standard_relation in ('=', '<', '>', '<=', '>=')
  and (lower(a.standard_type) in ({standard_type_placeholders}) or {standard_type_pattern_clause})
  and ({description_clause})
order by a.activity_id
limit ?
"""


def endpoint_query(rule: dict[str, Any]) -> tuple[str, list[Any]]:
    standard_types = [item.lower() for item in rule["standard_types"]]
    type_placeholders = ",".join("?" for _ in standard_types)
    standard_flags = list(rule.get("standard_flags", [1]))
    flag_placeholders = ",".join("?" for _ in standard_flags)
    standard_type_patterns = rule.get("standard_type_patterns", [])
    standard_type_pattern_clause = " or ".join("lower(a.standard_type) like ?" for _ in standard_type_patterns) or "0"
    description_patterns = rule["description_patterns"]
    description_clause = " or ".join("lower(assays.description) like ?" for _ in description_patterns)
    sql = BASE_QUERY.format(
        standard_flag_placeholders=flag_placeholders,
        standard_type_placeholders=type_placeholders,
        standard_type_pattern_clause=standard_type_pattern_clause,
        description_clause=description_clause,
    )
    return sql, [*standard_flags, *standard_types, *standard_type_patterns, *description_patterns]


def confidence(row: sqlite3.Row) -> float:
    score = 0.9
    if row["data_validity_comment"]:
        score -= 0.2
    if row["potential_duplicate"]:
        score -= 0.15
    if row["standard_relation"] != "=":
        score -= 0.1
    return round(max(0.3, score), 3)


def extract(args: argparse.Namespace) -> tuple[int, dict[str, Any]]:
    args.out.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(args.db)
    con.row_factory = sqlite3.Row
    counters: Counter[str] = Counter()
    endpoint_counts: Counter[str] = Counter()
    unit_counts: Counter[str] = Counter()
    seen_activity_ids: set[int] = set()
    written = 0

    with args.out.open("w", encoding="utf-8", newline="\n") as handle:
        for endpoint_name, rule in ENDPOINT_RULES.items():
            sql, params = endpoint_query(rule)
            normalizer: MeasurementNormalizer = rule["normalizer"]
            rows_seen_for_endpoint = 0
            rows_written_for_endpoint = 0

            for row in con.execute(sql, [*params, args.sql_limit_per_endpoint]):
                rows_seen_for_endpoint += 1
                counters[f"{endpoint_name}:sql_rows_seen"] += 1
                activity_id = int(row["activity_id"])
                if activity_id in seen_activity_ids:
                    counters["duplicate_activity_id_skipped"] += 1
                    continue
                normalized = normalizer(row)
                if normalized is None:
                    counters[f"{endpoint_name}:normalization_skipped"] += 1
                    continue
                smiles_canon = canon_smiles(row["canonical_smiles"])
                if smiles_canon is None:
                    counters["invalid_smiles_skipped"] += 1
                    continue

                value_canonical, unit_canonical = normalized
                record = {
                    "measurement_id": f"chembl36_act_{activity_id}",
                    "source": "ChEMBL",
                    "source_version": "36",
                    "endpoint_name": endpoint_name,
                    "label_type": "experimental",
                    "experimental_only_flag": True,
                    "molecule_chembl_id": row["molecule_chembl_id"],
                    "smiles_raw": row["canonical_smiles"],
                    "smiles_canon": smiles_canon,
                    "inchikey": row["standard_inchi_key"],
                    "connectivity_key": str(row["standard_inchi_key"]).split("-", 1)[0],
                    "value_raw": row["standard_value"],
                    "unit_raw": row["standard_units"],
                    "relation_raw": row["standard_relation"],
                    "type_raw": row["standard_type"],
                    "value_canonical": round(float(value_canonical), 6),
                    "unit_canonical": unit_canonical,
                    "confidence": confidence(row),
                    "assay": {
                        "assay_id": row["assay_id"],
                        "assay_chembl_id": row["assay_chembl_id"],
                        "description": row["assay_description"],
                        "assay_type": row["assay_type"],
                        "assay_organism": row["assay_organism"],
                        "assay_cell_type": row["assay_cell_type"],
                        "assay_category": row["assay_category"],
                    },
                    "document": {
                        "year": row["document_year"],
                        "doi": row["doi"],
                        "pubmed_id": row["pubmed_id"],
                        "title": row["document_title"],
                    },
                    "quality_flags": {
                        "data_validity_comment": row["data_validity_comment"],
                        "potential_duplicate": bool(row["potential_duplicate"]),
                        "standard_relation_exact": row["standard_relation"] == "=",
                    },
                }
                handle.write(json.dumps(record, ensure_ascii=False, separators=(",", ":")) + "\n")
                seen_activity_ids.add(activity_id)
                written += 1
                rows_written_for_endpoint += 1
                endpoint_counts[endpoint_name] += 1
                unit_counts[f"{endpoint_name}:{unit_canonical}"] += 1
                if rows_written_for_endpoint >= args.max_per_endpoint:
                    break

            counters[f"{endpoint_name}:rows_seen_for_endpoint"] = rows_seen_for_endpoint
            counters[f"{endpoint_name}:rows_written_for_endpoint"] = rows_written_for_endpoint

    con.close()
    stats = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": {"name": "ChEMBL", "version": "36", "db": str(args.db)},
        "output": str(args.out),
        "num_measurements": written,
        "endpoint_counts": dict(sorted(endpoint_counts.items())),
        "unit_counts": dict(sorted(unit_counts.items())),
        "counters": dict(sorted(counters.items())),
        "limitations": [
            "This table is an experimental measurement layer, not a query-answer benchmark by itself.",
            "Endpoint normalization is conservative; unsupported units are skipped.",
            "Downstream gold benchmark generation must still deduplicate replicates, bucket assay conditions, and apply leakage-safe splits.",
        ],
    }
    args.stats.write_text(json.dumps(stats, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return written, stats


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract normalized experimental ADMET-like measurements from ChEMBL SQLite.")
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("data/experimental/chembl_admet_measurements.jsonl"))
    parser.add_argument("--stats", type=Path, default=Path("data/experimental/chembl_admet_measurement_stats.json"))
    parser.add_argument("--max-per-endpoint", type=int, default=20000)
    parser.add_argument("--sql-limit-per-endpoint", type=int, default=200000)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.db.is_file():
        print(f"ERROR: ChEMBL SQLite DB not found: {args.db}", file=sys.stderr)
        return 1
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.stats.parent.mkdir(parents=True, exist_ok=True)
    written, stats = extract(args)
    if written == 0:
        print("ERROR: no measurements were extracted.", file=sys.stderr)
        return 1
    print(json.dumps(stats, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

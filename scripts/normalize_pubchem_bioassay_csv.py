from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import statistics
import zipfile
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

try:
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit.Contrib.SA_Score import sascorer
except Exception as exc:  # noqa: BLE001
    raise SystemExit(
        "ERROR: normalize_pubchem_bioassay_csv.py requires RDKit. "
        "Run it with .\\myenv311\\Scripts\\python.exe"
    ) from exc


FIELDS = [
    "observation_uid",
    "source",
    "source_version",
    "property_family",
    "endpoint_name",
    "label_type",
    "experimental_only_flag",
    "molecule_chembl_id",
    "smiles_canon",
    "inchikey",
    "connectivity_key",
    "murcko_scaffold",
    "condition_bucket",
    "value_raw",
    "unit_raw",
    "relation_raw",
    "type_raw",
    "value_canonical",
    "unit_canonical",
    "confidence",
    "target_chembl_id",
    "target_label",
    "assay_chembl_id",
    "assay_description",
    "assay_type",
    "assay_organism",
    "assay_cell_type",
    "document_year",
    "doi",
    "pubmed_id",
    "mw",
    "logp",
    "tpsa",
    "hba",
    "hbd",
    "rotatable_bonds",
    "heavy_atoms",
    "qed",
    "sa_score",
]

CONCENTRATION_TYPES = {
    "AC50",
    "CC50",
    "EC50",
    "GI50",
    "IC50",
    "IC90",
    "KD",
    "KI",
    "LC50",
    "POTENCY",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Normalize locally downloaded PubChem BioAssay CSV/Data zip shards into a separate observation table."
    )
    parser.add_argument(
        "--zip-dir",
        type=Path,
        default=Path("pos_way_admet_benchmark/raw/public/pubchem/bioassay_csv"),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/pubchem_bioassay"),
    )
    parser.add_argument("--source-version", default="pubchem_bioassay_csv_data_local")
    parser.add_argument("--max-rows-per-assay", type=int, default=0, help="0 means no cap.")
    parser.add_argument("--min-standard-value", type=float, default=1e-12)
    return parser.parse_args()


def as_float(value: Any) -> float | None:
    try:
        parsed = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def clean_text(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value or "").strip())


def aid_from_name(name: str) -> str:
    match = re.search(r"(\d+)\.csv\.gz$", name)
    return match.group(1) if match else ""


def canonical_type(raw: str) -> str:
    return clean_text(raw).upper().replace("-", "").replace("_", "")


def canonical_relation(raw: str) -> str:
    value = clean_text(raw)
    if value in {"=", "<", ">", "<=", ">="}:
        return value
    return value or ""


def canonical_value(std_type: str, value: float, unit: str) -> tuple[float, str] | None:
    unit_norm = clean_text(unit).lower()
    type_norm = canonical_type(std_type)
    if value <= 0:
        return None
    if type_norm in CONCENTRATION_TYPES and unit_norm in {"nm", "nanomolar"}:
        return round(9.0 - math.log10(value), 6), f"p{type_norm}"
    if type_norm in CONCENTRATION_TYPES and unit_norm in {"um", "micromolar", "Âµm"}:
        return round(6.0 - math.log10(value), 6), f"p{type_norm}"
    if type_norm in CONCENTRATION_TYPES and unit_norm in {"mm", "millimolar"}:
        return round(3.0 - math.log10(value), 6), f"p{type_norm}"
    return round(value, 6), unit or "unitless"


def descriptor_record(smiles: str, cache: dict[str, dict[str, Any] | None]) -> dict[str, Any] | None:
    if smiles in cache:
        return cache[smiles]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        cache[smiles] = None
        return None
    smiles_canon = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    inchikey = Chem.MolToInchiKey(mol)
    scaffold = ""
    scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold_mol is not None and scaffold_mol.GetNumHeavyAtoms() > 0:
        scaffold = Chem.MolToSmiles(scaffold_mol, canonical=True, isomericSmiles=False)
    record = {
        "smiles_canon": smiles_canon,
        "inchikey": inchikey,
        "connectivity_key": inchikey.split("-")[0] if inchikey else "",
        "murcko_scaffold": scaffold,
        "mw": round(float(Descriptors.MolWt(mol)), 4),
        "logp": round(float(Crippen.MolLogP(mol)), 4),
        "tpsa": round(float(rdMolDescriptors.CalcTPSA(mol)), 4),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "rotatable_bonds": int(Lipinski.NumRotatableBonds(mol)),
        "heavy_atoms": int(mol.GetNumHeavyAtoms()),
        "qed": round(float(QED.qed(mol)), 6),
        "sa_score": round(float(sascorer.calculateScore(mol)), 6),
    }
    cache[smiles] = record
    return record


def data_rows_from_csv_gz(payload: bytes) -> Iterable[dict[str, str]]:
    text = gzip.decompress(payload).decode("utf-8", errors="replace").splitlines()
    if not text:
        return
    reader = csv.DictReader(text)
    for row in reader:
        tag = clean_text(row.get("PUBCHEM_RESULT_TAG"))
        if not tag or tag.startswith("RESULT_"):
            continue
        yield row


def normalize_zip(path: Path, args: argparse.Namespace, desc_cache: dict[str, dict[str, Any] | None]) -> tuple[list[dict[str, Any]], Counter[str]]:
    output: list[dict[str, Any]] = []
    counters: Counter[str] = Counter()
    per_assay_counts: Counter[str] = Counter()
    with zipfile.ZipFile(path) as zf:
        for info in zf.infolist():
            if not info.filename.endswith(".csv.gz"):
                continue
            aid = aid_from_name(info.filename)
            counters["csv_gz_entries"] += 1
            payload = zf.read(info)
            for row in data_rows_from_csv_gz(payload):
                if args.max_rows_per_assay and per_assay_counts[aid] >= args.max_rows_per_assay:
                    counters["skipped:max_rows_per_assay"] += 1
                    continue
                cid = clean_text(row.get("PUBCHEM_CID"))
                sid = clean_text(row.get("PUBCHEM_SID"))
                smiles = clean_text(row.get("PUBCHEM_EXT_DATASOURCE_SMILES"))
                std_type_raw = clean_text(row.get("Standard Type"))
                std_type = canonical_type(std_type_raw)
                relation = canonical_relation(row.get("Standard Relation"))
                value = as_float(row.get("Standard Value"))
                unit = clean_text(row.get("Standard Units"))
                if value is None:
                    value = as_float(row.get("PubChem Standard Value"))
                    unit = unit or "uM"
                if not aid or not cid or not smiles:
                    counters["skipped:missing_aid_cid_or_smiles"] += 1
                    continue
                if not std_type or value is None:
                    counters["skipped:missing_standard_measurement"] += 1
                    continue
                if abs(value) < args.min_standard_value:
                    counters["skipped:near_zero_value"] += 1
                    continue
                canonical = canonical_value(std_type, value, unit)
                if canonical is None:
                    counters["skipped:invalid_canonical_value"] += 1
                    continue
                value_canonical, unit_canonical = canonical
                desc = descriptor_record(smiles, desc_cache)
                if desc is None:
                    counters["skipped:invalid_smiles"] += 1
                    continue
                relation_confidence = 0.8 if relation == "=" else 0.65 if relation in {"<", ">", "<=", ">="} else 0.55
                outcome = clean_text(row.get("PUBCHEM_ACTIVITY_OUTCOME"))
                endpoint_name = f"PubChem_AID{aid}_{unit_canonical}"
                condition_bucket = "|".join([endpoint_name, std_type or "type_unknown", f"AID:{aid}"])
                observation_uid = f"pubchem_aid{aid}_sid{sid or 'unknown'}_cid{cid}_{std_type or 'value'}"
                output.append(
                    {
                        "observation_uid": observation_uid,
                        "activity_outcome": outcome,
                        "pubchem_sid": sid,
                        "pubchem_cid": cid,
                        "source": "PubChem BioAssay",
                        "source_version": args.source_version,
                        "property_family": "pubchem_bioassay_activity",
                        "endpoint_name": endpoint_name,
                        "label_type": "experimental",
                        "experimental_only_flag": "true",
                        "molecule_chembl_id": "",
                        "smiles_canon": desc["smiles_canon"],
                        "inchikey": desc["inchikey"],
                        "connectivity_key": desc["connectivity_key"],
                        "murcko_scaffold": desc["murcko_scaffold"],
                        "condition_bucket": condition_bucket,
                        "value_raw": value,
                        "unit_raw": unit,
                        "relation_raw": relation,
                        "type_raw": std_type_raw or std_type,
                        "value_canonical": value_canonical,
                        "unit_canonical": unit_canonical,
                        "confidence": relation_confidence,
                        "target_chembl_id": "",
                        "target_label": "",
                        "assay_chembl_id": f"PubChem_AID:{aid}",
                        "assay_description": "",
                        "assay_type": "",
                        "assay_organism": "",
                        "assay_cell_type": "",
                        "document_year": "",
                        "doi": "",
                        "pubmed_id": "",
                        **{key: desc[key] for key in FIELDS if key in desc},
                    }
                )
                per_assay_counts[aid] += 1
                counters["kept"] += 1
                counters[f"type:{std_type or 'unknown'}"] += 1
                counters[f"outcome:{outcome or 'unknown'}"] += 1
    return output, counters


def write_csv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            count += 1
    return count


def write_property_summary(path: Path, rows: list[dict[str, Any]]) -> int:
    fields = [
        "endpoint_name",
        "observations",
        "unique_molecules",
        "unique_pubchem_cids",
        "unique_assays",
        "unique_condition_buckets",
        "unit_canonical",
        "top_type_raw",
        "exact_relation_rows",
        "inequality_relation_rows",
        "active_rows",
        "inactive_rows",
        "unspecified_rows",
        "value_min",
        "value_median",
        "value_max",
    ]
    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(str(row.get("endpoint_name", "")), []).append(row)
    summary_rows = []
    for endpoint, group in grouped.items():
        values = [as_float(row.get("value_canonical")) for row in group]
        values = [value for value in values if value is not None]
        type_counts = Counter(str(row.get("type_raw", "")) for row in group)
        outcomes = Counter(str(row.get("activity_outcome", "")) for row in group)
        cids = set()
        for row in group:
            match = re.search(r"_cid([^_]+)", str(row.get("observation_uid", "")))
            if match:
                cids.add(match.group(1))
        summary_rows.append(
            {
                "endpoint_name": endpoint,
                "observations": len(group),
                "unique_molecules": len({row.get("connectivity_key", "") for row in group if row.get("connectivity_key")}),
                "unique_pubchem_cids": len(cids),
                "unique_assays": len({row.get("assay_chembl_id", "") for row in group if row.get("assay_chembl_id")}),
                "unique_condition_buckets": len({row.get("condition_bucket", "") for row in group if row.get("condition_bucket")}),
                "unit_canonical": ";".join(sorted({str(row.get("unit_canonical", "")) for row in group if row.get("unit_canonical")})),
                "top_type_raw": type_counts.most_common(1)[0][0] if type_counts else "",
                "exact_relation_rows": sum(1 for row in group if row.get("relation_raw") == "="),
                "inequality_relation_rows": sum(1 for row in group if row.get("relation_raw") in {"<", ">", "<=", ">="}),
                "active_rows": outcomes.get("Active", 0),
                "inactive_rows": outcomes.get("Inactive", 0),
                "unspecified_rows": outcomes.get("Unspecified", 0),
                "value_min": min(values) if values else "",
                "value_median": statistics.median(values) if values else "",
                "value_max": max(values) if values else "",
            }
        )
    summary_rows.sort(key=lambda row: (-int(row["observations"]), row["endpoint_name"]))
    return write_csv(path, summary_rows, fields)


def main() -> int:
    args = parse_args()
    zip_paths = sorted(args.zip_dir.glob("*.zip"))
    if not zip_paths:
        raise SystemExit(f"ERROR: no PubChem BioAssay zip files found in {args.zip_dir}")

    rows: list[dict[str, Any]] = []
    counters: Counter[str] = Counter()
    desc_cache: dict[str, dict[str, Any] | None] = {}
    for path in zip_paths:
        counters["zip_files"] += 1
        parsed_rows, parsed_counters = normalize_zip(path, args, desc_cache)
        rows.extend(parsed_rows)
        counters.update(parsed_counters)

    out_csv = args.out_dir / "pubchem_bioassay_observations.csv"
    count = write_csv(out_csv, rows, FIELDS)
    summary_csv = args.out_dir / "pubchem_property_summary.csv"
    summary_count = write_property_summary(summary_csv, rows)
    unique_cids = {row["observation_uid"].split("_cid", 1)[-1].split("_", 1)[0] for row in rows}
    unique_molecules = {row["connectivity_key"] for row in rows if row.get("connectivity_key")}
    unique_endpoints = {row["endpoint_name"] for row in rows if row.get("endpoint_name")}
    condition_buckets = {row["condition_bucket"] for row in rows if row.get("condition_bucket")}
    stats = {
        "source": "PubChem BioAssay CSV/Data",
        "input_zip_dir": str(args.zip_dir),
        "output_csv": str(out_csv),
        "property_summary_csv": str(summary_csv),
        "rows": count,
        "summary_rows": summary_count,
        "unique_pubchem_cids": len(unique_cids),
        "unique_molecules_connectivity": len(unique_molecules),
        "unique_endpoints": len(unique_endpoints),
        "unique_condition_buckets": len(condition_buckets),
        "counters": dict(sorted(counters.items())),
        "notes": [
            "This is stored separately from the ChEMBL/BindingDB final dataset.",
            "Endpoint names are assay-local because PubChem CSV/Data shards do not include full target metadata.",
            "Exact concentration values are converted to p-scale when units are nM/uM/mM and type is concentration-like.",
        ],
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "pubchem_bioassay_stats.json").write_text(
        json.dumps(stats, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    (args.out_dir / "README.md").write_text(
        "# PubChem BioAssay Normalized Observations\n\n"
        "Separate PubChem processing output. These rows are not merged into the active ChEMBL-side dataset yet.\n\n"
        f"- Observations: `{out_csv.name}`\n"
        f"- Property summary: `{summary_csv.name}`\n"
        "- Schema: same columns as `data/normalized_csv_expanded/property_observations.csv`.\n"
        "- Source: locally downloaded PubChem BioAssay `CSV/Data/*.zip` shards.\n"
        "- Caveat: endpoints are assay-local until assay descriptions/targets are enriched.\n",
        encoding="utf-8",
    )
    print(json.dumps(stats, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


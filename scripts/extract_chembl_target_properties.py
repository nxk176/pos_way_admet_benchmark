from __future__ import annotations

import argparse
import json
import math
import sqlite3
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

try:
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors
    from rdkit.Contrib.SA_Score import sascorer
except Exception as exc:  # noqa: BLE001 - this extractor standardizes structures with RDKit.
    raise SystemExit(
        "ERROR: extract_chembl_target_properties.py requires RDKit. "
        "Run it with the project environment, for example .\\myenv311\\Scripts\\python.exe"
    ) from exc


TARGETS = {
    "DRD2_pChEMBL": {
        "target_chembl_id": "CHEMBL217",
        "target_label": "DRD2",
        "pref_name": "D(2) dopamine receptor",
        "preferred_direction": "decrease",
    },
    "GSK3B_pChEMBL": {
        "target_chembl_id": "CHEMBL262",
        "target_label": "GSK3B",
        "pref_name": "Glycogen synthase kinase-3 beta",
        "preferred_direction": "decrease",
    },
    "JNK3_pChEMBL": {
        "target_chembl_id": "CHEMBL2637",
        "target_label": "JNK3",
        "pref_name": "Mitogen-activated protein kinase 10",
        "preferred_direction": "decrease",
    },
}

SUPPORTED_TYPES = {"IC50", "EC50", "AC50", "Ki", "Kd", "Potency"}


QUERY = """
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
  td.tid,
  td.chembl_id as target_chembl_id,
  td.pref_name as target_pref_name,
  td.organism as target_organism,
  docs.year as document_year,
  docs.doi,
  docs.pubmed_id,
  docs.title as document_title
from target_dictionary td
join assays on assays.tid = td.tid
join activities a on a.assay_id = assays.assay_id
join molecule_dictionary md on a.molregno = md.molregno
join compound_structures cs on a.molregno = cs.molregno
left join docs on a.doc_id = docs.doc_id
where td.chembl_id = ?
  and a.standard_value is not null
  and a.pchembl_value is not null
  and a.standard_relation = '='
  and a.standard_flag = 1
  and a.standard_type in ({type_placeholders})
order by a.activity_id
limit ?
"""


def descriptor_dict(mol: Chem.Mol) -> dict[str, Any]:
    return {
        "mw": round(float(Descriptors.MolWt(mol)), 4),
        "logp": round(float(Crippen.MolLogP(mol)), 4),
        "tpsa": round(float(rdMolDescriptors.CalcTPSA(mol)), 4),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "rotatable_bonds": int(Lipinski.NumRotatableBonds(mol)),
        "heavy_atoms": int(mol.GetNumHeavyAtoms()),
        "qed": round(float(QED.qed(mol)), 4),
        "sa_score": round(float(sascorer.calculateScore(mol)), 4),
    }


def canon_smiles(smiles: str | None) -> tuple[str, Chem.Mol] | None:
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True), mol


def confidence(row: sqlite3.Row) -> float:
    score = 0.9
    if row["data_validity_comment"]:
        score -= 0.2
    if row["potential_duplicate"]:
        score -= 0.15
    return round(max(0.3, score), 3)


def query_for_types(activity_types: Iterable[str]) -> tuple[str, list[str]]:
    types = sorted(activity_types)
    placeholders = ",".join("?" for _ in types)
    return QUERY.format(type_placeholders=placeholders), types


def extract(args: argparse.Namespace) -> dict[str, Any]:
    con = sqlite3.connect(args.db)
    con.row_factory = sqlite3.Row
    sql, type_params = query_for_types(SUPPORTED_TYPES)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.stats.parent.mkdir(parents=True, exist_ok=True)

    counters: Counter[str] = Counter()
    endpoint_counts: Counter[str] = Counter()
    type_counts: Counter[str] = Counter()
    written = 0

    with args.out.open("w", encoding="utf-8", newline="\n") as handle:
        for endpoint_name, target in TARGETS.items():
            target_written = 0
            for row in con.execute(sql, [target["target_chembl_id"], *type_params, args.sql_limit_per_target]):
                counters[f"{endpoint_name}:sql_rows_seen"] += 1
                pchembl = row["pchembl_value"]
                if not isinstance(pchembl, (int, float)) or not math.isfinite(float(pchembl)):
                    counters[f"{endpoint_name}:invalid_pchembl"] += 1
                    continue
                parsed = canon_smiles(row["canonical_smiles"])
                if parsed is None:
                    counters[f"{endpoint_name}:invalid_smiles"] += 1
                    continue
                smiles_canon, mol = parsed
                record = {
                    "property_id": f"chembl36_target_act_{int(row['activity_id'])}",
                    "source": "ChEMBL",
                    "source_version": "36",
                    "endpoint_name": endpoint_name,
                    "family": "target_bioactivity",
                    "label_type": "experimental",
                    "experimental_only_flag": True,
                    "molecule_chembl_id": row["molecule_chembl_id"],
                    "smiles_raw": row["canonical_smiles"],
                    "smiles_canon": smiles_canon,
                    "inchikey": row["standard_inchi_key"],
                    "connectivity_key": str(row["standard_inchi_key"]).split("-", 1)[0],
                    "descriptors": descriptor_dict(mol),
                    "value_raw": row["standard_value"],
                    "unit_raw": row["standard_units"],
                    "relation_raw": row["standard_relation"],
                    "type_raw": row["standard_type"],
                    "value_canonical": round(float(pchembl), 6),
                    "unit_canonical": "pChEMBL",
                    "preferred_direction": target["preferred_direction"],
                    "confidence": confidence(row),
                    "target": {
                        "target_label": target["target_label"],
                        "target_chembl_id": row["target_chembl_id"],
                        "target_pref_name": row["target_pref_name"],
                        "target_organism": row["target_organism"],
                        "tid": row["tid"],
                    },
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
                        "standard_relation_exact": True,
                    },
                }
                handle.write(json.dumps(record, ensure_ascii=False, separators=(",", ":")) + "\n")
                written += 1
                target_written += 1
                endpoint_counts[endpoint_name] += 1
                type_counts[f"{endpoint_name}:{row['standard_type']}"] += 1
                if target_written >= args.max_per_target:
                    break

            counters[f"{endpoint_name}:rows_written"] = target_written

    con.close()
    stats = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": {"name": "ChEMBL", "version": "36", "db": str(args.db)},
        "output": str(args.out),
        "num_properties": written,
        "endpoint_counts": dict(sorted(endpoint_counts.items())),
        "activity_type_counts": dict(sorted(type_counts.items())),
        "targets": TARGETS,
        "counters": dict(sorted(counters.items())),
        "limitations": [
            "Target bioactivity properties are experimental ChEMBL pChEMBL values.",
            "Rows are not deduplicated into one label per molecule; downstream benchmark builders should aggregate by molecule, target, assay condition, and activity type.",
            "Lower pChEMBL is preferred for molecular editing tasks that reduce off-target activity.",
        ],
    }
    args.stats.write_text(json.dumps(stats, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return stats


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract DRD2, GSK3B, and JNK3 target bioactivity properties from ChEMBL.")
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path("data/properties/chembl_target_activity_properties.jsonl"))
    parser.add_argument("--stats", type=Path, default=Path("data/properties/chembl_target_activity_property_stats.json"))
    parser.add_argument("--max-per-target", type=int, default=20000)
    parser.add_argument("--sql-limit-per-target", type=int, default=100000)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.db.is_file():
        print(f"ERROR: ChEMBL SQLite DB not found: {args.db}", file=sys.stderr)
        return 1
    stats = extract(args)
    if stats["num_properties"] == 0:
        print("ERROR: no target properties were extracted.", file=sys.stderr)
        return 1
    print(json.dumps(stats, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

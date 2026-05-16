from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sqlite3
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit.Contrib.SA_Score import sascorer
except Exception as exc:  # noqa: BLE001
    raise SystemExit(
        "ERROR: extract_chembl_broad_pchembl_observations.py requires RDKit. "
        "Run it with the project environment, for example .\\myenv311\\Scripts\\python.exe"
    ) from exc


FIELDNAMES = [
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

SUPPORTED_TYPES = ("IC50", "Ki", "Kd", "EC50", "AC50", "Potency")
EXCLUDED_TARGET_NAMES = {"Unchecked", "NON-PROTEIN TARGET"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract a broad ChEMBL pChEMBL observation CSV.")
    parser.add_argument(
        "--db",
        type=Path,
        default=Path("pos_way_admet_benchmark/raw/public/chembl/chembl_36/chembl_36_sqlite/chembl_36.db"),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/broad/chembl_broad_pchembl_observations.csv"),
    )
    parser.add_argument(
        "--stats",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/broad/chembl_broad_pchembl_observation_stats.json"),
    )
    parser.add_argument("--max-targets", type=int, default=300)
    parser.add_argument("--min-rows-per-target", type=int, default=500)
    parser.add_argument("--max-rows-per-target", type=int, default=5000)
    parser.add_argument("--max-rows", type=int, default=750000)
    parser.add_argument("--min-target-confidence", type=int, default=5)
    parser.add_argument("--scan-activities", action="store_true", help="Fast path: scan indexed pChEMBL activities instead of preselecting top targets.")
    parser.add_argument("--skip-rdkit-descriptors", action="store_true")
    parser.add_argument("--compute-sa-score", action="store_true")
    return parser.parse_args()


def clean_text(value: Any) -> str:
    if value is None:
        return ""
    return str(value).replace("\r", " ").replace("\n", " ").strip()


def endpoint_name(target_chembl_id: str) -> str:
    return f"{target_chembl_id}_pChEMBL"


def connectivity_key(inchikey: str) -> str:
    return str(inchikey or "").split("-", 1)[0]


def confidence(row: sqlite3.Row) -> float:
    score = 0.9
    if row["data_validity_comment"]:
        score -= 0.2
    if row["potential_duplicate"]:
        score -= 0.15
    if row["assay_confidence_score"] is not None and int(row["assay_confidence_score"]) < 7:
        score -= 0.1
    return round(max(0.3, score), 3)


def descriptor_record(
    smiles: str,
    compute_sa_score: bool,
    skip_rdkit_descriptors: bool,
    cache: dict[str, dict[str, Any]],
) -> dict[str, Any] | None:
    if smiles in cache:
        return cache[smiles]
    if skip_rdkit_descriptors:
        record = {
            "smiles_canon": smiles,
            "murcko_scaffold": "",
            "mw": "",
            "logp": "",
            "tpsa": "",
            "hba": "",
            "hbd": "",
            "rotatable_bonds": "",
            "heavy_atoms": "",
            "qed": "",
            "sa_score": "",
        }
        cache[smiles] = record
        return record
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    scaffold = ""
    try:
        scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
        if scaffold_mol is not None and scaffold_mol.GetNumHeavyAtoms() > 0:
            scaffold = Chem.MolToSmiles(scaffold_mol, canonical=True, isomericSmiles=False)
    except Exception:  # noqa: BLE001
        scaffold = ""
    record = {
        "smiles_canon": Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
        "murcko_scaffold": scaffold,
        "mw": round(float(Descriptors.MolWt(mol)), 4),
        "logp": round(float(Crippen.MolLogP(mol)), 4),
        "tpsa": round(float(rdMolDescriptors.CalcTPSA(mol)), 4),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "rotatable_bonds": int(Lipinski.NumRotatableBonds(mol)),
        "heavy_atoms": int(mol.GetNumHeavyAtoms()),
        "qed": round(float(QED.qed(mol)), 4),
        "sa_score": round(float(sascorer.calculateScore(mol)), 4) if compute_sa_score else "",
    }
    cache[smiles] = record
    return record


def top_targets(con: sqlite3.Connection, args: argparse.Namespace) -> list[sqlite3.Row]:
    placeholders = ",".join("?" for _ in SUPPORTED_TYPES)
    excluded = ",".join("?" for _ in EXCLUDED_TARGET_NAMES)
    sql = f"""
    select
      td.tid,
      td.chembl_id as target_chembl_id,
      td.pref_name,
      td.target_type,
      td.organism,
      count(*) as row_count,
      count(distinct a.molregno) as molecule_count
    from activities a
    join assays ass on ass.assay_id = a.assay_id
    join target_dictionary td on td.tid = ass.tid
    join compound_structures cs on cs.molregno = a.molregno
    where a.pchembl_value is not null
      and a.standard_value is not null
      and a.standard_relation = '='
      and a.standard_flag = 1
      and a.standard_type in ({placeholders})
      and ass.confidence_score >= ?
      and td.pref_name not in ({excluded})
      and cs.canonical_smiles is not null
    group by td.tid, td.chembl_id, td.pref_name, td.target_type, td.organism
    having row_count >= ?
    order by row_count desc
    limit ?
    """
    params = [
        *SUPPORTED_TYPES,
        args.min_target_confidence,
        *sorted(EXCLUDED_TARGET_NAMES),
        args.min_rows_per_target,
        args.max_targets,
    ]
    return list(con.execute(sql, params))


def row_query(num_targets: int) -> str:
    type_placeholders = ",".join("?" for _ in SUPPORTED_TYPES)
    target_placeholders = ",".join("?" for _ in range(num_targets))
    return f"""
    select
      a.activity_id,
      a.standard_type,
      a.standard_relation,
      a.standard_value,
      a.standard_units,
      a.pchembl_value,
      a.data_validity_comment,
      a.potential_duplicate,
      md.chembl_id as molecule_chembl_id,
      cs.canonical_smiles,
      cs.standard_inchi_key,
      ass.chembl_id as assay_chembl_id,
      ass.description as assay_description,
      ass.assay_type,
      ass.assay_organism,
      ass.assay_cell_type,
      ass.confidence_score as assay_confidence_score,
      td.tid,
      td.chembl_id as target_chembl_id,
      td.pref_name as target_pref_name,
      docs.year as document_year,
      docs.doi,
      docs.pubmed_id
    from activities a
    join assays ass on ass.assay_id = a.assay_id
    join target_dictionary td on td.tid = ass.tid
    join molecule_dictionary md on md.molregno = a.molregno
    join compound_structures cs on cs.molregno = a.molregno
    left join docs on docs.doc_id = a.doc_id
    where a.assay_id in (
        select assay_id from assays
        where tid in ({target_placeholders})
          and confidence_score >= ?
    )
      and a.pchembl_value is not null
      and a.standard_value is not null
      and a.standard_relation = '='
      and a.standard_flag = 1
      and a.standard_type in ({type_placeholders})
      and cs.canonical_smiles is not null
    """


def scan_query() -> str:
    placeholders = ",".join("?" for _ in SUPPORTED_TYPES)
    excluded = ",".join("?" for _ in EXCLUDED_TARGET_NAMES)
    return f"""
    select
      a.activity_id,
      a.standard_type,
      a.standard_relation,
      a.standard_value,
      a.standard_units,
      a.pchembl_value,
      a.data_validity_comment,
      a.potential_duplicate,
      md.chembl_id as molecule_chembl_id,
      cs.canonical_smiles,
      cs.standard_inchi_key,
      ass.chembl_id as assay_chembl_id,
      ass.description as assay_description,
      ass.assay_type,
      ass.assay_organism,
      ass.assay_cell_type,
      ass.confidence_score as assay_confidence_score,
      td.tid,
      td.chembl_id as target_chembl_id,
      td.pref_name as target_pref_name,
      docs.year as document_year,
      docs.doi,
      docs.pubmed_id
    from activities a indexed by idx_act_pchembl
    join assays ass on ass.assay_id = a.assay_id
    join target_dictionary td on td.tid = ass.tid
    join molecule_dictionary md on md.molregno = a.molregno
    join compound_structures cs on cs.molregno = a.molregno
    left join docs on docs.doc_id = a.doc_id
    where a.pchembl_value is not null
      and a.standard_value is not null
      and a.standard_relation = '='
      and a.standard_flag = 1
      and a.standard_type in ({placeholders})
      and ass.confidence_score >= ?
      and td.pref_name not in ({excluded})
      and cs.canonical_smiles is not null
    limit ?
    """


def extract(args: argparse.Namespace) -> dict[str, Any]:
    con = sqlite3.connect(args.db)
    con.row_factory = sqlite3.Row
    targets = [] if args.scan_activities else top_targets(con, args)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.stats.parent.mkdir(parents=True, exist_ok=True)

    counters: Counter[str] = Counter()
    endpoint_counts: Counter[str] = Counter()
    descriptor_cache: dict[str, dict[str, Any]] = {}
    written = 0
    target_tids = [int(row["tid"]) for row in targets]
    target_by_tid = {int(row["tid"]): row for row in targets}

    with args.out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        if args.scan_activities:
            sql = scan_query()
            params = [*SUPPORTED_TYPES, args.min_target_confidence, *sorted(EXCLUDED_TARGET_NAMES), args.max_rows * 3]
            rows_iter = con.execute(sql, params)
        elif target_tids:
            sql = row_query(len(target_tids))
            params = [*target_tids, args.min_target_confidence, *SUPPORTED_TYPES]
            rows_iter = con.execute(sql, params)
        else:
            rows_iter = []
        for row in rows_iter:
            if written >= args.max_rows:
                break
            target = target_by_tid.get(int(row["tid"]))
            if target is None:
                endpoint = endpoint_name(row["target_chembl_id"])
            else:
                endpoint = endpoint_name(target["target_chembl_id"])
            if counters[f"{endpoint}:rows_written"] >= args.max_rows_per_target:
                continue
            pchembl = row["pchembl_value"]
            if not isinstance(pchembl, (int, float)) or not math.isfinite(float(pchembl)):
                counters[f"{endpoint}:invalid_pchembl"] += 1
                continue
            desc = descriptor_record(
                row["canonical_smiles"],
                args.compute_sa_score,
                args.skip_rdkit_descriptors,
                descriptor_cache,
            )
            if desc is None:
                counters[f"{endpoint}:invalid_smiles"] += 1
                continue
            condition_bucket = "|".join(
                [
                    endpoint,
                    clean_text(row["standard_type"]),
                    clean_text(row["assay_chembl_id"]) or "assay_unknown",
                ]
            )
            out = {
                "observation_uid": f"chembl36_act_{int(row['activity_id'])}",
                "source": "ChEMBL",
                "source_version": "36",
                "property_family": "broad_target_bioactivity",
                "endpoint_name": endpoint,
                "label_type": "experimental",
                "experimental_only_flag": "true",
                "molecule_chembl_id": row["molecule_chembl_id"],
                "smiles_canon": desc["smiles_canon"],
                "inchikey": row["standard_inchi_key"],
                "connectivity_key": connectivity_key(row["standard_inchi_key"]),
                "murcko_scaffold": desc["murcko_scaffold"],
                "condition_bucket": condition_bucket,
                "value_raw": row["standard_value"],
                "unit_raw": row["standard_units"],
                "relation_raw": row["standard_relation"],
                "type_raw": row["standard_type"],
                "value_canonical": round(float(pchembl), 6),
                "unit_canonical": "pChEMBL",
                "confidence": confidence(row),
                "target_chembl_id": row["target_chembl_id"],
                "target_label": clean_text(row["target_pref_name"]),
                "assay_chembl_id": row["assay_chembl_id"],
                "assay_description": clean_text(row["assay_description"]),
                "assay_type": row["assay_type"],
                "assay_organism": row["assay_organism"],
                "assay_cell_type": row["assay_cell_type"],
                "document_year": row["document_year"],
                "doi": row["doi"],
                "pubmed_id": row["pubmed_id"],
                **{key: desc.get(key, "") for key in FIELDNAMES if key in desc},
            }
            writer.writerow(out)
            written += 1
            counters[f"{endpoint}:rows_written"] += 1
            endpoint_counts[endpoint] += 1

    stats = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": {"name": "ChEMBL", "version": "36", "db": str(args.db)},
        "output": str(args.out),
        "num_observations": written,
        "num_targets_selected": len(targets),
        "num_endpoints_written": len(endpoint_counts),
        "parameters": {
            "max_targets": args.max_targets,
            "min_rows_per_target": args.min_rows_per_target,
            "max_rows_per_target": args.max_rows_per_target,
            "max_rows": args.max_rows,
            "min_target_confidence": args.min_target_confidence,
            "scan_activities": args.scan_activities,
            "skip_rdkit_descriptors": args.skip_rdkit_descriptors,
            "supported_types": list(SUPPORTED_TYPES),
        },
        "top_targets": [
            {
                "target_chembl_id": row["target_chembl_id"],
                "pref_name": row["pref_name"],
                "target_type": row["target_type"],
                "organism": row["organism"],
                "source_row_count": row["row_count"],
                "source_molecule_count": row["molecule_count"],
            }
            for row in targets[:50]
        ],
        "endpoint_counts_top50": dict(endpoint_counts.most_common(50)),
        "counters": dict(sorted(counters.items())),
    }
    args.stats.write_text(json.dumps(stats, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    con.close()
    return stats


def main() -> int:
    args = parse_args()
    stats = extract(args)
    print(json.dumps({k: stats[k] for k in ["num_observations", "num_targets_selected", "num_endpoints_written"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

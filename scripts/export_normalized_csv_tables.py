from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

try:
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
except Exception as exc:  # noqa: BLE001 - this script needs RDKit for scaffold statistics.
    raise SystemExit(
        "ERROR: export_normalized_csv_tables.py requires RDKit. "
        "Run it with the project environment, for example .\\myenv311\\Scripts\\python.exe"
    ) from exc


DESCRIPTOR_FIELDS = ["mw", "logp", "tpsa", "hba", "hbd", "rotatable_bonds", "heavy_atoms", "qed", "sa_score"]
RDKIT_PROPERTY_ROWS = [
    ("RDKit_MW", "mw", "Da", "proxy_physchem"),
    ("RDKit_MolLogP", "logp", "unitless", "proxy_physchem"),
    ("RDKit_TPSA", "tpsa", "A2", "proxy_physchem"),
    ("RDKit_HAcceptors", "hba", "count", "proxy_physchem"),
    ("RDKit_HDonors", "hbd", "count", "proxy_physchem"),
    ("RDKit_RotBonds", "rotatable_bonds", "count", "proxy_physchem"),
    ("RDKit_QED", "qed", "unitless", "proxy_druglikeness"),
    ("RDKit_SA_Score", "sa_score", "unitless", "proxy_synthesizability"),
]


def load_jsonl(path: Path) -> Iterable[dict[str, Any]]:
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line:
                yield json.loads(line)


def write_csv(path: Path, rows: Iterable[dict[str, Any]], fieldnames: list[str]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: clean_cell(row.get(key)) for key in fieldnames})
            count += 1
    return count


def clean_cell(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float, str)):
        return value
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def slug(value: Any, fallback: str = "unknown") -> str:
    cleaned = re.sub(r"[^a-z0-9]+", "_", str(value or "").strip().lower()).strip("_")
    return cleaned[:64] if cleaned else fallback


def scaffold_from_smiles(smiles: str, cache: dict[str, str]) -> str:
    if smiles in cache:
        return cache[smiles]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        cache[smiles] = ""
        return ""
    scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold_mol is None or scaffold_mol.GetNumHeavyAtoms() == 0:
        cache[smiles] = ""
    else:
        cache[smiles] = Chem.MolToSmiles(scaffold_mol, canonical=True, isomericSmiles=False)
    return cache[smiles]


def measurement_condition_bucket(row: dict[str, Any]) -> str:
    endpoint = row["endpoint_name"]
    assay = row.get("assay") or {}
    description = str(assay.get("description") or "").lower()
    organism = slug(assay.get("assay_organism"))
    cell_type = slug(assay.get("assay_cell_type"))
    if endpoint == "logS_mol_L":
        ph = re.search(r"ph\s*([0-9]+(?:\.[0-9]+)?)", description)
        return f"{endpoint}|ph_{ph.group(1)}" if ph else f"{endpoint}|ph_unknown"
    if endpoint == "Caco2_logPapp_cm_s":
        if "a to b" in description or "a-b" in description:
            direction = "a_to_b"
        elif "b to a" in description or "b-a" in description:
            direction = "b_to_a"
        else:
            direction = "bidirectional"
        return f"{endpoint}|{direction}"
    if endpoint == "hERG_pIC50":
        return f"{endpoint}|{organism}|{cell_type}"
    if endpoint == "microsomal_clearance_mL_min_kg":
        matrix = "microsome" if "microsom" in description else "clearance"
        return f"{endpoint}|{organism}|{matrix}"
    if endpoint == "half_life_min":
        matrix = "microsome" if "microsom" in description else "plasma_or_systemic"
        return f"{endpoint}|{organism}|{matrix}"
    return f"{endpoint}|generic"


def target_condition_bucket(row: dict[str, Any]) -> str:
    target = row.get("target") or {}
    assay = row.get("assay") or {}
    return "|".join(
        [
            str(row["endpoint_name"]),
            slug(row.get("type_raw")),
            slug(target.get("target_chembl_id")),
            slug(assay.get("assay_organism")),
            slug(assay.get("assay_cell_type")),
        ]
    )


def descriptor_columns(desc: dict[str, Any] | None) -> dict[str, Any]:
    desc = desc or {}
    return {field: desc.get(field, "") for field in DESCRIPTOR_FIELDS}


def molecule_rows(root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for row in load_jsonl(root / "data/derived/chembl_molecules.jsonl"):
        desc = row.get("descriptors") or {}
        rows.append(
            {
                "molecule_uid": row.get("inchikey") or row.get("connectivity_key") or row.get("chembl_id"),
                "chembl_id": row.get("chembl_id"),
                "smiles_raw": row.get("smiles_raw"),
                "smiles_canon": row.get("smiles_canon"),
                "inchikey": row.get("inchikey"),
                "connectivity_key": row.get("connectivity_key"),
                "murcko_scaffold": row.get("murcko_scaffold"),
                "source": row.get("source"),
                **descriptor_columns(desc),
            }
        )
    return rows


def observation_rows(root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    scaffold_cache: dict[str, str] = {}

    for row in load_jsonl(root / "data/experimental/chembl_admet_measurements.jsonl"):
        assay = row.get("assay") or {}
        document = row.get("document") or {}
        smiles = row.get("smiles_canon") or row.get("smiles_raw") or ""
        rows.append(
            {
                "observation_uid": row.get("measurement_id"),
                "source": row.get("source"),
                "source_version": row.get("source_version"),
                "property_family": "admet_experimental",
                "endpoint_name": row.get("endpoint_name"),
                "label_type": row.get("label_type"),
                "experimental_only_flag": row.get("experimental_only_flag"),
                "molecule_chembl_id": row.get("molecule_chembl_id"),
                "smiles_canon": smiles,
                "inchikey": row.get("inchikey"),
                "connectivity_key": row.get("connectivity_key"),
                "murcko_scaffold": scaffold_from_smiles(smiles, scaffold_cache),
                "condition_bucket": measurement_condition_bucket(row),
                "value_raw": row.get("value_raw"),
                "unit_raw": row.get("unit_raw"),
                "relation_raw": row.get("relation_raw"),
                "type_raw": row.get("type_raw"),
                "value_canonical": row.get("value_canonical"),
                "unit_canonical": row.get("unit_canonical"),
                "confidence": row.get("confidence"),
                "target_chembl_id": "",
                "target_label": "",
                "assay_chembl_id": assay.get("assay_chembl_id"),
                "assay_description": assay.get("description"),
                "assay_type": assay.get("assay_type"),
                "assay_organism": assay.get("assay_organism"),
                "assay_cell_type": assay.get("assay_cell_type"),
                "document_year": document.get("year"),
                "doi": document.get("doi"),
                "pubmed_id": document.get("pubmed_id"),
                **descriptor_columns(None),
            }
        )

    for row in load_jsonl(root / "data/properties/chembl_target_activity_properties.jsonl"):
        assay = row.get("assay") or {}
        target = row.get("target") or {}
        document = row.get("document") or {}
        smiles = row.get("smiles_canon") or row.get("smiles_raw") or ""
        rows.append(
            {
                "observation_uid": row.get("property_id"),
                "source": row.get("source"),
                "source_version": row.get("source_version"),
                "property_family": row.get("family") or "target_bioactivity",
                "endpoint_name": row.get("endpoint_name"),
                "label_type": row.get("label_type"),
                "experimental_only_flag": row.get("experimental_only_flag"),
                "molecule_chembl_id": row.get("molecule_chembl_id"),
                "smiles_canon": smiles,
                "inchikey": row.get("inchikey"),
                "connectivity_key": row.get("connectivity_key"),
                "murcko_scaffold": scaffold_from_smiles(smiles, scaffold_cache),
                "condition_bucket": target_condition_bucket(row),
                "value_raw": row.get("value_raw"),
                "unit_raw": row.get("unit_raw"),
                "relation_raw": row.get("relation_raw"),
                "type_raw": row.get("type_raw"),
                "value_canonical": row.get("value_canonical"),
                "unit_canonical": row.get("unit_canonical"),
                "confidence": row.get("confidence"),
                "target_chembl_id": target.get("target_chembl_id"),
                "target_label": target.get("target_label"),
                "assay_chembl_id": assay.get("assay_chembl_id"),
                "assay_description": assay.get("description"),
                "assay_type": assay.get("assay_type"),
                "assay_organism": assay.get("assay_organism"),
                "assay_cell_type": assay.get("assay_cell_type"),
                "document_year": document.get("year"),
                "doi": document.get("doi"),
                "pubmed_id": document.get("pubmed_id"),
                **descriptor_columns(row.get("descriptors")),
            }
        )
    return rows


def flat_query_answer_rows(root: Path, qpath: str, apath: str) -> list[dict[str, Any]]:
    queries = {row["query_id"]: row for row in load_jsonl(root / qpath)}
    rows: list[dict[str, Any]] = []
    for answer in load_jsonl(root / apath):
        query = queries[answer["query_id"]]
        endpoint = (query.get("target_endpoints") or [{}])[0]
        transform = answer.get("transform_detail") or {}
        prov = answer.get("provenance") or {}
        rows.append(
            {
                "row_uid": answer.get("answer_id"),
                "query_id": answer.get("query_id"),
                "answer_id": answer.get("answer_id"),
                "split": query.get("split"),
                "source_pool": query.get("source_pool"),
                "question_template": query.get("question_template"),
                "question_text": query.get("question_text"),
                "endpoint_name": answer.get("endpoint_name") or endpoint.get("endpoint_name"),
                "direction": endpoint.get("direction"),
                "input_chembl_id": query.get("source_chembl_id"),
                "input_smiles_canon": query.get("input_smiles_canon"),
                "input_inchikey": query.get("input_inchikey"),
                "input_connectivity_key": query.get("input_connectivity_key"),
                "input_murcko_scaffold": query.get("input_murcko_scaffold"),
                "target_chembl_id": transform.get("target_chembl_id"),
                "target_smiles_canon": answer.get("target_smiles_canon"),
                "target_inchikey": answer.get("target_inchikey"),
                "target_connectivity_key": answer.get("target_connectivity_key"),
                "value_before": answer.get("value_before"),
                "value_after": answer.get("value_after"),
                "delta_value": answer.get("delta_value"),
                "unit_canonical": answer.get("unit_canonical"),
                "label_type": answer.get("label_type"),
                "experimental_only_flag": answer.get("experimental_only_flag"),
                "confidence": answer.get("confidence"),
                "transform_class": answer.get("transform_class"),
                "condition_bucket": transform.get("condition_bucket") or prov.get("condition_bucket"),
                "source_measurement_ids": prov.get("source_measurement_ids"),
                "target_measurement_ids": prov.get("target_measurement_ids"),
                "constraint_flags": answer.get("constraint_flags"),
            }
        )
    return rows


def median(values: list[float]) -> float | None:
    values = [value for value in values if math.isfinite(value)]
    return float(statistics.median(values)) if values else None


def property_summary_rows(molecules: list[dict[str, Any]], observations: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    by_endpoint: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in observations:
        by_endpoint[str(row["endpoint_name"])].append(row)

    for endpoint, group in sorted(by_endpoint.items()):
        values = [float(row["value_canonical"]) for row in group if str(row.get("value_canonical", "")) not in {"", "nan"}]
        rows.append(
            {
                "endpoint_name": endpoint,
                "source_tier": "experimental",
                "property_family": Counter(str(row.get("property_family")) for row in group).most_common(1)[0][0],
                "sample_count": len(group),
                "unique_molecules": len({row.get("molecule_chembl_id") for row in group if row.get("molecule_chembl_id")}),
                "unique_connectivity_keys": len({row.get("connectivity_key") for row in group if row.get("connectivity_key")}),
                "unique_scaffolds": len({row.get("murcko_scaffold") for row in group if row.get("murcko_scaffold")}),
                "unique_condition_buckets": len({row.get("condition_bucket") for row in group if row.get("condition_bucket")}),
                "unit_canonical": "|".join(sorted({str(row.get("unit_canonical")) for row in group if row.get("unit_canonical")})),
                "median_value": median(values),
                "min_value": min(values) if values else "",
                "max_value": max(values) if values else "",
                "recommendation": recommendation(endpoint, len(group), "experimental"),
            }
        )

    for endpoint, desc_field, unit, family in RDKIT_PROPERTY_ROWS:
        values = [float(row[desc_field]) for row in molecules if str(row.get(desc_field, "")) not in {"", "nan"}]
        rows.append(
            {
                "endpoint_name": endpoint,
                "source_tier": "proxy",
                "property_family": family,
                "sample_count": len(values),
                "unique_molecules": len(molecules),
                "unique_connectivity_keys": len({row.get("connectivity_key") for row in molecules if row.get("connectivity_key")}),
                "unique_scaffolds": len({row.get("murcko_scaffold") for row in molecules if row.get("murcko_scaffold")}),
                "unique_condition_buckets": 0,
                "unit_canonical": unit,
                "median_value": median(values),
                "min_value": min(values) if values else "",
                "max_value": max(values) if values else "",
                "recommendation": recommendation(endpoint, len(values), "proxy"),
            }
        )
    return rows


def recommendation(endpoint: str, sample_count: int, tier: str) -> str:
    if tier == "proxy":
        return "Use as secondary constraint or silver training signal"
    if sample_count >= 10000:
        return "Strong candidate for benchmark property"
    if sample_count >= 3000:
        return "Usable candidate; inspect assay buckets"
    return "Lower coverage; use selectively or after expansion"


def write_property_selection_markdown(path: Path, summary_rows: list[dict[str, Any]]) -> None:
    ordered = sorted(summary_rows, key=lambda row: (row["source_tier"] != "experimental", -int(row["sample_count"]), row["endpoint_name"]))
    lines = [
        "# Property Selection Table",
        "",
        "This table is generated from normalized CSV exports and is intended for deciding which properties to target next.",
        "",
        "| Endpoint | Tier | Samples | Molecules | Scaffolds | Buckets | Unit | Recommendation |",
        "|---|---|---:|---:|---:|---:|---|---|",
    ]
    for row in ordered:
        lines.append(
            "| {endpoint_name} | {source_tier} | {sample_count} | {unique_molecules} | {unique_scaffolds} | "
            "{unique_condition_buckets} | {unit_canonical} | {recommendation} |".format(**row)
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export normalized CSV tables and property statistics.")
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--out-dir", type=Path, default=Path("data/normalized_csv"))
    parser.add_argument("--property-md", type=Path, default=Path("PROPERTY_SELECTION_TABLE.md"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.root
    out_dir = root / args.out_dir

    molecules = molecule_rows(root)
    observations = observation_rows(root)
    silver_pairs = flat_query_answer_rows(root, "data/queries.jsonl", "data/answers.jsonl")
    gold_pairs = flat_query_answer_rows(root, "data/gold/queries.jsonl", "data/gold/answers.jsonl")
    summary = property_summary_rows(molecules, observations)

    molecule_fields = [
        "molecule_uid",
        "chembl_id",
        "smiles_raw",
        "smiles_canon",
        "inchikey",
        "connectivity_key",
        "murcko_scaffold",
        "source",
        *DESCRIPTOR_FIELDS,
    ]
    observation_fields = [
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
        *DESCRIPTOR_FIELDS,
    ]
    pair_fields = [
        "row_uid",
        "query_id",
        "answer_id",
        "split",
        "source_pool",
        "question_template",
        "question_text",
        "endpoint_name",
        "direction",
        "input_chembl_id",
        "input_smiles_canon",
        "input_inchikey",
        "input_connectivity_key",
        "input_murcko_scaffold",
        "target_chembl_id",
        "target_smiles_canon",
        "target_inchikey",
        "target_connectivity_key",
        "value_before",
        "value_after",
        "delta_value",
        "unit_canonical",
        "label_type",
        "experimental_only_flag",
        "confidence",
        "transform_class",
        "condition_bucket",
        "source_measurement_ids",
        "target_measurement_ids",
        "constraint_flags",
    ]
    summary_fields = [
        "endpoint_name",
        "source_tier",
        "property_family",
        "sample_count",
        "unique_molecules",
        "unique_connectivity_keys",
        "unique_scaffolds",
        "unique_condition_buckets",
        "unit_canonical",
        "median_value",
        "min_value",
        "max_value",
        "recommendation",
    ]

    counts = {
        "molecules": write_csv(out_dir / "molecules.csv", molecules, molecule_fields),
        "property_observations": write_csv(out_dir / "property_observations.csv", observations, observation_fields),
        "silver_query_answer_pairs": write_csv(out_dir / "silver_query_answer_pairs.csv", silver_pairs, pair_fields),
        "gold_query_answer_pairs": write_csv(out_dir / "gold_query_answer_pairs.csv", gold_pairs, pair_fields),
        "property_summary": write_csv(out_dir / "property_summary.csv", summary, summary_fields),
    }
    write_property_selection_markdown(root / args.property_md, summary)
    print(json.dumps({"out_dir": str(out_dir), "counts": counts}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

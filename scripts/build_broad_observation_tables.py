from __future__ import annotations

import argparse
import csv
import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build broad normalized observation CSV and property summary.")
    parser.add_argument(
        "--base-observations",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/normalized_csv/property_observations.csv"),
    )
    parser.add_argument(
        "--broad-observations",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/broad/chembl_broad_pchembl_observations.csv"),
    )
    parser.add_argument(
        "--extra-observations",
        type=Path,
        action="append",
        default=[],
        help="Additional normalized observation CSVs to merge into the output table.",
    )
    parser.add_argument("--out-dir", type=Path, default=Path("pos_way_admet_benchmark/data/normalized_csv_broad"))
    return parser.parse_args()


def clean_cell(value: Any) -> Any:
    if value is None:
        return ""
    return value


def as_float(value: str) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed


def stream_csv(path: Path):
    if not path.is_file():
        return
    with path.open("r", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle)


def recommendation(row: dict[str, Any]) -> str:
    tier = row["source_tier"]
    samples = row["sample_count"]
    molecules = row["unique_molecules"]
    buckets = row["unique_condition_buckets"]
    if tier == "proxy":
        return "Use as secondary constraint or silver training signal"
    if samples >= 5000 and molecules >= 2500 and buckets >= 50:
        return "Strong broad candidate; inspect assay families"
    if samples >= 1000 and molecules >= 700:
        return "Usable broad candidate; condition filtering recommended"
    return "Lower coverage; use selectively"


def build(args: argparse.Namespace) -> dict[str, Any]:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_observations = args.out_dir / "property_observations.csv"
    out_summary = args.out_dir / "property_summary.csv"
    out_status = args.out_dir / "broad_build_stats.json"

    seen: set[str] = set()
    endpoint_counts: Counter[str] = Counter()
    endpoint_family: dict[str, str] = {}
    endpoint_tier: dict[str, str] = {}
    endpoint_unit: dict[str, str] = {}
    molecules: dict[str, set[str]] = defaultdict(set)
    connectivity: dict[str, set[str]] = defaultdict(set)
    scaffolds: dict[str, set[str]] = defaultdict(set)
    buckets: dict[str, set[str]] = defaultdict(set)
    values: dict[str, list[float]] = defaultdict(list)
    rows_written = 0
    source_counts: Counter[str] = Counter()

    with out_observations.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        for source_path in [args.base_observations, args.broad_observations, *args.extra_observations]:
            for row in stream_csv(source_path):
                uid = row.get("observation_uid", "")
                if not uid or uid in seen:
                    continue
                seen.add(uid)
                normalized = {field: clean_cell(row.get(field, "")) for field in FIELDNAMES}
                writer.writerow(normalized)
                rows_written += 1

                endpoint = normalized["endpoint_name"]
                endpoint_counts[endpoint] += 1
                source_counts[str(source_path)] += 1
                endpoint_family.setdefault(endpoint, normalized["property_family"])
                endpoint_tier.setdefault(endpoint, "proxy" if normalized["label_type"] == "proxy" else "experimental")
                endpoint_unit.setdefault(endpoint, normalized["unit_canonical"])
                if normalized["molecule_chembl_id"]:
                    molecules[endpoint].add(normalized["molecule_chembl_id"])
                if normalized["connectivity_key"]:
                    connectivity[endpoint].add(normalized["connectivity_key"])
                if normalized["murcko_scaffold"]:
                    scaffolds[endpoint].add(normalized["murcko_scaffold"])
                if normalized["condition_bucket"]:
                    buckets[endpoint].add(normalized["condition_bucket"])
                value = as_float(normalized["value_canonical"])
                if value is not None:
                    values[endpoint].append(value)

    summary_rows = []
    for endpoint, count in endpoint_counts.items():
        endpoint_values = values.get(endpoint, [])
        summary_rows.append(
            {
                "endpoint_name": endpoint,
                "source_tier": endpoint_tier.get(endpoint, "experimental"),
                "property_family": endpoint_family.get(endpoint, ""),
                "sample_count": count,
                "unique_molecules": len(molecules[endpoint]),
                "unique_connectivity_keys": len(connectivity[endpoint]),
                "unique_scaffolds": len(scaffolds[endpoint]),
                "unique_condition_buckets": len(buckets[endpoint]),
                "unit_canonical": endpoint_unit.get(endpoint, ""),
                "median_value": round(float(statistics.median(endpoint_values)), 6) if endpoint_values else "",
                "min_value": round(float(min(endpoint_values)), 6) if endpoint_values else "",
                "max_value": round(float(max(endpoint_values)), 6) if endpoint_values else "",
            }
        )
    for row in summary_rows:
        row["recommendation"] = recommendation(row)
    summary_rows.sort(key=lambda row: (0 if row["source_tier"] == "experimental" else 1, -row["sample_count"], row["endpoint_name"]))

    with out_summary.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
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
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)

    stats = {
        "property_observations": str(out_observations),
        "property_summary": str(out_summary),
        "rows_written": rows_written,
        "unique_observation_uids": len(seen),
        "num_endpoints": len(endpoint_counts),
        "source_counts": dict(source_counts),
        "top_endpoints": summary_rows[:30],
    }
    out_status.write_text(json.dumps(stats, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return stats


def main() -> int:
    args = parse_args()
    stats = build(args)
    print(json.dumps({k: stats[k] for k in ["rows_written", "num_endpoints"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Merge PubChem supported-observation CSV parts with observation_uid dedupe.")
    parser.add_argument(
        "--inputs",
        nargs="+",
        type=Path,
        default=[
            Path("pos_way_admet_benchmark/data/pubchem_bioassay/pubchem_final_supported_observations.csv"),
            Path("pos_way_admet_benchmark/data/pubchem_bioassay/pubchem_final_supported_observations_part3.csv"),
        ],
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/pubchem_bioassay/pubchem_final_supported_observations_merged.csv"),
    )
    parser.add_argument(
        "--stats",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/pubchem_bioassay/pubchem_final_supported_merged_stats.json"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    seen: set[str] = set()
    fieldnames: list[str] | None = None
    input_stats = []
    total_rows = 0
    duplicate_rows = 0
    unique_endpoints: set[str] = set()
    unique_molecules: set[str] = set()

    with args.out_csv.open("w", encoding="utf-8", newline="") as out_handle:
        writer = None
        for path in args.inputs:
            rows_in_file = 0
            kept_in_file = 0
            dup_in_file = 0
            with path.open("r", encoding="utf-8", newline="") as in_handle:
                reader = csv.DictReader(in_handle)
                if fieldnames is None:
                    fieldnames = list(reader.fieldnames or [])
                    writer = csv.DictWriter(out_handle, fieldnames=fieldnames, extrasaction="ignore")
                    writer.writeheader()
                elif list(reader.fieldnames or []) != fieldnames:
                    raise SystemExit(f"ERROR: schema mismatch in {path}")
                assert writer is not None
                for row in reader:
                    rows_in_file += 1
                    total_rows += 1
                    uid = row.get("observation_uid") or ""
                    if not uid:
                        uid = "|".join(
                            [
                                row.get("endpoint_name", ""),
                                row.get("connectivity_key", ""),
                                row.get("value_canonical", ""),
                                row.get("condition_bucket", ""),
                            ]
                        )
                    if uid in seen:
                        duplicate_rows += 1
                        dup_in_file += 1
                        continue
                    seen.add(uid)
                    writer.writerow(row)
                    kept_in_file += 1
                    if row.get("endpoint_name"):
                        unique_endpoints.add(row["endpoint_name"])
                    if row.get("connectivity_key"):
                        unique_molecules.add(row["connectivity_key"])
            input_stats.append(
                {
                    "path": str(path),
                    "rows": rows_in_file,
                    "kept": kept_in_file,
                    "duplicates": dup_in_file,
                }
            )

    stats = {
        "output_csv": str(args.out_csv),
        "input_files": input_stats,
        "input_rows": total_rows,
        "deduplicated_rows": len(seen),
        "duplicate_rows": duplicate_rows,
        "unique_molecules": len(unique_molecules),
        "unique_endpoints": len(unique_endpoints),
        "dedupe_key": "observation_uid",
    }
    args.stats.write_text(json.dumps(stats, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(stats, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any


FIELDS = [
    "query_id",
    "split",
    "instruction",
    "input_smiles_canon",
    "input_connectivity_key",
    "input_chembl_id",
    "primary_endpoint",
    "condition_bucket",
    "primary_direction",
    "num_property_objectives",
    "primary_objective_json",
    "preserved_property_json",
    "local_constraints_json",
    "num_positive_answers",
    "positive_answer_smiles_json",
    "positive_answers_json",
    "mmp_cores_json",
    "source_positive_sample_ids_json",
    "selection_rule",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build ChEMBL-side rows with exactly 3 property objective groups and exactly 2 positive answers."
    )
    parser.add_argument("--input-dir", type=Path, default=Path("pos_way_admet_benchmark/data/final_multianswer_multiprop"))
    parser.add_argument("--out-dir", type=Path, default=Path("pos_way_admet_benchmark/data/chembl_3prop_2pos"))
    return parser.parse_args()


def safe_json(raw: str, fallback: Any) -> Any:
    try:
        return json.loads(raw)
    except (TypeError, json.JSONDecodeError):
        return fallback


def compact_json(payload: Any) -> str:
    return json.dumps(payload, ensure_ascii=False, separators=(",", ":"))


def as_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("-inf")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def normalize_instruction(row: dict[str, str], primary: dict[str, Any], preserved: dict[str, Any]) -> str:
    direction = row.get("primary_direction") or "increase"
    verb = "Decrease" if direction == "decrease" else "Increase"
    endpoint = primary.get("endpoint_name") or row.get("primary_endpoint", "")
    preserved_name = preserved.get("endpoint_name", "")
    tolerance_text = preserved.get("tolerance_text", "the specified tolerance")
    return (
        f"{verb} {endpoint} while preserving {preserved_name} within {tolerance_text}; "
        "also keep MW, LogP, QED, and synthetic accessibility within local edit constraints."
    )


def choose_rows(source_rows: list[dict[str, str]], prefix: str) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    output: list[dict[str, Any]] = []
    skipped = Counter()
    for row in source_rows:
        positives = safe_json(row.get("positive_answers_json", "[]"), [])
        objectives = safe_json(row.get("property_objectives_json", "{}"), {})
        secondary = objectives.get("secondary", []) if isinstance(objectives, dict) else []
        primary = objectives.get("primary", {}) if isinstance(objectives, dict) else {}
        if len(positives) < 2:
            skipped["less_than_2_positive"] += 1
            continue
        if not secondary:
            skipped["no_secondary_property"] += 1
            continue
        preserved = secondary[0]
        selected_positives = sorted(
            positives,
            key=lambda item: (-as_float(item.get("delta_value")), item.get("target_smiles_canon", "")),
        )[:2]
        local_constraints = {
            "objective_group": "local_constraints",
            "properties": ["MW", "LogP", "QED", "SA"],
            "rule": "Keep within local edit constraints recorded during MMP candidate generation.",
        }
        output.append(
            {
                "query_id": f"{prefix}_{len(output) + 1:08d}",
                "split": row.get("split", ""),
                "instruction": normalize_instruction(row, primary, preserved),
                "input_smiles_canon": row.get("input_smiles_canon", ""),
                "input_connectivity_key": row.get("input_connectivity_key", ""),
                "input_chembl_id": row.get("input_chembl_id", ""),
                "primary_endpoint": row.get("primary_endpoint", ""),
                "condition_bucket": row.get("condition_bucket", ""),
                "primary_direction": row.get("primary_direction", ""),
                "num_property_objectives": 3,
                "primary_objective_json": compact_json(primary),
                "preserved_property_json": compact_json(preserved),
                "local_constraints_json": compact_json(local_constraints),
                "num_positive_answers": 2,
                "positive_answer_smiles_json": compact_json(
                    [item.get("target_smiles_canon", "") for item in selected_positives]
                ),
                "positive_answers_json": compact_json(selected_positives),
                "mmp_cores_json": row.get("mmp_cores_json", "[]"),
                "source_positive_sample_ids_json": compact_json(
                    [item.get("sample_id", "") for item in selected_positives]
                ),
                "selection_rule": "Top 2 positive answers by primary delta. Negative candidates are omitted by policy.",
            }
        )
    return output, {"kept": len(output), "skipped": dict(skipped)}


def stats(rows: list[dict[str, Any]]) -> dict[str, Any]:
    primary = {row["primary_endpoint"] for row in rows}
    secondary = set()
    inputs = {row["input_connectivity_key"] for row in rows}
    pos_outputs = set()
    for row in rows:
        preserved = safe_json(row.get("preserved_property_json", "{}"), {})
        if isinstance(preserved, dict) and preserved.get("endpoint_name"):
            secondary.add(preserved["endpoint_name"])
        for item in safe_json(row.get("positive_answers_json", "[]"), []):
            pos_outputs.add(item.get("target_connectivity_key") or item.get("target_smiles_canon"))
    return {
        "queries": len(rows),
        "unique_inputs": len(inputs),
        "positive_answers": len(rows) * 2,
        "answers_per_query": 2,
        "property_objectives_per_query": 3,
        "primary_endpoints": len(primary),
        "preserved_secondary_endpoints": len(secondary),
        "unique_positive_outputs": len(pos_outputs),
    }


def molecule_keys(rows: list[dict[str, Any]]) -> set[str]:
    keys = {row.get("input_connectivity_key", "") for row in rows if row.get("input_connectivity_key")}
    for row in rows:
        for item in safe_json(row.get("positive_answers_json", "[]"), []):
            key = item.get("target_connectivity_key") or ""
            if key:
                keys.add(key)
    return keys


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for stale_name in ["pretrain.csv", "strict_all.csv"]:
        stale_path = args.out_dir / stale_name
        if stale_path.exists():
            stale_path.unlink()
    summary: dict[str, Any] = {
        "schema": "chembl_one_input_3property_instruction_to_2positive_answers_v1",
        "row_policy": "Keep rows with at least 2 positive ChEMBL-derived outputs. Emit exactly 2 positive answers and no negative candidates.",
        "property_policy": "Every instruction has exactly 3 objective groups: primary experimental endpoint, one preserved experimental secondary endpoint, and local MW/LogP/QED/SA constraints.",
        "naming_policy": "ChEMBL follows the same source-specific file layout as PubChem: all/train/val/test. The all file is the union of the leakage-safe train/val/test splits.",
        "leakage_check": {},
        "files": {},
    }
    split_rows_by_name: dict[str, list[dict[str, Any]]] = {}
    split_build_by_name: dict[str, dict[str, Any]] = {}
    for name in ["train", "val", "test"]:
        rows, build_stats = choose_rows(read_rows(args.input_dir / f"{name}.csv"), f"q3p2pos_{name}")
        split_rows_by_name[name] = rows
        split_build_by_name[name] = build_stats
        write_rows(args.out_dir / f"{name}.csv", rows)

    all_rows = split_rows_by_name["train"] + split_rows_by_name["val"] + split_rows_by_name["test"]
    write_rows(args.out_dir / "all.csv", all_rows)
    summary["files"]["all"] = {
        "path": str(args.out_dir / "all.csv"),
        "build": {"source": "union of train, val, and test"},
        "stats": stats(all_rows),
    }
    for name in ["train", "val", "test"]:
        summary["files"][name] = {
            "path": str(args.out_dir / f"{name}.csv"),
            "build": split_build_by_name[name],
            "stats": stats(split_rows_by_name[name]),
        }
    split_keys = {name: molecule_keys(rows) for name, rows in split_rows_by_name.items()}
    overlap_counts = {
        "train_val": len(split_keys["train"] & split_keys["val"]),
        "train_test": len(split_keys["train"] & split_keys["test"]),
        "val_test": len(split_keys["val"] & split_keys["test"]),
    }
    summary["leakage_check"] = {
        "molecule_overlap_counts": overlap_counts,
        "passes_no_molecule_overlap": all(value == 0 for value in overlap_counts.values()),
    }
    (args.out_dir / "dataset_stats.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())



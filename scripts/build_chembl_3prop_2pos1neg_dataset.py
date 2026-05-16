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
    "num_candidates",
    "num_positive_answers",
    "positive_answer_smiles_json",
    "positive_answers_json",
    "num_negative_candidates",
    "negative_candidate_smiles_json",
    "negative_candidates_json",
    "candidate_answers_json",
    "mmp_cores_json",
    "source_positive_sample_ids_json",
    "source_negative_sample_ids_json",
    "selection_rule",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build final rows with exactly 3 property objective groups and exactly 2 positive + 1 negative candidates."
    )
    parser.add_argument("--input-dir", type=Path, default=Path("pos_way_admet_benchmark/data/final_multianswer_multiprop"))
    parser.add_argument("--out-dir", type=Path, default=Path("pos_way_admet_benchmark/data/chembl_3prop_2pos1neg"))
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
    tolerance_text = preserved.get("tolerance_text", "within the specified tolerance")
    return (
        f"{verb} {endpoint} while preserving {preserved_name} within {tolerance_text}; "
        "also keep MW, LogP, QED, and synthetic accessibility within local edit constraints."
    )


def choose_rows(source_rows: list[dict[str, str]], prefix: str) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    output: list[dict[str, Any]] = []
    skipped = Counter()
    for row in source_rows:
        positives = safe_json(row.get("positive_answers_json", "[]"), [])
        negatives = safe_json(row.get("negative_candidates_json", "[]"), [])
        objectives = safe_json(row.get("property_objectives_json", "{}"), {})
        secondary = objectives.get("secondary", []) if isinstance(objectives, dict) else []
        primary = objectives.get("primary", {}) if isinstance(objectives, dict) else {}
        if len(positives) < 2:
            skipped["less_than_2_positive"] += 1
            continue
        if len(negatives) < 1:
            skipped["no_negative"] += 1
            continue
        if not secondary:
            skipped["no_secondary_property"] += 1
            continue
        preserved = secondary[0]
        selected_positives = sorted(
            positives,
            key=lambda item: (-as_float(item.get("delta_value")), item.get("target_smiles_canon", "")),
        )[:2]
        selected_negatives = sorted(
            negatives,
            key=lambda item: (-as_float(item.get("tanimoto_similarity")), -as_float(item.get("delta_value"))),
        )[:1]
        local_constraints = {
            "objective_group": "local_constraints",
            "properties": ["MW", "LogP", "QED", "SA"],
            "rule": "Keep within local edit constraints recorded during MMP candidate generation.",
        }
        candidate_answers = [
            {**item, "candidate_label": "positive"} for item in selected_positives
        ] + [
            {**item, "candidate_label": "negative"} for item in selected_negatives
        ]
        query_id = f"{prefix}_{len(output) + 1:08d}"
        output.append(
            {
                "query_id": query_id,
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
                "num_candidates": 3,
                "num_positive_answers": 2,
                "positive_answer_smiles_json": compact_json(
                    [item.get("target_smiles_canon", "") for item in selected_positives]
                ),
                "positive_answers_json": compact_json(selected_positives),
                "num_negative_candidates": 1,
                "negative_candidate_smiles_json": compact_json(
                    [item.get("target_smiles_canon", "") for item in selected_negatives]
                ),
                "negative_candidates_json": compact_json(selected_negatives),
                "candidate_answers_json": compact_json(candidate_answers),
                "mmp_cores_json": row.get("mmp_cores_json", "[]"),
                "source_positive_sample_ids_json": compact_json(
                    [item.get("sample_id", "") for item in selected_positives]
                ),
                "source_negative_sample_ids_json": compact_json(
                    [item.get("sample_id", "") for item in selected_negatives]
                ),
                "selection_rule": "Top 2 positive answers by primary delta; top 1 negative candidate by target-input Tanimoto then primary delta.",
            }
        )
    return output, {"kept": len(output), "skipped": dict(skipped)}


def stats(rows: list[dict[str, Any]]) -> dict[str, Any]:
    primary = {row["primary_endpoint"] for row in rows}
    secondary = set()
    inputs = {row["input_connectivity_key"] for row in rows}
    pos_outputs = set()
    neg_outputs = set()
    for row in rows:
        preserved = safe_json(row.get("preserved_property_json", "{}"), {})
        if isinstance(preserved, dict) and preserved.get("endpoint_name"):
            secondary.add(preserved["endpoint_name"])
        for item in safe_json(row.get("positive_answers_json", "[]"), []):
            pos_outputs.add(item.get("target_connectivity_key") or item.get("target_smiles_canon"))
        for item in safe_json(row.get("negative_candidates_json", "[]"), []):
            neg_outputs.add(item.get("target_connectivity_key") or item.get("target_smiles_canon"))
    return {
        "queries": len(rows),
        "unique_inputs": len(inputs),
        "positive_answers": len(rows) * 2,
        "negative_candidates": len(rows),
        "candidate_molecules_per_query": 3,
        "positive_answers_per_query": 2,
        "negative_candidates_per_query": 1,
        "property_objectives_per_query": 3,
        "primary_endpoints": len(primary),
        "preserved_secondary_endpoints": len(secondary),
        "unique_positive_outputs": len(pos_outputs),
        "unique_negative_outputs": len(neg_outputs),
    }


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary: dict[str, Any] = {
        "schema": "one_input_3property_instruction_to_2positive_1negative_v1",
        "row_policy": "Keep only rows with at least 2 positive answers and at least 1 negative candidate. Emit exactly 2 positive answers and 1 negative candidate.",
        "property_policy": "Every instruction has exactly 3 objective groups: primary experimental endpoint, one preserved experimental secondary endpoint, and local MW/LogP/QED/SA constraints.",
        "files": {},
    }
    for name in ["pretrain", "train", "val", "test", "strict_all"]:
        rows, build_stats = choose_rows(read_rows(args.input_dir / f"{name}.csv"), f"q3p_{name}")
        write_rows(args.out_dir / f"{name}.csv", rows)
        summary["files"][name] = {
            "path": str(args.out_dir / f"{name}.csv"),
            "build": build_stats,
            "stats": stats(rows),
        }
    (args.out_dir / "dataset_stats.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


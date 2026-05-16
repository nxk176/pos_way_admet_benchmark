from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
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
    "num_positive_answers",
    "positive_answer_smiles_json",
    "positive_answers_json",
    "num_negative_candidates",
    "negative_candidate_smiles_json",
    "negative_candidates_json",
    "num_secondary_objectives",
    "secondary_objectives_json",
    "property_objectives_json",
    "mmp_cores_json",
    "source_positive_sample_ids_json",
    "source_negative_sample_ids_json",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build the canonical one-input multi-property instruction to multi-answer dataset."
    )
    parser.add_argument(
        "--large-pair",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/multiproperty_large_multiprop/fragment_multiproperty_samples.csv"),
    )
    parser.add_argument(
        "--strict-splits-dir",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/multiproperty_strict_multiprop/splits"),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/final_multianswer_multiprop"),
    )
    parser.add_argument("--min-positive-answers", type=int, default=2)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def safe_json(raw: str, fallback: Any) -> Any:
    try:
        return json.loads(raw)
    except (TypeError, json.JSONDecodeError):
        return fallback


def as_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("-inf")


def direction_from_instruction(instruction: str) -> str:
    return "decrease" if instruction.strip().lower().startswith("decrease") else "increase"


def group_key(row: dict[str, str]) -> tuple[str, str, str, str, str, str]:
    return (
        row.get("split", ""),
        row.get("instruction", ""),
        row.get("input_smiles_canon", ""),
        row.get("input_connectivity_key", ""),
        row.get("primary_endpoint", ""),
        row.get("condition_bucket", ""),
    )


def answer_item(row: dict[str, str]) -> dict[str, Any]:
    return {
        "sample_id": row.get("sample_id", ""),
        "target_smiles_canon": row.get("target_smiles_canon", ""),
        "target_connectivity_key": row.get("target_connectivity_key", ""),
        "target_chembl_id": row.get("target_chembl_id", ""),
        "value_before": row.get("value_before", ""),
        "value_after": row.get("value_after", ""),
        "delta_value": row.get("delta_value", ""),
        "relative_delta": row.get("relative_delta", ""),
        "tanimoto_similarity": row.get("tanimoto_similarity", ""),
        "mmp_core": row.get("mmp_core", ""),
        "input_variable_fragment": row.get("input_variable_fragment", ""),
        "target_variable_fragment": row.get("target_variable_fragment", ""),
        "secondary_objectives": safe_json(row.get("secondary_objectives", "[]"), []),
        "input_observation_ids": safe_json(row.get("input_observation_ids", "[]"), []),
        "target_observation_ids": safe_json(row.get("target_observation_ids", "[]"), []),
        "input_assay_ids": safe_json(row.get("input_assay_ids", "[]"), []),
        "target_assay_ids": safe_json(row.get("target_assay_ids", "[]"), []),
    }


def negative_item(row: dict[str, str]) -> dict[str, Any]:
    item = answer_item(row)
    item["negative_failure_reason"] = row.get("negative_failure_reason", "")
    return item


def objective_names(objectives: list[dict[str, Any]]) -> list[str]:
    names = []
    for objective in objectives:
        if isinstance(objective, dict) and objective.get("endpoint_name"):
            names.append(str(objective["endpoint_name"]))
    return names


def build_rows(rows: list[dict[str, str]], prefix: str, min_positive_answers: int) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    groups: dict[tuple[str, str, str, str, str, str], dict[str, Any]] = {}
    seen_positive: dict[tuple[str, str, str, str, str, str], set[str]] = defaultdict(set)
    seen_negative: dict[tuple[str, str, str, str, str, str], set[str]] = defaultdict(set)

    for row in rows:
        key = group_key(row)
        if key not in groups:
            split, instruction, input_smiles, input_key, primary_endpoint, condition_bucket = key
            secondary_objectives = safe_json(row.get("secondary_objectives", "[]"), [])
            groups[key] = {
                "split": split,
                "instruction": instruction,
                "input_smiles_canon": input_smiles,
                "input_connectivity_key": input_key,
                "input_chembl_id": row.get("input_chembl_id", ""),
                "primary_endpoint": primary_endpoint,
                "condition_bucket": condition_bucket,
                "primary_direction": direction_from_instruction(instruction),
                "secondary_objectives": secondary_objectives,
                "property_objectives": {
                    "primary": safe_json(row.get("primary_objective", "{}"), {}),
                    "secondary": secondary_objectives,
                    "local_constraints": ["MW", "LogP", "QED", "SA"],
                },
                "positive_answers": [],
                "negative_candidates": [],
                "mmp_cores": set(),
            }

        target_key = row.get("target_connectivity_key", "") or row.get("target_smiles_canon", "")
        if row.get("mmp_core"):
            groups[key]["mmp_cores"].add(row["mmp_core"])

        if row.get("sample_type") == "positive":
            if target_key and target_key not in seen_positive[key]:
                seen_positive[key].add(target_key)
                groups[key]["positive_answers"].append(answer_item(row))
        elif row.get("sample_type") == "negative":
            if target_key and target_key not in seen_negative[key]:
                seen_negative[key].add(target_key)
                groups[key]["negative_candidates"].append(negative_item(row))

    output: list[dict[str, Any]] = []
    for item in groups.values():
        positives = item["positive_answers"]
        if len(positives) < min_positive_answers:
            continue
        positives = sorted(
            positives,
            key=lambda answer: (-as_float(answer.get("delta_value")), answer.get("target_smiles_canon", "")),
        )
        negatives = sorted(
            item["negative_candidates"],
            key=lambda answer: (-as_float(answer.get("tanimoto_similarity")), -as_float(answer.get("delta_value"))),
        )
        output.append(
            {
                "query_id": f"{prefix}_{len(output) + 1:08d}",
                "split": item["split"],
                "instruction": item["instruction"],
                "input_smiles_canon": item["input_smiles_canon"],
                "input_connectivity_key": item["input_connectivity_key"],
                "input_chembl_id": item["input_chembl_id"],
                "primary_endpoint": item["primary_endpoint"],
                "condition_bucket": item["condition_bucket"],
                "primary_direction": item["primary_direction"],
                "num_positive_answers": len(positives),
                "positive_answer_smiles_json": compact_json([answer["target_smiles_canon"] for answer in positives]),
                "positive_answers_json": compact_json(positives),
                "num_negative_candidates": len(negatives),
                "negative_candidate_smiles_json": compact_json([answer["target_smiles_canon"] for answer in negatives]),
                "negative_candidates_json": compact_json(negatives),
                "num_secondary_objectives": len(item["secondary_objectives"]),
                "secondary_objectives_json": compact_json(item["secondary_objectives"]),
                "property_objectives_json": compact_json(item["property_objectives"]),
                "mmp_cores_json": compact_json(sorted(item["mmp_cores"])),
                "source_positive_sample_ids_json": compact_json([answer["sample_id"] for answer in positives]),
                "source_negative_sample_ids_json": compact_json([answer["sample_id"] for answer in negatives]),
            }
        )

    return output, stats(output)


def compact_json(payload: Any) -> str:
    return json.dumps(payload, ensure_ascii=False, separators=(",", ":"))


def stats(rows: list[dict[str, Any]]) -> dict[str, Any]:
    positive_counts = [int(row["num_positive_answers"]) for row in rows]
    negative_counts = [int(row["num_negative_candidates"]) for row in rows]
    secondary_names = set()
    for row in rows:
        for objective in safe_json(row.get("secondary_objectives_json", "[]"), []):
            if isinstance(objective, dict) and objective.get("endpoint_name"):
                secondary_names.add(objective["endpoint_name"])
    return {
        "queries": len(rows),
        "positive_answers": sum(positive_counts),
        "negative_candidates": sum(negative_counts),
        "queries_with_negative_candidates": sum(1 for count in negative_counts if count > 0),
        "primary_endpoints": len({row["primary_endpoint"] for row in rows}),
        "secondary_endpoints": len(secondary_names),
        "max_positive_answers_per_query": max(positive_counts) if positive_counts else 0,
        "max_negative_candidates_per_query": max(negative_counts) if negative_counts else 0,
        "positive_answers_per_query": dict(sorted(Counter(positive_counts).items())),
    }


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    large_rows, large_stats = build_rows(read_csv(args.large_pair), "final_pretrain", args.min_positive_answers)
    write_csv(args.out_dir / "pretrain.csv", large_rows)

    strict_stats: dict[str, Any] = {}
    strict_all: list[dict[str, Any]] = []
    for split in ["train", "val", "test"]:
        split_rows, split_stats = build_rows(
            read_csv(args.strict_splits_dir / f"{split}.csv"),
            f"final_{split}",
            args.min_positive_answers,
        )
        write_csv(args.out_dir / f"{split}.csv", split_rows)
        strict_all.extend(split_rows)
        strict_stats[split] = split_stats

    write_csv(args.out_dir / "strict_all.csv", strict_all)
    stats_payload = {
        "schema": "one_input_multiproperty_instruction_to_multianswer_v1",
        "min_positive_answers": args.min_positive_answers,
        "inputs": {
            "large_pair": str(args.large_pair),
            "strict_splits_dir": str(args.strict_splits_dir),
        },
        "outputs": {
            "pretrain": str(args.out_dir / "pretrain.csv"),
            "train": str(args.out_dir / "train.csv"),
            "val": str(args.out_dir / "val.csv"),
            "test": str(args.out_dir / "test.csv"),
            "strict_all": str(args.out_dir / "strict_all.csv"),
        },
        "pretrain": large_stats,
        "strict": {
            **strict_stats,
            "all": stats(strict_all),
        },
    }
    (args.out_dir / "dataset_stats.json").write_text(
        json.dumps(stats_payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(stats_payload, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

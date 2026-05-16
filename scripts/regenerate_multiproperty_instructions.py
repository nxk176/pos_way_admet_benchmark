from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any


PAIR_EXTRA_FIELDS = [
    "original_instruction",
    "primary_objective",
    "secondary_objectives",
    "num_secondary_objectives",
    "multiproperty_success",
    "negative_failure_reason",
]

SFT_FIELDS = [
    "example_id",
    "split",
    "source_sample_id",
    "instruction",
    "input_smiles_canon",
    "output_smiles_canon",
    "input_connectivity_key",
    "output_connectivity_key",
    "primary_endpoint",
    "condition_bucket",
    "value_before",
    "value_after",
    "delta_value",
    "secondary_objectives",
]

TRIPLET_EXTRA_FIELDS = [
    "original_instruction",
    "secondary_objectives",
    "num_secondary_objectives",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Regenerate instructions so every retained row has explicit secondary property objectives."
    )
    parser.add_argument(
        "--strict-splits-dir",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/multiproperty_expanded_strict/splits"),
    )
    parser.add_argument(
        "--strict-full",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/multiproperty_expanded_strict/fragment_multiproperty_samples.csv"),
    )
    parser.add_argument(
        "--large-pair",
        type=Path,
        default=None,
        help="Optional broad training pair CSV to regenerate with the same multi-property instruction policy.",
    )
    parser.add_argument(
        "--triplet-full",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/ranking_strict/strict_hard_negative_triplets.csv"),
    )
    parser.add_argument(
        "--triplet-same-input",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/ranking_strict/strict_hard_negative_triplets_same_input.csv"),
    )
    parser.add_argument(
        "--out-pair-dir",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/multiproperty_strict_multiprop"),
    )
    parser.add_argument(
        "--out-sft-dir",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/sft_multiprop_strict"),
    )
    parser.add_argument(
        "--out-ranking-dir",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/ranking_strict_multiprop"),
    )
    parser.add_argument(
        "--out-large-pair-dir",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/multiproperty_large_multiprop"),
    )
    parser.add_argument(
        "--out-large-sft-dir",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/sft_multiprop_large"),
    )
    parser.add_argument("--max-secondary-objectives", type=int, default=2)
    parser.add_argument("--min-secondary-objectives", type=int, default=1)
    parser.add_argument("--p-scale-preserve-tolerance", type=float, default=0.3)
    parser.add_argument("--log-scale-preserve-tolerance", type=float, default=0.3)
    parser.add_argument("--relative-preserve-tolerance", type=float, default=0.25)
    return parser.parse_args()


def read_csv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader), list(reader.fieldnames or [])


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: clean_cell(row.get(key)) for key in fieldnames})
    return len(rows)


def write_jsonl(path: Path, rows: list[dict[str, Any]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False, separators=(",", ":")) + "\n")
    return len(rows)


def clean_cell(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float, str)):
        return value
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def as_float(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def bool_text(value: str) -> bool:
    return str(value or "").strip().lower() == "true"


def parse_shared_secondary(row: dict[str, str]) -> list[dict[str, Any]]:
    raw = row.get("shared_secondary_experimental_endpoints", "")
    if not raw:
        return []
    try:
        payload = json.loads(raw)
    except json.JSONDecodeError:
        return []
    items = payload.get("items", []) if isinstance(payload, dict) else []
    return [item for item in items if isinstance(item, dict)]


def tolerance_for(endpoint: str, input_value: float, args: argparse.Namespace) -> tuple[float, str]:
    lower = endpoint.lower()
    if any(token in endpoint for token in ["pChEMBL", "pIC50", "pKi", "pKd", "pEC50"]):
        return args.p_scale_preserve_tolerance, f"+/- {args.p_scale_preserve_tolerance:g} log unit"
    if "logs" in lower or "logpapp" in lower:
        return args.log_scale_preserve_tolerance, f"+/- {args.log_scale_preserve_tolerance:g} log unit"
    if "half_life" in lower:
        abs_tol = max(30.0, abs(input_value) * args.relative_preserve_tolerance)
        return abs_tol, f"+/- {args.relative_preserve_tolerance:g} relative or 30 min"
    if "clearance" in lower:
        abs_tol = max(5.0, abs(input_value) * args.relative_preserve_tolerance)
        return abs_tol, f"+/- {args.relative_preserve_tolerance:g} relative or 5 unit"
    abs_tol = max(0.3, abs(input_value) * args.relative_preserve_tolerance)
    return abs_tol, f"+/- {args.relative_preserve_tolerance:g} relative"


def preserved_secondary_objectives(row: dict[str, str], args: argparse.Namespace) -> list[dict[str, Any]]:
    objectives: list[dict[str, Any]] = []
    for item in parse_shared_secondary(row):
        endpoint = str(item.get("endpoint_name", ""))
        input_value = as_float(item.get("input_value"))
        target_value = as_float(item.get("target_value"))
        delta = as_float(item.get("delta"))
        if not endpoint or input_value is None or target_value is None:
            continue
        if delta is None:
            delta = target_value - input_value
        tolerance, tolerance_text = tolerance_for(endpoint, input_value, args)
        abs_delta = abs(delta)
        if abs_delta > tolerance:
            continue
        objectives.append(
            {
                "endpoint_name": endpoint,
                "input_value": round(input_value, 6),
                "target_value": round(target_value, 6),
                "delta": round(delta, 6),
                "abs_delta": round(abs_delta, 6),
                "tolerance": round(tolerance, 6),
                "tolerance_text": tolerance_text,
            }
        )
    objectives.sort(key=lambda obj: (obj["abs_delta"] / max(obj["tolerance"], 1e-12), obj["endpoint_name"]))
    return objectives[: args.max_secondary_objectives]


def primary_direction(row: dict[str, str]) -> str:
    instruction = row.get("instruction", "").strip().lower()
    if instruction.startswith("decrease"):
        return "Decrease"
    return "Increase"


def build_instruction(row: dict[str, str], objectives: list[dict[str, Any]]) -> str:
    direction = primary_direction(row)
    primary = row.get("primary_endpoint", "")
    if len(objectives) == 1:
        secondary_text = f"{objectives[0]['endpoint_name']} within {objectives[0]['tolerance_text']}"
    else:
        parts = [f"{obj['endpoint_name']} within {obj['tolerance_text']}" for obj in objectives]
        secondary_text = ", ".join(parts[:-1]) + f", and {parts[-1]}"
    return (
        f"{direction} {primary} while preserving {secondary_text}; "
        "also keep MW, LogP, QED, and synthetic accessibility within local edit constraints."
    )


def row_failure_reason(row: dict[str, str]) -> str:
    if bool_text(row.get("primary_success", "")):
        return "not_negative_primary_failure"
    if not bool_text(row.get("secondary_success", "")):
        return "failed_rdkit_secondary_constraints"
    return "primary_not_improved_enough"


def transform_pair_row(row: dict[str, str], args: argparse.Namespace) -> dict[str, Any] | None:
    objectives = preserved_secondary_objectives(row, args)
    if len(objectives) < args.min_secondary_objectives:
        return None

    sample_type = row.get("sample_type", "")
    primary_ok = bool_text(row.get("primary_success", ""))
    rdkit_secondary_ok = bool_text(row.get("secondary_success", ""))
    if sample_type == "positive":
        if not primary_ok or not rdkit_secondary_ok:
            return None
        multiproperty_success = True
        negative_reason = ""
    elif sample_type == "negative":
        if primary_ok or not rdkit_secondary_ok:
            return None
        multiproperty_success = False
        negative_reason = row_failure_reason(row)
    else:
        return None

    out = dict(row)
    out["original_instruction"] = row.get("instruction", "")
    out["instruction"] = build_instruction(row, objectives)
    out["primary_objective"] = json.dumps(
        {
            "endpoint_name": row.get("primary_endpoint", ""),
            "direction": primary_direction(row).lower(),
            "value_before": as_float(row.get("value_before")),
            "value_after": as_float(row.get("value_after")),
            "delta_value": as_float(row.get("delta_value")),
        },
        ensure_ascii=False,
        separators=(",", ":"),
    )
    out["secondary_objectives"] = json.dumps(objectives, ensure_ascii=False, separators=(",", ":"))
    out["num_secondary_objectives"] = len(objectives)
    out["multiproperty_success"] = multiproperty_success
    out["negative_failure_reason"] = negative_reason
    return out


def transform_pair_rows(rows: list[dict[str, str]], args: argparse.Namespace) -> tuple[list[dict[str, Any]], Counter[str]]:
    counters: Counter[str] = Counter()
    output = []
    for row in rows:
        transformed = transform_pair_row(row, args)
        if transformed is None:
            counters[f"filtered:{row.get('sample_type', 'unknown')}"] += 1
            continue
        counters[f"kept:{row.get('sample_type', 'unknown')}"] += 1
        output.append(transformed)
    return output, counters


def sft_rows_from_pairs(rows: list[dict[str, Any]], split: str) -> list[dict[str, Any]]:
    output = []
    for row in rows:
        if row.get("sample_type") != "positive":
            continue
        output.append(
            {
                "example_id": f"sft_{split}_{len(output) + 1:08d}",
                "split": split,
                "source_sample_id": row.get("sample_id", ""),
                "instruction": row.get("instruction", ""),
                "input_smiles_canon": row.get("input_smiles_canon", ""),
                "output_smiles_canon": row.get("target_smiles_canon", ""),
                "input_connectivity_key": row.get("input_connectivity_key", ""),
                "output_connectivity_key": row.get("target_connectivity_key", ""),
                "primary_endpoint": row.get("primary_endpoint", ""),
                "condition_bucket": row.get("condition_bucket", ""),
                "value_before": row.get("value_before", ""),
                "value_after": row.get("value_after", ""),
                "delta_value": row.get("delta_value", ""),
                "secondary_objectives": row.get("secondary_objectives", ""),
            }
        )
    return output


def objective_map(row: dict[str, Any]) -> dict[str, dict[str, Any]]:
    try:
        objectives = json.loads(row.get("secondary_objectives", "[]"))
    except json.JSONDecodeError:
        return {}
    return {obj["endpoint_name"]: obj for obj in objectives if isinstance(obj, dict) and obj.get("endpoint_name")}


def triplet_instruction(base_row: dict[str, str], objectives: list[dict[str, Any]]) -> str:
    pseudo_pair = {
        "instruction": base_row.get("instruction", ""),
        "primary_endpoint": base_row.get("primary_endpoint", ""),
    }
    return build_instruction(pseudo_pair, objectives)


def transform_triplet_rows(
    rows: list[dict[str, str]],
    pair_by_sample_id: dict[str, dict[str, Any]],
) -> tuple[list[dict[str, Any]], Counter[str]]:
    output = []
    counters: Counter[str] = Counter()
    for row in rows:
        positive = pair_by_sample_id.get(row.get("positive_sample_id", ""))
        negative = pair_by_sample_id.get(row.get("negative_sample_id", ""))
        if positive is None or negative is None:
            counters["filtered:missing_filtered_pair"] += 1
            continue
        pos_objectives = objective_map(positive)
        neg_objectives = objective_map(negative)
        shared_names = sorted(set(pos_objectives) & set(neg_objectives))
        if not shared_names:
            counters["filtered:no_shared_preserved_secondary"] += 1
            continue
        objectives = []
        for name in shared_names[:2]:
            pos_obj = pos_objectives[name]
            neg_obj = neg_objectives[name]
            objectives.append(
                {
                    "endpoint_name": name,
                    "positive_delta": pos_obj.get("delta"),
                    "negative_delta": neg_obj.get("delta"),
                    "tolerance": pos_obj.get("tolerance"),
                    "tolerance_text": pos_obj.get("tolerance_text"),
                }
            )
        out = dict(row)
        out["original_instruction"] = row.get("instruction", "")
        out["instruction"] = triplet_instruction(row, objectives)
        out["secondary_objectives"] = json.dumps(objectives, ensure_ascii=False, separators=(",", ":"))
        out["num_secondary_objectives"] = len(objectives)
        output.append(out)
        counters[f"kept:{row.get('split', 'unknown')}"] += 1
    return output, counters


def add_fields(fieldnames: list[str], extra: list[str]) -> list[str]:
    return [*fieldnames, *[field for field in extra if field not in fieldnames]]


def process_pairs(args: argparse.Namespace) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    split_stats: dict[str, Any] = {}
    pair_by_sample_id: dict[str, dict[str, Any]] = {}
    all_rows: list[dict[str, Any]] = []

    for split in ["train", "val", "test"]:
        rows, fieldnames = read_csv(args.strict_splits_dir / f"{split}.csv")
        transformed, counters = transform_pair_rows(rows, args)
        for row in transformed:
            row["split"] = split
            pair_by_sample_id[row.get("sample_id", "")] = row
        all_rows.extend(transformed)
        out_fields = add_fields(add_fields(fieldnames, ["split"]), PAIR_EXTRA_FIELDS)
        write_csv(args.out_pair_dir / "splits" / f"{split}.csv", transformed, out_fields)
        write_jsonl(args.out_pair_dir / "splits" / f"{split}.jsonl", transformed)
        sft_rows = sft_rows_from_pairs(transformed, split)
        write_csv(args.out_sft_dir / f"{split}.csv", sft_rows, SFT_FIELDS)
        write_jsonl(args.out_sft_dir / f"{split}.jsonl", sft_rows)
        split_stats[split] = {
            "rows": len(transformed),
            "positive": sum(1 for row in transformed if row.get("sample_type") == "positive"),
            "negative": sum(1 for row in transformed if row.get("sample_type") == "negative"),
            "sft_positive_examples": len(sft_rows),
            "counters": dict(sorted(counters.items())),
        }

    _, full_fieldnames = read_csv(args.strict_full)
    out_fields = add_fields(add_fields(full_fieldnames, ["split"]), PAIR_EXTRA_FIELDS)
    write_csv(args.out_pair_dir / "fragment_multiproperty_samples.csv", all_rows, out_fields)
    write_jsonl(args.out_pair_dir / "fragment_multiproperty_samples.jsonl", all_rows)
    stats = {
        "source_strict_splits_dir": str(args.strict_splits_dir),
        "source_strict_full": str(args.strict_full),
        "row_policy": (
            "Keep positives that satisfy primary improvement, RDKit secondary constraints, and at least one "
            "preserved shared experimental secondary endpoint. Keep negatives only when they preserve the "
            "secondary endpoint and RDKit constraints but fail primary improvement."
        ),
        "instruction_policy": "Explicitly name preserved shared experimental secondary endpoint(s) in every instruction.",
        "thresholds": {
            "p_scale_preserve_tolerance": args.p_scale_preserve_tolerance,
            "log_scale_preserve_tolerance": args.log_scale_preserve_tolerance,
            "relative_preserve_tolerance": args.relative_preserve_tolerance,
            "min_secondary_objectives": args.min_secondary_objectives,
            "max_secondary_objectives": args.max_secondary_objectives,
        },
        "splits": split_stats,
        "total_rows": len(all_rows),
        "total_positive": sum(1 for row in all_rows if row.get("sample_type") == "positive"),
        "total_negative": sum(1 for row in all_rows if row.get("sample_type") == "negative"),
    }
    args.out_pair_dir.mkdir(parents=True, exist_ok=True)
    (args.out_pair_dir / "multiproperty_instruction_stats.json").write_text(
        json.dumps(stats, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return pair_by_sample_id, stats


def process_large_pairs(args: argparse.Namespace) -> dict[str, Any] | None:
    if args.large_pair is None:
        return None

    rows, fieldnames = read_csv(args.large_pair)
    transformed, counters = transform_pair_rows(rows, args)
    for row in transformed:
        row["split"] = "pretrain"

    out_fields = add_fields(add_fields(fieldnames, ["split"]), PAIR_EXTRA_FIELDS)
    write_csv(args.out_large_pair_dir / "fragment_multiproperty_samples.csv", transformed, out_fields)
    write_jsonl(args.out_large_pair_dir / "fragment_multiproperty_samples.jsonl", transformed)

    sft_rows = sft_rows_from_pairs(transformed, "pretrain")
    write_csv(args.out_large_sft_dir / "all.csv", sft_rows, SFT_FIELDS)
    write_jsonl(args.out_large_sft_dir / "all.jsonl", sft_rows)

    stats = {
        "source_large_pair": str(args.large_pair),
        "row_policy": (
            "Keep positives that satisfy primary improvement, RDKit secondary constraints, and at least one "
            "preserved shared experimental secondary endpoint. Keep negatives only when they preserve the "
            "secondary endpoint and RDKit constraints but fail primary improvement."
        ),
        "instruction_policy": "Explicitly name preserved shared experimental secondary endpoint(s) in every instruction.",
        "thresholds": {
            "p_scale_preserve_tolerance": args.p_scale_preserve_tolerance,
            "log_scale_preserve_tolerance": args.log_scale_preserve_tolerance,
            "relative_preserve_tolerance": args.relative_preserve_tolerance,
            "min_secondary_objectives": args.min_secondary_objectives,
            "max_secondary_objectives": args.max_secondary_objectives,
        },
        "total_rows": len(transformed),
        "total_positive": sum(1 for row in transformed if row.get("sample_type") == "positive"),
        "total_negative": sum(1 for row in transformed if row.get("sample_type") == "negative"),
        "sft_positive_examples": len(sft_rows),
        "counters": dict(sorted(counters.items())),
    }
    args.out_large_pair_dir.mkdir(parents=True, exist_ok=True)
    (args.out_large_pair_dir / "multiproperty_instruction_stats.json").write_text(
        json.dumps(stats, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return stats


def process_triplets(
    args: argparse.Namespace,
    pair_by_sample_id: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    stats: dict[str, Any] = {}
    for name, source in [
        ("strict_hard_negative_triplets", args.triplet_full),
        ("strict_hard_negative_triplets_same_input", args.triplet_same_input),
    ]:
        rows, fieldnames = read_csv(source)
        transformed, counters = transform_triplet_rows(rows, pair_by_sample_id)
        out_fields = add_fields(fieldnames, TRIPLET_EXTRA_FIELDS)
        write_csv(args.out_ranking_dir / f"{name}.csv", transformed, out_fields)
        write_jsonl(args.out_ranking_dir / f"{name}.jsonl", transformed)
        stats[name] = {
            "source": str(source),
            "rows": len(transformed),
            "counts_by_split": dict(sorted(Counter(row.get("split", "") for row in transformed).items())),
            "counters": dict(sorted(counters.items())),
        }
    args.out_ranking_dir.mkdir(parents=True, exist_ok=True)
    (args.out_ranking_dir / "multiproperty_ranking_stats.json").write_text(
        json.dumps(stats, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return stats


def main() -> int:
    args = parse_args()
    large_pair_stats = process_large_pairs(args)
    pair_by_sample_id, pair_stats = process_pairs(args)
    ranking_stats = process_triplets(args, pair_by_sample_id)
    print(
        json.dumps(
            {
                "large_multiproperty_pair_rows": large_pair_stats["total_rows"] if large_pair_stats else None,
                "large_multiproperty_positive_rows": large_pair_stats["total_positive"] if large_pair_stats else None,
                "large_multiproperty_negative_rows": large_pair_stats["total_negative"] if large_pair_stats else None,
                "multiproperty_pair_rows": pair_stats["total_rows"],
                "multiproperty_positive_rows": pair_stats["total_positive"],
                "multiproperty_negative_rows": pair_stats["total_negative"],
                "ranking_rows": {key: value["rows"] for key, value in ranking_stats.items()},
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

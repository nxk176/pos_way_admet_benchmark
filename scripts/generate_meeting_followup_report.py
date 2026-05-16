from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from statistics import median
from typing import Any, Iterable


MULTIPROPERTY_DATASETS = [
    ("ChEMBL", "chembl_3prop_2pos"),
    ("PubChem", "pubchem_3prop_2pos"),
    ("Papyrus", "papyrus_3prop_2pos"),
]

BINDINGDB_SPLITS = ["train", "val", "test_seen_target", "test_unseen_target"]


def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description="Generate a meeting follow-up report: dataset inventory, quality audit, and first BindingDB baselines."
    )
    parser.add_argument("--root", type=Path, default=repo_root)
    parser.add_argument("--out-dir", type=Path, default=repo_root / "reports")
    return parser.parse_args()


def read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path: Path) -> Iterable[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle)


def safe_json(raw: Any, fallback: Any) -> Any:
    if raw is None:
        return fallback
    try:
        return json.loads(raw)
    except (TypeError, json.JSONDecodeError):
        return fallback


def as_float(value: Any) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def as_int(value: Any) -> int | None:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def compact_number(value: float | int | None, digits: int = 4) -> float | None:
    if value is None:
        return None
    return round(float(value), digits)


def tie_aware_accuracy(score_positive: float | None, score_negative: float | None) -> float | None:
    if score_positive is None or score_negative is None:
        return None
    if score_positive > score_negative:
        return 1.0
    if score_positive < score_negative:
        return 0.0
    return 0.5


def average(values: list[float]) -> float | None:
    return sum(values) / len(values) if values else None


def collect_file_stats(dataset_stats: dict[str, Any], folder: Path, split: str) -> dict[str, Any]:
    files = dataset_stats.get("files", {})
    if split in files and isinstance(files[split], dict):
        stats = files[split].get("stats", {})
        if isinstance(stats, dict):
            return dict(stats)
    if split == "all" and "final_queries" in dataset_stats:
        return {
            "queries": dataset_stats.get("final_queries"),
            "positive_answers": dataset_stats.get("positive_answers"),
            "answers_per_query": dataset_stats.get("answers_per_query"),
            "property_objectives_per_query": dataset_stats.get("property_objectives_per_query"),
            "primary_endpoints": dataset_stats.get("primary_endpoint_count"),
            "preserved_secondary_endpoints": dataset_stats.get("preserved_secondary_endpoint_count"),
        }
    path = folder / f"{split}.csv"
    if not path.is_file():
        return {}
    count = sum(1 for _ in read_csv(path))
    return {"queries": count}


def audit_multiproperty_dataset(root: Path, display_name: str, folder_name: str) -> dict[str, Any]:
    folder = root / "data" / folder_name
    stats = read_json(folder / "dataset_stats.json")
    all_path = folder / "all.csv"
    split_counter: Counter[str] = Counter()
    issue_counter: Counter[str] = Counter()
    primary_improvement_values: list[float] = []
    tanimoto_values: list[float] = []
    descriptor_success_values: list[float] = []
    answer_payload_shapes: Counter[str] = Counter()
    rows = 0

    if all_path.is_file():
        for row in read_csv(all_path):
            rows += 1
            split_counter[row.get("split", "")] += 1

            if row.get("num_property_objectives") != "3":
                issue_counter["num_property_objectives_not_3"] += 1
            if row.get("num_positive_answers") != "2":
                issue_counter["num_positive_answers_not_2"] += 1
            if not row.get("instruction"):
                issue_counter["empty_instruction"] += 1
            if not row.get("input_smiles_canon"):
                issue_counter["empty_input_smiles"] += 1

            answer_smiles = safe_json(row.get("positive_answer_smiles_json"), [])
            if not isinstance(answer_smiles, list) or len(answer_smiles) != 2:
                issue_counter["positive_answer_smiles_json_not_len_2"] += 1

            answers = safe_json(row.get("positive_answers_json"), [])
            direction = (row.get("primary_direction") or "increase").strip().lower()
            if isinstance(answers, list) and answers and isinstance(answers[0], dict):
                answer_payload_shapes["dict_answers"] += 1
                for answer in answers:
                    delta = as_float(answer.get("delta_value"))
                    if delta is not None:
                        primary_improvement_values.append(-delta if direction == "decrease" else delta)
                    tanimoto = as_float(answer.get("tanimoto_similarity"))
                    if tanimoto is not None:
                        tanimoto_values.append(tanimoto)
                    descriptor_constraints = answer.get("descriptor_constraints")
                    if isinstance(descriptor_constraints, dict) and isinstance(descriptor_constraints.get("success"), bool):
                        descriptor_success_values.append(1.0 if descriptor_constraints["success"] else 0.0)
            elif isinstance(answers, list):
                answer_payload_shapes["non_dict_answer_list"] += 1
            else:
                answer_payload_shapes["unparsed_answers"] += 1

    leakage = stats.get("leakage_check", {})
    if isinstance(leakage, str):
        leakage_summary: Any = leakage
    elif isinstance(leakage, dict):
        leakage_summary = {
            "passes_no_molecule_overlap": leakage.get("passes_no_molecule_overlap"),
            "molecule_overlap_counts": leakage.get("molecule_overlap_counts"),
        }
    else:
        leakage_summary = None

    return {
        "name": display_name,
        "folder": str(folder.relative_to(root)),
        "rows_read_from_all_csv": rows,
        "file_stats": {
            split: collect_file_stats(stats, folder, split) for split in ["all", "train", "val", "test"]
        },
        "split_rows_seen_in_all_csv": dict(split_counter),
        "schema_audit": {
            "issue_counts": dict(sorted(issue_counter.items())),
            "passes_basic_shape_audit": sum(issue_counter.values()) == 0,
            "answer_payload_shapes": dict(answer_payload_shapes),
        },
        "label_quality_signals": {
            "median_answer_primary_improvement": compact_number(median(primary_improvement_values), 4)
            if primary_improvement_values
            else None,
            "median_answer_tanimoto": compact_number(median(tanimoto_values), 4) if tanimoto_values else None,
            "descriptor_success_rate": compact_number(average(descriptor_success_values), 4)
            if descriptor_success_values
            else None,
            "answer_rows_with_improvement_values": len(primary_improvement_values),
        },
        "leakage": leakage_summary,
    }


def audit_bindingdb_split(path: Path) -> dict[str, Any]:
    rows = 0
    issue_counter: Counter[str] = Counter()
    input_similarity_scores: list[float] = []
    evidence_scores: list[float] = []
    positive_deltas: list[float] = []
    negative_deltas: list[float] = []
    input_positive_tanimoto_values: list[float] = []
    input_negative_tanimoto_values: list[float] = []
    positive_negative_tanimoto_values: list[float] = []
    negative_more_input_similar = 0

    for row in read_csv(path):
        rows += 1
        pos_delta = as_float(row.get("positive_delta"))
        neg_delta = as_float(row.get("negative_delta"))
        pos_sim = as_float(row.get("input_positive_tanimoto"))
        neg_sim = as_float(row.get("input_negative_tanimoto"))
        pos_neg_sim = as_float(row.get("positive_negative_tanimoto"))
        pos_evidence = as_int(row.get("positive_evidence_count"))
        neg_evidence = as_int(row.get("negative_evidence_count"))

        if pos_delta is None or pos_delta <= 0:
            issue_counter["positive_delta_not_positive"] += 1
        else:
            positive_deltas.append(pos_delta)
        if neg_delta is None or neg_delta >= 0:
            issue_counter["negative_delta_not_negative"] += 1
        else:
            negative_deltas.append(neg_delta)
        if row.get("input_connectivity_key") == row.get("positive_connectivity_key"):
            issue_counter["input_equals_positive"] += 1
        if row.get("input_connectivity_key") == row.get("negative_connectivity_key"):
            issue_counter["input_equals_negative"] += 1
        if row.get("positive_connectivity_key") == row.get("negative_connectivity_key"):
            issue_counter["positive_equals_negative"] += 1

        acc = tie_aware_accuracy(pos_sim, neg_sim)
        if acc is not None:
            input_similarity_scores.append(acc)
        ev_acc = tie_aware_accuracy(float(pos_evidence) if pos_evidence is not None else None, float(neg_evidence) if neg_evidence is not None else None)
        if ev_acc is not None:
            evidence_scores.append(ev_acc)

        if pos_sim is not None:
            input_positive_tanimoto_values.append(pos_sim)
        if neg_sim is not None:
            input_negative_tanimoto_values.append(neg_sim)
        if pos_neg_sim is not None:
            positive_negative_tanimoto_values.append(pos_neg_sim)
        if pos_sim is not None and neg_sim is not None and neg_sim >= pos_sim:
            negative_more_input_similar += 1

    return {
        "rows": rows,
        "quality_audit": {
            "issue_counts": dict(sorted(issue_counter.items())),
            "passes_basic_triplet_audit": sum(issue_counter.values()) == 0,
        },
        "baseline_accuracy": {
            "random_expected": 0.5 if rows else None,
            "choose_candidate_more_similar_to_input": compact_number(average(input_similarity_scores), 4),
            "choose_candidate_with_more_evidence": compact_number(average(evidence_scores), 4),
        },
        "difficulty_signals": {
            "negative_at_least_as_similar_to_input_rate": compact_number(negative_more_input_similar / rows, 4)
            if rows
            else None,
            "median_positive_delta": compact_number(median(positive_deltas), 4) if positive_deltas else None,
            "median_negative_delta": compact_number(median(negative_deltas), 4) if negative_deltas else None,
            "median_input_positive_tanimoto": compact_number(median(input_positive_tanimoto_values), 4)
            if input_positive_tanimoto_values
            else None,
            "median_input_negative_tanimoto": compact_number(median(input_negative_tanimoto_values), 4)
            if input_negative_tanimoto_values
            else None,
            "median_positive_negative_tanimoto": compact_number(median(positive_negative_tanimoto_values), 4)
            if positive_negative_tanimoto_values
            else None,
        },
    }


def audit_bindingdb(root: Path) -> dict[str, Any]:
    folder = root / "data" / "bindingdb_target_conditioned"
    stats = read_json(folder / "dataset_stats.json")
    split_reports = {}
    for split in BINDINGDB_SPLITS:
        path = folder / f"{split}.csv"
        if path.is_file():
            split_reports[split] = audit_bindingdb_split(path)

    return {
        "folder": str(folder.relative_to(root)),
        "stats_summary": {
            "final_triplets": stats.get("ranking", {}).get("triplets"),
            "rank_ready_rows": stats.get("rank_ready", {}).get("rank_ready_rows"),
            "aggregated_rows": stats.get("aggregation", {}).get("aggregated_rows"),
            "parsed_observations": stats.get("normalization", {}).get("parsed_observations"),
            "split_rows": stats.get("split", {}).get("rows"),
            "leakage": stats.get("split", {}).get("leakage"),
        },
        "split_reports": split_reports,
    }


def render_table(rows: list[list[Any]]) -> str:
    if not rows:
        return ""
    normalized_rows = [[format_cell(cell) for cell in row] for row in rows]
    widths = [0] * len(rows[0])
    for row in normalized_rows:
        for idx, cell in enumerate(row):
            widths[idx] = max(widths[idx], len(str(cell)))
    rendered = []
    for ridx, row in enumerate(normalized_rows):
        rendered.append("| " + " | ".join(str(cell).ljust(widths[idx]) for idx, cell in enumerate(row)) + " |")
        if ridx == 0:
            rendered.append("| " + " | ".join("-" * widths[idx] for idx in range(len(widths))) + " |")
    return "\n".join(rendered)


def format_cell(value: Any) -> str:
    if value is None:
        return "n/a"
    if value is True:
        return "pass"
    if value is False:
        return "fail"
    return str(value)


def value_at(mapping: dict[str, Any], *keys: str, default: Any = "") -> Any:
    current: Any = mapping
    for key in keys:
        if not isinstance(current, dict):
            return default
        current = current.get(key, default)
    return current


def render_report(payload: dict[str, Any]) -> str:
    generated = payload["generated_at_utc"]
    dataset_rows = [["Dataset", "Queries", "Positive answers", "Train", "Val", "Test", "Audit"]]
    for dataset in payload["multiproperty_datasets"]:
        file_stats = dataset["file_stats"]
        dataset_rows.append(
            [
                dataset["name"],
                value_at(file_stats, "all", "queries"),
                value_at(file_stats, "all", "positive_answers"),
                value_at(file_stats, "train", "queries"),
                value_at(file_stats, "val", "queries"),
                value_at(file_stats, "test", "queries"),
                "pass" if dataset["schema_audit"]["passes_basic_shape_audit"] else "check",
            ]
        )

    baseline_rows = [["Split", "Rows", "Random", "Input-sim", "Evidence", "Hard-neg rate"]]
    for split, report in payload["bindingdb"]["split_reports"].items():
        baseline_rows.append(
            [
                split,
                report["rows"],
                report["baseline_accuracy"]["random_expected"],
                report["baseline_accuracy"]["choose_candidate_more_similar_to_input"],
                report["baseline_accuracy"]["choose_candidate_with_more_evidence"],
                report["difficulty_signals"]["negative_at_least_as_similar_to_input_rate"],
            ]
        )

    quality_rows = [["Dataset", "Median improvement", "Median Tanimoto", "Descriptor success", "Leakage"]]
    for dataset in payload["multiproperty_datasets"]:
        leakage = dataset["leakage"]
        if isinstance(leakage, dict):
            leakage_text = leakage.get("passes_no_molecule_overlap")
        else:
            leakage_text = leakage
        if isinstance(leakage_text, str) and leakage_text.lower() == "passed":
            leakage_text = True
        quality_rows.append(
            [
                dataset["name"],
                dataset["label_quality_signals"]["median_answer_primary_improvement"],
                dataset["label_quality_signals"]["median_answer_tanimoto"],
                dataset["label_quality_signals"]["descriptor_success_rate"],
                leakage_text,
            ]
        )

    return f"""# Meeting Follow-up Report

Generated at UTC: {generated}

## Boss Requirements Parsed

- Build the datasets quickly, but keep source-specific quality boundaries clear.
- Install and run baselines, then report where each baseline works and fails.
- Treat rare/experimental data as the main contribution; filtering must be accurate and defensible.
- Keep the BindingDB target-conditioned direction because it supports protein context, ranking, and hard negatives.
- Do not overclaim wet-lab validation yet; prepare candidate-selection outputs that can later be handed to a lab.

## Current Dataset Inventory

{render_table(dataset_rows)}

BindingDB target-conditioned dataset:

- Parsed observations: {value_at(payload, "bindingdb", "stats_summary", "parsed_observations")}
- Aggregated ligand-target-measurement rows: {value_at(payload, "bindingdb", "stats_summary", "aggregated_rows")}
- Rank-ready rows: {value_at(payload, "bindingdb", "stats_summary", "rank_ready_rows")}
- Final triplets: {value_at(payload, "bindingdb", "stats_summary", "final_triplets")}

## Data Quality Signals

{render_table(quality_rows)}

Notes:

- ChEMBL, PubChem, and Papyrus all expose the requested `1 input + 3 property groups -> 2 positive answers` shape.
- PubChem and Papyrus are intentionally not merged into ChEMBL. They should stay separate until target/assay enrichment and source-tier weighting are finalized.
- Papyrus currently stores positive answer IDs/connectivity keys rather than full answer evidence objects, so its per-answer delta/Tanimoto audit is limited compared with ChEMBL/PubChem.

## First BindingDB Baselines

The BindingDB task already has a positive and hard-negative candidate in each row, so simple candidate-ranking baselines can be measured immediately.

{render_table(baseline_rows)}

Interpretation:

- `Input-sim` chooses the candidate with higher Tanimoto similarity to the input ligand.
- `Evidence` chooses the candidate with more source evidence records.
- `Hard-neg rate` is the fraction of rows where the negative is at least as similar to the input as the positive. Higher values mean harder structural distractors.
- These baselines do not use activity labels. The p-value fields are only used for audit checks, not for candidate selection.

## Immediate Next Steps

1. Use this report as the next meeting status snapshot.
2. Add a non-cheating MMP/retrieval baseline for ChEMBL/PubChem/Papyrus that predicts from train-only molecules.
3. Add target/assay enrichment for PubChem and Papyrus before any combined benchmark claim.
4. Keep BindingDB as the direction-2 baseline track: compare random, similarity, evidence, and then protein-conditioned learned ranking.
5. Prepare a short candidate-selection interface for future wet-lab discussion only after baseline failure modes are clear.
"""


def main() -> int:
    args = parse_args()
    root = args.root.resolve()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    payload = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "root": str(root),
        "multiproperty_datasets": [
            audit_multiproperty_dataset(root, display_name, folder_name)
            for display_name, folder_name in MULTIPROPERTY_DATASETS
        ],
        "bindingdb": audit_bindingdb(root),
    }

    json_path = out_dir / "meeting_followup_report.json"
    md_path = out_dir / "meeting_followup_report.md"
    json_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    md_path.write_text(render_report(payload), encoding="utf-8")

    print(json.dumps({"json": str(json_path), "markdown": str(md_path)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

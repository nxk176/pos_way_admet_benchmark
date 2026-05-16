from __future__ import annotations

import argparse
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


DEFAULT_DATASETS = {
    "admet_edit_merged": {
        "dir": Path("data/llm_sft/admet_edit_3prop_merged"),
        "splits": ["train", "val", "test"],
    },
    "bindingdb_target_conditioned": {
        "dir": Path("data/llm_sft_bindingdb_full/bindingdb_target_conditioned"),
        "splits": ["train", "val", "test_seen_target", "test_unseen_target"],
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit POS-WAY SFT datasets for merge fairness and leakage indicators.")
    parser.add_argument("--out", type=Path, default=Path("reports/server_training/dataset_audit.json"))
    return parser.parse_args()


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        raise SystemExit(f"ERROR: missing dataset split: {path}")
    rows = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def row_molecules(row: dict[str, Any]) -> set[str]:
    out = set()
    input_smiles = str(row.get("input_smiles") or "").strip()
    if input_smiles:
        out.add(input_smiles)
    gold = row.get("gold_smiles") or []
    if isinstance(gold, str):
        gold = [gold]
    for value in gold:
        text = str(value or "").strip()
        if text:
            out.add(text)
    hard_negative = str(row.get("hard_negative_smiles") or "").strip()
    if hard_negative:
        out.add(hard_negative)
    return out


def prompt_chars(row: dict[str, Any]) -> int:
    messages = row.get("messages") or []
    if not isinstance(messages, list):
        return 0
    return sum(len(str(msg.get("content") or "")) for msg in messages if isinstance(msg, dict))


def split_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    source_counts = Counter(str(row.get("source_dataset") or "unknown") for row in rows)
    task_counts = Counter(str(row.get("task") or "unknown") for row in rows)
    prompt_lengths = [prompt_chars(row) for row in rows]
    return {
        "rows": len(rows),
        "source_counts": dict(sorted(source_counts.items())),
        "task_counts": dict(sorted(task_counts.items())),
        "mean_prompt_chars": round(sum(prompt_lengths) / len(prompt_lengths), 2) if prompt_lengths else 0.0,
        "max_prompt_chars": max(prompt_lengths) if prompt_lengths else 0,
    }


def overlap_count(split_sets: dict[str, set[str]], left: str, right: str) -> int:
    return len(split_sets.get(left, set()) & split_sets.get(right, set()))


def audit_dataset(name: str, spec: dict[str, Any]) -> dict[str, Any]:
    dataset_dir = spec["dir"]
    splits = spec["splits"]
    rows_by_split = {split: read_jsonl(dataset_dir / f"{split}.jsonl") for split in splits}
    molecule_sets = {
        split: set().union(*(row_molecules(row) for row in rows)) if rows else set()
        for split, rows in rows_by_split.items()
    }
    target_sets = {
        split: {str(row.get("target_id") or "").strip() for row in rows if str(row.get("target_id") or "").strip()}
        for split, rows in rows_by_split.items()
    }

    pairwise_molecule_overlap = {}
    pairwise_target_overlap = {}
    for idx, left in enumerate(splits):
        for right in splits[idx + 1 :]:
            pairwise_molecule_overlap[f"{left}__{right}"] = overlap_count(molecule_sets, left, right)
            if target_sets[left] or target_sets[right]:
                pairwise_target_overlap[f"{left}__{right}"] = overlap_count(target_sets, left, right)

    warnings = []
    train_overlaps = {
        key: value for key, value in pairwise_molecule_overlap.items() if key.startswith("train__") and value
    }
    if train_overlaps:
        warnings.append({"level": "high", "message": "Train/held-out molecule overlap detected.", "details": train_overlaps})
    if name == "admet_edit_merged":
        test_counts = split_summary(rows_by_split["test"])["source_counts"]
        if len(test_counts) > 1:
            total = sum(test_counts.values())
            max_share = max(test_counts.values()) / total if total else 0.0
            if max_share > 0.6:
                warnings.append(
                    {
                        "level": "medium",
                        "message": "Merged ADMET test set is source-imbalanced; report per-source and source-macro metrics.",
                        "details": test_counts,
                    }
                )
    if name == "bindingdb_target_conditioned":
        unseen_target_overlap = {
            key: value
            for key, value in pairwise_target_overlap.items()
            if "test_unseen_target" in key and value
        }
        if unseen_target_overlap:
            warnings.append(
                {
                    "level": "high",
                    "message": "Unseen-target split has target overlap with another split.",
                    "details": unseen_target_overlap,
                }
            )

    return {
        "dataset_dir": str(dataset_dir),
        "splits": {split: split_summary(rows) for split, rows in rows_by_split.items()},
        "pairwise_molecule_overlap_counts": pairwise_molecule_overlap,
        "pairwise_target_overlap_counts": pairwise_target_overlap,
        "warnings": warnings,
    }


def main() -> int:
    args = parse_args()
    payload = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "recommendation": (
            "Merged ADMET training is reasonable because the task schema is shared, but reporting must include "
            "per-source and source-macro metrics. BindingDB remains a separate task and should be reported with "
            "seen-target and unseen-target splits."
        ),
        "datasets": {name: audit_dataset(name, spec) for name, spec in DEFAULT_DATASETS.items()},
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

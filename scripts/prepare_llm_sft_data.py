from __future__ import annotations

import argparse
import csv
import json
import random
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


EDIT_SOURCES = [
    ("chembl", Path("data/chembl_3prop_2pos")),
    ("pubchem", Path("data/pubchem_3prop_2pos")),
    ("papyrus", Path("data/papyrus_3prop_2pos")),
]

EDIT_SPLITS = ["train", "val", "test"]
BINDINGDB_SPLITS = ["train", "val", "test_seen_target", "test_unseen_target"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare chat-style JSONL files for LLM SFT/LoRA. "
            "The default merge policy keeps ChEMBL/PubChem/Papyrus together because they share "
            "the same 3-property/2-positive edit schema, and keeps BindingDB separate because "
            "it is target-conditioned ranking/triplet data."
        )
    )
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--out-dir", type=Path, default=Path("data/llm_sft"))
    parser.add_argument(
        "--preset",
        choices=["edit-merged", "bindingdb", "all"],
        default="all",
        help="Which SFT datasets to write.",
    )
    parser.add_argument("--max-train-per-source", type=int, default=0)
    parser.add_argument("--max-val-per-source", type=int, default=0)
    parser.add_argument("--max-test-per-source", type=int, default=0)
    parser.add_argument("--seed", type=int, default=29)
    return parser.parse_args()


def clean_text(value: Any) -> str:
    return str(value or "").replace("\r", " ").replace("\n", " ").strip()


def safe_json(value: Any, fallback: Any) -> Any:
    try:
        return json.loads(value)
    except (TypeError, json.JSONDecodeError):
        return fallback


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise SystemExit(f"ERROR: missing input file: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_jsonl(path: Path, rows: Iterable[dict[str, Any]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False, separators=(",", ":")) + "\n")
            count += 1
    return count


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def sample_rows(rows: list[dict[str, str]], limit: int, seed: int, label: str) -> list[dict[str, str]]:
    if not limit or len(rows) <= limit:
        return rows
    rng = random.Random(f"{seed}:{label}")
    selected = rows[:]
    rng.shuffle(selected)
    selected = selected[:limit]
    selected.sort(key=lambda row: clean_text(row.get("query_id") or row.get("sample_id")))
    return selected


def positive_smiles(row: dict[str, str]) -> list[str]:
    values = safe_json(row.get("positive_answer_smiles_json"), [])
    if not isinstance(values, list):
        return []
    return [clean_text(value) for value in values if clean_text(value)]


def edit_molecule_keys(row: dict[str, str]) -> set[str]:
    keys = {clean_text(row.get("input_smiles_canon")), clean_text(row.get("input_connectivity_key"))}
    keys.update(positive_smiles(row))
    return {key for key in keys if key}


def edit_system_prompt() -> str:
    return (
        "You are a careful medicinal chemistry molecule editor. "
        "Return only valid JSON and no explanation."
    )


def edit_user_prompt(row: dict[str, str], source: str) -> str:
    return (
        "Task: edit the input molecule according to the instruction.\n"
        "Return exactly this JSON shape:\n"
        '{"edited_smiles":["SMILES_1","SMILES_2"]}\n\n'
        f"Input SMILES: {clean_text(row.get('input_smiles_canon'))}\n"
        f"Instruction: {clean_text(row.get('instruction'))}\n\n"
        "Rules:\n"
        "- Return exactly 2 edited SMILES strings.\n"
        "- Do not return the unchanged input molecule.\n"
        "- Prefer local, chemically plausible edits.\n"
        "- Preserve the requested secondary property and local constraints."
    )


def edit_assistant(row: dict[str, str]) -> str:
    return json.dumps({"edited_smiles": positive_smiles(row)[:2]}, ensure_ascii=False, separators=(",", ":"))


def edit_sft_row(row: dict[str, str], source: str) -> dict[str, Any]:
    return {
        "id": f"{source}:{clean_text(row.get('query_id'))}",
        "task": "admet_3property_2positive_edit",
        "source_dataset": source,
        "split": clean_text(row.get("split")),
        "input_smiles": clean_text(row.get("input_smiles_canon")),
        "gold_smiles": positive_smiles(row)[:2],
        "primary_endpoint": clean_text(row.get("primary_endpoint")),
        "primary_objective": safe_json(row.get("primary_objective_json"), {}),
        "preserved_property": safe_json(row.get("preserved_property_json"), {}),
        "local_constraints": safe_json(row.get("local_constraints_json"), {}),
        "positive_answers": safe_json(row.get("positive_answers_json"), []),
        "source_positive_sample_ids": safe_json(row.get("source_positive_sample_ids_json"), []),
        "messages": [
            {"role": "system", "content": edit_system_prompt()},
            {"role": "user", "content": edit_user_prompt(row, source)},
            {"role": "assistant", "content": edit_assistant(row)},
        ],
    }


def load_edit_rows(repo_root: Path, seed: int, args: argparse.Namespace) -> dict[str, list[tuple[str, dict[str, str]]]]:
    by_split: dict[str, list[tuple[str, dict[str, str]]]] = {split: [] for split in EDIT_SPLITS}
    for source, rel_dir in EDIT_SOURCES:
        data_dir = repo_root / rel_dir
        for split in EDIT_SPLITS:
            rows = [
                row
                for row in read_csv(data_dir / f"{split}.csv")
                if clean_text(row.get("query_id")) and len(positive_smiles(row)) >= 2
            ]
            limit = {
                "train": args.max_train_per_source,
                "val": args.max_val_per_source,
                "test": args.max_test_per_source,
            }[split]
            rows = sample_rows(rows, limit, seed, f"{source}:{split}")
            by_split[split].extend((source, row) for row in rows)
    return by_split


def remove_edit_val_test_overlap(
    by_split: dict[str, list[tuple[str, dict[str, str]]]]
) -> tuple[dict[str, list[tuple[str, dict[str, str]]]], dict[str, Any]]:
    val_keys: set[str] = set()
    for _, row in by_split["val"]:
        val_keys.update(edit_molecule_keys(row))

    kept_test: list[tuple[str, dict[str, str]]] = []
    dropped_test: list[dict[str, str]] = []
    for source, row in by_split["test"]:
        overlap = edit_molecule_keys(row) & val_keys
        if overlap:
            dropped_test.append(
                {
                    "source": source,
                    "query_id": clean_text(row.get("query_id")),
                    "overlap_count": str(len(overlap)),
                }
            )
        else:
            kept_test.append((source, row))

    cleaned = dict(by_split)
    cleaned["test"] = kept_test
    audit = {
        "policy": "Kept merged validation rows and dropped merged test rows whose input/gold molecules overlap validation.",
        "dropped_test_rows": len(dropped_test),
        "remaining_test_rows": len(kept_test),
        "dropped_test_examples_preview": dropped_test[:25],
    }
    return cleaned, audit


def remove_edit_train_leakage(
    by_split: dict[str, list[tuple[str, dict[str, str]]]]
) -> tuple[dict[str, list[tuple[str, dict[str, str]]]], dict[str, Any]]:
    heldout_keys: set[str] = set()
    for split in ["val", "test"]:
        for _, row in by_split[split]:
            heldout_keys.update(edit_molecule_keys(row))

    kept_train: list[tuple[str, dict[str, str]]] = []
    dropped_train: list[dict[str, str]] = []
    for source, row in by_split["train"]:
        overlap = edit_molecule_keys(row) & heldout_keys
        if overlap:
            dropped_train.append(
                {
                    "source": source,
                    "query_id": clean_text(row.get("query_id")),
                    "overlap_count": str(len(overlap)),
                }
            )
        else:
            kept_train.append((source, row))

    cleaned = dict(by_split)
    cleaned["train"] = kept_train

    split_key_sets: dict[str, set[str]] = defaultdict(set)
    for split, items in cleaned.items():
        for _, row in items:
            split_key_sets[split].update(edit_molecule_keys(row))

    leakage = {
        "policy": "Dropped merged train rows whose input or gold-output molecule identifiers overlap any merged val/test row.",
        "dropped_train_rows": len(dropped_train),
        "remaining_train_rows": len(kept_train),
        "heldout_molecule_keys": len(heldout_keys),
        "post_filter_overlap_counts": {
            "train_val": len(split_key_sets["train"] & split_key_sets["val"]),
            "train_test": len(split_key_sets["train"] & split_key_sets["test"]),
            "val_test": len(split_key_sets["val"] & split_key_sets["test"]),
        },
        "dropped_train_examples_preview": dropped_train[:25],
    }
    return cleaned, leakage


def write_edit_merged(repo_root: Path, out_dir: Path, args: argparse.Namespace) -> dict[str, Any]:
    by_split = load_edit_rows(repo_root, args.seed, args)
    by_split, heldout_overlap = remove_edit_val_test_overlap(by_split)
    by_split, leakage = remove_edit_train_leakage(by_split)

    dataset_dir = out_dir / "admet_edit_3prop_merged"
    stats: dict[str, Any] = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "dataset": "admet_edit_3prop_merged",
        "merge_policy": (
            "ChEMBL, PubChem, and Papyrus are merged because they share the same "
            "input SMILES + 3-property instruction -> exactly 2 edited SMILES schema. "
            "Source labels are retained in metadata but not injected into the model prompt."
        ),
        "heldout_overlap": heldout_overlap,
        "leakage": leakage,
        "splits": {},
    }
    for split in EDIT_SPLITS:
        rows = [edit_sft_row(row, source) for source, row in by_split[split]]
        out_path = dataset_dir / f"{split}.jsonl"
        count = write_jsonl(out_path, rows)
        source_counts = Counter(source for source, _ in by_split[split])
        stats["splits"][split] = {
            "path": str(out_path),
            "rows": count,
            "source_counts": dict(sorted(source_counts.items())),
        }
    write_json(dataset_dir / "dataset_stats.json", stats)
    return stats


def bindingdb_system_prompt() -> str:
    return (
        "You are a careful medicinal chemistry assistant for target-conditioned ligand optimization. "
        "Return only valid JSON and no explanation."
    )


def bindingdb_user_prompt(row: dict[str, str]) -> str:
    target_name = clean_text(row.get("target_name"))
    organism = clean_text(row.get("target_organism"))
    measurement_type = clean_text(row.get("measurement_type"))
    measurement_group = clean_text(row.get("measurement_group"))
    sequence = clean_text(row.get("target_sequence"))
    if len(sequence) > 800:
        sequence = sequence[:800] + "...[truncated]"
    return (
        "Task: propose a structurally related ligand with stronger measured activity for the target.\n"
        "Return exactly this JSON shape:\n"
        '{"stronger_smiles":"SMILES","hard_negative_smiles":"SMILES"}\n\n'
        f"Input SMILES: {clean_text(row.get('input_smiles'))}\n"
        f"Target: {target_name}\n"
        f"Organism: {organism}\n"
        f"Measurement: {measurement_type} ({measurement_group}); higher p-scale means stronger.\n"
        f"Instruction: {clean_text(row.get('instruction'))}\n"
        f"Target sequence prefix: {sequence}\n\n"
        "Rules:\n"
        "- The stronger SMILES should be structurally related to the input ligand.\n"
        "- The hard negative is a related ligand with weaker measured activity.\n"
        "- Do not explain."
    )


def bindingdb_assistant(row: dict[str, str]) -> str:
    return json.dumps(
        {
            "stronger_smiles": clean_text(row.get("positive_smiles")),
            "hard_negative_smiles": clean_text(row.get("negative_smiles")),
        },
        ensure_ascii=False,
        separators=(",", ":"),
    )


def bindingdb_sft_row(row: dict[str, str], split: str) -> dict[str, Any]:
    return {
        "id": clean_text(row.get("sample_id")),
        "task": "bindingdb_target_conditioned_triplet",
        "source_dataset": "bindingdb",
        "split": split,
        "input_smiles": clean_text(row.get("input_smiles")),
        "target_id": clean_text(row.get("target_id")),
        "target_name": clean_text(row.get("target_name")),
        "measurement_type": clean_text(row.get("measurement_type")),
        "gold_smiles": [clean_text(row.get("positive_smiles"))],
        "hard_negative_smiles": clean_text(row.get("negative_smiles")),
        "messages": [
            {"role": "system", "content": bindingdb_system_prompt()},
            {"role": "user", "content": bindingdb_user_prompt(row)},
            {"role": "assistant", "content": bindingdb_assistant(row)},
        ],
    }


def write_bindingdb(repo_root: Path, out_dir: Path, args: argparse.Namespace) -> dict[str, Any]:
    source_dir = repo_root / "data/bindingdb_target_conditioned"
    dataset_dir = out_dir / "bindingdb_target_conditioned"
    stats: dict[str, Any] = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "dataset": "bindingdb_target_conditioned",
        "merge_policy": (
            "Kept separate from ADMET edit data because the input includes protein target context "
            "and the supervised answer is a stronger ligand plus hard negative."
        ),
        "splits": {},
    }
    for split in BINDINGDB_SPLITS:
        rows = read_csv(source_dir / f"{split}.csv")
        limit = args.max_train_per_source if split == "train" else args.max_val_per_source
        if split.startswith("test"):
            limit = args.max_test_per_source
        rows = sample_rows(rows, limit, args.seed, f"bindingdb:{split}")
        sft_rows = [bindingdb_sft_row(row, split) for row in rows if clean_text(row.get("positive_smiles"))]
        out_path = dataset_dir / f"{split}.jsonl"
        count = write_jsonl(out_path, sft_rows)
        stats["splits"][split] = {"path": str(out_path), "rows": count}
    write_json(dataset_dir / "dataset_stats.json", stats)
    return stats


def main() -> int:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    out_dir = (repo_root / args.out_dir).resolve()

    summary: dict[str, Any] = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "repo_root": str(repo_root),
        "out_dir": str(out_dir),
        "preset": args.preset,
        "outputs": {},
    }
    if args.preset in {"edit-merged", "all"}:
        summary["outputs"]["edit_merged"] = write_edit_merged(repo_root, out_dir, args)
    if args.preset in {"bindingdb", "all"}:
        summary["outputs"]["bindingdb"] = write_bindingdb(repo_root, out_dir, args)

    write_json(out_dir / "prepare_llm_sft_summary.json", summary)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

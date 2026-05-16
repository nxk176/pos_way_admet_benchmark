from __future__ import annotations

import argparse
import csv
import json
import random
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Any

try:
    from rdkit import Chem, DataStructs
except Exception as exc:  # noqa: BLE001
    raise SystemExit(
        "ERROR: build_strict_splits_and_triplets.py requires RDKit. "
        "Run it with the project environment, for example .\\myenv311\\Scripts\\python.exe"
    ) from exc


TRIPLET_FIELDNAMES = [
    "triplet_id",
    "split",
    "hard_negative_scope",
    "primary_endpoint",
    "condition_bucket",
    "mmp_core",
    "instruction",
    "input_chembl_id",
    "input_smiles_canon",
    "input_connectivity_key",
    "positive_sample_id",
    "positive_target_chembl_id",
    "positive_target_smiles_canon",
    "positive_target_connectivity_key",
    "positive_value_before",
    "positive_value_after",
    "positive_delta_value",
    "positive_improvement_score",
    "positive_input_target_tanimoto",
    "negative_sample_id",
    "negative_target_chembl_id",
    "negative_target_smiles_canon",
    "negative_target_connectivity_key",
    "negative_value_before",
    "negative_value_after",
    "negative_delta_value",
    "negative_improvement_score",
    "negative_input_target_tanimoto",
    "positive_negative_tanimoto",
    "improvement_margin",
    "positive_primary_success",
    "positive_secondary_success",
    "negative_primary_success",
    "negative_secondary_success",
]


@dataclass
class Component:
    component_id: int
    nodes: set[str]
    row_indices: list[int]
    row_count: int
    positive_count: int
    negative_count: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create leakage-aware strict splits and positive-negative triplets.")
    parser.add_argument(
        "--samples",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/multiproperty_expanded_strict/fragment_multiproperty_samples.csv"),
    )
    parser.add_argument(
        "--splits-dir",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/multiproperty_expanded_strict/splits"),
    )
    parser.add_argument(
        "--triplets-dir",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/ranking_strict"),
    )
    parser.add_argument("--train-ratio", type=float, default=0.8)
    parser.add_argument("--val-ratio", type=float, default=0.1)
    parser.add_argument("--test-ratio", type=float, default=0.1)
    parser.add_argument("--max-triplets-per-positive", type=int, default=1)
    parser.add_argument("--min-positive-negative-tanimoto", type=float, default=0.35)
    parser.add_argument("--allow-core-fallback", action="store_true")
    parser.add_argument("--seed", type=int, default=29)
    return parser.parse_args()


def read_rows(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    return rows, fieldnames


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
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


def as_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def bool_text(value: str) -> bool:
    return str(value).strip().lower() == "true"


def improvement_score(row: dict[str, str]) -> float:
    delta = as_float(row.get("delta_value"))
    instruction = row.get("instruction", "").lower()
    if instruction.startswith("decrease"):
        return -delta
    return delta


def fingerprint(smiles: str, cache: dict[str, Any]) -> Any:
    if smiles in cache:
        return cache[smiles]
    mol = Chem.MolFromSmiles(smiles)
    fp = Chem.RDKFingerprint(mol) if mol is not None else None
    cache[smiles] = fp
    return fp


def tanimoto(smiles_a: str, smiles_b: str, cache: dict[str, Any]) -> float:
    fp_a = fingerprint(smiles_a, cache)
    fp_b = fingerprint(smiles_b, cache)
    if fp_a is None or fp_b is None:
        return 0.0
    return float(DataStructs.FingerprintSimilarity(fp_a, fp_b))


def build_components(rows: list[dict[str, str]]) -> tuple[list[Component], dict[str, int]]:
    adjacency: dict[str, set[str]] = defaultdict(set)
    node_rows: dict[str, list[int]] = defaultdict(list)
    for idx, row in enumerate(rows):
        source = row.get("input_connectivity_key", "")
        target = row.get("target_connectivity_key", "")
        if not source or not target:
            continue
        adjacency[source].add(target)
        adjacency[target].add(source)
        node_rows[source].append(idx)
        node_rows[target].append(idx)

    components: list[Component] = []
    node_to_component: dict[str, int] = {}
    seen: set[str] = set()
    for node in sorted(adjacency):
        if node in seen:
            continue
        component_id = len(components)
        queue: deque[str] = deque([node])
        seen.add(node)
        nodes: set[str] = set()
        row_indices: set[int] = set()
        while queue:
            current = queue.popleft()
            nodes.add(current)
            node_to_component[current] = component_id
            row_indices.update(node_rows[current])
            for neighbor in adjacency[current]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)
        indices = sorted(row_indices)
        counts = Counter(rows[idx].get("sample_type", "") for idx in indices)
        components.append(
            Component(
                component_id=component_id,
                nodes=nodes,
                row_indices=indices,
                row_count=len(indices),
                positive_count=counts["positive"],
                negative_count=counts["negative"],
            )
        )
    return components, node_to_component


def assign_splits(components: list[Component], args: argparse.Namespace) -> dict[int, str]:
    ratios = {"train": args.train_ratio, "val": args.val_ratio, "test": args.test_ratio}
    total_ratio = sum(ratios.values())
    ratios = {key: value / total_ratio for key, value in ratios.items()}
    total_rows = sum(component.row_count for component in components)
    targets = {split: total_rows * ratio for split, ratio in ratios.items()}
    split_counts = {split: 0 for split in ratios}
    assignments: dict[int, str] = {}

    rng = random.Random(args.seed)
    shuffled = list(components)
    rng.shuffle(shuffled)
    shuffled.sort(key=lambda component: component.row_count, reverse=True)
    for component in shuffled:
        split = min(
            targets,
            key=lambda candidate: (split_counts[candidate] + component.row_count) / max(targets[candidate], 1.0),
        )
        assignments[component.component_id] = split
        split_counts[split] += component.row_count
    return assignments


def split_rows(
    rows: list[dict[str, str]],
    components: list[Component],
    assignments: dict[int, str],
) -> tuple[dict[str, list[dict[str, str]]], dict[int, str]]:
    row_to_split: dict[int, str] = {}
    split_payloads: dict[str, list[dict[str, str]]] = {"train": [], "val": [], "test": []}
    for component in components:
        split = assignments[component.component_id]
        for idx in component.row_indices:
            row_to_split[idx] = split
    for idx, row in enumerate(rows):
        split = row_to_split.get(idx, "train")
        out = dict(row)
        out["split"] = split
        split_payloads[split].append(out)
    return split_payloads, row_to_split


def choose_hard_negative(
    positive: dict[str, str],
    candidates: list[dict[str, str]],
    fp_cache: dict[str, Any],
    min_tanimoto: float,
) -> tuple[dict[str, str], float] | None:
    scored: list[tuple[float, float, dict[str, str]]] = []
    positive_smiles = positive.get("target_smiles_canon", "")
    positive_improvement = improvement_score(positive)
    for candidate in candidates:
        if candidate.get("sample_id") == positive.get("sample_id"):
            continue
        if candidate.get("target_connectivity_key") == positive.get("target_connectivity_key"):
            continue
        sim = tanimoto(positive_smiles, candidate.get("target_smiles_canon", ""), fp_cache)
        if sim < min_tanimoto:
            continue
        margin = positive_improvement - improvement_score(candidate)
        scored.append((sim, margin, candidate))
    if not scored:
        return None
    scored.sort(key=lambda item: (item[0], item[1]), reverse=True)
    sim, _margin, candidate = scored[0]
    return candidate, sim


def build_triplets(
    split_payloads: dict[str, list[dict[str, str]]],
    args: argparse.Namespace,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    fp_cache: dict[str, Any] = {}
    triplets: list[dict[str, Any]] = []
    counters: Counter[str] = Counter()

    for split, rows in split_payloads.items():
        negatives_by_input: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
        negatives_by_core: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
        positives = []
        for row in rows:
            if row.get("sample_type") == "positive":
                positives.append(row)
                continue
            if row.get("sample_type") != "negative":
                continue
            key_input = (
                row.get("primary_endpoint", ""),
                row.get("condition_bucket", ""),
                row.get("input_connectivity_key", ""),
                row.get("mmp_core", ""),
            )
            key_core = (
                row.get("primary_endpoint", ""),
                row.get("condition_bucket", ""),
                row.get("mmp_core", ""),
            )
            negatives_by_input[key_input].append(row)
            negatives_by_core[key_core].append(row)

        for positive in positives:
            input_key = (
                positive.get("primary_endpoint", ""),
                positive.get("condition_bucket", ""),
                positive.get("input_connectivity_key", ""),
                positive.get("mmp_core", ""),
            )
            core_key = (
                positive.get("primary_endpoint", ""),
                positive.get("condition_bucket", ""),
                positive.get("mmp_core", ""),
            )
            picked = choose_hard_negative(
                positive,
                negatives_by_input.get(input_key, []),
                fp_cache,
                args.min_positive_negative_tanimoto,
            )
            scope = "same_input_condition_core"
            if picked is None and args.allow_core_fallback:
                picked = choose_hard_negative(
                    positive,
                    negatives_by_core.get(core_key, []),
                    fp_cache,
                    args.min_positive_negative_tanimoto,
                )
                scope = "same_condition_core"
            if picked is None:
                counters[f"{split}:unmatched_positive"] += 1
                continue

            negative, positive_negative_similarity = picked
            pos_improvement = improvement_score(positive)
            neg_improvement = improvement_score(negative)
            triplet_idx = len(triplets) + 1
            triplets.append(
                {
                    "triplet_id": f"strict_triplet_{triplet_idx:08d}",
                    "split": split,
                    "hard_negative_scope": scope,
                    "primary_endpoint": positive.get("primary_endpoint", ""),
                    "condition_bucket": positive.get("condition_bucket", ""),
                    "mmp_core": positive.get("mmp_core", ""),
                    "instruction": positive.get("instruction", ""),
                    "input_chembl_id": positive.get("input_chembl_id", ""),
                    "input_smiles_canon": positive.get("input_smiles_canon", ""),
                    "input_connectivity_key": positive.get("input_connectivity_key", ""),
                    "positive_sample_id": positive.get("sample_id", ""),
                    "positive_target_chembl_id": positive.get("target_chembl_id", ""),
                    "positive_target_smiles_canon": positive.get("target_smiles_canon", ""),
                    "positive_target_connectivity_key": positive.get("target_connectivity_key", ""),
                    "positive_value_before": positive.get("value_before", ""),
                    "positive_value_after": positive.get("value_after", ""),
                    "positive_delta_value": positive.get("delta_value", ""),
                    "positive_improvement_score": round(pos_improvement, 6),
                    "positive_input_target_tanimoto": positive.get("tanimoto_similarity", ""),
                    "negative_sample_id": negative.get("sample_id", ""),
                    "negative_target_chembl_id": negative.get("target_chembl_id", ""),
                    "negative_target_smiles_canon": negative.get("target_smiles_canon", ""),
                    "negative_target_connectivity_key": negative.get("target_connectivity_key", ""),
                    "negative_value_before": negative.get("value_before", ""),
                    "negative_value_after": negative.get("value_after", ""),
                    "negative_delta_value": negative.get("delta_value", ""),
                    "negative_improvement_score": round(neg_improvement, 6),
                    "negative_input_target_tanimoto": negative.get("tanimoto_similarity", ""),
                    "positive_negative_tanimoto": round(positive_negative_similarity, 6),
                    "improvement_margin": round(pos_improvement - neg_improvement, 6),
                    "positive_primary_success": positive.get("primary_success", ""),
                    "positive_secondary_success": positive.get("secondary_success", ""),
                    "negative_primary_success": negative.get("primary_success", ""),
                    "negative_secondary_success": negative.get("secondary_success", ""),
                }
            )
            counters[f"{split}:{scope}"] += 1
    stats = {
        "num_triplets": len(triplets),
        "counts": dict(sorted(counters.items())),
        "parameters": {
            "min_positive_negative_tanimoto": args.min_positive_negative_tanimoto,
            "allow_core_fallback": args.allow_core_fallback,
            "max_triplets_per_positive": args.max_triplets_per_positive,
        },
    }
    return triplets, stats


def split_stats(
    split_payloads: dict[str, list[dict[str, str]]],
    components: list[Component],
    assignments: dict[int, str],
) -> dict[str, Any]:
    payload: dict[str, Any] = {}
    molecule_sets: dict[str, set[str]] = {}
    for split, rows in split_payloads.items():
        molecules = {
            value
            for row in rows
            for value in [row.get("input_connectivity_key", ""), row.get("target_connectivity_key", "")]
            if value
        }
        molecule_sets[split] = molecules
        counts = Counter(row.get("sample_type", "") for row in rows)
        payload[split] = {
            "rows": len(rows),
            "positive": counts["positive"],
            "negative": counts["negative"],
            "unique_molecules": len(molecules),
            "unique_primary_endpoints": len({row.get("primary_endpoint", "") for row in rows}),
            "components": sum(1 for component in components if assignments[component.component_id] == split),
        }

    overlaps = {}
    for left, right in [("train", "val"), ("train", "test"), ("val", "test")]:
        overlaps[f"{left}_{right}"] = len(molecule_sets[left] & molecule_sets[right])
    payload["leakage_check"] = {
        "molecule_overlap_counts": overlaps,
        "passes_no_molecule_overlap": all(value == 0 for value in overlaps.values()),
    }
    payload["component_stats"] = {
        "num_components": len(components),
        "largest_components_by_molecule_count": sorted((len(component.nodes) for component in components), reverse=True)[:20],
        "largest_components_by_row_count": sorted((component.row_count for component in components), reverse=True)[:20],
    }
    return payload


def main() -> int:
    args = parse_args()
    if not args.samples.is_file():
        raise SystemExit(f"ERROR: samples file not found: {args.samples}")

    rows, fieldnames = read_rows(args.samples)
    components, _node_to_component = build_components(rows)
    assignments = assign_splits(components, args)
    split_payloads, _row_to_split = split_rows(rows, components, assignments)

    split_fieldnames = [*fieldnames, "split"] if "split" not in fieldnames else fieldnames
    for split in ["train", "val", "test"]:
        write_csv(args.splits_dir / f"{split}.csv", split_payloads[split], split_fieldnames)
        write_jsonl(args.splits_dir / f"{split}.jsonl", split_payloads[split])

    split_summary = {
        "source_samples": str(args.samples),
        "split_policy": "connected_components_over_input_and_target_connectivity_keys",
        "ratios": {
            "train": args.train_ratio,
            "val": args.val_ratio,
            "test": args.test_ratio,
        },
        **split_stats(split_payloads, components, assignments),
    }
    args.splits_dir.mkdir(parents=True, exist_ok=True)
    (args.splits_dir / "split_stats.json").write_text(
        json.dumps(split_summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    triplets, triplet_summary = build_triplets(split_payloads, args)
    triplet_summary.update(
        {
            "source_samples": str(args.samples),
            "source_split_stats": str(args.splits_dir / "split_stats.json"),
            "triplet_rule": (
                "For each positive row, choose the nearest negative target by target-target RDKit "
                "Tanimoto, first within the same input+endpoint+condition+MMP core; optionally "
                "fallback to same endpoint+condition+MMP core."
            ),
        }
    )
    args.triplets_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.triplets_dir / "strict_hard_negative_triplets.csv", triplets, TRIPLET_FIELDNAMES)
    write_jsonl(args.triplets_dir / "strict_hard_negative_triplets.jsonl", triplets)
    (args.triplets_dir / "strict_hard_negative_triplet_stats.json").write_text(
        json.dumps(triplet_summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "split_rows": {split: len(split_payloads[split]) for split in ["train", "val", "test"]},
                "triplets": len(triplets),
                "leakage_check": split_summary["leakage_check"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

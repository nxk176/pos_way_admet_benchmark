from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

try:
    from rdkit import Chem
    from rdkit.Chem import BRICS, Crippen, Descriptors, Lipinski, QED, rdMolDescriptors
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit.Contrib.SA_Score import sascorer
except Exception as exc:  # noqa: BLE001 - RDKit is a hard requirement here.
    raise SystemExit(
        "ERROR: build_fragment_multiproperty_samples_fast.py requires RDKit. "
        "Run it with the project Python environment, for example .\\myenv311\\Scripts\\python.exe"
    ) from exc


TASKS: dict[str, dict[str, Any]] = {
    "logS_mol_L": {"direction": "increase", "min_abs_delta": 0.5},
    "Caco2_logPapp_cm_s": {"direction": "increase", "min_abs_delta": 0.35},
    "hERG_pIC50": {"direction": "decrease", "min_abs_delta": 0.5},
    "microsomal_clearance_mL_min_kg": {
        "direction": "decrease",
        "min_abs_delta": 5.0,
        "min_relative_delta": 0.25,
    },
    "half_life_min": {"direction": "increase", "min_abs_delta": 30.0, "min_relative_delta": 0.25},
    "DRD2_pChEMBL": {"direction": "decrease", "min_abs_delta": 0.5},
    "GSK3B_pChEMBL": {"direction": "decrease", "min_abs_delta": 0.5},
    "JNK3_pChEMBL": {"direction": "decrease", "min_abs_delta": 0.5},
}

DESCRIPTOR_CONSTRAINTS = {
    "max_abs_delta_mw": 120.0,
    "max_abs_delta_logp": 2.0,
    "max_sa_increase": 1.0,
    "max_qed_drop": 0.15,
}


@dataclass
class Aggregate:
    endpoint_name: str
    condition_bucket: str
    molecule_chembl_id: str
    smiles_canon: str
    inchikey: str
    connectivity_key: str
    unit_canonical: str
    values: list[float] = field(default_factory=list)
    observation_ids: list[str] = field(default_factory=list)
    assay_ids: set[str] = field(default_factory=set)
    value: float = 0.0
    descriptors: dict[str, float] = field(default_factory=dict)
    murcko_scaffold: str = ""
    fragments: set[str] = field(default_factory=set)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a fast fragment-based multi-property sample table from normalized observations. "
            "This is intended as a bounded decision artifact, not the exhaustive full benchmark builder."
        )
    )
    parser.add_argument(
        "--observations",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/normalized_csv/property_observations.csv"),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/multiproperty"),
    )
    parser.add_argument("--seed", type=int, default=13)
    parser.add_argument("--max-aggregates-per-endpoint", type=int, default=3000)
    parser.add_argument("--max-positive-per-endpoint", type=int, default=80)
    parser.add_argument("--max-negative-per-endpoint", type=int, default=80)
    parser.add_argument("--max-fragments-per-molecule", type=int, default=4)
    parser.add_argument("--max-bucket-size", type=int, default=90)
    parser.add_argument("--max-targets-per-source", type=int, default=40)
    parser.add_argument("--min-shared-properties", type=int, default=2)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def clean_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    if not math.isfinite(parsed):
        return None
    return parsed


def median(values: list[float]) -> float:
    return float(statistics.median(values))


def clean_cell(value: Any) -> Any:
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, ensure_ascii=False, separators=(",", ":"))
    if isinstance(value, set):
        return json.dumps(sorted(value), ensure_ascii=False, separators=(",", ":"))
    return value


def write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False, separators=(",", ":")) + "\n")


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: clean_cell(row.get(key, "")) for key in fieldnames})


def descriptors_from_smiles(smiles: str) -> tuple[dict[str, float], str, set[str]]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}, "", set()

    descriptors = {
        "mw": float(Descriptors.MolWt(mol)),
        "logp": float(Crippen.MolLogP(mol)),
        "tpsa": float(rdMolDescriptors.CalcTPSA(mol)),
        "hba": float(Lipinski.NumHAcceptors(mol)),
        "hbd": float(Lipinski.NumHDonors(mol)),
        "rotatable_bonds": float(Lipinski.NumRotatableBonds(mol)),
        "qed": float(QED.qed(mol)),
        "sa_score": float(sascorer.calculateScore(mol)),
    }

    scaffold = ""
    try:
        scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
        scaffold = Chem.MolToSmiles(scaffold_mol, isomericSmiles=False) if scaffold_mol else ""
    except Exception:  # noqa: BLE001 - scaffold extraction should not kill the batch.
        scaffold = ""

    fragments: set[str] = set()
    try:
        fragments = {frag for frag in BRICS.BRICSDecompose(mol) if fragment_heavy_atoms(frag) >= 4}
    except Exception:  # noqa: BLE001 - BRICS may fail on uncommon structures.
        fragments = set()
    if not fragments and scaffold:
        fragments = {scaffold}
    return descriptors, scaffold, fragments


def fragment_heavy_atoms(fragment_smiles: str) -> int:
    mol = Chem.MolFromSmiles(fragment_smiles)
    if mol is None:
        return 0
    return int(mol.GetNumHeavyAtoms())


def top_fragments(fragments: set[str], limit: int) -> list[str]:
    return sorted(fragments, key=lambda frag: (fragment_heavy_atoms(frag), len(frag), frag), reverse=True)[:limit]


def aggregate_observations(
    rows: list[dict[str, str]],
    min_shared_properties: int,
) -> tuple[list[Aggregate], dict[str, dict[str, float]], dict[str, dict[str, list[str]]], Counter[str]]:
    aggregates: dict[tuple[str, str, str], Aggregate] = {}
    profile_values: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    profile_observations: dict[str, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    raw_counts: Counter[str] = Counter()

    for row in rows:
        endpoint = row.get("endpoint_name", "")
        if endpoint not in TASKS:
            continue
        value = clean_float(row.get("value_canonical"))
        connectivity_key = row.get("connectivity_key", "")
        smiles = row.get("smiles_canon", "")
        if value is None or not connectivity_key or not smiles:
            continue

        raw_counts[endpoint] += 1
        profile_values[connectivity_key][endpoint].append(value)
        obs_id = row.get("observation_uid", "")
        if obs_id:
            profile_observations[connectivity_key][endpoint].append(obs_id)

        condition_bucket = row.get("condition_bucket") or endpoint
        key = (endpoint, condition_bucket, connectivity_key)
        if key not in aggregates:
            aggregates[key] = Aggregate(
                endpoint_name=endpoint,
                condition_bucket=condition_bucket,
                molecule_chembl_id=row.get("molecule_chembl_id", ""),
                smiles_canon=smiles,
                inchikey=row.get("inchikey", ""),
                connectivity_key=connectivity_key,
                unit_canonical=row.get("unit_canonical", ""),
            )
        aggregate = aggregates[key]
        aggregate.values.append(value)
        if obs_id:
            aggregate.observation_ids.append(obs_id)
        assay_id = row.get("assay_chembl_id", "")
        if assay_id:
            aggregate.assay_ids.add(assay_id)

    profiles = {
        connectivity_key: {endpoint: median(values) for endpoint, values in endpoints.items()}
        for connectivity_key, endpoints in profile_values.items()
        if len(endpoints) >= min_shared_properties
    }
    multi_property_keys = set(profiles)

    kept: list[Aggregate] = []
    for aggregate in aggregates.values():
        if aggregate.connectivity_key not in multi_property_keys:
            continue
        aggregate.value = median(aggregate.values)
        aggregate.observation_ids = sorted(set(aggregate.observation_ids))
        kept.append(aggregate)
    return kept, profiles, profile_observations, raw_counts


def prepare_structures(aggregates: list[Aggregate]) -> None:
    cache: dict[str, tuple[dict[str, float], str, set[str]]] = {}
    for aggregate in aggregates:
        if aggregate.smiles_canon not in cache:
            cache[aggregate.smiles_canon] = descriptors_from_smiles(aggregate.smiles_canon)
        descriptors, scaffold, fragments = cache[aggregate.smiles_canon]
        aggregate.descriptors = descriptors
        aggregate.murcko_scaffold = scaffold
        aggregate.fragments = fragments


def select_aggregates(
    aggregates: list[Aggregate],
    max_per_endpoint: int,
    rng: random.Random,
) -> list[Aggregate]:
    by_endpoint: dict[str, list[Aggregate]] = defaultdict(list)
    for aggregate in aggregates:
        by_endpoint[aggregate.endpoint_name].append(aggregate)

    selected: list[Aggregate] = []
    for endpoint, values in sorted(by_endpoint.items()):
        values = list(values)
        rng.shuffle(values)
        selected.extend(values[:max_per_endpoint])
    return selected


def primary_delta(endpoint: str, before: float, after: float) -> float:
    direction = TASKS[endpoint]["direction"]
    return after - before if direction == "increase" else before - after


def improvement_pass(endpoint: str, before: float, after: float, multiplier: float = 1.0) -> bool:
    spec = TASKS[endpoint]
    delta = primary_delta(endpoint, before, after)
    min_abs = float(spec["min_abs_delta"]) * multiplier
    if delta < min_abs:
        return False
    min_relative = spec.get("min_relative_delta")
    if min_relative is not None and abs(before) > 1e-9:
        if delta / abs(before) < float(min_relative) * multiplier:
            return False
    return True


def significant_worse(endpoint: str, before: float, after: float) -> bool:
    spec = TASKS[endpoint]
    min_abs = float(spec["min_abs_delta"]) * 0.5
    unfavorable_delta = -primary_delta(endpoint, before, after)
    return unfavorable_delta >= min_abs


def descriptor_constraints(before: dict[str, float], after: dict[str, float]) -> dict[str, Any]:
    if not before or not after:
        return {"pass": False}

    abs_delta_mw = abs(after["mw"] - before["mw"])
    abs_delta_logp = abs(after["logp"] - before["logp"])
    delta_qed = after["qed"] - before["qed"]
    delta_sa = after["sa_score"] - before["sa_score"]
    passed = (
        abs_delta_mw <= DESCRIPTOR_CONSTRAINTS["max_abs_delta_mw"]
        and abs_delta_logp <= DESCRIPTOR_CONSTRAINTS["max_abs_delta_logp"]
        and delta_sa <= DESCRIPTOR_CONSTRAINTS["max_sa_increase"]
        and delta_qed >= -DESCRIPTOR_CONSTRAINTS["max_qed_drop"]
    )
    return {
        "pass": passed,
        "abs_delta_mw": round(abs_delta_mw, 6),
        "abs_delta_logp": round(abs_delta_logp, 6),
        "delta_qed": round(delta_qed, 6),
        "delta_sa_score": round(delta_sa, 6),
    }


def secondary_outcomes(
    source: Aggregate,
    target: Aggregate,
    profiles: dict[str, dict[str, float]],
    profile_observations: dict[str, dict[str, list[str]]],
) -> list[dict[str, Any]]:
    source_profile = profiles[source.connectivity_key]
    target_profile = profiles[target.connectivity_key]
    shared = sorted((set(source_profile) & set(target_profile)) - {source.endpoint_name})

    outcomes: list[dict[str, Any]] = []
    for endpoint in shared:
        before = source_profile[endpoint]
        after = target_profile[endpoint]
        outcomes.append(
            {
                "endpoint": endpoint,
                "source_value": round(before, 6),
                "target_value": round(after, 6),
                "direction": TASKS[endpoint]["direction"],
                "favorable_delta": round(primary_delta(endpoint, before, after), 6),
                "improved": improvement_pass(endpoint, before, after, multiplier=0.5),
                "materially_worse": significant_worse(endpoint, before, after),
                "source_observation_ids": profile_observations[source.connectivity_key].get(endpoint, [])[:5],
                "target_observation_ids": profile_observations[target.connectivity_key].get(endpoint, [])[:5],
            }
        )
    return outcomes


def build_instruction(endpoint: str, secondary: list[dict[str, Any]]) -> str:
    direction = TASKS[endpoint]["direction"]
    action = "increase" if direction == "increase" else "decrease"
    secondary_names = ", ".join(item["endpoint"] for item in secondary[:3])
    if secondary_names:
        return (
            f"Given the input molecule, propose a structurally related candidate that shares a BRICS fragment, "
            f"{action}s {endpoint}, preserves or improves measured secondary properties ({secondary_names}), "
            "and keeps MW, logP, QED, and synthetic accessibility within the benchmark constraints."
        )
    return (
        f"Given the input molecule, propose a structurally related candidate that shares a BRICS fragment, "
        f"{action}s {endpoint}, and keeps MW, logP, QED, and synthetic accessibility within the benchmark constraints."
    )


def sample_row(
    sample_id: str,
    sample_type: str,
    fragment_key: str,
    source: Aggregate,
    target: Aggregate,
    profiles: dict[str, dict[str, float]],
    profile_observations: dict[str, dict[str, list[str]]],
) -> dict[str, Any] | None:
    secondary = secondary_outcomes(source, target, profiles, profile_observations)
    if not secondary:
        return None

    descriptor_result = descriptor_constraints(source.descriptors, target.descriptors)
    primary_success = improvement_pass(source.endpoint_name, source.value, target.value)
    secondary_pass = not any(item["materially_worse"] for item in secondary)
    secondary_improved_count = sum(1 for item in secondary if item["improved"])
    is_positive = primary_success and secondary_pass and bool(descriptor_result["pass"])
    if sample_type == "positive" and not is_positive:
        return None
    if sample_type == "negative" and is_positive:
        return None

    direction = TASKS[source.endpoint_name]["direction"]
    row = {
        "sample_id": sample_id,
        "sample_type": sample_type,
        "primary_endpoint": source.endpoint_name,
        "primary_direction": direction,
        "condition_bucket": source.condition_bucket,
        "fragment_key": fragment_key,
        "same_murcko_scaffold": bool(source.murcko_scaffold and source.murcko_scaffold == target.murcko_scaffold),
        "instruction": build_instruction(source.endpoint_name, secondary),
        "input_chembl_id": source.molecule_chembl_id,
        "target_chembl_id": target.molecule_chembl_id,
        "input_connectivity_key": source.connectivity_key,
        "target_connectivity_key": target.connectivity_key,
        "input_smiles_canon": source.smiles_canon,
        "target_smiles_canon": target.smiles_canon,
        "value_before": round(source.value, 6),
        "value_after": round(target.value, 6),
        "favorable_delta": round(primary_delta(source.endpoint_name, source.value, target.value), 6),
        "primary_success": primary_success,
        "secondary_pass": secondary_pass,
        "secondary_improved_count": secondary_improved_count,
        "secondary_outcomes": secondary,
        "descriptor_constraints_pass": bool(descriptor_result["pass"]),
        "abs_delta_mw": descriptor_result.get("abs_delta_mw", ""),
        "abs_delta_logp": descriptor_result.get("abs_delta_logp", ""),
        "delta_qed": descriptor_result.get("delta_qed", ""),
        "delta_sa_score": descriptor_result.get("delta_sa_score", ""),
        "input_primary_observation_ids": source.observation_ids[:5],
        "target_primary_observation_ids": target.observation_ids[:5],
        "input_assay_ids": sorted(source.assay_ids)[:5],
        "target_assay_ids": sorted(target.assay_ids)[:5],
        "multi_property_evidence_count": 1 + len(secondary),
    }
    return row


def build_samples(
    aggregates: list[Aggregate],
    profiles: dict[str, dict[str, float]],
    profile_observations: dict[str, dict[str, list[str]]],
    args: argparse.Namespace,
    rng: random.Random,
) -> list[dict[str, Any]]:
    buckets: dict[tuple[str, str, str], list[Aggregate]] = defaultdict(list)
    for aggregate in aggregates:
        for fragment in top_fragments(aggregate.fragments, args.max_fragments_per_molecule):
            buckets[(aggregate.endpoint_name, aggregate.condition_bucket, fragment)].append(aggregate)

    endpoint_counts: dict[str, Counter[str]] = defaultdict(Counter)
    samples: list[dict[str, Any]] = []
    seen_pairs: set[tuple[str, str, str, str, str]] = set()

    for (endpoint, _condition_bucket, fragment_key), members in sorted(buckets.items()):
        if len(members) < 2:
            continue
        members = list(members)
        rng.shuffle(members)
        members = members[: args.max_bucket_size]

        for source in members:
            target_candidates = [member for member in members if member.connectivity_key != source.connectivity_key]
            rng.shuffle(target_candidates)
            for target in target_candidates[: args.max_targets_per_source]:
                pair_key = (
                    endpoint,
                    source.condition_bucket,
                    fragment_key,
                    source.connectivity_key,
                    target.connectivity_key,
                )
                if pair_key in seen_pairs:
                    continue
                seen_pairs.add(pair_key)

                tentative = sample_row(
                    sample_id="",
                    sample_type="positive",
                    fragment_key=fragment_key,
                    source=source,
                    target=target,
                    profiles=profiles,
                    profile_observations=profile_observations,
                )
                if tentative is not None:
                    if endpoint_counts[endpoint]["positive"] < args.max_positive_per_endpoint:
                        endpoint_counts[endpoint]["positive"] += 1
                        tentative["sample_id"] = (
                            f"fragmp_{endpoint}_pos_{endpoint_counts[endpoint]['positive']:05d}"
                        )
                        samples.append(tentative)
                    continue

                tentative = sample_row(
                    sample_id="",
                    sample_type="negative",
                    fragment_key=fragment_key,
                    source=source,
                    target=target,
                    profiles=profiles,
                    profile_observations=profile_observations,
                )
                if tentative is not None and endpoint_counts[endpoint]["negative"] < args.max_negative_per_endpoint:
                    endpoint_counts[endpoint]["negative"] += 1
                    tentative["sample_id"] = f"fragmp_{endpoint}_neg_{endpoint_counts[endpoint]['negative']:05d}"
                    samples.append(tentative)

                if (
                    endpoint_counts[endpoint]["positive"] >= args.max_positive_per_endpoint
                    and endpoint_counts[endpoint]["negative"] >= args.max_negative_per_endpoint
                ):
                    break
            if (
                endpoint_counts[endpoint]["positive"] >= args.max_positive_per_endpoint
                and endpoint_counts[endpoint]["negative"] >= args.max_negative_per_endpoint
            ):
                break

    samples.sort(key=lambda row: (row["primary_endpoint"], row["sample_type"], row["sample_id"]))
    return samples


def write_stats(
    path: Path,
    samples: list[dict[str, Any]],
    raw_counts: Counter[str],
    kept_aggregates: int,
    selected_aggregates: int,
    profiles: dict[str, dict[str, float]],
) -> None:
    by_endpoint: dict[str, Counter[str]] = defaultdict(Counter)
    secondary_counter: Counter[str] = Counter()
    same_scaffold = 0
    for row in samples:
        by_endpoint[row["primary_endpoint"]][row["sample_type"]] += 1
        same_scaffold += int(bool(row["same_murcko_scaffold"]))
        for outcome in row["secondary_outcomes"]:
            secondary_counter[outcome["endpoint"]] += 1

    stats = {
        "description": (
            "Fast fragment-based multi-property decision artifact. Positive and negative rows share a BRICS "
            "fragment and include at least one measured secondary property in addition to the primary endpoint."
        ),
        "raw_observation_counts": dict(sorted(raw_counts.items())),
        "molecules_with_at_least_two_task_properties": len(profiles),
        "kept_primary_aggregates": kept_aggregates,
        "selected_primary_aggregates": selected_aggregates,
        "sample_count": len(samples),
        "positive_count": sum(1 for row in samples if row["sample_type"] == "positive"),
        "negative_count": sum(1 for row in samples if row["sample_type"] == "negative"),
        "same_murcko_scaffold_count": same_scaffold,
        "fragment_only_count": len(samples) - same_scaffold,
        "by_primary_endpoint": {endpoint: dict(counter) for endpoint, counter in sorted(by_endpoint.items())},
        "secondary_endpoint_evidence_counts": dict(sorted(secondary_counter.items())),
        "descriptor_constraints": DESCRIPTOR_CONSTRAINTS,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(stats, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    rng = random.Random(args.seed)

    rows = read_csv(args.observations)
    aggregates, profiles, profile_observations, raw_counts = aggregate_observations(
        rows,
        min_shared_properties=args.min_shared_properties,
    )
    selected = select_aggregates(aggregates, args.max_aggregates_per_endpoint, rng)
    prepare_structures(selected)
    selected = [aggregate for aggregate in selected if aggregate.descriptors and aggregate.fragments]
    samples = build_samples(selected, profiles, profile_observations, args, rng)

    jsonl_path = args.out_dir / "fragment_multiproperty_samples.fast.jsonl"
    csv_path = args.out_dir / "fragment_multiproperty_samples.fast.csv"
    stats_path = args.out_dir / "fragment_multiproperty_stats.fast.json"
    write_jsonl(jsonl_path, samples)

    fieldnames = [
        "sample_id",
        "sample_type",
        "primary_endpoint",
        "primary_direction",
        "condition_bucket",
        "fragment_key",
        "same_murcko_scaffold",
        "instruction",
        "input_chembl_id",
        "target_chembl_id",
        "input_connectivity_key",
        "target_connectivity_key",
        "input_smiles_canon",
        "target_smiles_canon",
        "value_before",
        "value_after",
        "favorable_delta",
        "primary_success",
        "secondary_pass",
        "secondary_improved_count",
        "secondary_outcomes",
        "descriptor_constraints_pass",
        "abs_delta_mw",
        "abs_delta_logp",
        "delta_qed",
        "delta_sa_score",
        "input_primary_observation_ids",
        "target_primary_observation_ids",
        "input_assay_ids",
        "target_assay_ids",
        "multi_property_evidence_count",
    ]
    write_csv(csv_path, samples, fieldnames)
    write_stats(
        stats_path,
        samples,
        raw_counts=raw_counts,
        kept_aggregates=len(aggregates),
        selected_aggregates=len(selected),
        profiles=profiles,
    )

    print(
        json.dumps(
            {
                "jsonl": str(jsonl_path),
                "csv": str(csv_path),
                "stats": str(stats_path),
                "samples": len(samples),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

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
from typing import Any, Iterable

try:
    from rdkit import Chem
    from rdkit import DataStructs
    from rdkit.Chem import BRICS, Crippen, Descriptors, Lipinski, QED, rdMolDescriptors
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit.Contrib.SA_Score import sascorer
except Exception as exc:  # noqa: BLE001 - this builder is RDKit-backed.
    raise SystemExit(
        "ERROR: build_fragment_multiproperty_samples.py requires RDKit. "
        "Run it with the project environment, for example .\\myenv311\\Scripts\\python.exe"
    ) from exc


TASKS: dict[str, dict[str, Any]] = {
    "logS_mol_L": {"direction": "increase", "min_abs_delta": 0.5},
    "Caco2_logPapp_cm_s": {"direction": "increase", "min_abs_delta": 0.35},
    "hERG_pIC50": {"direction": "decrease", "min_abs_delta": 0.5},
    "microsomal_clearance_mL_min_kg": {"direction": "decrease", "min_abs_delta": 5.0, "min_relative_delta": 0.25},
    "half_life_min": {"direction": "increase", "min_abs_delta": 30.0, "min_relative_delta": 0.25},
    "DRD2_pChEMBL": {"direction": "decrease", "min_abs_delta": 0.5},
    "GSK3B_pChEMBL": {"direction": "decrease", "min_abs_delta": 0.5},
    "JNK3_pChEMBL": {"direction": "decrease", "min_abs_delta": 0.5},
}


def is_dynamic_pchembl(endpoint: str, args: argparse.Namespace) -> bool:
    dynamic_suffixes = ("_pChEMBL", "_pKi", "_pIC50", "_pKd", "_pEC50")
    return bool(args.allow_dynamic_pchembl_endpoints and endpoint.endswith(dynamic_suffixes))


def task_for_endpoint(endpoint: str, args: argparse.Namespace) -> dict[str, Any] | None:
    if endpoint in TASKS:
        return TASKS[endpoint]
    if is_dynamic_pchembl(endpoint, args):
        return {
            "direction": args.dynamic_pchembl_direction,
            "min_abs_delta": args.dynamic_pchembl_min_abs_delta,
        }
    return None


@dataclass(frozen=True)
class MMPDecomposition:
    core: str
    variable: str
    core_heavy_atoms: int
    variable_heavy_atoms: int


@dataclass
class Aggregate:
    endpoint_name: str
    condition_bucket: str
    molecule_chembl_id: str
    smiles_canon: str
    inchikey: str
    connectivity_key: str
    condition_buckets: set[str] = field(default_factory=set)
    values: list[float] = field(default_factory=list)
    observation_ids: list[str] = field(default_factory=list)
    assay_ids: set[str] = field(default_factory=set)
    value: float = 0.0
    descriptors: dict[str, float] = field(default_factory=dict)
    murcko_scaffold: str = ""
    fragments: set[str] = field(default_factory=set)
    mmp_decompositions: list[MMPDecomposition] = field(default_factory=list)
    fingerprint: Any = None


def read_csv(path: Path) -> Iterable[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle)


def write_csv(path: Path, rows: Iterable[dict[str, Any]], fieldnames: list[str]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: clean_cell(row.get(key)) for key in fieldnames})
            count += 1
    return count


def write_jsonl(path: Path, rows: Iterable[dict[str, Any]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False, separators=(",", ":")) + "\n")
            count += 1
    return count


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
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def descriptors(mol: Chem.Mol) -> dict[str, float]:
    return {
        "mw": float(Descriptors.MolWt(mol)),
        "logp": float(Crippen.MolLogP(mol)),
        "tpsa": float(rdMolDescriptors.CalcTPSA(mol)),
        "hba": float(Lipinski.NumHAcceptors(mol)),
        "hbd": float(Lipinski.NumHDonors(mol)),
        "rotatable_bonds": float(Lipinski.NumRotatableBonds(mol)),
        "heavy_atoms": float(mol.GetNumHeavyAtoms()),
        "qed": float(QED.qed(mol)),
    }


def fragment_heavy_atoms(fragment: str) -> int:
    mol = Chem.MolFromSmiles(fragment)
    return int(mol.GetNumHeavyAtoms()) if mol is not None else 0


def mmp_decompositions_from_brics(mol: Chem.Mol) -> list[MMPDecomposition]:
    decompositions: dict[tuple[str, str], MMPDecomposition] = {}
    for atom_pair, _labels in BRICS.FindBRICSBonds(mol):
        bond = mol.GetBondBetweenAtoms(int(atom_pair[0]), int(atom_pair[1]))
        if bond is None:
            continue
        fragmented = Chem.FragmentOnBonds(mol, [bond.GetIdx()], addDummies=True)
        try:
            parts = Chem.GetMolFrags(fragmented, asMols=True, sanitizeFrags=True)
        except Exception:  # noqa: BLE001 - uncommon sanitization failures should not kill the build.
            continue
        if len(parts) != 2:
            continue

        scored_parts = []
        for part in parts:
            heavy_atoms = int(part.GetNumHeavyAtoms())
            smiles = Chem.MolToSmiles(part, canonical=True, isomericSmiles=False)
            scored_parts.append((heavy_atoms, len(smiles), smiles))
        scored_parts.sort(reverse=True)
        core_heavy, _core_len, core = scored_parts[0]
        variable_heavy, _variable_len, variable = scored_parts[1]
        if core == variable:
            continue
        decompositions[(core, variable)] = MMPDecomposition(
            core=core,
            variable=variable,
            core_heavy_atoms=core_heavy,
            variable_heavy_atoms=variable_heavy,
        )
    return list(decompositions.values())


def fragment_signature(
    smiles: str,
    cache: dict[str, tuple[str, dict[str, float], set[str], list[MMPDecomposition], Any]],
) -> tuple[str, dict[str, float], set[str], list[MMPDecomposition], Any] | None:
    if smiles in cache:
        return cache[smiles]
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    scaffold = ""
    brics_bond_ids: list[int] = []
    for atom_pair, _labels in BRICS.FindBRICSBonds(mol):
        bond = mol.GetBondBetweenAtoms(int(atom_pair[0]), int(atom_pair[1]))
        if bond is not None:
            brics_bond_ids.append(bond.GetIdx())

    fragments: set[str] = set()
    if brics_bond_ids:
        fragmented = Chem.FragmentOnBonds(mol, brics_bond_ids, addDummies=True)
        for fragment_mol in Chem.GetMolFrags(fragmented, asMols=True, sanitizeFrags=False):
            if fragment_mol.GetNumHeavyAtoms() < 4:
                continue
            fragment = Chem.MolToSmiles(fragment_mol, canonical=True, isomericSmiles=False)
            if 6 <= len(fragment) <= 180:
                fragments.add(fragment)
    if not fragments:
        fallback = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
        fragments.add(f"fallback:{fallback}")

    result = (scaffold, descriptors(mol), fragments, mmp_decompositions_from_brics(mol), Chem.RDKFingerprint(mol))
    cache[smiles] = result
    return result


def discover_supported_endpoints(path: Path, args: argparse.Namespace) -> set[str]:
    counts: Counter[str] = Counter()
    for row in read_csv(path):
        endpoint = row.get("endpoint_name", "")
        if endpoint in TASKS or is_dynamic_pchembl(endpoint, args):
            if as_float(row.get("value_canonical")) is not None:
                counts[endpoint] += 1

    supported = {endpoint for endpoint in TASKS if counts.get(endpoint, 0) > 0}
    dynamic = [
        (endpoint, count)
        for endpoint, count in counts.items()
        if endpoint not in TASKS and is_dynamic_pchembl(endpoint, args) and count >= args.min_endpoint_observations
    ]
    dynamic.sort(key=lambda item: (-item[1], item[0]))
    if args.max_dynamic_endpoints:
        dynamic = dynamic[: args.max_dynamic_endpoints]
    supported.update(endpoint for endpoint, _count in dynamic)
    return supported


def load_aggregates(
    path: Path,
    min_shared_secondary_endpoints: int,
    supported_endpoints: set[str],
    strict_condition_bucket: bool,
) -> tuple[list[Aggregate], dict[str, dict[str, float]]]:
    grouped: dict[tuple[str, str, str], Aggregate] = {}
    endpoint_values_by_molecule: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))

    for row in read_csv(path):
        endpoint = row.get("endpoint_name", "")
        if endpoint not in supported_endpoints:
            continue
        value = as_float(row.get("value_canonical"))
        if value is None:
            continue
        smiles = row.get("smiles_canon", "")
        if not smiles:
            continue
        connectivity_key = row.get("connectivity_key") or row.get("inchikey") or row.get("molecule_chembl_id")
        if not connectivity_key:
            continue
        condition_bucket = row.get("condition_bucket", "") if strict_condition_bucket else "endpoint_level_median"
        key = (endpoint, condition_bucket, connectivity_key)
        if key not in grouped:
            grouped[key] = Aggregate(
                endpoint_name=endpoint,
                condition_bucket=condition_bucket or "condition_unknown",
                molecule_chembl_id=row.get("molecule_chembl_id", ""),
                smiles_canon=smiles,
                inchikey=row.get("inchikey", ""),
                connectivity_key=connectivity_key,
            )
        agg = grouped[key]
        agg.values.append(value)
        if row.get("condition_bucket"):
            agg.condition_buckets.add(row["condition_bucket"])
        if row.get("observation_uid"):
            agg.observation_ids.append(row["observation_uid"])
        if row.get("assay_chembl_id"):
            agg.assay_ids.add(row["assay_chembl_id"])
        endpoint_values_by_molecule[connectivity_key][endpoint].append(value)

    endpoint_values = {
        mol_key: {endpoint: float(statistics.median(values)) for endpoint, values in endpoint_map.items()}
        for mol_key, endpoint_map in endpoint_values_by_molecule.items()
    }

    min_endpoint_count = min_shared_secondary_endpoints + 1
    eligible_molecules = {
        mol_key for mol_key, endpoint_map in endpoint_values.items() if len(endpoint_map) >= min_endpoint_count
    }

    aggregates = []
    mol_cache: dict[str, tuple[str, dict[str, float], set[str], list[MMPDecomposition], Any]] = {}
    for agg in grouped.values():
        if agg.connectivity_key not in eligible_molecules:
            continue
        parsed = fragment_signature(agg.smiles_canon, mol_cache)
        if parsed is None:
            continue
        scaffold, desc, fragments, decompositions, fingerprint = parsed
        if not decompositions:
            continue
        agg.value = float(statistics.median(agg.values))
        agg.descriptors = desc
        agg.murcko_scaffold = scaffold
        agg.fragments = fragments
        agg.mmp_decompositions = decompositions
        agg.fingerprint = fingerprint
        aggregates.append(agg)
    return aggregates, endpoint_values


def primary_improvement(args: argparse.Namespace, source: Aggregate, target: Aggregate) -> tuple[bool, float, float | None]:
    task = task_for_endpoint(source.endpoint_name, args)
    if task is None:
        return False, 0.0, None
    delta = target.value - source.value
    direction = task["direction"]
    if direction == "increase":
        ok = delta >= float(task["min_abs_delta"])
        rel = delta / abs(source.value) if source.value else None
    else:
        ok = delta <= -float(task["min_abs_delta"])
        rel = (source.value - target.value) / abs(source.value) if source.value else None
    if "min_relative_delta" in task:
        ok = ok and rel is not None and rel >= float(task["min_relative_delta"])
    return ok, delta, rel


def secondary_descriptor_constraints(source: Aggregate, target: Aggregate) -> dict[str, Any]:
    sd = source.descriptors
    td = target.descriptors
    ensure_sa_score(source)
    ensure_sa_score(target)
    checks = {
        "abs_delta_mw_lte_100": abs(td["mw"] - sd["mw"]) <= 100.0,
        "delta_sa_lte_1": (td["sa_score"] - sd["sa_score"]) <= 1.0,
        "qed_drop_lte_0_1": (sd["qed"] - td["qed"]) <= 0.1,
        "abs_delta_logp_lte_2": abs(td["logp"] - sd["logp"]) <= 2.0,
    }
    return {
        "success": all(checks.values()),
        "checks": checks,
        "delta_mw": round(td["mw"] - sd["mw"], 4),
        "delta_logp": round(td["logp"] - sd["logp"], 4),
        "delta_qed": round(td["qed"] - sd["qed"], 4),
        "delta_sa_score": round(td["sa_score"] - sd["sa_score"], 4),
    }


def ensure_murcko_scaffold(aggregate: Aggregate) -> None:
    if aggregate.murcko_scaffold:
        return
    mol = Chem.MolFromSmiles(aggregate.smiles_canon)
    if mol is None:
        return
    scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold_mol is not None and scaffold_mol.GetNumHeavyAtoms() > 0:
        aggregate.murcko_scaffold = Chem.MolToSmiles(scaffold_mol, canonical=True, isomericSmiles=False)


def ensure_sa_score(aggregate: Aggregate) -> None:
    if "sa_score" in aggregate.descriptors:
        return
    mol = Chem.MolFromSmiles(aggregate.smiles_canon)
    aggregate.descriptors["sa_score"] = float(sascorer.calculateScore(mol)) if mol is not None else 99.0


def shared_secondary_endpoints(
    source: Aggregate,
    target: Aggregate,
    endpoint_values_by_molecule: dict[str, dict[str, float]],
) -> dict[str, Any]:
    source_values = endpoint_values_by_molecule.get(source.connectivity_key, {})
    target_values = endpoint_values_by_molecule.get(target.connectivity_key, {})
    shared = sorted((set(source_values) & set(target_values)) - {source.endpoint_name})
    payload = []
    for endpoint in shared:
        payload.append(
            {
                "endpoint_name": endpoint,
                "input_value": round(source_values[endpoint], 6),
                "target_value": round(target_values[endpoint], 6),
                "delta": round(target_values[endpoint] - source_values[endpoint], 6),
            }
        )
    return {"count": len(payload), "items": payload[:5]}


def fragment_similarity(source: Aggregate, target: Aggregate) -> tuple[int, float, list[str]]:
    shared = source.fragments & target.fragments
    union = source.fragments | target.fragments
    sim = len(shared) / len(union) if union else 0.0
    return len(shared), sim, sorted(shared)[:10]


def tanimoto_similarity(source: Aggregate, target: Aggregate) -> float:
    if source.fingerprint is None or target.fingerprint is None:
        return 0.0
    return float(DataStructs.FingerprintSimilarity(source.fingerprint, target.fingerprint))


def candidate_groups(args: argparse.Namespace, aggregates: list[Aggregate]) -> dict[tuple[str, str, str], list[tuple[Aggregate, MMPDecomposition]]]:
    groups: dict[tuple[str, str, str], list[tuple[Aggregate, MMPDecomposition]]] = defaultdict(list)
    for agg in aggregates:
        heavy_atoms = max(float(agg.descriptors.get("heavy_atoms", 0.0)), 1.0)
        for decomposition in agg.mmp_decompositions:
            if decomposition.core_heavy_atoms < args.min_mmp_core_heavy_atoms:
                continue
            if decomposition.variable_heavy_atoms > args.max_mmp_variable_heavy_atoms:
                continue
            if decomposition.core_heavy_atoms / heavy_atoms < args.min_mmp_core_ratio:
                continue
            condition_key = agg.condition_bucket if args.strict_condition_bucket else "brics_mmp_endpoint_pool"
            groups[(agg.endpoint_name, condition_key, decomposition.core)].append((agg, decomposition))
    return groups


def build_samples(args: argparse.Namespace, aggregates: list[Aggregate], endpoint_values: dict[str, dict[str, float]]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rng = random.Random(args.seed)
    groups = candidate_groups(args, aggregates)
    positive_counts: Counter[str] = Counter()
    negative_counts: Counter[str] = Counter()
    stats: Counter[str] = Counter()
    samples: list[dict[str, Any]] = []
    seen_pairs: set[tuple[str, str, str]] = set()
    sample_idx = 1

    for (endpoint, condition_bucket, mmp_core), rows in sorted(groups.items(), key=lambda item: (-len(item[1]), item[0])):
        if positive_counts[endpoint] >= args.max_positive_per_endpoint and negative_counts[endpoint] >= args.max_negative_per_endpoint:
            continue
        if len(rows) < 2:
            continue
        source_indices = list(range(len(rows)))
        rng.shuffle(source_indices)
        source_indices = source_indices[: args.max_sources_per_group]

        for source_idx in source_indices:
            source, source_decomposition = rows[source_idx]
            ordered_candidates = [
                idx
                for idx, (candidate, candidate_decomposition) in enumerate(rows)
                if idx != source_idx
                and candidate.connectivity_key != source.connectivity_key
                and candidate_decomposition.variable != source_decomposition.variable
                and abs(candidate_decomposition.variable_heavy_atoms - source_decomposition.variable_heavy_atoms)
                <= args.max_mmp_variable_heavy_delta
            ]
            rng.shuffle(ordered_candidates)
            ordered_candidates = ordered_candidates[: args.max_candidates_per_source]

            for target_idx in ordered_candidates:
                target, target_decomposition = rows[target_idx]
                pair_key = (endpoint, source.connectivity_key, target.connectivity_key)
                if pair_key in seen_pairs:
                    continue
                seen_pairs.add(pair_key)
                shared_count, frag_sim, shared_fragments = fragment_similarity(source, target)
                if shared_count < args.min_shared_fragments:
                    stats[f"{endpoint}:skipped_no_shared_fragment"] += 1
                    continue
                tanimoto = tanimoto_similarity(source, target)
                if tanimoto < args.min_tanimoto_similarity:
                    stats[f"{endpoint}:skipped_low_tanimoto"] += 1
                    continue

                primary_ok, delta, rel = primary_improvement(args, source, target)
                secondary = secondary_descriptor_constraints(source, target)
                shared_secondary = shared_secondary_endpoints(source, target, endpoint_values)
                if shared_secondary["count"] < args.min_shared_secondary_endpoints:
                    stats[f"{endpoint}:skipped_no_shared_secondary"] += 1
                    continue
                is_positive = primary_ok and secondary["success"]

                if is_positive:
                    if positive_counts[endpoint] >= args.max_positive_per_endpoint:
                        continue
                    sample_type = "positive"
                    positive_counts[endpoint] += 1
                else:
                    if negative_counts[endpoint] >= args.max_negative_per_endpoint:
                        continue
                    sample_type = "negative"
                    negative_counts[endpoint] += 1

                ensure_murcko_scaffold(source)
                ensure_murcko_scaffold(target)
                sample = {
                    "sample_id": f"fragmp_{sample_idx:08d}",
                    "sample_type": sample_type,
                    "primary_endpoint": endpoint,
                    "condition_bucket": condition_bucket,
                    "instruction": instruction(args, endpoint),
                    "input_chembl_id": source.molecule_chembl_id,
                    "input_smiles_canon": source.smiles_canon,
                    "input_inchikey": source.inchikey,
                    "input_connectivity_key": source.connectivity_key,
                    "target_chembl_id": target.molecule_chembl_id,
                    "target_smiles_canon": target.smiles_canon,
                    "target_inchikey": target.inchikey,
                    "target_connectivity_key": target.connectivity_key,
                    "value_before": round(source.value, 6),
                    "value_after": round(target.value, 6),
                    "delta_value": round(delta, 6),
                    "relative_delta": round(rel, 6) if rel is not None else "",
                    "primary_success": primary_ok,
                    "secondary_success": secondary["success"],
                    "pairing_rule": "brics_mmp_local_replacement",
                    "mmp_core": mmp_core,
                    "input_variable_fragment": source_decomposition.variable,
                    "target_variable_fragment": target_decomposition.variable,
                    "mmp_core_heavy_atoms": source_decomposition.core_heavy_atoms,
                    "input_variable_heavy_atoms": source_decomposition.variable_heavy_atoms,
                    "target_variable_heavy_atoms": target_decomposition.variable_heavy_atoms,
                    "tanimoto_similarity": round(tanimoto, 6),
                    "fragment_shared_count": shared_count,
                    "fragment_jaccard": round(frag_sim, 6),
                    "shared_fragments": shared_fragments,
                    "same_murcko_scaffold": source.murcko_scaffold == target.murcko_scaffold and bool(source.murcko_scaffold),
                    "input_murcko_scaffold": source.murcko_scaffold,
                    "target_murcko_scaffold": target.murcko_scaffold,
                    "input_observation_ids": source.observation_ids[:10],
                    "target_observation_ids": target.observation_ids[:10],
                    "input_assay_ids": sorted(source.assay_ids)[:10],
                    "target_assay_ids": sorted(target.assay_ids)[:10],
                    "input_condition_buckets": sorted(source.condition_buckets)[:10],
                    "target_condition_buckets": sorted(target.condition_buckets)[:10],
                    "secondary_constraints": secondary,
                    "shared_secondary_experimental_endpoints": shared_secondary,
                }
                samples.append(sample)
                sample_idx += 1

                if positive_counts[endpoint] >= args.max_positive_per_endpoint and negative_counts[endpoint] >= args.max_negative_per_endpoint:
                    break
            if positive_counts[endpoint] >= args.max_positive_per_endpoint and negative_counts[endpoint] >= args.max_negative_per_endpoint:
                break

        stats[f"{endpoint}:groups_seen"] += 1

    stats_payload = {
        "num_samples": len(samples),
        "positive_counts": dict(sorted(positive_counts.items())),
        "negative_counts": dict(sorted(negative_counts.items())),
        "generation_counters": dict(sorted(stats.items())),
        "parameters": {
            "pairing_rule": "brics_mmp_local_replacement",
            "strict_condition_bucket": args.strict_condition_bucket,
            "allow_dynamic_pchembl_endpoints": args.allow_dynamic_pchembl_endpoints,
            "min_endpoint_observations": args.min_endpoint_observations,
            "max_dynamic_endpoints": args.max_dynamic_endpoints,
            "dynamic_pchembl_direction": args.dynamic_pchembl_direction,
            "dynamic_pchembl_min_abs_delta": args.dynamic_pchembl_min_abs_delta,
            "min_mmp_core_heavy_atoms": args.min_mmp_core_heavy_atoms,
            "min_mmp_core_ratio": args.min_mmp_core_ratio,
            "max_mmp_variable_heavy_atoms": args.max_mmp_variable_heavy_atoms,
            "max_mmp_variable_heavy_delta": args.max_mmp_variable_heavy_delta,
            "min_tanimoto_similarity": args.min_tanimoto_similarity,
            "min_shared_fragments": args.min_shared_fragments,
            "max_positive_per_endpoint": args.max_positive_per_endpoint,
            "max_negative_per_endpoint": args.max_negative_per_endpoint,
            "max_sources_per_group": args.max_sources_per_group,
            "max_candidates_per_source": args.max_candidates_per_source,
            "min_shared_secondary_endpoints": args.min_shared_secondary_endpoints,
            "seed": args.seed,
        },
        "notes": [
            "Samples are grouped by endpoint-level BRICS-derived matched molecular pair cores, not exact Murcko scaffold.",
            "Each pair represents a local replacement: same MMP core, different variable fragment.",
            (
                "Primary endpoint values are medians within the same condition bucket."
                if args.strict_condition_bucket
                else "Primary endpoint values are endpoint-level medians across available condition buckets; source condition buckets remain in provenance fields."
            ),
            "Positive samples require primary experimental improvement plus RDKit descriptor constraints.",
            "Every sample requires at least the configured number of shared secondary experimental endpoints.",
            "Negative samples share the MMP core but fail primary improvement or secondary property constraints.",
        ],
    }
    print(
        "[build] final_counts="
        + json.dumps(
            {
                "positive": dict(sorted(positive_counts.items())),
                "negative": dict(sorted(negative_counts.items())),
                "total_samples": len(samples),
            },
            ensure_ascii=False,
            sort_keys=True,
        ),
        flush=True,
    )
    return samples, stats_payload


def instruction(args: argparse.Namespace, endpoint: str) -> str:
    task = task_for_endpoint(endpoint, args) or {"direction": "decrease"}
    direction = "increase" if task["direction"] == "increase" else "decrease"
    return (
        f"{direction.capitalize()} {endpoint} while preserving relevant fragments and satisfying "
        "secondary MW/LogP/QED/SA constraints."
    )


def flat_sample(row: dict[str, Any]) -> dict[str, Any]:
    flattened = dict(row)
    flattened["shared_fragments"] = row.get("shared_fragments")
    flattened["input_observation_ids"] = row.get("input_observation_ids")
    flattened["target_observation_ids"] = row.get("target_observation_ids")
    flattened["input_assay_ids"] = row.get("input_assay_ids")
    flattened["target_assay_ids"] = row.get("target_assay_ids")
    flattened["input_condition_buckets"] = row.get("input_condition_buckets")
    flattened["target_condition_buckets"] = row.get("target_condition_buckets")
    flattened["secondary_constraints"] = row.get("secondary_constraints")
    flattened["shared_secondary_experimental_endpoints"] = row.get("shared_secondary_experimental_endpoints")
    return flattened


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build BRICS/MMP-style fragment multi-property positive/negative samples."
    )
    parser.add_argument("--observations", type=Path, default=Path("data/normalized_csv/property_observations.csv"))
    parser.add_argument("--out-jsonl", type=Path, default=Path("data/multiproperty/fragment_multiproperty_samples.jsonl"))
    parser.add_argument("--out-csv", type=Path, default=Path("data/multiproperty/fragment_multiproperty_samples.csv"))
    parser.add_argument("--stats", type=Path, default=Path("data/multiproperty/fragment_multiproperty_stats.json"))
    parser.add_argument("--max-positive-per-endpoint", type=int, default=1200)
    parser.add_argument("--max-negative-per-endpoint", type=int, default=1200)
    parser.add_argument("--max-sources-per-group", type=int, default=1200)
    parser.add_argument("--max-candidates-per-source", type=int, default=200)
    parser.add_argument("--strict-condition-bucket", action="store_true")
    parser.add_argument("--allow-dynamic-pchembl-endpoints", action="store_true")
    parser.add_argument("--min-endpoint-observations", type=int, default=1000)
    parser.add_argument("--max-dynamic-endpoints", type=int, default=0, help="0 means no explicit dynamic endpoint cap.")
    parser.add_argument("--dynamic-pchembl-direction", choices=["increase", "decrease"], default="decrease")
    parser.add_argument("--dynamic-pchembl-min-abs-delta", type=float, default=0.5)
    parser.add_argument("--min-mmp-core-heavy-atoms", type=int, default=8)
    parser.add_argument("--min-mmp-core-ratio", type=float, default=0.35)
    parser.add_argument("--max-mmp-variable-heavy-atoms", type=int, default=24)
    parser.add_argument("--max-mmp-variable-heavy-delta", type=int, default=12)
    parser.add_argument("--min-tanimoto-similarity", type=float, default=0.45)
    parser.add_argument("--min-shared-fragments", type=int, default=1)
    parser.add_argument("--min-shared-secondary-endpoints", type=int, default=1)
    parser.add_argument("--seed", type=int, default=29)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.observations.is_file():
        raise SystemExit(f"ERROR: normalized observation table not found: {args.observations}")
    print(f"[load] observations={args.observations}", flush=True)
    supported_endpoints = discover_supported_endpoints(args.observations, args)
    print(f"[load] supported_endpoints={len(supported_endpoints)}", flush=True)
    aggregates, endpoint_values = load_aggregates(
        args.observations,
        args.min_shared_secondary_endpoints,
        supported_endpoints,
        args.strict_condition_bucket,
    )
    print(
        f"[load] aggregates={len(aggregates)} molecules_with_endpoint_values={len(endpoint_values)}",
        flush=True,
    )
    samples, stats = build_samples(args, aggregates, endpoint_values)
    if not samples:
        raise SystemExit("ERROR: no fragment multi-property samples were generated.")

    fieldnames = [
        "sample_id",
        "sample_type",
        "primary_endpoint",
        "condition_bucket",
        "instruction",
        "input_chembl_id",
        "input_smiles_canon",
        "input_inchikey",
        "input_connectivity_key",
        "target_chembl_id",
        "target_smiles_canon",
        "target_inchikey",
        "target_connectivity_key",
        "value_before",
        "value_after",
        "delta_value",
        "relative_delta",
        "primary_success",
        "secondary_success",
        "pairing_rule",
        "mmp_core",
        "input_variable_fragment",
        "target_variable_fragment",
        "mmp_core_heavy_atoms",
        "input_variable_heavy_atoms",
        "target_variable_heavy_atoms",
        "tanimoto_similarity",
        "fragment_shared_count",
        "fragment_jaccard",
        "same_murcko_scaffold",
        "input_murcko_scaffold",
        "target_murcko_scaffold",
        "shared_fragments",
        "input_observation_ids",
        "target_observation_ids",
        "input_assay_ids",
        "target_assay_ids",
        "input_condition_buckets",
        "target_condition_buckets",
        "secondary_constraints",
        "shared_secondary_experimental_endpoints",
    ]

    write_jsonl(args.out_jsonl, samples)
    write_csv(args.out_csv, (flat_sample(row) for row in samples), fieldnames)
    args.stats.parent.mkdir(parents=True, exist_ok=True)
    args.stats.write_text(json.dumps(stats, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(stats, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

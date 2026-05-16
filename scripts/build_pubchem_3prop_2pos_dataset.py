from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import re
import statistics
from collections import Counter, defaultdict, deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import BRICS
except Exception as exc:  # noqa: BLE001
    raise SystemExit(
        "ERROR: build_pubchem_3prop_2pos_dataset.py requires RDKit. "
        "Run it with .\\myenv311\\Scripts\\python.exe"
    ) from exc


FIELDS = [
    "query_id",
    "split",
    "instruction",
    "input_smiles_canon",
    "input_connectivity_key",
    "input_pubchem_cid",
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

PRIMARY_ENDPOINT_RE = re.compile(r"_p(?:IC50|KI|KD|EC50|AC50|CC50|GI50|LC50)$")


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
    pubchem_cid: str
    smiles_canon: str
    inchikey: str
    connectivity_key: str
    values: list[float] = field(default_factory=list)
    value: float = 0.0
    observation_ids: list[str] = field(default_factory=list)
    descriptors: dict[str, float] = field(default_factory=dict)
    mmp_decompositions: list[MMPDecomposition] = field(default_factory=list)
    fingerprint: Any = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a PubChem-only source-specific dataset: one input + 3-property instruction -> 2 positive answers."
    )
    parser.add_argument(
        "--observations",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/pubchem_bioassay/pubchem_bioassay_observations.csv"),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/pubchem_3prop_2pos"),
    )
    parser.add_argument("--min-primary-molecules", type=int, default=5)
    parser.add_argument("--max-primary-endpoints", type=int, default=5000, help="Keep largest primary endpoints before MMP generation; 0 means no cap.")
    parser.add_argument("--max-aggregates-per-endpoint", type=int, default=600, help="Deterministic cap per endpoint before MMP generation; 0 means no cap.")
    parser.add_argument("--min-primary-delta", type=float, default=0.5)
    parser.add_argument("--min-tanimoto", type=float, default=0.35)
    parser.add_argument("--min-mmp-core-heavy-atoms", type=int, default=6)
    parser.add_argument("--min-mmp-core-ratio", type=float, default=0.25)
    parser.add_argument("--max-mmp-variable-heavy-atoms", type=int, default=30)
    parser.add_argument("--max-mmp-variable-heavy-delta", type=int, default=16)
    parser.add_argument("--max-sources-per-group", type=int, default=120)
    parser.add_argument("--max-candidates-per-source", type=int, default=120)
    parser.add_argument("--max-queries", type=int, default=0, help="0 means no cap.")
    parser.add_argument("--seed", type=int, default=29)
    return parser.parse_args()


def read_csv(path: Path) -> Iterable[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle)


def write_csv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: clean_cell(row.get(key)) for key in fields})
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


def compact_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def as_float(value: Any) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def pubchem_cid_from_uid(uid: str) -> str:
    match = re.search(r"_cid([^_]+)", uid or "")
    return match.group(1) if match else ""


def descriptor_payload(row: dict[str, str]) -> dict[str, float] | None:
    keys = ["mw", "logp", "qed", "sa_score", "heavy_atoms"]
    payload: dict[str, float] = {}
    for key in keys:
        value = as_float(row.get(key))
        if value is None:
            return None
        payload[key] = value
    return payload


def mmp_decompositions_from_smiles(smiles: str) -> tuple[list[MMPDecomposition], Any] | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    decompositions: dict[tuple[str, str], MMPDecomposition] = {}
    for atom_pair, _labels in BRICS.FindBRICSBonds(mol):
        bond = mol.GetBondBetweenAtoms(int(atom_pair[0]), int(atom_pair[1]))
        if bond is None:
            continue
        fragmented = Chem.FragmentOnBonds(mol, [bond.GetIdx()], addDummies=True)
        try:
            parts = Chem.GetMolFrags(fragmented, asMols=True, sanitizeFrags=True)
        except Exception:  # noqa: BLE001
            continue
        if len(parts) != 2:
            continue
        scored = []
        for part in parts:
            smiles_part = Chem.MolToSmiles(part, canonical=True, isomericSmiles=False)
            scored.append((int(part.GetNumHeavyAtoms()), len(smiles_part), smiles_part))
        scored.sort(reverse=True)
        core_heavy, _core_len, core = scored[0]
        variable_heavy, _var_len, variable = scored[1]
        if core != variable:
            decompositions[(core, variable)] = MMPDecomposition(core, variable, core_heavy, variable_heavy)
    return list(decompositions.values()), Chem.RDKFingerprint(mol)


def supported_primary_endpoint(endpoint: str) -> bool:
    return bool(PRIMARY_ENDPOINT_RE.search(endpoint or ""))


def load_pubchem_aggregates(
    path: Path,
    min_primary_molecules: int,
    max_primary_endpoints: int,
    max_aggregates_per_endpoint: int,
    seed: int,
) -> tuple[list[Aggregate], dict[str, dict[str, float]], dict[str, Any]]:
    primary_counts: Counter[str] = Counter()
    all_endpoint_values_raw: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    rows = []
    for row in read_csv(path):
        value = as_float(row.get("value_canonical"))
        if value is None:
            continue
        if row.get("relation_raw") != "=":
            continue
        key = row.get("connectivity_key", "")
        endpoint = row.get("endpoint_name", "")
        if not key or not endpoint:
            continue
        all_endpoint_values_raw[key][endpoint].append(value)
        if supported_primary_endpoint(endpoint):
            primary_counts[endpoint] += 1
        rows.append(row)

    supported_ranked = [
        endpoint for endpoint, count in primary_counts.most_common() if count >= min_primary_molecules
    ]
    if max_primary_endpoints:
        supported_ranked = supported_ranked[:max_primary_endpoints]
    supported = set(supported_ranked)
    grouped: dict[tuple[str, str, str], Aggregate] = {}
    for row in rows:
        endpoint = row.get("endpoint_name", "")
        if endpoint not in supported:
            continue
        value = as_float(row.get("value_canonical"))
        desc = descriptor_payload(row)
        if value is None or desc is None:
            continue
        key = row.get("connectivity_key", "")
        group_key = (endpoint, row.get("condition_bucket", ""), key)
        if group_key not in grouped:
            grouped[group_key] = Aggregate(
                endpoint_name=endpoint,
                condition_bucket=row.get("condition_bucket", ""),
                pubchem_cid=pubchem_cid_from_uid(row.get("observation_uid", "")),
                smiles_canon=row.get("smiles_canon", ""),
                inchikey=row.get("inchikey", ""),
                connectivity_key=key,
                descriptors=desc,
            )
        grouped[group_key].values.append(value)
        if row.get("observation_uid"):
            grouped[group_key].observation_ids.append(row["observation_uid"])

    endpoint_values = {
        key: {endpoint: float(statistics.median(values)) for endpoint, values in endpoint_map.items()}
        for key, endpoint_map in all_endpoint_values_raw.items()
    }

    mol_cache: dict[str, tuple[list[MMPDecomposition], Any] | None] = {}
    aggregates = []
    counters: Counter[str] = Counter()
    grouped_values = list(grouped.values())
    if max_aggregates_per_endpoint:
        rng = random.Random(seed)
        by_endpoint: dict[str, list[Aggregate]] = defaultdict(list)
        for agg in grouped_values:
            by_endpoint[agg.endpoint_name].append(agg)
        grouped_values = []
        for endpoint in sorted(by_endpoint):
            members = by_endpoint[endpoint]
            rng.shuffle(members)
            grouped_values.extend(members[:max_aggregates_per_endpoint])

    for agg in grouped_values:
        agg.value = float(statistics.median(agg.values))
        parsed = mol_cache.get(agg.smiles_canon)
        if agg.smiles_canon not in mol_cache:
            parsed = mmp_decompositions_from_smiles(agg.smiles_canon)
            mol_cache[agg.smiles_canon] = parsed
        if parsed is None:
            counters["skipped_invalid_smiles"] += 1
            continue
        decompositions, fp = parsed
        if not decompositions:
            counters["skipped_no_brics_mmp_decomposition"] += 1
            continue
        agg.mmp_decompositions = decompositions
        agg.fingerprint = fp
        aggregates.append(agg)

    return aggregates, endpoint_values, {
        "supported_primary_endpoints": len(supported),
        "max_primary_endpoints": max_primary_endpoints,
        "max_aggregates_per_endpoint": max_aggregates_per_endpoint,
        "raw_exact_rows_seen": len(rows),
        "primary_aggregate_rows": len(grouped),
        "selected_primary_aggregate_rows": len(grouped_values),
        "usable_primary_aggregates": len(aggregates),
        "counters": dict(counters),
    }


def tanimoto(source: Aggregate, target: Aggregate) -> float:
    if source.fingerprint is None or target.fingerprint is None:
        return 0.0
    return float(DataStructs.FingerprintSimilarity(source.fingerprint, target.fingerprint))


def descriptor_constraints(source: Aggregate, target: Aggregate) -> dict[str, Any]:
    sd = source.descriptors
    td = target.descriptors
    checks = {
        "abs_delta_mw_lte_100": abs(td["mw"] - sd["mw"]) <= 100.0,
        "abs_delta_logp_lte_2": abs(td["logp"] - sd["logp"]) <= 2.0,
        "qed_drop_lte_0_1": (sd["qed"] - td["qed"]) <= 0.1,
        "delta_sa_lte_1": (td["sa_score"] - sd["sa_score"]) <= 1.0,
    }
    return {
        "success": all(checks.values()),
        "checks": checks,
        "delta_mw": round(td["mw"] - sd["mw"], 4),
        "delta_logp": round(td["logp"] - sd["logp"], 4),
        "delta_qed": round(td["qed"] - sd["qed"], 6),
        "delta_sa_score": round(td["sa_score"] - sd["sa_score"], 6),
    }


def tolerance_for(endpoint: str, input_value: float) -> tuple[float, str]:
    if "_p" in endpoint or endpoint.endswith(("pIC50", "pKI", "pKD", "pEC50", "pAC50")):
        return 0.3, "+/- 0.3 log unit"
    return max(0.3, abs(input_value) * 0.25), "+/- 0.25 relative or 0.3 absolute"


def preserved_secondary(
    source: Aggregate,
    target: Aggregate,
    endpoint_values: dict[str, dict[str, float]],
) -> dict[str, Any] | None:
    source_values = endpoint_values.get(source.connectivity_key, {})
    target_values = endpoint_values.get(target.connectivity_key, {})
    shared = sorted((set(source_values) & set(target_values)) - {source.endpoint_name})
    candidates = []
    for endpoint in shared:
        input_value = source_values[endpoint]
        target_value = target_values[endpoint]
        delta = target_value - input_value
        tolerance, tolerance_text = tolerance_for(endpoint, input_value)
        abs_delta = abs(delta)
        if abs_delta <= tolerance:
            candidates.append(
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
    candidates.sort(key=lambda item: (item["abs_delta"] / max(item["tolerance"], 1e-12), item["endpoint_name"]))
    return candidates[0] if candidates else None


def candidate_groups(
    aggregates: list[Aggregate],
    args: argparse.Namespace,
) -> dict[tuple[str, str, str], list[tuple[Aggregate, MMPDecomposition]]]:
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
            groups[(agg.endpoint_name, agg.condition_bucket, decomposition.core)].append((agg, decomposition))
    return groups


def build_pair_candidates(
    aggregates: list[Aggregate],
    endpoint_values: dict[str, dict[str, float]],
    args: argparse.Namespace,
) -> tuple[dict[tuple[str, str, str, str], dict[str, Any]], dict[str, Any]]:
    rng = random.Random(args.seed)
    groups = candidate_groups(aggregates, args)
    query_map: dict[tuple[str, str, str, str], dict[str, Any]] = {}
    counters: Counter[str] = Counter()
    sample_idx = 1
    seen_pairs: set[tuple[str, str, str]] = set()

    for (_endpoint, _condition, _core), members in sorted(groups.items(), key=lambda item: (-len(item[1]), item[0])):
        if len(members) < 3:
            counters["skipped_group_lt_3"] += 1
            continue
        source_indices = list(range(len(members)))
        rng.shuffle(source_indices)
        source_indices = source_indices[: args.max_sources_per_group]
        for source_idx in source_indices:
            source, source_decomp = members[source_idx]
            candidates = [
                idx
                for idx, (target, target_decomp) in enumerate(members)
                if idx != source_idx
                and target.connectivity_key != source.connectivity_key
                and target_decomp.variable != source_decomp.variable
                and abs(target_decomp.variable_heavy_atoms - source_decomp.variable_heavy_atoms)
                <= args.max_mmp_variable_heavy_delta
            ]
            rng.shuffle(candidates)
            candidates = candidates[: args.max_candidates_per_source]
            for target_idx in candidates:
                target, target_decomp = members[target_idx]
                pair_key = (source.endpoint_name, source.connectivity_key, target.connectivity_key)
                if pair_key in seen_pairs:
                    continue
                seen_pairs.add(pair_key)
                delta = target.value - source.value
                if delta < args.min_primary_delta:
                    counters["skipped_primary_delta"] += 1
                    continue
                sim = tanimoto(source, target)
                if sim < args.min_tanimoto:
                    counters["skipped_low_tanimoto"] += 1
                    continue
                descriptor_result = descriptor_constraints(source, target)
                if not descriptor_result["success"]:
                    counters["skipped_descriptor_constraints"] += 1
                    continue
                preserved = preserved_secondary(source, target, endpoint_values)
                if preserved is None:
                    counters["skipped_no_preserved_secondary"] += 1
                    continue

                qkey = (source.endpoint_name, source.condition_bucket, source.connectivity_key, preserved["endpoint_name"])
                if qkey not in query_map:
                    query_map[qkey] = {
                        "source": source,
                        "preserved": preserved,
                        "answers": {},
                        "mmp_cores": set(),
                    }
                answer_key = target.connectivity_key
                answer = {
                    "sample_id": f"pubchem_fragmp_{sample_idx:08d}",
                    "target_smiles_canon": target.smiles_canon,
                    "target_connectivity_key": target.connectivity_key,
                    "target_pubchem_cid": target.pubchem_cid,
                    "value_before": round(source.value, 6),
                    "value_after": round(target.value, 6),
                    "delta_value": round(delta, 6),
                    "relative_delta": round(delta / abs(source.value), 6) if source.value else "",
                    "tanimoto_similarity": round(sim, 6),
                    "mmp_core": source_decomp.core,
                    "input_variable_fragment": source_decomp.variable,
                    "target_variable_fragment": target_decomp.variable,
                    "secondary_objectives": [preserved],
                    "descriptor_constraints": descriptor_result,
                    "input_observation_ids": source.observation_ids[:10],
                    "target_observation_ids": target.observation_ids[:10],
                    "input_assay_ids": [source.condition_bucket],
                    "target_assay_ids": [target.condition_bucket],
                }
                old = query_map[qkey]["answers"].get(answer_key)
                if old is None or answer["delta_value"] > old["delta_value"]:
                    query_map[qkey]["answers"][answer_key] = answer
                query_map[qkey]["mmp_cores"].add(source_decomp.core)
                sample_idx += 1
                counters["candidate_positive"] += 1

    return query_map, {
        "groups": len(groups),
        "candidate_queries": len(query_map),
        "counters": dict(counters),
    }


def instruction(primary_endpoint: str, preserved: dict[str, Any]) -> str:
    return (
        f"Increase {primary_endpoint} while preserving {preserved['endpoint_name']} within "
        f"{preserved['tolerance_text']}; also keep MW, LogP, QED, and synthetic accessibility "
        "within local edit constraints."
    )


def materialize_queries(query_map: dict[tuple[str, str, str, str], dict[str, Any]], args: argparse.Namespace) -> list[dict[str, Any]]:
    rows = []
    for item in query_map.values():
        answers = sorted(
            item["answers"].values(),
            key=lambda answer: (-float(answer["delta_value"]), -float(answer["tanimoto_similarity"]), answer["target_smiles_canon"]),
        )
        if len(answers) < 2:
            continue
        selected = answers[:2]
        source: Aggregate = item["source"]
        preserved = item["preserved"]
        primary = {
            "endpoint_name": source.endpoint_name,
            "direction": "increase",
            "value_before": round(source.value, 6),
            "min_delta": args.min_primary_delta,
        }
        local_constraints = {
            "objective_group": "local_constraints",
            "properties": ["MW", "LogP", "QED", "SA"],
            "rule": "Keep within PubChem BRICS/MMP local edit constraints.",
        }
        rows.append(
            {
                "query_id": "",
                "split": "",
                "instruction": instruction(source.endpoint_name, preserved),
                "input_smiles_canon": source.smiles_canon,
                "input_connectivity_key": source.connectivity_key,
                "input_pubchem_cid": source.pubchem_cid,
                "primary_endpoint": source.endpoint_name,
                "condition_bucket": source.condition_bucket,
                "primary_direction": "increase",
                "num_property_objectives": 3,
                "primary_objective_json": compact_json(primary),
                "preserved_property_json": compact_json(preserved),
                "local_constraints_json": compact_json(local_constraints),
                "num_positive_answers": 2,
                "positive_answer_smiles_json": compact_json([answer["target_smiles_canon"] for answer in selected]),
                "positive_answers_json": compact_json(selected),
                "mmp_cores_json": compact_json(sorted(item["mmp_cores"])),
                "source_positive_sample_ids_json": compact_json([answer["sample_id"] for answer in selected]),
                "selection_rule": "PubChem-only: same assay-local primary endpoint, same BRICS/MMP core, top 2 positive targets by primary delta; negative candidates omitted.",
            }
        )
    rows.sort(key=lambda row: (row["primary_endpoint"], row["input_connectivity_key"], row["condition_bucket"]))
    if args.max_queries:
        rows = rows[: args.max_queries]
    return rows


def query_nodes(row: dict[str, Any]) -> set[str]:
    nodes = {row["input_connectivity_key"]}
    for answer in json.loads(row["positive_answers_json"]):
        if answer.get("target_connectivity_key"):
            nodes.add(answer["target_connectivity_key"])
    return nodes


def split_rows(rows: list[dict[str, Any]], seed: int) -> dict[str, list[dict[str, Any]]]:
    adjacency: dict[str, set[str]] = defaultdict(set)
    node_to_rows: dict[str, list[int]] = defaultdict(list)
    for idx, row in enumerate(rows):
        nodes = sorted(query_nodes(row))
        for node in nodes:
            node_to_rows[node].append(idx)
        for left in nodes:
            adjacency[left].update(node for node in nodes if node != left)

    components = []
    seen = set()
    for node in sorted(adjacency):
        if node in seen:
            continue
        queue = deque([node])
        seen.add(node)
        nodes = set()
        indices = set()
        while queue:
            current = queue.popleft()
            nodes.add(current)
            indices.update(node_to_rows[current])
            for neighbor in adjacency[current]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)
        components.append(sorted(indices))

    rng = random.Random(seed)
    rng.shuffle(components)
    components.sort(key=len, reverse=True)
    targets = {"train": len(rows) * 0.8, "val": len(rows) * 0.1, "test": len(rows) * 0.1}
    split_counts = {split: 0 for split in targets}
    out = {"train": [], "val": [], "test": []}
    for component in components:
        split = min(targets, key=lambda name: (split_counts[name] + len(component)) / max(targets[name], 1.0))
        split_counts[split] += len(component)
        for idx in component:
            row = dict(rows[idx])
            row["split"] = split
            out[split].append(row)
    return out


def add_ids(rows_by_split: dict[str, list[dict[str, Any]]]) -> list[dict[str, Any]]:
    all_rows = []
    for split in ["train", "val", "test"]:
        rows_by_split[split].sort(key=lambda row: (row["primary_endpoint"], row["input_connectivity_key"]))
        for row in rows_by_split[split]:
            payload = f"{split}|{row['primary_endpoint']}|{row['condition_bucket']}|{row['input_connectivity_key']}|{row['positive_answer_smiles_json']}"
            digest = hashlib.sha1(payload.encode("utf-8")).hexdigest()[:12]
            row["query_id"] = f"pubchem3p2pos_{split}_{digest}"
            all_rows.append(row)
    all_rows.sort(key=lambda row: (row["split"], row["query_id"]))
    return all_rows


def stats(rows: list[dict[str, Any]]) -> dict[str, Any]:
    outputs = set()
    secondary = set()
    for row in rows:
        for answer in json.loads(row["positive_answers_json"]):
            outputs.add(answer.get("target_connectivity_key") or answer.get("target_smiles_canon"))
        preserved = json.loads(row["preserved_property_json"])
        if preserved.get("endpoint_name"):
            secondary.add(preserved["endpoint_name"])
    return {
        "queries": len(rows),
        "positive_answers": len(rows) * 2,
        "answers_per_query": 2,
        "property_objectives_per_query": 3,
        "unique_inputs": len({row["input_connectivity_key"] for row in rows}),
        "unique_positive_outputs": len(outputs),
        "primary_endpoints": len({row["primary_endpoint"] for row in rows}),
        "preserved_secondary_endpoints": len(secondary),
    }


def leakage_check(rows_by_split: dict[str, list[dict[str, Any]]]) -> dict[str, Any]:
    molecule_sets = {split: set() for split in rows_by_split}
    for split, rows in rows_by_split.items():
        for row in rows:
            molecule_sets[split].update(query_nodes(row))
    overlaps = {}
    for left, right in [("train", "val"), ("train", "test"), ("val", "test")]:
        overlaps[f"{left}_{right}"] = len(molecule_sets[left] & molecule_sets[right])
    return {
        "molecule_overlap_counts": overlaps,
        "passes_no_molecule_overlap": all(value == 0 for value in overlaps.values()),
    }


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    aggregates, endpoint_values, load_stats = load_pubchem_aggregates(
        args.observations,
        args.min_primary_molecules,
        args.max_primary_endpoints,
        args.max_aggregates_per_endpoint,
        args.seed,
    )
    query_map, candidate_stats = build_pair_candidates(aggregates, endpoint_values, args)
    rows = materialize_queries(query_map, args)
    rows_by_split = split_rows(rows, args.seed)
    all_rows = add_ids(rows_by_split)

    write_csv(args.out_dir / "all.csv", all_rows, FIELDS)
    for split in ["train", "val", "test"]:
        write_csv(args.out_dir / f"{split}.csv", rows_by_split[split], FIELDS)

    parameters = {key: (str(value) if isinstance(value, Path) else value) for key, value in vars(args).items()}
    summary = {
        "schema": "pubchem_one_input_3property_instruction_to_2positive_answers_v1",
        "source_observations": str(args.observations),
        "row_policy": "Keep rows with at least 2 positive PubChem outputs. Emit exactly 2 positive answers and no negative candidates.",
        "property_policy": "Every instruction has exactly 3 objective groups: primary PubChem assay-local endpoint, one preserved PubChem experimental secondary endpoint, and local MW/LogP/QED/SA constraints.",
        "selection_policy": {
            "primary": "Only exact PubChem measurements with p-scale concentration-like endpoints are used as primary targets.",
            "candidate": "Same assay-local primary endpoint and same BRICS/MMP core; target must improve primary p-scale value.",
            "secondary": "Input and target must share at least one secondary PubChem endpoint preserved within tolerance.",
            "splitting": "Connected components over input/output connectivity keys to reduce molecule leakage.",
        },
        "parameters": parameters,
        "load": load_stats,
        "candidate_generation": candidate_stats,
        "leakage_check": leakage_check(rows_by_split),
        "files": {
            "all": {"path": str(args.out_dir / "all.csv"), "stats": stats(all_rows)},
            "train": {"path": str(args.out_dir / "train.csv"), "stats": stats(rows_by_split["train"])},
            "val": {"path": str(args.out_dir / "val.csv"), "stats": stats(rows_by_split["val"])},
            "test": {"path": str(args.out_dir / "test.csv"), "stats": stats(rows_by_split["test"])},
        },
    }
    (args.out_dir / "dataset_stats.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    (args.out_dir / "README.md").write_text(
        "# PubChem 3-Property 2-Positive Dataset\n\n"
        "Separate PubChem-only 3-property/2-positive dataset. It is not merged with `data/chembl_3prop_2pos`.\n\n"
        "Format:\n\n"
        "```text\n"
        "1 input molecule + 3-property English instruction -> exactly 2 positive PubChem output molecules\n"
        "```\n\n"
        "Files: `all.csv`, `train.csv`, `val.csv`, `test.csv`, `dataset_stats.json`.\n\n"
        "Caveat: primary endpoints are PubChem assay-local endpoint names such as `PubChem_AID4411_pKI`; "
        "target/assay metadata enrichment has not yet been applied.\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())



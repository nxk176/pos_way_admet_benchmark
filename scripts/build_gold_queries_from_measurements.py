from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
import re
import statistics
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdFingerprintGenerator, rdMolDescriptors
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit.Contrib.SA_Score import sascorer
except Exception as exc:  # noqa: BLE001 - this build stage is intentionally RDKit-backed.
    raise SystemExit(
        "ERROR: build_gold_queries_from_measurements.py requires RDKit. "
        "Run it with the project environment, for example .\\myenv311\\Scripts\\python.exe"
    ) from exc


SCHEMA_VERSION = "posway-admet-edit-v0.3-chembl-gold"
ORGANIC_ATOMS = {6, 7, 8, 9, 15, 16, 17, 35, 53}


ENDPOINT_TASKS: dict[str, dict[str, Any]] = {
    "logS_mol_L": {
        "direction": "increase",
        "min_abs_delta": 0.5,
        "unit_canonical": "log10_mol_L",
        "question_template": "increase_experimental_solubility_keep_scaffold",
        "question_text": (
            "Increase experimentally measured aqueous solubility by at least 0.5 log10 mol/L "
            "while preserving the Murcko scaffold."
        ),
        "hard_constraints": {"scaffold_retained": True, "min_tanimoto_similarity": 0.3},
    },
    "Caco2_logPapp_cm_s": {
        "direction": "increase",
        "min_abs_delta": 0.35,
        "unit_canonical": "log10_cm_s",
        "question_template": "increase_experimental_caco2_papp_keep_scaffold",
        "question_text": (
            "Increase experimentally measured apparent permeability by at least 0.35 log10 cm/s "
            "while preserving the Murcko scaffold."
        ),
        "hard_constraints": {"scaffold_retained": True, "min_tanimoto_similarity": 0.3},
    },
    "hERG_pIC50": {
        "direction": "decrease",
        "min_abs_delta": 0.5,
        "unit_canonical": "pIC50",
        "question_template": "decrease_experimental_herg_inhibition_keep_scaffold",
        "question_text": (
            "Decrease experimentally measured hERG inhibition potency by at least 0.5 pIC50 "
            "while preserving the Murcko scaffold."
        ),
        "hard_constraints": {"scaffold_retained": True, "min_tanimoto_similarity": 0.3},
    },
    "microsomal_clearance_mL_min_kg": {
        "direction": "decrease",
        "min_abs_delta": 5.0,
        "min_relative_delta": 0.25,
        "unit_canonical": "mL_min_kg",
        "question_template": "decrease_experimental_microsomal_clearance_keep_scaffold",
        "question_text": (
            "Decrease experimentally measured microsomal clearance by at least 25 percent "
            "and at least 5 mL/min/kg while preserving the Murcko scaffold."
        ),
        "hard_constraints": {
            "scaffold_retained": True,
            "min_tanimoto_similarity": 0.3,
            "min_relative_improvement": 0.25,
        },
    },
    "half_life_min": {
        "direction": "increase",
        "min_abs_delta": 30.0,
        "min_relative_delta": 0.25,
        "unit_canonical": "min",
        "question_template": "increase_experimental_half_life_keep_scaffold",
        "question_text": (
            "Increase experimentally measured half-life by at least 25 percent and at least "
            "30 minutes while preserving the Murcko scaffold."
        ),
        "hard_constraints": {
            "scaffold_retained": True,
            "min_tanimoto_similarity": 0.3,
            "min_relative_improvement": 0.25,
        },
    },
}


@dataclass
class Aggregate:
    endpoint_name: str
    condition_bucket: str
    molecule_chembl_id: str
    smiles_raw: str
    smiles_canon: str
    inchikey: str
    connectivity_key: str
    scaffold: str
    values: list[float] = field(default_factory=list)
    confidences: list[float] = field(default_factory=list)
    measurement_ids: list[str] = field(default_factory=list)
    assay_ids: set[str] = field(default_factory=set)
    assay_chembl_ids: set[str] = field(default_factory=set)
    document_years: set[int] = field(default_factory=set)
    descriptions: list[str] = field(default_factory=list)
    mol: Any = None
    fp: Any = None

    @property
    def value(self) -> float:
        return float(statistics.median(self.values))

    @property
    def confidence(self) -> float:
        return float(statistics.mean(self.confidences)) if self.confidences else 0.7

    @property
    def value_range(self) -> float:
        return float(max(self.values) - min(self.values)) if self.values else 0.0


def load_jsonl(path: Path) -> Iterable[dict[str, Any]]:
    with path.open("r", encoding="utf-8") as handle:
        for lineno, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            try:
                row = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"{path}:{lineno}: invalid JSON: {exc}") from exc
            if not isinstance(row, dict):
                raise ValueError(f"{path}:{lineno}: row must be a JSON object")
            yield row


def write_jsonl(path: Path, rows: Iterable[dict[str, Any]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False, separators=(",", ":")) + "\n")
            count += 1
    return count


def stable_split(key: str, train: float, valid: float) -> str:
    value = int(hashlib.sha1(key.encode("utf-8")).hexdigest()[:8], 16) / 0xFFFFFFFF
    if value < train:
        return "train"
    if value < train + valid:
        return "valid"
    return "test"


def slug(text: Any, fallback: str = "unknown") -> str:
    cleaned = re.sub(r"[^a-z0-9]+", "_", str(text or "").strip().lower()).strip("_")
    return cleaned[:48] if cleaned else fallback


def condition_bucket(row: dict[str, Any]) -> str:
    endpoint = row["endpoint_name"]
    assay = row.get("assay") or {}
    description = str(assay.get("description") or "").lower()
    organism = slug(assay.get("assay_organism"))
    cell_type = slug(assay.get("assay_cell_type"))

    if endpoint == "logS_mol_L":
        ph = re.search(r"ph\s*([0-9]+(?:\.[0-9]+)?)", description)
        ph_bucket = f"ph_{ph.group(1)}" if ph else "ph_unknown"
        return f"{endpoint}|{ph_bucket}"
    if endpoint == "Caco2_logPapp_cm_s":
        direction = "bidirectional"
        if "a to b" in description or "a-b" in description:
            direction = "a_to_b"
        elif "b to a" in description or "b-a" in description:
            direction = "b_to_a"
        return f"{endpoint}|{direction}"
    if endpoint == "hERG_pIC50":
        return f"{endpoint}|{organism}|{cell_type}"
    if endpoint == "microsomal_clearance_mL_min_kg":
        matrix = "microsome" if "microsom" in description else "clearance"
        return f"{endpoint}|{organism}|{matrix}"
    if endpoint == "half_life_min":
        matrix = "microsome" if "microsom" in description else "plasma_or_systemic"
        return f"{endpoint}|{organism}|{matrix}"
    return f"{endpoint}|generic"


def passes_structure_filters(mol: Chem.Mol, min_heavy_atoms: int, max_heavy_atoms: int) -> bool:
    heavy_atoms = mol.GetNumHeavyAtoms()
    if heavy_atoms < min_heavy_atoms or heavy_atoms > max_heavy_atoms:
        return False
    if any(atom.GetAtomicNum() not in ORGANIC_ATOMS for atom in mol.GetAtoms()):
        return False
    if abs(Chem.GetFormalCharge(mol)) > 1:
        return False
    return True


def descriptors(mol: Chem.Mol) -> dict[str, Any]:
    return {
        "mw": round(float(Descriptors.MolWt(mol)), 4),
        "logp": round(float(Crippen.MolLogP(mol)), 4),
        "tpsa": round(float(rdMolDescriptors.CalcTPSA(mol)), 4),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "rotatable_bonds": int(Lipinski.NumRotatableBonds(mol)),
        "heavy_atoms": int(mol.GetNumHeavyAtoms()),
        "qed": round(float(QED.qed(mol)), 4),
        "sa_score": round(float(sascorer.calculateScore(mol)), 4),
    }


def aggregate_measurements(args: argparse.Namespace) -> tuple[list[Aggregate], dict[str, Any]]:
    fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    aggregates: dict[tuple[str, str, str], Aggregate] = {}
    counters: Counter[str] = Counter()

    for row in load_jsonl(args.measurements):
        counters["measurement_rows_seen"] += 1
        endpoint = row.get("endpoint_name")
        if endpoint not in ENDPOINT_TASKS:
            counters["unsupported_endpoint"] += 1
            continue
        if row.get("label_type") != "experimental" or row.get("experimental_only_flag") is not True:
            counters["not_experimental"] += 1
            continue
        if args.exact_relation_only and not (row.get("quality_flags") or {}).get("standard_relation_exact", False):
            counters["non_exact_relation"] += 1
            continue
        confidence = float(row.get("confidence") or 0)
        if confidence < args.min_confidence:
            counters["low_confidence"] += 1
            continue
        value = row.get("value_canonical")
        if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
            counters["invalid_value"] += 1
            continue

        mol = Chem.MolFromSmiles(str(row.get("smiles_canon") or row.get("smiles_raw") or ""))
        if mol is None:
            counters["invalid_smiles"] += 1
            continue
        if not passes_structure_filters(mol, args.min_heavy_atoms, args.max_heavy_atoms):
            counters["structure_filtered"] += 1
            continue
        scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
        if scaffold_mol is None or scaffold_mol.GetNumHeavyAtoms() < args.min_scaffold_heavy_atoms:
            counters["missing_or_small_scaffold"] += 1
            continue

        scaffold = Chem.MolToSmiles(scaffold_mol, canonical=True, isomericSmiles=False)
        bucket = condition_bucket(row)
        connectivity_key = str(row.get("connectivity_key") or str(row.get("inchikey", "")).split("-", 1)[0])
        key = (endpoint, bucket, connectivity_key)
        assay = row.get("assay") or {}
        document = row.get("document") or {}

        if key not in aggregates:
            aggregates[key] = Aggregate(
                endpoint_name=endpoint,
                condition_bucket=bucket,
                molecule_chembl_id=str(row.get("molecule_chembl_id") or ""),
                smiles_raw=str(row.get("smiles_raw") or row.get("smiles_canon") or ""),
                smiles_canon=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
                inchikey=str(row.get("inchikey") or ""),
                connectivity_key=connectivity_key,
                scaffold=scaffold,
                mol=mol,
                fp=fpgen.GetFingerprint(mol),
            )

        agg = aggregates[key]
        agg.values.append(float(value))
        agg.confidences.append(confidence)
        agg.measurement_ids.append(str(row.get("measurement_id") or ""))
        if assay.get("assay_id") is not None:
            agg.assay_ids.add(str(assay["assay_id"]))
        if assay.get("assay_chembl_id"):
            agg.assay_chembl_ids.add(str(assay["assay_chembl_id"]))
        if isinstance(document.get("year"), int):
            agg.document_years.add(int(document["year"]))
        if assay.get("description") and len(agg.descriptions) < 3:
            agg.descriptions.append(str(assay["description"]))

    kept: list[Aggregate] = []
    for agg in aggregates.values():
        if len(agg.values) < args.min_replicates:
            counters["below_min_replicates"] += 1
            continue
        kept.append(agg)
    counters["aggregates_kept"] = len(kept)
    return kept, dict(counters)


def improved(source: Aggregate, target: Aggregate) -> tuple[bool, float, float | None]:
    task = ENDPOINT_TASKS[source.endpoint_name]
    before = source.value
    after = target.value
    delta = after - before
    min_abs_delta = float(task["min_abs_delta"])
    min_rel = task.get("min_relative_delta")

    if task["direction"] == "increase":
        ok = delta >= min_abs_delta
        rel = (delta / abs(before)) if before else None
    else:
        ok = delta <= -min_abs_delta
        rel = ((before - after) / abs(before)) if before else None

    if min_rel is not None:
        ok = ok and rel is not None and rel >= float(min_rel)
    return ok, delta, rel


def make_query(query_id: str, source: Aggregate, split: str, num_answers: int) -> dict[str, Any]:
    task = ENDPOINT_TASKS[source.endpoint_name]
    query = {
        "query_id": query_id,
        "input_smiles_raw": source.smiles_raw,
        "input_smiles_canon": source.smiles_canon,
        "input_inchikey": source.inchikey,
        "input_connectivity_key": source.connectivity_key,
        "input_murcko_scaffold": source.scaffold,
        "input_descriptors": descriptors(source.mol),
        "input_experimental_endpoint": {
            "endpoint_name": source.endpoint_name,
            "value": round(source.value, 6),
            "unit_canonical": task["unit_canonical"],
            "condition_bucket": source.condition_bucket,
            "n_measurements": len(source.values),
            "value_range": round(source.value_range, 6),
            "confidence": round(source.confidence, 3),
        },
        "question_text": task["question_text"],
        "question_template": task["question_template"],
        "target_endpoints": [
            {
                "endpoint_name": source.endpoint_name,
                "direction": task["direction"],
                "min_abs_delta": task["min_abs_delta"],
                "unit_canonical": task["unit_canonical"],
            }
        ],
        "hard_constraints": task["hard_constraints"],
        "source_pool": "gold",
        "split": split,
        "num_answers": num_answers,
        "schema_version": SCHEMA_VERSION,
        "source_chembl_id": source.molecule_chembl_id,
    }
    if "min_relative_delta" in task:
        query["target_endpoints"][0]["min_relative_delta"] = task["min_relative_delta"]
    return query


def make_answer(
    answer_id: str,
    query_id: str,
    source: Aggregate,
    target: Aggregate,
    delta: float,
    rel: float | None,
    similarity: float,
) -> dict[str, Any]:
    task = ENDPOINT_TASKS[source.endpoint_name]
    constraint_flags: dict[str, Any] = {
        "success": True,
        "scaffold_retained": source.scaffold == target.scaffold,
        "similarity_input_target": round(similarity, 4),
    }
    if rel is not None:
        constraint_flags["relative_improvement"] = round(rel, 6)

    return {
        "answer_id": answer_id,
        "query_id": query_id,
        "target_smiles_canon": target.smiles_canon,
        "target_inchikey": target.inchikey,
        "target_connectivity_key": target.connectivity_key,
        "target_descriptors": descriptors(target.mol),
        "transform_class": "same_murcko_experimental_analog",
        "transform_detail": {
            "source_chembl_id": source.molecule_chembl_id,
            "target_chembl_id": target.molecule_chembl_id,
            "murcko_scaffold": source.scaffold,
            "condition_bucket": source.condition_bucket,
        },
        "label_type": "experimental",
        "endpoint_name": source.endpoint_name,
        "value_before": round(source.value, 6),
        "value_after": round(target.value, 6),
        "delta_value": round(delta, 6),
        "unit_canonical": task["unit_canonical"],
        "confidence": round(min(source.confidence, target.confidence), 3),
        "experimental_only_flag": True,
        "provenance": {
            "source": "ChEMBL 36 SQLite experimental measurements",
            "source_pool": "gold_experimental",
            "source_measurement_ids": source.measurement_ids[:20],
            "target_measurement_ids": target.measurement_ids[:20],
            "source_n_measurements": len(source.values),
            "target_n_measurements": len(target.values),
            "source_assay_chembl_ids": sorted(source.assay_chembl_ids)[:20],
            "target_assay_chembl_ids": sorted(target.assay_chembl_ids)[:20],
            "source_document_years": sorted(source.document_years),
            "target_document_years": sorted(target.document_years),
            "condition_bucket": source.condition_bucket,
        },
        "constraint_flags": constraint_flags,
    }


def build_gold(args: argparse.Namespace, aggregates: list[Aggregate]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    rng = random.Random(args.seed)
    groups: dict[tuple[str, str, str], list[Aggregate]] = defaultdict(list)
    for agg in aggregates:
        groups[(agg.endpoint_name, agg.condition_bucket, agg.scaffold)].append(agg)

    group_items = [(key, rows) for key, rows in groups.items() if len(rows) >= args.answers_per_query + 1]
    group_items.sort(key=lambda item: (-len(item[1]), item[0]))

    queries: list[dict[str, Any]] = []
    answers: list[dict[str, Any]] = []
    query_index = 1
    answer_index = 1
    per_endpoint_counts: Counter[str] = Counter()
    group_counters: Counter[str] = Counter()

    for (endpoint, _bucket, scaffold), rows in group_items:
        if len(queries) >= args.max_queries:
            break
        if per_endpoint_counts[endpoint] >= args.max_queries_per_endpoint:
            continue

        shuffled = list(rows)
        rng.shuffle(shuffled)
        per_group_queries = 0
        split = stable_split(scaffold, args.train_fraction, args.valid_fraction)

        for source in shuffled:
            if len(queries) >= args.max_queries:
                break
            if per_endpoint_counts[endpoint] >= args.max_queries_per_endpoint:
                break
            if per_group_queries >= args.max_queries_per_group:
                break

            candidates: list[tuple[Aggregate, float, float | None, float]] = []
            for target in rows:
                if target.connectivity_key == source.connectivity_key:
                    continue
                similarity = float(DataStructs.TanimotoSimilarity(source.fp, target.fp))
                if similarity < args.min_tanimoto_similarity:
                    continue
                ok, delta, rel = improved(source, target)
                if not ok:
                    continue
                candidates.append((target, delta, rel, similarity))

            if len(candidates) < args.answers_per_query:
                group_counters["sources_without_enough_answers"] += 1
                continue

            direction = ENDPOINT_TASKS[endpoint]["direction"]
            if direction == "increase":
                candidates.sort(key=lambda item: (item[1], item[3], item[0].confidence), reverse=True)
            else:
                candidates.sort(key=lambda item: (-item[1], item[3], item[0].confidence), reverse=True)
            selected = candidates[: args.answers_per_query]

            query_id = f"qgold_{query_index:06d}_{endpoint}"
            queries.append(make_query(query_id, source, split, len(selected)))
            for target, delta, rel, similarity in selected:
                answers.append(
                    make_answer(f"agold_{answer_index:08d}", query_id, source, target, delta, rel, similarity)
                )
                answer_index += 1

            query_index += 1
            per_endpoint_counts[endpoint] += 1
            per_group_queries += 1

    stats = {
        "groups_considered": len(group_items),
        "group_counters": dict(group_counters),
        "queries_per_endpoint": dict(sorted(per_endpoint_counts.items())),
    }
    return queries, answers, stats


def prediction_rows(answers: list[dict[str, Any]]) -> Iterable[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for answer in answers:
        grouped[answer["query_id"]].append(answer)
    for query_id in sorted(grouped):
        rows = sorted(grouped[query_id], key=lambda row: (-float(row["confidence"]), row["answer_id"]))
        yield {
            "query_id": query_id,
            "ranked_smiles": [
                {
                    "target_smiles_canon": row["target_smiles_canon"],
                    "confidence": row["confidence"],
                    "model_metadata": {"baseline": "gold_oracle_answer_set"},
                }
                for row in rows
            ],
        }


def build_stats(
    args: argparse.Namespace,
    aggregate_counters: dict[str, Any],
    aggregates: list[Aggregate],
    queries: list[dict[str, Any]],
    answers: list[dict[str, Any]],
    build_extra: dict[str, Any],
) -> dict[str, Any]:
    split_counts = Counter(query["split"] for query in queries)
    endpoint_counts = Counter(answer["endpoint_name"] for answer in answers)
    aggregate_endpoint_counts = Counter(agg.endpoint_name for agg in aggregates)
    scaffold_counts = Counter(agg.scaffold for agg in aggregates)
    condition_counts = Counter(agg.condition_bucket for agg in aggregates)
    return {
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": {
            "measurements": str(args.measurements),
            "source_name": "ChEMBL",
            "source_version": "36",
            "label_status": "gold_experimental",
        },
        "parameters": {
            "max_queries": args.max_queries,
            "max_queries_per_endpoint": args.max_queries_per_endpoint,
            "answers_per_query": args.answers_per_query,
            "min_tanimoto_similarity": args.min_tanimoto_similarity,
            "min_confidence": args.min_confidence,
            "exact_relation_only": args.exact_relation_only,
            "min_replicates": args.min_replicates,
            "train_fraction": args.train_fraction,
            "valid_fraction": args.valid_fraction,
            "seed": args.seed,
        },
        "aggregate_counters": aggregate_counters,
        "num_aggregated_molecule_endpoint_records": len(aggregates),
        "aggregate_endpoint_counts": dict(sorted(aggregate_endpoint_counts.items())),
        "num_scaffolds": len(scaffold_counts),
        "num_condition_buckets": len(condition_counts),
        "num_queries": len(queries),
        "num_answers": len(answers),
        "splits": dict(sorted(split_counts.items())),
        "answer_endpoints": dict(sorted(endpoint_counts.items())),
        "label_types": {"experimental": len(answers)},
        "experimental_only_answers": len(answers),
        "build_extra": build_extra,
        "limitations": [
            "Gold records are generated from ChEMBL experimental measurements only.",
            "ToxCast/Tox21 experimental query-answer integration is still a separate stage.",
            "Condition buckets are conservative but not a substitute for manual assay harmonization.",
            "Oracle predictions are for evaluator smoke tests, not model performance claims.",
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build gold experimental query-answer sets from normalized ChEMBL measurements.")
    parser.add_argument("--measurements", type=Path, default=Path("data/experimental/chembl_admet_measurements.jsonl"))
    parser.add_argument("--out-dir", type=Path, default=Path("data/gold"))
    parser.add_argument("--predictions-out", type=Path, default=Path("predictions/gold_oracle_predictions.jsonl"))
    parser.add_argument("--max-queries", type=int, default=3000)
    parser.add_argument("--max-queries-per-endpoint", type=int, default=900)
    parser.add_argument("--max-queries-per-group", type=int, default=20)
    parser.add_argument("--answers-per-query", type=int, default=2)
    parser.add_argument("--min-tanimoto-similarity", type=float, default=0.3)
    parser.add_argument("--min-confidence", type=float, default=0.65)
    parser.add_argument("--min-replicates", type=int, default=1)
    parser.add_argument("--min-heavy-atoms", type=int, default=8)
    parser.add_argument("--max-heavy-atoms", type=int, default=60)
    parser.add_argument("--min-scaffold-heavy-atoms", type=int, default=6)
    parser.add_argument("--train-fraction", type=float, default=0.8)
    parser.add_argument("--valid-fraction", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--allow-non-exact-relation", action="store_true")
    args = parser.parse_args()
    args.exact_relation_only = not args.allow_non_exact_relation
    return args


def main() -> int:
    args = parse_args()
    if not args.measurements.is_file():
        print(f"ERROR: measurement file not found: {args.measurements}", file=sys.stderr)
        return 1
    if args.answers_per_query < 2:
        print("ERROR: --answers-per-query must be at least 2.", file=sys.stderr)
        return 1
    if args.train_fraction <= 0 or args.valid_fraction < 0 or args.train_fraction + args.valid_fraction >= 1:
        print("ERROR: split fractions must leave a positive test fraction.", file=sys.stderr)
        return 1

    aggregates, aggregate_counters = aggregate_measurements(args)
    queries, answers, build_extra = build_gold(args, aggregates)
    if not queries:
        print("ERROR: no gold query-answer records could be generated with the selected thresholds.", file=sys.stderr)
        return 1

    write_jsonl(args.out_dir / "queries.jsonl", queries)
    write_jsonl(args.out_dir / "answers.jsonl", answers)
    write_jsonl(args.predictions_out, prediction_rows(answers))

    stats = build_stats(args, aggregate_counters, aggregates, queries, answers, build_extra)
    stats_path = args.out_dir / "dataset_stats.json"
    stats_path.parent.mkdir(parents=True, exist_ok=True)
    stats_path.write_text(json.dumps(stats, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(stats, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

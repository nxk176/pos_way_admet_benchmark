from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import random
import statistics
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdFingerprintGenerator, rdMolDescriptors
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit.Contrib.SA_Score import sascorer
except Exception as exc:  # noqa: BLE001 - this script is intentionally RDKit-backed.
    raise SystemExit(
        "ERROR: build_chembl_proxy_dataset.py requires RDKit. "
        "Run it with the project environment, for example .\\myenv311\\Scripts\\python.exe"
    ) from exc


SCHEMA_VERSION = "posway-admet-edit-v0.2-chembl-proxy"
ORGANIC_ATOMS = {6, 7, 8, 9, 15, 16, 17, 35, 53}


@dataclass(frozen=True)
class MoleculeRecord:
    chembl_id: str
    smiles_raw: str
    smiles_canon: str
    inchikey: str
    connectivity_key: str
    scaffold: str
    mw: float
    logp: float
    tpsa: float
    hba: int
    hbd: int
    rotatable_bonds: int
    heavy_atoms: int
    qed: float
    sa_score: float
    fp: Any


TASKS = {
    "decrease_logp": {
        "endpoint_name": "RDKit_MolLogP",
        "direction": "decrease",
        "unit_canonical": "unitless",
        "min_abs_delta": 0.5,
        "question_template": "decrease_calculated_logp_keep_scaffold_mw_window",
        "question_text": (
            "Decrease calculated logP by at least 0.5 while preserving the Murcko scaffold "
            "and keeping molecular-weight change within 80 Da."
        ),
        "hard_constraints": {
            "scaffold_retained": True,
            "max_abs_delta_mw": 80.0,
            "min_tanimoto_similarity": 0.35,
        },
    },
    "decrease_tpsa": {
        "endpoint_name": "RDKit_TPSA",
        "direction": "decrease",
        "unit_canonical": "A2",
        "min_abs_delta": 10.0,
        "question_template": "decrease_tpsa_keep_scaffold_logp_window",
        "question_text": (
            "Decrease topological polar surface area by at least 10 A2 while preserving the "
            "Murcko scaffold and keeping calculated logP within 1.0."
        ),
        "hard_constraints": {
            "scaffold_retained": True,
            "max_abs_delta_logp": 1.0,
            "min_tanimoto_similarity": 0.35,
        },
    },
}


def stable_split(key: str, train: float, valid: float) -> str:
    value = int(hashlib.sha1(key.encode("utf-8")).hexdigest()[:8], 16) / 0xFFFFFFFF
    if value < train:
        return "train"
    if value < train + valid:
        return "valid"
    return "test"


def rounded(value: float) -> float:
    return round(float(value), 4)


def descriptor_dict(row: MoleculeRecord) -> dict[str, Any]:
    return {
        "mw": rounded(row.mw),
        "logp": rounded(row.logp),
        "tpsa": rounded(row.tpsa),
        "hba": row.hba,
        "hbd": row.hbd,
        "rotatable_bonds": row.rotatable_bonds,
        "heavy_atoms": row.heavy_atoms,
        "qed": rounded(row.qed),
        "sa_score": rounded(row.sa_score),
    }


def passes_basic_filters(mol: Chem.Mol, min_heavy_atoms: int, max_heavy_atoms: int) -> bool:
    heavy_atoms = mol.GetNumHeavyAtoms()
    if heavy_atoms < min_heavy_atoms or heavy_atoms > max_heavy_atoms:
        return False
    if any(atom.GetAtomicNum() not in ORGANIC_ATOMS for atom in mol.GetAtoms()):
        return False
    if abs(Chem.GetFormalCharge(mol)) > 1:
        return False
    return True


def iter_chembl_chemreps(path: Path) -> Iterable[tuple[str, str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        expected = ["chembl_id", "canonical_smiles", "standard_inchi", "standard_inchi_key"]
        if header[:4] != expected:
            raise ValueError(f"Unexpected ChEMBL chemreps header: {header!r}")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            yield parts[0], parts[1], parts[3]


def load_molecules(args: argparse.Namespace) -> tuple[list[MoleculeRecord], dict[str, Any]]:
    fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    molecules: list[MoleculeRecord] = []
    seen_connectivity: set[str] = set()
    counters: Counter[str] = Counter()

    for chembl_id, smiles, inchikey in iter_chembl_chemreps(args.chemreps):
        counters["raw_rows_seen"] += 1
        if args.max_raw_rows and counters["raw_rows_seen"] > args.max_raw_rows:
            break
        if not smiles or not inchikey:
            counters["missing_smiles_or_inchikey"] += 1
            continue

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            counters["invalid_smiles"] += 1
            continue
        if not passes_basic_filters(mol, args.min_heavy_atoms, args.max_heavy_atoms):
            counters["filtered_basic_structure"] += 1
            continue

        connectivity_key = inchikey.split("-", 1)[0]
        if connectivity_key in seen_connectivity:
            counters["duplicate_connectivity"] += 1
            continue

        scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
        if scaffold_mol is None or scaffold_mol.GetNumHeavyAtoms() < args.min_scaffold_heavy_atoms:
            counters["missing_or_small_scaffold"] += 1
            continue
        scaffold = Chem.MolToSmiles(scaffold_mol, canonical=True, isomericSmiles=False)
        if not scaffold:
            counters["empty_scaffold"] += 1
            continue

        smiles_canon = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        seen_connectivity.add(connectivity_key)
        molecules.append(
            MoleculeRecord(
                chembl_id=chembl_id,
                smiles_raw=smiles,
                smiles_canon=smiles_canon,
                inchikey=inchikey,
                connectivity_key=connectivity_key,
                scaffold=scaffold,
                mw=Descriptors.MolWt(mol),
                logp=Crippen.MolLogP(mol),
                tpsa=rdMolDescriptors.CalcTPSA(mol),
                hba=Lipinski.NumHAcceptors(mol),
                hbd=Lipinski.NumHDonors(mol),
                rotatable_bonds=Lipinski.NumRotatableBonds(mol),
                heavy_atoms=mol.GetNumHeavyAtoms(),
                qed=QED.qed(mol),
                sa_score=sascorer.calculateScore(mol),
                fp=fpgen.GetFingerprint(mol),
            )
        )
        counters["kept_molecules"] += 1
        if len(molecules) >= args.max_molecules:
            break

    return molecules, dict(counters)


def answer_candidate(
    source: MoleculeRecord,
    target: MoleculeRecord,
    task_name: str,
    min_similarity: float,
) -> dict[str, Any] | None:
    if source.connectivity_key == target.connectivity_key:
        return None
    similarity = float(DataStructs.TanimotoSimilarity(source.fp, target.fp))
    if similarity < min_similarity:
        return None

    if task_name == "decrease_logp":
        before = source.logp
        after = target.logp
        delta = after - before
        if delta > -TASKS[task_name]["min_abs_delta"]:
            return None
        mw_delta = target.mw - source.mw
        if abs(mw_delta) > TASKS[task_name]["hard_constraints"]["max_abs_delta_mw"]:
            return None
        extra_flags = {"mw_delta": rounded(mw_delta), "similarity_input_target": rounded(similarity)}
    elif task_name == "decrease_tpsa":
        before = source.tpsa
        after = target.tpsa
        delta = after - before
        if delta > -TASKS[task_name]["min_abs_delta"]:
            return None
        logp_delta = target.logp - source.logp
        if abs(logp_delta) > TASKS[task_name]["hard_constraints"]["max_abs_delta_logp"]:
            return None
        extra_flags = {"delta_logp": rounded(logp_delta), "similarity_input_target": rounded(similarity)}
    else:
        raise ValueError(f"Unknown task: {task_name}")

    return {
        "target": target,
        "value_before": rounded(before),
        "value_after": rounded(after),
        "delta_value": rounded(delta),
        "similarity": similarity,
        "constraint_flags": {
            "success": True,
            "scaffold_retained": source.scaffold == target.scaffold,
            **extra_flags,
        },
    }


def make_query(
    query_id: str,
    source: MoleculeRecord,
    task_name: str,
    split: str,
    num_answers: int,
) -> dict[str, Any]:
    task = TASKS[task_name]
    return {
        "query_id": query_id,
        "input_smiles_raw": source.smiles_raw,
        "input_smiles_canon": source.smiles_canon,
        "input_inchikey": source.inchikey,
        "input_connectivity_key": source.connectivity_key,
        "input_murcko_scaffold": source.scaffold,
        "input_descriptors": descriptor_dict(source),
        "question_text": task["question_text"],
        "question_template": task["question_template"],
        "target_endpoints": [
            {
                "endpoint_name": task["endpoint_name"],
                "direction": task["direction"],
                "min_abs_delta": task["min_abs_delta"],
                "unit_canonical": task["unit_canonical"],
            }
        ],
        "hard_constraints": task["hard_constraints"],
        "source_pool": "silver",
        "split": split,
        "num_answers": num_answers,
        "schema_version": SCHEMA_VERSION,
        "source_chembl_id": source.chembl_id,
    }


def make_answer(
    answer_id: str,
    query_id: str,
    source: MoleculeRecord,
    candidate: dict[str, Any],
    task_name: str,
) -> dict[str, Any]:
    task = TASKS[task_name]
    target: MoleculeRecord = candidate["target"]
    confidence = min(0.55, 0.25 + 0.25 * candidate["similarity"])
    return {
        "answer_id": answer_id,
        "query_id": query_id,
        "target_smiles_canon": target.smiles_canon,
        "target_inchikey": target.inchikey,
        "target_connectivity_key": target.connectivity_key,
        "target_descriptors": descriptor_dict(target),
        "transform_class": "same_murcko_analog_retrieval",
        "transform_detail": {
            "source_chembl_id": source.chembl_id,
            "target_chembl_id": target.chembl_id,
            "murcko_scaffold": source.scaffold,
            "selection_rule": task_name,
        },
        "label_type": "proxy",
        "endpoint_name": task["endpoint_name"],
        "value_before": candidate["value_before"],
        "value_after": candidate["value_after"],
        "delta_value": candidate["delta_value"],
        "unit_canonical": task["unit_canonical"],
        "confidence": rounded(confidence),
        "experimental_only_flag": False,
        "provenance": {
            "source": "ChEMBL 36 chemreps + RDKit descriptor proxy",
            "source_pool": "silver_proxy",
            "source_document": "raw/public/chembl/chembl_36_chemreps.txt.gz",
            "source_chembl_id": source.chembl_id,
            "target_chembl_id": target.chembl_id,
            "rdkit_proxy": True,
        },
        "constraint_flags": candidate["constraint_flags"],
    }


def build_dataset(args: argparse.Namespace, molecules: list[MoleculeRecord]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    groups: dict[str, list[MoleculeRecord]] = defaultdict(list)
    for row in molecules:
        groups[row.scaffold].append(row)

    eligible_groups = [rows for rows in groups.values() if len(rows) >= args.answers_per_query + 1]
    eligible_groups.sort(key=lambda rows: (-len(rows), rows[0].scaffold))

    queries: list[dict[str, Any]] = []
    answers: list[dict[str, Any]] = []
    query_index = 1
    answer_index = 1
    rng = random.Random(args.seed)

    for group in eligible_groups:
        if len(queries) >= args.max_queries:
            break

        split = stable_split(group[0].scaffold, args.train_fraction, args.valid_fraction)
        shuffled_group = list(group)
        rng.shuffle(shuffled_group)
        per_scaffold_queries = 0

        for source in shuffled_group:
            if len(queries) >= args.max_queries or per_scaffold_queries >= args.max_queries_per_scaffold:
                break
            for task_name in sorted(TASKS):
                if len(queries) >= args.max_queries or per_scaffold_queries >= args.max_queries_per_scaffold:
                    break

                candidates = [
                    candidate
                    for target in group
                    if (candidate := answer_candidate(source, target, task_name, args.min_tanimoto_similarity)) is not None
                ]
                if len(candidates) < args.answers_per_query:
                    continue

                candidates.sort(
                    key=lambda item: (
                        abs(float(item["delta_value"])),
                        float(item["similarity"]),
                        item["target"].chembl_id,
                    ),
                    reverse=True,
                )
                selected = candidates[: args.answers_per_query]
                query_id = f"qchembl_{query_index:06d}_{task_name}"
                queries.append(make_query(query_id, source, task_name, split, len(selected)))

                for candidate in selected:
                    answers.append(make_answer(f"achembl_{answer_index:08d}", query_id, source, candidate, task_name))
                    answer_index += 1

                query_index += 1
                per_scaffold_queries += 1

    return queries, answers


def write_jsonl(path: Path, rows: Iterable[dict[str, Any]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False, separators=(",", ":")) + "\n")
            count += 1
    return count


def molecule_json_rows(molecules: Iterable[MoleculeRecord]) -> Iterable[dict[str, Any]]:
    for row in molecules:
        yield {
            "chembl_id": row.chembl_id,
            "smiles_raw": row.smiles_raw,
            "smiles_canon": row.smiles_canon,
            "inchikey": row.inchikey,
            "connectivity_key": row.connectivity_key,
            "murcko_scaffold": row.scaffold,
            "descriptors": descriptor_dict(row),
            "source": "ChEMBL 36 chemreps",
        }


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
                    "model_metadata": {"baseline": "oracle_same_murcko_proxy"},
                }
                for row in rows
            ],
        }


def build_stats(
    args: argparse.Namespace,
    counters: dict[str, Any],
    molecules: list[MoleculeRecord],
    queries: list[dict[str, Any]],
    answers: list[dict[str, Any]],
) -> dict[str, Any]:
    scaffold_counts = Counter(row.scaffold for row in molecules)
    split_counts = Counter(query["split"] for query in queries)
    endpoint_counts = Counter(answer["endpoint_name"] for answer in answers)
    label_counts = Counter(answer["label_type"] for answer in answers)
    group_sizes = list(scaffold_counts.values())
    return {
        "schema_version": SCHEMA_VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": {
            "chemreps": str(args.chemreps),
            "chembl_version": "36",
            "label_status": "silver_proxy_rdkit_descriptors",
        },
        "parameters": {
            "max_molecules": args.max_molecules,
            "max_queries": args.max_queries,
            "answers_per_query": args.answers_per_query,
            "min_tanimoto_similarity": args.min_tanimoto_similarity,
            "min_heavy_atoms": args.min_heavy_atoms,
            "max_heavy_atoms": args.max_heavy_atoms,
            "min_scaffold_heavy_atoms": args.min_scaffold_heavy_atoms,
            "train_fraction": args.train_fraction,
            "valid_fraction": args.valid_fraction,
            "seed": args.seed,
        },
        "raw_filter_counters": counters,
        "num_molecules": len(molecules),
        "num_scaffolds": len(scaffold_counts),
        "scaffold_group_size": {
            "min": min(group_sizes) if group_sizes else 0,
            "median": statistics.median(group_sizes) if group_sizes else 0,
            "max": max(group_sizes) if group_sizes else 0,
        },
        "num_queries": len(queries),
        "num_answers": len(answers),
        "splits": dict(sorted(split_counts.items())),
        "endpoints": dict(sorted(endpoint_counts.items())),
        "label_types": dict(sorted(label_counts.items())),
        "answers_per_query": args.answers_per_query,
        "limitations": [
            "This build is a ChEMBL-derived silver/proxy molecular editing benchmark.",
            "Endpoint values are RDKit descriptors, not experimental ADMET measurements.",
            "Use the ChEMBL SQLite/ToxCast parser stage before claiming gold experimental benchmark status.",
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a ChEMBL-derived silver/proxy POS-WAY ADMET editing dataset.")
    parser.add_argument("--chemreps", type=Path, default=Path("raw/public/chembl/chembl_36_chemreps.txt.gz"))
    parser.add_argument("--out-root", type=Path, default=Path("."))
    parser.add_argument("--max-molecules", type=int, default=50000)
    parser.add_argument("--max-raw-rows", type=int, default=0, help="0 means no explicit raw-row cap.")
    parser.add_argument("--max-queries", type=int, default=5000)
    parser.add_argument("--max-queries-per-scaffold", type=int, default=12)
    parser.add_argument("--answers-per-query", type=int, default=3)
    parser.add_argument("--min-tanimoto-similarity", type=float, default=0.35)
    parser.add_argument("--min-heavy-atoms", type=int, default=8)
    parser.add_argument("--max-heavy-atoms", type=int, default=60)
    parser.add_argument("--min-scaffold-heavy-atoms", type=int, default=6)
    parser.add_argument("--train-fraction", type=float, default=0.8)
    parser.add_argument("--valid-fraction", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=13)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.chemreps.is_file():
        print(f"ERROR: ChEMBL chemreps file not found: {args.chemreps}", file=sys.stderr)
        return 1
    if args.answers_per_query < 2:
        print("ERROR: --answers-per-query must be at least 2.", file=sys.stderr)
        return 1
    if args.train_fraction <= 0 or args.valid_fraction < 0 or args.train_fraction + args.valid_fraction >= 1:
        print("ERROR: split fractions must leave a positive test fraction.", file=sys.stderr)
        return 1

    molecules, counters = load_molecules(args)
    if not molecules:
        print("ERROR: no molecules passed filters.", file=sys.stderr)
        return 1

    queries, answers = build_dataset(args, molecules)
    if not queries:
        print("ERROR: no query-answer records could be generated with the selected thresholds.", file=sys.stderr)
        return 1

    out_root = args.out_root
    write_jsonl(out_root / "data" / "derived" / "chembl_molecules.jsonl", molecule_json_rows(molecules))
    write_jsonl(out_root / "data" / "queries.jsonl", queries)
    write_jsonl(out_root / "data" / "answers.jsonl", answers)
    write_jsonl(out_root / "predictions" / "baseline_predictions.jsonl", prediction_rows(answers))

    stats = build_stats(args, counters, molecules, queries, answers)
    stats_path = out_root / "data" / "dataset_stats.json"
    stats_path.parent.mkdir(parents=True, exist_ok=True)
    stats_path.write_text(json.dumps(stats, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(stats, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

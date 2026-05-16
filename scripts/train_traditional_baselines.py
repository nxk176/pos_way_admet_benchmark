from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

try:
    import numpy as np
    from joblib import dump
    from scipy.sparse import csr_matrix
    from sklearn.linear_model import SGDClassifier
    from sklearn.metrics import average_precision_score, roc_auc_score
except Exception:  # noqa: BLE001 - reported cleanly at runtime.
    np = None
    dump = None
    csr_matrix = None
    SGDClassifier = None
    average_precision_score = None
    roc_auc_score = None

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, Crippen, Descriptors, QED
except Exception:  # noqa: BLE001 - reported cleanly at runtime.
    Chem = None
    RDLogger = None
    AllChem = None
    Crippen = None
    Descriptors = None
    QED = None


TOKEN_RE = re.compile(r"[A-Za-z0-9_+\-.]+")


@dataclass
class CandidateExample:
    query_id: str
    split: str
    input_smiles: str
    candidate_smiles: str
    instruction: str
    label: int
    source: str = ""
    primary_endpoint: str = ""
    primary_direction: str = ""
    target_id: str = ""
    target_name: str = ""
    measurement_type: str = ""
    measurement_group: str = ""
    evidence_count: float = 0.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Train traditional supervised baselines for molecule edit/ranking tasks. "
            "These are candidate rankers, not de novo generators."
        )
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    edit = subparsers.add_parser(
        "edit-ranker",
        help="Train a candidate ranker for ChEMBL/PubChem/Papyrus 3-property edit datasets.",
    )
    edit.add_argument("--data-dir", type=Path, default=Path("data/chembl_3prop_2pos"))
    edit.add_argument("--source-name", default="")
    edit.add_argument("--eval-splits", nargs="+", default=["val", "test"], choices=["val", "test"])
    edit.add_argument("--decoys-per-query", type=int, default=20)
    edit.add_argument("--max-train-queries", type=int, default=0)
    edit.add_argument("--max-eval-queries", type=int, default=0)
    add_common_args(edit)

    bindingdb = subparsers.add_parser(
        "bindingdb-ranker",
        help="Train a pairwise candidate ranker for BindingDB target-conditioned triplets.",
    )
    bindingdb.add_argument("--data-dir", type=Path, default=Path("data/bindingdb_target_conditioned"))
    bindingdb.add_argument(
        "--eval-splits",
        nargs="+",
        default=["val", "test_seen_target", "test_unseen_target"],
        choices=["val", "test_seen_target", "test_unseen_target"],
    )
    bindingdb.add_argument("--max-train-rows", type=int, default=0)
    bindingdb.add_argument("--max-eval-rows", type=int, default=0)
    add_common_args(bindingdb)

    return parser.parse_args()


def add_common_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--out-dir", type=Path, default=Path("reports/traditional_baselines"))
    parser.add_argument("--fp-bits", type=int, default=2048)
    parser.add_argument("--fp-radius", type=int, default=2)
    parser.add_argument("--text-features", type=int, default=1024)
    parser.add_argument("--categorical-features", type=int, default=1024)
    parser.add_argument("--max-instruction-tokens", type=int, default=96)
    parser.add_argument("--seed", type=int, default=29)
    parser.add_argument("--max-iter", type=int, default=25)
    parser.add_argument("--alpha", type=float, default=1e-5)
    parser.add_argument("--save-model", action="store_true")


def require_dependencies() -> None:
    missing = []
    if np is None or csr_matrix is None or SGDClassifier is None:
        missing.append("numpy scipy scikit-learn joblib")
    if Chem is None or AllChem is None:
        missing.append("rdkit")
    if missing:
        raise SystemExit(
            "ERROR: missing baseline dependencies: "
            + "; ".join(missing)
            + ". Use the project RDKit environment, for example ..\\myenv311\\Scripts\\python.exe."
        )
    if RDLogger is not None:
        RDLogger.DisableLog("rdApp.warning")
        RDLogger.DisableLog("rdApp.error")


def stable_hash(value: str) -> int:
    return int(hashlib.sha1(value.encode("utf-8", errors="ignore")).hexdigest()[:12], 16)


def read_csv(path: Path, limit: int = 0) -> list[dict[str, str]]:
    if not path.is_file():
        raise SystemExit(f"ERROR: input file not found: {path}")
    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append(row)
            if limit and len(rows) >= limit:
                break
    return rows


def safe_json(raw: Any, fallback: Any) -> Any:
    try:
        return json.loads(raw)
    except (TypeError, json.JSONDecodeError):
        return fallback


def clean_text(value: Any) -> str:
    return str(value or "").replace("\r", " ").replace("\n", " ").strip()


def positive_smiles_from_row(row: dict[str, str]) -> list[str]:
    smiles = safe_json(row.get("positive_answer_smiles_json"), [])
    if not isinstance(smiles, list):
        return []
    return [clean_text(item) for item in smiles if clean_text(item)]


def infer_source_name(data_dir: Path, explicit: str) -> str:
    if explicit:
        return explicit
    name = data_dir.name.lower()
    if "chembl" in name:
        return "ChEMBL"
    if "pubchem" in name:
        return "PubChem"
    if "papyrus" in name:
        return "Papyrus"
    return data_dir.name


class MoleculeFeaturizer:
    def __init__(
        self,
        fp_bits: int,
        fp_radius: int,
        text_features: int,
        categorical_features: int,
        max_instruction_tokens: int,
    ) -> None:
        self.fp_bits = fp_bits
        self.fp_radius = fp_radius
        self.text_features = text_features
        self.categorical_features = categorical_features
        self.max_instruction_tokens = max_instruction_tokens
        self.fp_cache: dict[str, tuple[int, ...]] = {}
        self.desc_cache: dict[str, dict[str, float]] = {}
        self.num_features = 9
        self.total_features = (3 * fp_bits) + text_features + categorical_features + self.num_features

    @property
    def text_offset(self) -> int:
        return 3 * self.fp_bits

    @property
    def categorical_offset(self) -> int:
        return self.text_offset + self.text_features

    @property
    def numeric_offset(self) -> int:
        return self.categorical_offset + self.categorical_features

    def fp(self, smiles: str) -> tuple[int, ...]:
        smiles = clean_text(smiles)
        if smiles in self.fp_cache:
            return self.fp_cache[smiles]
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        if mol is None:
            bits: tuple[int, ...] = ()
        else:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, self.fp_radius, nBits=self.fp_bits)
            bits = tuple(int(bit) for bit in fp.GetOnBits())
        self.fp_cache[smiles] = bits
        return bits

    def descriptors(self, smiles: str) -> dict[str, float]:
        smiles = clean_text(smiles)
        if smiles in self.desc_cache:
            return self.desc_cache[smiles]
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        if mol is None:
            desc = {"valid": 0.0, "mw": 0.0, "logp": 0.0, "qed": 0.0, "heavy_atoms": 0.0}
        else:
            desc = {
                "valid": 1.0,
                "mw": float(Descriptors.MolWt(mol)),
                "logp": float(Crippen.MolLogP(mol)),
                "qed": float(QED.qed(mol)),
                "heavy_atoms": float(mol.GetNumHeavyAtoms()),
            }
        self.desc_cache[smiles] = desc
        return desc

    def tanimoto(self, smiles_a: str, smiles_b: str) -> float:
        a = set(self.fp(smiles_a))
        b = set(self.fp(smiles_b))
        if not a and not b:
            return 0.0
        union = len(a | b)
        return float(len(a & b) / union) if union else 0.0

    def tokenize(self, text: str) -> list[str]:
        return TOKEN_RE.findall(text.lower())[: self.max_instruction_tokens]

    def add(self, cols: list[int], vals: list[float], col: int, value: float = 1.0) -> None:
        if value != 0.0:
            cols.append(col)
            vals.append(float(value))

    def row_features(self, example: CandidateExample) -> tuple[list[int], list[float]]:
        cols: list[int] = []
        vals: list[float] = []

        input_bits = set(self.fp(example.input_smiles))
        candidate_bits = set(self.fp(example.candidate_smiles))
        xor_bits = input_bits ^ candidate_bits
        for bit in candidate_bits:
            self.add(cols, vals, bit)
        for bit in input_bits:
            self.add(cols, vals, self.fp_bits + bit)
        for bit in xor_bits:
            self.add(cols, vals, (2 * self.fp_bits) + bit)

        for token in self.tokenize(example.instruction):
            col = self.text_offset + (stable_hash(f"tok:{token}") % self.text_features)
            self.add(cols, vals, col)

        categorical_values = {
            "source": example.source,
            "primary_endpoint": example.primary_endpoint,
            "primary_direction": example.primary_direction,
            "target_id": example.target_id,
            "measurement_type": example.measurement_type,
            "measurement_group": example.measurement_group,
        }
        for name, value in categorical_values.items():
            if value:
                col = self.categorical_offset + (stable_hash(f"{name}:{value}") % self.categorical_features)
                self.add(cols, vals, col)

        input_desc = self.descriptors(example.input_smiles)
        cand_desc = self.descriptors(example.candidate_smiles)
        sim = self.tanimoto(example.input_smiles, example.candidate_smiles)
        numeric = [
            sim,
            cand_desc["valid"],
            (cand_desc["mw"] - input_desc["mw"]) / 100.0,
            (cand_desc["logp"] - input_desc["logp"]) / 5.0,
            cand_desc["qed"] - input_desc["qed"],
            (cand_desc["heavy_atoms"] - input_desc["heavy_atoms"]) / 50.0,
            cand_desc["mw"] / 1000.0,
            cand_desc["logp"] / 10.0,
            math.log1p(max(example.evidence_count, 0.0)) / 5.0,
        ]
        for idx, value in enumerate(numeric):
            self.add(cols, vals, self.numeric_offset + idx, value)
        return cols, vals

    def matrix(self, examples: list[CandidateExample]) -> Any:
        data: list[float] = []
        indices: list[int] = []
        indptr = [0]
        for example in examples:
            cols, vals = self.row_features(example)
            indices.extend(cols)
            data.extend(vals)
            indptr.append(len(data))
        return csr_matrix((data, indices, indptr), shape=(len(examples), self.total_features), dtype=np.float32)

    def similarity_scores(self, examples: list[CandidateExample]) -> list[float]:
        return [self.tanimoto(example.input_smiles, example.candidate_smiles) for example in examples]


def unique_preserve_order(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        value = clean_text(value)
        if value and value not in seen:
            seen.add(value)
            out.append(value)
    return out


def build_edit_candidate_universe(train_rows: list[dict[str, str]]) -> tuple[dict[str, list[str]], list[str]]:
    by_endpoint: dict[str, list[str]] = defaultdict(list)
    global_candidates: list[str] = []
    for row in train_rows:
        endpoint = clean_text(row.get("primary_endpoint"))
        positives = positive_smiles_from_row(row)
        by_endpoint[endpoint].extend(positives)
        global_candidates.extend(positives)
    return (
        {endpoint: unique_preserve_order(values) for endpoint, values in by_endpoint.items()},
        unique_preserve_order(global_candidates),
    )


def sample_decoys(
    row: dict[str, str],
    positives: list[str],
    by_endpoint: dict[str, list[str]],
    global_candidates: list[str],
    decoys_per_query: int,
    rng: random.Random,
) -> list[str]:
    blocked = set(positives)
    blocked.add(clean_text(row.get("input_smiles_canon")))
    endpoint = clean_text(row.get("primary_endpoint"))
    candidates = [item for item in by_endpoint.get(endpoint, []) if item not in blocked]
    if len(candidates) < decoys_per_query:
        candidates.extend(item for item in global_candidates if item not in blocked and item not in candidates)
    rng.shuffle(candidates)
    return unique_preserve_order(candidates[:decoys_per_query])


def make_edit_examples(
    rows: list[dict[str, str]],
    by_endpoint: dict[str, list[str]],
    global_candidates: list[str],
    decoys_per_query: int,
    max_queries: int,
    source_name: str,
    rng: random.Random,
) -> list[CandidateExample]:
    selected_rows = list(rows)
    if max_queries:
        rng.shuffle(selected_rows)
        selected_rows = selected_rows[:max_queries]
    examples: list[CandidateExample] = []
    for row in selected_rows:
        positives = positive_smiles_from_row(row)
        if not positives:
            continue
        query_id = clean_text(row.get("query_id"))
        common = {
            "query_id": query_id,
            "split": clean_text(row.get("split")),
            "input_smiles": clean_text(row.get("input_smiles_canon")),
            "instruction": clean_text(row.get("instruction")),
            "source": source_name,
            "primary_endpoint": clean_text(row.get("primary_endpoint")),
            "primary_direction": clean_text(row.get("primary_direction")),
        }
        for smiles in positives:
            examples.append(CandidateExample(candidate_smiles=smiles, label=1, **common))
        for smiles in sample_decoys(row, positives, by_endpoint, global_candidates, decoys_per_query, rng):
            examples.append(CandidateExample(candidate_smiles=smiles, label=0, **common))
    return examples


def make_bindingdb_examples(rows: list[dict[str, str]]) -> list[CandidateExample]:
    examples: list[CandidateExample] = []
    for row in rows:
        common = {
            "query_id": clean_text(row.get("sample_id")),
            "split": "",
            "input_smiles": clean_text(row.get("input_smiles")),
            "instruction": clean_text(row.get("instruction")),
            "source": "BindingDB",
            "primary_endpoint": clean_text(row.get("measurement_type")),
            "primary_direction": "increase",
            "target_id": clean_text(row.get("target_id")),
            "target_name": clean_text(row.get("target_name")),
            "measurement_type": clean_text(row.get("measurement_type")),
            "measurement_group": clean_text(row.get("measurement_group")),
        }
        examples.append(
            CandidateExample(
                candidate_smiles=clean_text(row.get("positive_smiles")),
                label=1,
                evidence_count=float(row.get("positive_evidence_count") or 0.0),
                **common,
            )
        )
        examples.append(
            CandidateExample(
                candidate_smiles=clean_text(row.get("negative_smiles")),
                label=0,
                evidence_count=float(row.get("negative_evidence_count") or 0.0),
                **common,
            )
        )
    return examples


def fit_classifier(train_x: Any, train_y: Any, args: argparse.Namespace) -> Any:
    if len(set(int(value) for value in train_y.tolist())) < 2:
        raise SystemExit("ERROR: training data must contain both positive and negative candidates.")
    clf = SGDClassifier(
        loss="log_loss",
        penalty="l2",
        alpha=args.alpha,
        max_iter=args.max_iter,
        class_weight="balanced",
        random_state=args.seed,
        n_jobs=-1,
    )
    clf.fit(train_x, train_y)
    return clf


def model_scores(model: Any, matrix: Any) -> list[float]:
    if hasattr(model, "decision_function"):
        return [float(value) for value in model.decision_function(matrix)]
    return [float(value) for value in model.predict_proba(matrix)[:, 1]]


def evaluate_ranked_candidates(examples: list[CandidateExample], scores: list[float], top_k: int = 2) -> dict[str, Any]:
    grouped: dict[str, list[tuple[CandidateExample, float]]] = defaultdict(list)
    for example, score in zip(examples, scores):
        grouped[example.query_id].append((example, float(score)))

    hit_at_1 = []
    recall_at_k = []
    all_gold_at_k = []
    mrr_values = []
    candidate_counts = []
    positive_counts = []
    for rows in grouped.values():
        rows.sort(key=lambda item: item[1], reverse=True)
        labels = [item[0].label for item in rows]
        positives = sum(labels)
        if positives == 0:
            continue
        candidate_counts.append(len(rows))
        positive_counts.append(positives)
        top_labels = labels[:top_k]
        hit_at_1.append(1.0 if labels[0] == 1 else 0.0)
        recall_at_k.append(sum(top_labels) / positives)
        all_gold_at_k.append(1.0 if sum(top_labels) == min(top_k, positives) else 0.0)
        rank = next((idx + 1 for idx, label in enumerate(labels) if label == 1), None)
        mrr_values.append(1.0 / rank if rank else 0.0)

    return {
        "queries": len(hit_at_1),
        "mean_candidates_per_query": round(sum(candidate_counts) / len(candidate_counts), 4) if candidate_counts else 0.0,
        "mean_positives_per_query": round(sum(positive_counts) / len(positive_counts), 4) if positive_counts else 0.0,
        "hit_at_1": mean(hit_at_1),
        f"recall_at_{top_k}": mean(recall_at_k),
        f"all_gold_at_{top_k}": mean(all_gold_at_k),
        "mrr": mean(mrr_values),
    }


def evaluate_pairwise_candidates(examples: list[CandidateExample], scores: list[float]) -> dict[str, Any]:
    grouped: dict[str, dict[int, float]] = defaultdict(dict)
    labels = []
    raw_scores = []
    for example, score in zip(examples, scores):
        grouped[example.query_id][example.label] = float(score)
        labels.append(example.label)
        raw_scores.append(float(score))

    accuracy = []
    for pair in grouped.values():
        if 1 not in pair or 0 not in pair:
            continue
        if pair[1] > pair[0]:
            accuracy.append(1.0)
        elif pair[1] < pair[0]:
            accuracy.append(0.0)
        else:
            accuracy.append(0.5)

    metrics = {
        "pairs": len(accuracy),
        "pairwise_accuracy": mean(accuracy),
    }
    if len(set(labels)) == 2 and roc_auc_score is not None and average_precision_score is not None:
        metrics["candidate_auc"] = round(float(roc_auc_score(labels, raw_scores)), 6)
        metrics["candidate_average_precision"] = round(float(average_precision_score(labels, raw_scores)), 6)
    return metrics


def mean(values: list[float]) -> float:
    return round(float(sum(values) / len(values)), 6) if values else 0.0


def labels_array(examples: list[CandidateExample]) -> Any:
    return np.asarray([example.label for example in examples], dtype=np.int8)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def output_dir(base_dir: Path, task_name: str) -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return base_dir / f"{task_name}_{stamp}"


def run_edit_ranker(args: argparse.Namespace) -> int:
    require_dependencies()
    rng = random.Random(args.seed)
    data_dir = args.data_dir
    source_name = infer_source_name(data_dir, args.source_name)
    train_rows = read_csv(data_dir / "train.csv", args.max_train_queries)
    by_endpoint, global_candidates = build_edit_candidate_universe(train_rows)
    train_examples = make_edit_examples(
        train_rows,
        by_endpoint,
        global_candidates,
        args.decoys_per_query,
        0,
        source_name,
        rng,
    )
    featurizer = MoleculeFeaturizer(
        args.fp_bits,
        args.fp_radius,
        args.text_features,
        args.categorical_features,
        args.max_instruction_tokens,
    )
    train_x = featurizer.matrix(train_examples)
    train_y = labels_array(train_examples)
    model = fit_classifier(train_x, train_y, args)

    metrics: dict[str, Any] = {
        "task": "edit_candidate_ranker",
        "data_dir": str(data_dir),
        "source_name": source_name,
        "train_candidates": len(train_examples),
        "train_queries": len({example.query_id for example in train_examples}),
        "decoys_per_query": args.decoys_per_query,
        "splits": {},
    }
    for split in args.eval_splits:
        split_rows = read_csv(data_dir / f"{split}.csv", args.max_eval_queries)
        eval_examples = make_edit_examples(
            split_rows,
            by_endpoint,
            global_candidates,
            args.decoys_per_query,
            0,
            source_name,
            rng,
        )
        eval_x = featurizer.matrix(eval_examples)
        learned_scores = model_scores(model, eval_x)
        similarity_scores = featurizer.similarity_scores(eval_examples)
        metrics["splits"][split] = {
            "candidate_count": len(eval_examples),
            "query_count": len({example.query_id for example in eval_examples}),
            "similarity_baseline": evaluate_ranked_candidates(eval_examples, similarity_scores, top_k=2),
            "supervised_sgd_ranker": evaluate_ranked_candidates(eval_examples, learned_scores, top_k=2),
        }

    out_dir = output_dir(args.out_dir, f"edit_ranker_{source_name.lower()}")
    write_json(out_dir / "metrics.json", metrics)
    write_json(out_dir / "config.json", serializable_args(args))
    if args.save_model and dump is not None:
        dump({"model": model, "featurizer_config": featurizer_config(args)}, out_dir / "model.joblib")
    print(json.dumps({"out_dir": str(out_dir), "metrics": metrics}, indent=2, ensure_ascii=False))
    return 0


def run_bindingdb_ranker(args: argparse.Namespace) -> int:
    require_dependencies()
    train_rows = read_csv(args.data_dir / "train.csv", args.max_train_rows)
    train_examples = make_bindingdb_examples(train_rows)
    featurizer = MoleculeFeaturizer(
        args.fp_bits,
        args.fp_radius,
        args.text_features,
        args.categorical_features,
        args.max_instruction_tokens,
    )
    train_x = featurizer.matrix(train_examples)
    train_y = labels_array(train_examples)
    model = fit_classifier(train_x, train_y, args)

    metrics: dict[str, Any] = {
        "task": "bindingdb_pairwise_ranker",
        "data_dir": str(args.data_dir),
        "train_candidates": len(train_examples),
        "train_pairs": len(train_examples) // 2,
        "splits": {},
    }
    for split in args.eval_splits:
        rows = read_csv(args.data_dir / f"{split}.csv", args.max_eval_rows)
        eval_examples = make_bindingdb_examples(rows)
        eval_x = featurizer.matrix(eval_examples)
        learned_scores = model_scores(model, eval_x)
        similarity_scores = featurizer.similarity_scores(eval_examples)
        evidence_scores = [example.evidence_count for example in eval_examples]
        metrics["splits"][split] = {
            "candidate_count": len(eval_examples),
            "pair_count": len(eval_examples) // 2,
            "similarity_baseline": evaluate_pairwise_candidates(eval_examples, similarity_scores),
            "evidence_baseline": evaluate_pairwise_candidates(eval_examples, evidence_scores),
            "supervised_sgd_ranker": evaluate_pairwise_candidates(eval_examples, learned_scores),
        }

    out_dir = output_dir(args.out_dir, "bindingdb_ranker")
    write_json(out_dir / "metrics.json", metrics)
    write_json(out_dir / "config.json", serializable_args(args))
    if args.save_model and dump is not None:
        dump({"model": model, "featurizer_config": featurizer_config(args)}, out_dir / "model.joblib")
    print(json.dumps({"out_dir": str(out_dir), "metrics": metrics}, indent=2, ensure_ascii=False))
    return 0


def serializable_args(args: argparse.Namespace) -> dict[str, Any]:
    out = {}
    for key, value in vars(args).items():
        out[key] = str(value) if isinstance(value, Path) else value
    return out


def featurizer_config(args: argparse.Namespace) -> dict[str, Any]:
    return {
        "fp_bits": args.fp_bits,
        "fp_radius": args.fp_radius,
        "text_features": args.text_features,
        "categorical_features": args.categorical_features,
        "max_instruction_tokens": args.max_instruction_tokens,
    }


def main() -> int:
    args = parse_args()
    if args.command == "edit-ranker":
        return run_edit_ranker(args)
    if args.command == "bindingdb-ranker":
        return run_bindingdb_ranker(args)
    raise SystemExit(f"unknown command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())

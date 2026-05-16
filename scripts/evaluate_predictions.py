from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

try:
    from rdkit import Chem  # type: ignore
except Exception:  # noqa: BLE001 - RDKit is optional for this script.
    Chem = None


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
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
            rows.append(row)
    return rows


def normalize_smiles(smiles: str) -> str | None:
    smiles = smiles.strip()
    if not smiles:
        return None
    if Chem is None:
        return smiles
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def is_valid_smiles(smiles: str) -> bool:
    if not smiles or not isinstance(smiles, str):
        return False
    if Chem is None:
        return any(ch.isalpha() for ch in smiles)
    return Chem.MolFromSmiles(smiles) is not None


def dcg(relevances: list[float]) -> float:
    return sum(rel / math.log2(rank + 2) for rank, rel in enumerate(relevances))


def ece_score(confidence_outcomes: list[tuple[float, int]], bins: int = 10) -> float | None:
    if not confidence_outcomes:
        return None
    total = len(confidence_outcomes)
    ece = 0.0
    for bin_idx in range(bins):
        lo = bin_idx / bins
        hi = (bin_idx + 1) / bins
        if bin_idx == bins - 1:
            bucket = [(conf, y) for conf, y in confidence_outcomes if lo <= conf <= hi]
        else:
            bucket = [(conf, y) for conf, y in confidence_outcomes if lo <= conf < hi]
        if not bucket:
            continue
        avg_conf = sum(conf for conf, _ in bucket) / len(bucket)
        avg_acc = sum(y for _, y in bucket) / len(bucket)
        ece += (len(bucket) / total) * abs(avg_acc - avg_conf)
    return ece


def evaluate(
    queries: list[dict[str, Any]],
    answers: list[dict[str, Any]],
    predictions: list[dict[str, Any]],
    k: int,
) -> dict[str, Any]:
    query_by_id = {query["query_id"]: query for query in queries}
    accepted: dict[str, dict[str, dict[str, Any]]] = defaultdict(dict)
    answer_confidences: dict[str, list[float]] = defaultdict(list)

    for answer in answers:
        norm = normalize_smiles(answer["target_smiles_canon"])
        if norm is None:
            continue
        accepted[answer["query_id"]][norm] = answer
        answer_confidences[answer["query_id"]].append(float(answer.get("confidence", 1.0)))

    pred_by_query = {pred["query_id"]: pred for pred in predictions}
    missing_predictions = sorted(set(query_by_id) - set(pred_by_query))
    extra_predictions = sorted(set(pred_by_query) - set(query_by_id))

    success_values: list[int] = []
    experimental_values: list[int] = []
    scaffold_values: list[int] = []
    diversity_values: list[float] = []
    mrr_values: list[float] = []
    ndcg_values: list[float] = []
    validity_flags: list[int] = []
    all_topk_norms: list[str] = []
    confidence_outcomes: list[tuple[float, int]] = []

    for query_id in sorted(query_by_id):
        pred = pred_by_query.get(query_id, {"ranked_smiles": []})
        ranked = pred.get("ranked_smiles", [])[:k]
        if not isinstance(ranked, list):
            ranked = []

        topk_norms: list[str | None] = []
        matched_answers: list[dict[str, Any] | None] = []
        relevance: list[float] = []

        for item in ranked:
            if isinstance(item, str):
                smiles = item
                confidence = None
            else:
                smiles = str(item.get("target_smiles_canon", ""))
                confidence = item.get("confidence")

            validity_flags.append(1 if is_valid_smiles(smiles) else 0)
            norm = normalize_smiles(smiles)
            topk_norms.append(norm)
            if norm is not None:
                all_topk_norms.append(norm)

            match = accepted.get(query_id, {}).get(norm or "")
            matched_answers.append(match)
            is_match = int(match is not None and match.get("constraint_flags", {}).get("success", True) is True)
            relevance.append(float(match.get("confidence", 1.0)) if is_match and match is not None else 0.0)

            if isinstance(confidence, (int, float)):
                confidence_outcomes.append((float(confidence), is_match))

        success_values.append(1 if any(rel > 0 for rel in relevance) else 0)
        experimental_values.append(
            1
            if any(
                match is not None
                and (match.get("experimental_only_flag") is True or match.get("label_type") == "experimental")
                for match in matched_answers
            )
            else 0
        )
        scaffold_values.append(
            1
            if any(match is not None and match.get("constraint_flags", {}).get("scaffold_retained") is True for match in matched_answers)
            else 0
        )

        unique_topk = {norm for norm in topk_norms if norm is not None}
        diversity_values.append(len(unique_topk) / len(topk_norms) if topk_norms else 0.0)

        rr = 0.0
        for idx, rel in enumerate(relevance, 1):
            if rel > 0:
                rr = 1.0 / idx
                break
        mrr_values.append(rr)

        ideal = sorted(answer_confidences.get(query_id, []), reverse=True)[:k]
        ideal_dcg = dcg(ideal)
        ndcg_values.append(dcg(relevance) / ideal_dcg if ideal_dcg > 0 else 0.0)

    def avg(values: list[float | int]) -> float:
        return float(sum(values) / len(values)) if values else 0.0

    metrics = {
        "k": k,
        "num_queries": len(queries),
        "num_answers": len(answers),
        "num_predictions": len(predictions),
        "validity": avg(validity_flags),
        "validity_mode": "rdkit" if Chem is not None else "syntax_only",
        "uniqueness_at_k_global": len(set(all_topk_norms)) / len(all_topk_norms) if all_topk_norms else 0.0,
        "query_success_at_k": avg(success_values),
        "experimental_fidelity_at_k": avg(experimental_values),
        "scaffold_retention_at_k": avg(scaffold_values),
        "diversity_at_k": avg(diversity_values),
        "mrr_at_k": avg(mrr_values),
        "ndcg_at_k": avg(ndcg_values),
        "ece": ece_score(confidence_outcomes),
        "missing_predictions": missing_predictions,
        "extra_predictions": extra_predictions,
    }
    return metrics


def main() -> int:
    parser = argparse.ArgumentParser(description="Evaluate predictions for POS-WAY ADMET editing benchmark.")
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument("--answers", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    if args.k < 1:
        print("ERROR: --k must be >= 1", file=sys.stderr)
        return 1

    try:
        queries = load_jsonl(args.queries)
        answers = load_jsonl(args.answers)
        predictions = load_jsonl(args.predictions)
    except Exception as exc:  # noqa: BLE001 - CLI should report load failures cleanly.
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    metrics = evaluate(queries, answers, predictions, args.k)
    rendered = json.dumps(metrics, ensure_ascii=False, indent=2)
    print(rendered)

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(rendered + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


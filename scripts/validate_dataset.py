from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


QUERY_REQUIRED = {
    "query_id",
    "input_smiles_raw",
    "input_smiles_canon",
    "input_inchikey",
    "input_connectivity_key",
    "question_text",
    "question_template",
    "target_endpoints",
    "hard_constraints",
    "source_pool",
    "split",
    "num_answers",
    "schema_version",
}

ANSWER_REQUIRED = {
    "answer_id",
    "query_id",
    "target_smiles_canon",
    "target_inchikey",
    "transform_class",
    "label_type",
    "endpoint_name",
    "value_before",
    "value_after",
    "delta_value",
    "unit_canonical",
    "confidence",
    "experimental_only_flag",
    "provenance",
    "constraint_flags",
}

ALLOWED_SPLITS = {"train", "valid", "test"}
ALLOWED_SOURCE_POOLS = {"gold", "silver", "mixed"}
ALLOWED_LABEL_TYPES = {"experimental", "predicted", "proxy"}


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


def ensure_object(value: Any, row_id: str, field: str, errors: list[str]) -> None:
    if isinstance(value, dict):
        return
    if isinstance(value, str):
        try:
            parsed = json.loads(value)
        except json.JSONDecodeError:
            errors.append(f"{row_id}: {field} is neither object nor JSON string")
            return
        if isinstance(parsed, dict):
            return
    errors.append(f"{row_id}: {field} must be an object")


def validate(queries: list[dict[str, Any]], answers: list[dict[str, Any]]) -> tuple[list[str], list[str], dict[str, Any]]:
    errors: list[str] = []
    warnings: list[str] = []

    query_ids: set[str] = set()
    query_by_id: dict[str, dict[str, Any]] = {}
    split_by_input: dict[str, set[str]] = defaultdict(set)
    split_by_target: dict[str, set[str]] = defaultdict(set)
    answers_by_query: dict[str, list[dict[str, Any]]] = defaultdict(list)

    for idx, query in enumerate(queries, 1):
        row_id = str(query.get("query_id", f"query_row_{idx}"))
        missing = QUERY_REQUIRED - set(query)
        if missing:
            errors.append(f"{row_id}: missing query fields: {sorted(missing)}")

        if row_id in query_ids:
            errors.append(f"{row_id}: duplicate query_id")
        query_ids.add(row_id)
        query_by_id[row_id] = query

        if query.get("split") not in ALLOWED_SPLITS:
            errors.append(f"{row_id}: invalid split {query.get('split')!r}")
        if query.get("source_pool") not in ALLOWED_SOURCE_POOLS:
            errors.append(f"{row_id}: invalid source_pool {query.get('source_pool')!r}")
        if not isinstance(query.get("target_endpoints"), list) or not query.get("target_endpoints"):
            errors.append(f"{row_id}: target_endpoints must be a non-empty list")
        ensure_object(query.get("hard_constraints"), row_id, "hard_constraints", errors)
        if not isinstance(query.get("num_answers"), int) or query.get("num_answers", 0) < 1:
            errors.append(f"{row_id}: num_answers must be a positive integer")

        split_by_input[str(query.get("input_inchikey"))].add(str(query.get("split")))

    answer_ids: set[str] = set()
    label_counts: Counter[str] = Counter()

    for idx, answer in enumerate(answers, 1):
        row_id = str(answer.get("answer_id", f"answer_row_{idx}"))
        missing = ANSWER_REQUIRED - set(answer)
        if missing:
            errors.append(f"{row_id}: missing answer fields: {sorted(missing)}")

        if row_id in answer_ids:
            errors.append(f"{row_id}: duplicate answer_id")
        answer_ids.add(row_id)

        query_id = str(answer.get("query_id"))
        if query_id not in query_by_id:
            errors.append(f"{row_id}: unknown query_id {query_id!r}")
        else:
            answers_by_query[query_id].append(answer)
            source_pool = query_by_id[query_id].get("source_pool")
            if source_pool == "gold" and answer.get("label_type") != "experimental":
                errors.append(f"{row_id}: gold query cannot contain non-experimental answer")

        if answer.get("label_type") not in ALLOWED_LABEL_TYPES:
            errors.append(f"{row_id}: invalid label_type {answer.get('label_type')!r}")
        else:
            label_counts[str(answer["label_type"])] += 1

        confidence = answer.get("confidence")
        if not isinstance(confidence, (int, float)) or not (0 <= float(confidence) <= 1):
            errors.append(f"{row_id}: confidence must be in [0, 1]")

        if not isinstance(answer.get("experimental_only_flag"), bool):
            errors.append(f"{row_id}: experimental_only_flag must be boolean")
        ensure_object(answer.get("provenance"), row_id, "provenance", errors)
        ensure_object(answer.get("constraint_flags"), row_id, "constraint_flags", errors)

        target_key = str(answer.get("target_inchikey"))
        if query_id in query_by_id:
            split_by_target[target_key].add(str(query_by_id[query_id].get("split")))

    for query_id, query in query_by_id.items():
        actual = len(answers_by_query.get(query_id, []))
        expected = query.get("num_answers")
        if actual != expected:
            errors.append(f"{query_id}: num_answers={expected}, actual answers={actual}")
        if actual < 2:
            errors.append(f"{query_id}: multi-answer benchmark requires at least 2 answers")

    for inchikey, splits in split_by_input.items():
        if len(splits) > 1:
            errors.append(f"input_inchikey leakage across splits: {inchikey} -> {sorted(splits)}")

    for inchikey, splits in split_by_target.items():
        if len(splits) > 1:
            errors.append(f"target_inchikey leakage across splits: {inchikey} -> {sorted(splits)}")

    test_answers = [
        answer
        for answer in answers
        if answer.get("query_id") in query_by_id and query_by_id[answer["query_id"]].get("split") == "test"
    ]
    if test_answers and not any(answer.get("experimental_only_flag") for answer in test_answers):
        warnings.append("test split has no experimental_only answers; this is acceptable for silver/proxy builds but not for a strict public benchmark")

    answer_counts = [len(rows) for rows in answers_by_query.values()]
    summary = {
        "num_queries": len(queries),
        "num_answers": len(answers),
        "splits": Counter(str(query.get("split")) for query in queries),
        "source_pools": Counter(str(query.get("source_pool")) for query in queries),
        "label_types": label_counts,
        "answers_per_query": {
            "min": min(answer_counts) if answer_counts else 0,
            "max": max(answer_counts) if answer_counts else 0,
            "mean": sum(answer_counts) / len(answer_counts) if answer_counts else 0,
        },
    }
    return errors, warnings, summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate POS-WAY ADMET editing dataset files.")
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument("--answers", type=Path, required=True)
    args = parser.parse_args()

    try:
        queries = load_jsonl(args.queries)
        answers = load_jsonl(args.answers)
    except Exception as exc:  # noqa: BLE001 - CLI should report all load failures cleanly.
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    errors, warnings, summary = validate(queries, answers)

    if errors:
        print("Validation failed:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 1

    print("Validation passed.")
    print(json.dumps(summary, ensure_ascii=False, indent=2, default=dict))
    if warnings:
        print("Warnings:")
        for warning in warnings:
            print(f"- {warning}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

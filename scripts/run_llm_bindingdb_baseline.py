from __future__ import annotations

import argparse
import csv
import json
import os
import random
import re
import time
import urllib.error
import urllib.request
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

try:
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import AllChem
except Exception:  # noqa: BLE001
    Chem = None
    DataStructs = None
    RDLogger = None
    AllChem = None


JSON_OBJECT_RE = re.compile(r"\{.*\}", re.DOTALL)
SMILES_LINE_RE = re.compile(r"^[A-Za-z0-9@+\-\[\]\(\)=#$\\/%.:]+$")


@dataclass
class BindingRecord:
    sample_id: str
    split: str
    input_smiles: str
    instruction: str
    target_id: str
    target_name: str
    target_organism: str
    target_sequence: str
    measurement_type: str
    measurement_group: str
    positive_smiles: str
    negative_smiles: str
    raw_row: dict[str, str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run or evaluate local/open LLM baselines for BindingDB target-conditioned ligand optimization."
        )
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    export = subparsers.add_parser("export-prompts")
    add_dataset_args(export)
    add_prompt_args(export)
    export.add_argument("--out", type=Path, default=Path("reports/llm_bindingdb_baselines/prompts.jsonl"))

    run = subparsers.add_parser("run")
    add_dataset_args(run)
    add_prompt_args(run)
    run.add_argument("--provider", choices=["openai-compatible"], required=True)
    run.add_argument("--model", required=True)
    run.add_argument("--out", type=Path, default=Path("reports/llm_bindingdb_baselines/predictions.jsonl"))
    run.add_argument("--base-url", default="http://127.0.0.1:8080/v1")
    run.add_argument("--api-key-env", default="LOCAL_LLM_API_KEY")
    run.add_argument("--temperature", type=float, default=0.0)
    run.add_argument("--max-tokens", type=int, default=160)
    run.add_argument("--timeout-seconds", type=int, default=90)
    run.add_argument("--sleep-seconds", type=float, default=0.0)
    run.add_argument("--resume", action="store_true")

    evaluate = subparsers.add_parser("evaluate")
    add_dataset_args(evaluate)
    evaluate.add_argument("--predictions", type=Path, required=True)
    evaluate.add_argument("--out", type=Path, default=Path("reports/llm_bindingdb_baselines/metrics.json"))

    return parser.parse_args()


def add_dataset_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--data-dir", type=Path, default=Path("data/bindingdb_target_conditioned"))
    parser.add_argument(
        "--split",
        default="test_seen_target",
        choices=["train", "val", "test_seen_target", "test_unseen_target"],
    )
    parser.add_argument("--max-rows", type=int, default=0)
    parser.add_argument("--seed", type=int, default=29)


def add_prompt_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--few-shot-k", type=int, default=0)
    parser.add_argument("--few-shot-train-split", default="train")
    parser.add_argument("--system-prompt", default="")
    parser.add_argument("--max-sequence-chars", type=int, default=800)


def init_rdkit() -> None:
    if RDLogger is not None:
        RDLogger.DisableLog("rdApp.warning")
        RDLogger.DisableLog("rdApp.error")


def clean_text(value: Any) -> str:
    return str(value or "").replace("\r", " ").replace("\n", " ").strip()


def read_csv(path: Path, max_rows: int = 0) -> list[dict[str, str]]:
    if not path.is_file():
        raise SystemExit(f"ERROR: file not found: {path}")
    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append(row)
            if max_rows and len(rows) >= max_rows:
                break
    return rows


def row_to_record(row: dict[str, str], split: str) -> BindingRecord:
    return BindingRecord(
        sample_id=clean_text(row.get("sample_id")),
        split=split,
        input_smiles=clean_text(row.get("input_smiles")),
        instruction=clean_text(row.get("instruction")),
        target_id=clean_text(row.get("target_id")),
        target_name=clean_text(row.get("target_name")),
        target_organism=clean_text(row.get("target_organism")),
        target_sequence=clean_text(row.get("target_sequence")),
        measurement_type=clean_text(row.get("measurement_type")),
        measurement_group=clean_text(row.get("measurement_group")),
        positive_smiles=clean_text(row.get("positive_smiles")),
        negative_smiles=clean_text(row.get("negative_smiles")),
        raw_row=row,
    )


def load_records(data_dir: Path, split: str, max_rows: int, seed: int) -> list[BindingRecord]:
    rows = read_csv(data_dir / f"{split}.csv")
    records = [row_to_record(row, split) for row in rows if row.get("sample_id") and row.get("positive_smiles")]
    if max_rows:
        rng = random.Random(seed)
        rng.shuffle(records)
        records = records[:max_rows]
    records.sort(key=lambda item: item.sample_id)
    return records


def default_system_prompt() -> str:
    return (
        "You are a careful medicinal chemistry assistant for target-conditioned ligand optimization. "
        "Return only valid JSON. Do not explain. Do not include markdown."
    )


def make_user_prompt(record: BindingRecord, max_sequence_chars: int) -> str:
    sequence = record.target_sequence
    if max_sequence_chars and len(sequence) > max_sequence_chars:
        sequence = sequence[:max_sequence_chars] + "...[truncated]"
    return (
        "Task: propose one structurally related ligand with stronger measured activity for the target.\n"
        "Return exactly this JSON shape:\n"
        '{"stronger_smiles":"SMILES"}\n\n'
        f"Input SMILES: {record.input_smiles}\n"
        f"Target: {record.target_name}\n"
        f"Organism: {record.target_organism}\n"
        f"Measurement: {record.measurement_type} ({record.measurement_group}); higher p-scale means stronger.\n"
        f"Instruction: {record.instruction}\n"
        f"Target sequence prefix: {sequence}\n\n"
        "Rules:\n"
        "- Return one valid SMILES string.\n"
        "- Do not return the unchanged input molecule.\n"
        "- Prefer a local, chemically plausible edit.\n"
        "- Output JSON only."
    )


def select_few_shot_examples(
    target: BindingRecord,
    train_records: list[BindingRecord],
    k: int,
    rng: random.Random,
) -> list[BindingRecord]:
    if k <= 0:
        return []
    same_target_metric = [
        item
        for item in train_records
        if item.target_id == target.target_id
        and item.measurement_type == target.measurement_type
        and item.sample_id != target.sample_id
    ]
    same_metric = [
        item
        for item in train_records
        if item.measurement_type == target.measurement_type and item.sample_id != target.sample_id
    ]
    fallback = [item for item in train_records if item.sample_id != target.sample_id]
    rng.shuffle(same_target_metric)
    rng.shuffle(same_metric)
    rng.shuffle(fallback)
    selected: list[BindingRecord] = []
    seen: set[str] = set()
    for item in same_target_metric + same_metric + fallback:
        if item.sample_id in seen:
            continue
        selected.append(item)
        seen.add(item.sample_id)
        if len(selected) >= k:
            break
    return selected


def build_messages(
    record: BindingRecord,
    train_records: list[BindingRecord],
    few_shot_k: int,
    seed: int,
    system_prompt: str,
    max_sequence_chars: int,
) -> list[dict[str, str]]:
    rng = random.Random(f"{seed}:{record.sample_id}")
    examples = select_few_shot_examples(record, train_records, few_shot_k, rng)
    messages = [{"role": "system", "content": system_prompt or default_system_prompt()}]
    for example in examples:
        messages.append({"role": "user", "content": make_user_prompt(example, max_sequence_chars)})
        messages.append(
            {
                "role": "assistant",
                "content": json.dumps(
                    {"stronger_smiles": example.positive_smiles},
                    separators=(",", ":"),
                    ensure_ascii=False,
                ),
            }
        )
    messages.append({"role": "user", "content": make_user_prompt(record, max_sequence_chars)})
    return messages


def prompt_record(
    record: BindingRecord,
    train_records: list[BindingRecord],
    args: argparse.Namespace,
) -> dict[str, Any]:
    return {
        "sample_id": record.sample_id,
        "split": record.split,
        "input_smiles": record.input_smiles,
        "target_id": record.target_id,
        "target_name": record.target_name,
        "measurement_type": record.measurement_type,
        "gold_smiles": record.positive_smiles,
        "negative_smiles": record.negative_smiles,
        "messages": build_messages(
            record,
            train_records,
            args.few_shot_k,
            args.seed,
            args.system_prompt,
            args.max_sequence_chars,
        ),
    }


def write_jsonl(path: Path, rows: Iterable[dict[str, Any]], append: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    mode = "a" if append else "w"
    with path.open(mode, encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False, separators=(",", ":")) + "\n")


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        return []
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def export_prompts(args: argparse.Namespace) -> int:
    records = load_records(args.data_dir, args.split, args.max_rows, args.seed)
    train_records = load_records(args.data_dir, args.few_shot_train_split, 0, args.seed) if args.few_shot_k else []
    rows = [prompt_record(record, train_records, args) for record in records]
    write_jsonl(args.out, rows)
    print(json.dumps({"out": str(args.out), "records": len(rows)}, indent=2))
    return 0


def call_openai_compatible(messages: list[dict[str, str]], args: argparse.Namespace) -> str:
    api_key = os.environ.get(args.api_key_env, "")
    if not api_key:
        raise SystemExit(f"ERROR: missing API key environment variable: {args.api_key_env}")
    url = args.base_url.rstrip("/") + "/chat/completions"
    payload = {
        "model": args.model,
        "messages": messages,
        "temperature": args.temperature,
        "max_tokens": args.max_tokens,
    }
    request = urllib.request.Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json", "Authorization": f"Bearer {api_key}"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=args.timeout_seconds) as response:
            data = json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"HTTP {exc.code}: {body}") from exc
    return clean_text(data["choices"][0]["message"]["content"])


def existing_sample_ids(path: Path) -> set[str]:
    return {clean_text(row.get("sample_id")) for row in read_jsonl(path)}


def run_model(args: argparse.Namespace) -> int:
    records = load_records(args.data_dir, args.split, args.max_rows, args.seed)
    train_records = load_records(args.data_dir, args.few_shot_train_split, 0, args.seed) if args.few_shot_k else []
    done = existing_sample_ids(args.out) if args.resume else set()
    wrote = 0
    for record in records:
        if record.sample_id in done:
            continue
        messages = build_messages(
            record,
            train_records,
            args.few_shot_k,
            args.seed,
            args.system_prompt,
            args.max_sequence_chars,
        )
        started = time.time()
        output = ""
        error = ""
        try:
            output = call_openai_compatible(messages, args)
        except Exception as exc:  # noqa: BLE001
            error = str(exc)
        row = {
            "sample_id": record.sample_id,
            "split": record.split,
            "provider": args.provider,
            "model": args.model,
            "few_shot_k": args.few_shot_k,
            "input_smiles": record.input_smiles,
            "target_id": record.target_id,
            "target_name": record.target_name,
            "measurement_type": record.measurement_type,
            "gold_smiles": record.positive_smiles,
            "negative_smiles": record.negative_smiles,
            "raw_output": output,
            "parsed_smiles": parse_smiles_from_output(output),
            "error": error,
            "latency_seconds": round(time.time() - started, 3),
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
        }
        write_jsonl(args.out, [row], append=True)
        wrote += 1
        if args.sleep_seconds:
            time.sleep(args.sleep_seconds)
    print(json.dumps({"out": str(args.out), "written": wrote, "skipped_existing": len(done)}, indent=2))
    return 0


def parse_smiles_from_output(output: str) -> str:
    output = clean_text(output)
    if not output:
        return ""
    match = JSON_OBJECT_RE.search(output)
    if match:
        try:
            payload = json.loads(match.group(0))
        except json.JSONDecodeError:
            payload = {}
        if isinstance(payload, dict):
            for key in ["stronger_smiles", "edited_smiles", "smiles", "molecule", "prediction"]:
                value = payload.get(key)
                if isinstance(value, list) and value:
                    return clean_text(value[0])
                if isinstance(value, str):
                    return clean_text(value)
    for raw_line in output.replace(",", "\n").splitlines():
        line = raw_line.strip().strip("-*`'\" ")
        if len(line) >= 2 and SMILES_LINE_RE.match(line):
            return line
    tokens = re.findall(r"[A-Za-z0-9@+\-\[\]\(\)=#$\\/%.:]{3,}", output)
    return clean_text(tokens[0]) if tokens else ""


def canonical_smiles(smiles: str) -> str | None:
    smiles = clean_text(smiles)
    if not smiles:
        return None
    if Chem is None:
        return smiles
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def mol_fp(smiles: str) -> Any | None:
    if Chem is None or AllChem is None:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def tanimoto(smiles_a: str, smiles_b: str) -> float | None:
    if DataStructs is None:
        return None
    fp_a = mol_fp(smiles_a)
    fp_b = mol_fp(smiles_b)
    if fp_a is None or fp_b is None:
        return None
    return float(DataStructs.FingerprintSimilarity(fp_a, fp_b))


def mean(values: list[float]) -> float:
    return round(sum(values) / len(values), 6) if values else 0.0


def evaluate_predictions(args: argparse.Namespace) -> int:
    init_rdkit()
    records = {record.sample_id: record for record in load_records(args.data_dir, args.split, args.max_rows, args.seed)}
    predictions = read_jsonl(args.predictions)
    counters: Counter[str] = Counter()
    exact_hit: list[float] = []
    negative_exact: list[float] = []
    valid_rate: list[float] = []
    unchanged_rate: list[float] = []
    tanimoto_to_positive: list[float] = []
    tanimoto_to_input: list[float] = []
    per_sample: list[dict[str, Any]] = []

    for pred in predictions:
        sample_id = clean_text(pred.get("sample_id"))
        record = records.get(sample_id)
        if record is None:
            counters["prediction_without_matching_sample"] += 1
            continue
        parsed = clean_text(pred.get("parsed_smiles")) or parse_smiles_from_output(clean_text(pred.get("raw_output")))
        pred_canon = canonical_smiles(parsed)
        gold_canon = canonical_smiles(record.positive_smiles)
        negative_canon = canonical_smiles(record.negative_smiles)
        input_canon = canonical_smiles(record.input_smiles)

        valid_rate.append(1.0 if pred_canon else 0.0)
        exact_hit.append(1.0 if pred_canon and gold_canon and pred_canon == gold_canon else 0.0)
        negative_exact.append(1.0 if pred_canon and negative_canon and pred_canon == negative_canon else 0.0)
        unchanged_rate.append(1.0 if pred_canon and input_canon and pred_canon == input_canon else 0.0)

        if pred_canon and gold_canon:
            sim = tanimoto(pred_canon, gold_canon)
            if sim is not None:
                tanimoto_to_positive.append(sim)
        if pred_canon and input_canon:
            sim = tanimoto(pred_canon, input_canon)
            if sim is not None:
                tanimoto_to_input.append(sim)
        if pred.get("error"):
            counters["provider_error"] += 1
        if not parsed:
            counters["empty_prediction"] += 1

        per_sample.append(
            {
                "sample_id": sample_id,
                "parsed_smiles": parsed,
                "gold_smiles": record.positive_smiles,
                "negative_smiles": record.negative_smiles,
                "exact_hit": bool(exact_hit[-1]),
            }
        )

    metrics = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "data_dir": str(args.data_dir),
        "split": args.split,
        "predictions": str(args.predictions),
        "rdkit_available": Chem is not None,
        "num_samples_in_split": len(records),
        "num_predictions": len(predictions),
        "num_evaluated": len(per_sample),
        "coverage": round(len(per_sample) / len(records), 6) if records else 0.0,
        "valid_smiles_rate": mean(valid_rate),
        "exact_positive_hit_rate": mean(exact_hit),
        "exact_negative_return_rate": mean(negative_exact),
        "unchanged_input_rate": mean(unchanged_rate),
        "mean_tanimoto_to_positive": mean(tanimoto_to_positive),
        "mean_tanimoto_to_input": mean(tanimoto_to_input),
        "counters": dict(sorted(counters.items())),
        "per_sample": per_sample,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(metrics, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"out": str(args.out), "summary": {k: v for k, v in metrics.items() if k != "per_sample"}}, indent=2))
    return 0


def main() -> int:
    init_rdkit()
    args = parse_args()
    if args.command == "export-prompts":
        return export_prompts(args)
    if args.command == "run":
        return run_model(args)
    if args.command == "evaluate":
        return evaluate_predictions(args)
    raise SystemExit(f"unknown command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())

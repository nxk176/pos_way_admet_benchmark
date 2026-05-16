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
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

try:
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import AllChem, Crippen, Descriptors, QED
except Exception:  # noqa: BLE001 - evaluation can still export prompts without RDKit.
    Chem = None
    DataStructs = None
    RDLogger = None
    AllChem = None
    Crippen = None
    Descriptors = None
    QED = None


JSON_OBJECT_RE = re.compile(r"\{.*\}", re.DOTALL)
JSON_ARRAY_RE = re.compile(r"\[.*\]", re.DOTALL)
SMILES_LINE_RE = re.compile(r"^[A-Za-z0-9@+\-\[\]\(\)=#$\\/%.:]+$")


@dataclass
class QueryRecord:
    query_id: str
    split: str
    input_smiles: str
    instruction: str
    primary_endpoint: str
    primary_direction: str
    gold_smiles: list[str]
    raw_row: dict[str, str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run or evaluate LLM direct-edit baselines: input molecule + instruction -> edited SMILES. "
            "This is separate from traditional candidate-ranker baselines."
        )
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    export = subparsers.add_parser("export-prompts", help="Export prompt records without calling any model.")
    add_dataset_args(export)
    add_prompt_args(export)
    export.add_argument("--out", type=Path, default=Path("reports/llm_edit_baselines/prompts.jsonl"))

    run = subparsers.add_parser("run", help="Call an LLM provider and write predictions JSONL.")
    add_dataset_args(run)
    add_prompt_args(run)
    run.add_argument("--provider", choices=["openai-compatible", "hf-local"], required=True)
    run.add_argument("--model", required=True)
    run.add_argument("--out", type=Path, default=Path("reports/llm_edit_baselines/predictions.jsonl"))
    run.add_argument("--base-url", default="https://api.openai.com/v1")
    run.add_argument("--api-key-env", default="OPENAI_API_KEY")
    run.add_argument("--temperature", type=float, default=0.0)
    run.add_argument("--max-tokens", type=int, default=256)
    run.add_argument("--timeout-seconds", type=int, default=90)
    run.add_argument("--sleep-seconds", type=float, default=0.0)
    run.add_argument("--resume", action="store_true", help="Skip query IDs already present in --out.")
    run.add_argument("--trust-remote-code", action="store_true", help="For local HF/Qwen-style models.")
    run.add_argument("--device-map", default="auto")
    run.add_argument("--torch-dtype", default="auto")

    evaluate = subparsers.add_parser("evaluate", help="Evaluate a predictions JSONL file.")
    add_dataset_args(evaluate)
    evaluate.add_argument("--predictions", type=Path, required=True)
    evaluate.add_argument("--out", type=Path, default=Path("reports/llm_edit_baselines/metrics.json"))
    evaluate.add_argument("--top-k", type=int, default=2)
    evaluate.add_argument("--mw-delta", type=float, default=100.0)
    evaluate.add_argument("--logp-delta", type=float, default=2.0)
    evaluate.add_argument("--qed-drop", type=float, default=0.1)

    return parser.parse_args()


def add_dataset_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--data-dir", type=Path, default=Path("data/chembl_3prop_2pos"))
    parser.add_argument("--split", default="test", choices=["train", "val", "test"])
    parser.add_argument("--max-rows", type=int, default=0)
    parser.add_argument("--seed", type=int, default=29)


def add_prompt_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--few-shot-k", type=int, default=0)
    parser.add_argument("--few-shot-train-split", default="train")
    parser.add_argument("--system-prompt", default="")


def init_rdkit() -> None:
    if RDLogger is not None:
        RDLogger.DisableLog("rdApp.warning")
        RDLogger.DisableLog("rdApp.error")


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


def safe_json(raw: Any, fallback: Any) -> Any:
    try:
        return json.loads(raw)
    except (TypeError, json.JSONDecodeError):
        return fallback


def clean_text(value: Any) -> str:
    return str(value or "").replace("\r", " ").replace("\n", " ").strip()


def positive_smiles_from_row(row: dict[str, str]) -> list[str]:
    values = safe_json(row.get("positive_answer_smiles_json"), [])
    if not isinstance(values, list):
        return []
    return [clean_text(value) for value in values if clean_text(value)]


def row_to_record(row: dict[str, str]) -> QueryRecord:
    return QueryRecord(
        query_id=clean_text(row.get("query_id")),
        split=clean_text(row.get("split")),
        input_smiles=clean_text(row.get("input_smiles_canon")),
        instruction=clean_text(row.get("instruction")),
        primary_endpoint=clean_text(row.get("primary_endpoint")),
        primary_direction=clean_text(row.get("primary_direction")),
        gold_smiles=positive_smiles_from_row(row),
        raw_row=row,
    )


def load_records(data_dir: Path, split: str, max_rows: int, seed: int) -> list[QueryRecord]:
    rows = read_csv(data_dir / f"{split}.csv")
    records = [row_to_record(row) for row in rows if row.get("query_id") and positive_smiles_from_row(row)]
    if max_rows:
        rng = random.Random(seed)
        rng.shuffle(records)
        records = records[:max_rows]
    records.sort(key=lambda item: item.query_id)
    return records


def default_system_prompt() -> str:
    return (
        "You are a careful medicinal chemistry molecule editor. "
        "Given one input SMILES and one editing instruction, propose edited molecules that satisfy the instruction. "
        "Return only valid JSON. Do not explain. Do not include markdown."
    )


def make_user_prompt(record: QueryRecord) -> str:
    return (
        "Task: edit the input molecule according to the instruction.\n"
        "Return exactly this JSON shape:\n"
        '{"edited_smiles":["SMILES_1","SMILES_2"]}\n\n'
        f"Input SMILES: {record.input_smiles}\n"
        f"Instruction: {record.instruction}\n\n"
        "Rules:\n"
        "- Return exactly 2 edited SMILES strings.\n"
        "- Do not return the unchanged input molecule.\n"
        "- Prefer local, chemically plausible edits.\n"
        "- Keep the requested preserved property and local constraints in mind.\n"
        "- Output JSON only."
    )


def few_shot_messages(record: QueryRecord, examples: list[QueryRecord]) -> list[dict[str, str]]:
    messages: list[dict[str, str]] = []
    for example in examples:
        messages.append({"role": "user", "content": make_user_prompt(example)})
        messages.append(
            {
                "role": "assistant",
                "content": json.dumps({"edited_smiles": example.gold_smiles[:2]}, separators=(",", ":")),
            }
        )
    messages.append({"role": "user", "content": make_user_prompt(record)})
    return messages


def select_few_shot_examples(
    target: QueryRecord,
    train_records: list[QueryRecord],
    k: int,
    rng: random.Random,
) -> list[QueryRecord]:
    if k <= 0:
        return []
    same_endpoint = [
        item for item in train_records if item.primary_endpoint == target.primary_endpoint and item.query_id != target.query_id
    ]
    fallback = [item for item in train_records if item.query_id != target.query_id]
    rng.shuffle(same_endpoint)
    rng.shuffle(fallback)
    selected: list[QueryRecord] = []
    seen: set[str] = set()
    for item in same_endpoint + fallback:
        if item.query_id in seen:
            continue
        selected.append(item)
        seen.add(item.query_id)
        if len(selected) >= k:
            break
    return selected


def build_messages(
    record: QueryRecord,
    train_records: list[QueryRecord],
    few_shot_k: int,
    seed: int,
    system_prompt: str,
) -> list[dict[str, str]]:
    rng = random.Random(f"{seed}:{record.query_id}")
    examples = select_few_shot_examples(record, train_records, few_shot_k, rng)
    messages = [{"role": "system", "content": system_prompt or default_system_prompt()}]
    messages.extend(few_shot_messages(record, examples))
    return messages


def prompt_record(
    record: QueryRecord,
    train_records: list[QueryRecord],
    few_shot_k: int,
    seed: int,
    system_prompt: str,
) -> dict[str, Any]:
    return {
        "query_id": record.query_id,
        "split": record.split,
        "input_smiles": record.input_smiles,
        "instruction": record.instruction,
        "primary_endpoint": record.primary_endpoint,
        "primary_direction": record.primary_direction,
        "gold_smiles": record.gold_smiles,
        "messages": build_messages(record, train_records, few_shot_k, seed, system_prompt),
    }


def write_jsonl(path: Path, rows: Iterable[dict[str, Any]], append: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    mode = "a" if append else "w"
    with path.open(mode, encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False, separators=(",", ":")) + "\n")


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if not path.is_file():
        return rows
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            rows.append(json.loads(line))
    return rows


def run_export_prompts(args: argparse.Namespace) -> int:
    records = load_records(args.data_dir, args.split, args.max_rows, args.seed)
    train_records = load_records(args.data_dir, args.few_shot_train_split, 0, args.seed) if args.few_shot_k else []
    rows = [
        prompt_record(record, train_records, args.few_shot_k, args.seed, args.system_prompt)
        for record in records
    ]
    write_jsonl(args.out, rows)
    print(json.dumps({"out": str(args.out), "records": len(rows)}, indent=2))
    return 0


def call_openai_compatible(
    messages: list[dict[str, str]],
    args: argparse.Namespace,
) -> str:
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
        headers={
            "Content-Type": "application/json",
            "Authorization": f"Bearer {api_key}",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=args.timeout_seconds) as response:
            data = json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"HTTP {exc.code}: {body}") from exc
    return clean_text(data["choices"][0]["message"]["content"])


def load_hf_model(args: argparse.Namespace) -> tuple[Any, Any]:
    try:
        import torch
        from transformers import AutoModelForCausalLM, AutoTokenizer
    except Exception as exc:  # noqa: BLE001
        raise SystemExit(
            "ERROR: hf-local requires torch and transformers. Install them in the selected environment."
        ) from exc

    dtype = args.torch_dtype
    torch_dtype = "auto"
    if dtype != "auto":
        torch_dtype = getattr(torch, dtype)
    tokenizer = AutoTokenizer.from_pretrained(args.model, trust_remote_code=args.trust_remote_code)
    model = AutoModelForCausalLM.from_pretrained(
        args.model,
        device_map=args.device_map,
        torch_dtype=torch_dtype,
        trust_remote_code=args.trust_remote_code,
    )
    model.eval()
    return tokenizer, model


def call_hf_local(messages: list[dict[str, str]], tokenizer: Any, model: Any, args: argparse.Namespace) -> str:
    try:
        import torch
    except Exception as exc:  # noqa: BLE001
        raise SystemExit("ERROR: hf-local requires torch.") from exc
    if hasattr(tokenizer, "apply_chat_template"):
        input_ids = tokenizer.apply_chat_template(
            messages,
            add_generation_prompt=True,
            return_tensors="pt",
        )
    else:
        prompt = "\n\n".join(f"{msg['role'].upper()}: {msg['content']}" for msg in messages) + "\n\nASSISTANT:"
        input_ids = tokenizer(prompt, return_tensors="pt").input_ids
    input_ids = input_ids.to(model.device)
    with torch.no_grad():
        output = model.generate(
            input_ids,
            max_new_tokens=args.max_tokens,
            do_sample=args.temperature > 0.0,
            temperature=max(args.temperature, 1e-6),
            pad_token_id=tokenizer.eos_token_id,
        )
    new_tokens = output[0][input_ids.shape[-1] :]
    return clean_text(tokenizer.decode(new_tokens, skip_special_tokens=True))


def existing_query_ids(path: Path) -> set[str]:
    return {clean_text(row.get("query_id")) for row in read_jsonl(path)}


def run_model(args: argparse.Namespace) -> int:
    records = load_records(args.data_dir, args.split, args.max_rows, args.seed)
    train_records = load_records(args.data_dir, args.few_shot_train_split, 0, args.seed) if args.few_shot_k else []
    done = existing_query_ids(args.out) if args.resume else set()
    tokenizer = model = None
    if args.provider == "hf-local":
        tokenizer, model = load_hf_model(args)

    wrote = 0
    for record in records:
        if record.query_id in done:
            continue
        messages = build_messages(record, train_records, args.few_shot_k, args.seed, args.system_prompt)
        started = time.time()
        error = ""
        output = ""
        try:
            if args.provider == "openai-compatible":
                output = call_openai_compatible(messages, args)
            elif args.provider == "hf-local":
                output = call_hf_local(messages, tokenizer, model, args)
            else:
                raise RuntimeError(f"unknown provider: {args.provider}")
        except Exception as exc:  # noqa: BLE001 - keep long runs resumable.
            error = str(exc)
        parsed = parse_smiles_from_output(output)
        row = {
            "query_id": record.query_id,
            "split": record.split,
            "provider": args.provider,
            "model": args.model,
            "few_shot_k": args.few_shot_k,
            "input_smiles": record.input_smiles,
            "instruction": record.instruction,
            "primary_endpoint": record.primary_endpoint,
            "gold_smiles": record.gold_smiles,
            "raw_output": output,
            "parsed_smiles": parsed,
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


def parse_smiles_from_output(output: str) -> list[str]:
    output = clean_text(output)
    if not output:
        return []
    parsed = parse_json_smiles(output)
    if parsed:
        return parsed[:2]

    candidates: list[str] = []
    for raw_line in output.replace(",", "\n").splitlines():
        line = raw_line.strip().strip("-*`'\" ")
        if not line or len(line) < 2:
            continue
        if SMILES_LINE_RE.match(line):
            candidates.append(line)
    if candidates:
        return unique_preserve_order(candidates)[:2]

    tokens = re.findall(r"[A-Za-z0-9@+\-\[\]\(\)=#$\\/%.:]{3,}", output)
    return unique_preserve_order(tokens)[:2]


def parse_json_smiles(output: str) -> list[str]:
    for pattern in [JSON_OBJECT_RE, JSON_ARRAY_RE]:
        match = pattern.search(output)
        if not match:
            continue
        try:
            payload = json.loads(match.group(0))
        except json.JSONDecodeError:
            continue
        smiles = smiles_from_payload(payload)
        if smiles:
            return smiles
    return []


def smiles_from_payload(payload: Any) -> list[str]:
    if isinstance(payload, list):
        values = payload
    elif isinstance(payload, dict):
        values = []
        for key in ["edited_smiles", "smiles", "molecules", "answers", "predictions"]:
            if key in payload:
                values = payload[key]
                break
    else:
        values = []
    if isinstance(values, str):
        values = [values]
    out: list[str] = []
    if isinstance(values, list):
        for value in values:
            if isinstance(value, dict):
                value = value.get("smiles") or value.get("SMILES") or value.get("target_smiles")
            text = clean_text(value)
            if text:
                out.append(text)
    return unique_preserve_order(out)


def unique_preserve_order(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        value = clean_text(value)
        if value and value not in seen:
            seen.add(value)
            out.append(value)
    return out


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


def descriptor_payload(smiles: str) -> dict[str, float] | None:
    if Chem is None or Descriptors is None or Crippen is None or QED is None:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "mw": float(Descriptors.MolWt(mol)),
        "logp": float(Crippen.MolLogP(mol)),
        "qed": float(QED.qed(mol)),
    }


def local_constraint_pass(input_smiles: str, pred_smiles: str, args: argparse.Namespace) -> bool | None:
    input_desc = descriptor_payload(input_smiles)
    pred_desc = descriptor_payload(pred_smiles)
    if input_desc is None or pred_desc is None:
        return None
    return (
        abs(pred_desc["mw"] - input_desc["mw"]) <= args.mw_delta
        and abs(pred_desc["logp"] - input_desc["logp"]) <= args.logp_delta
        and (input_desc["qed"] - pred_desc["qed"]) <= args.qed_drop
    )


def evaluate_predictions(args: argparse.Namespace) -> int:
    init_rdkit()
    records = {record.query_id: record for record in load_records(args.data_dir, args.split, args.max_rows, args.seed)}
    predictions = read_jsonl(args.predictions)
    counters: Counter[str] = Counter()
    hit_at_1: list[float] = []
    hit_at_k: list[float] = []
    recall_at_k: list[float] = []
    all_gold_at_k: list[float] = []
    valid_rates: list[float] = []
    duplicate_flags: list[float] = []
    best_pair_tanimoto: list[float] = []
    gold_coverage_tanimoto: list[float] = []
    local_pass_rates: list[float] = []
    per_query: list[dict[str, Any]] = []

    for pred in predictions:
        query_id = clean_text(pred.get("query_id"))
        record = records.get(query_id)
        if record is None:
            counters["prediction_without_matching_query"] += 1
            continue
        parsed = pred.get("parsed_smiles")
        if not isinstance(parsed, list):
            parsed = parse_smiles_from_output(clean_text(pred.get("raw_output")))
        parsed = [clean_text(item) for item in parsed[: args.top_k] if clean_text(item)]
        gold_canon = [canonical_smiles(item) for item in record.gold_smiles]
        gold_canon = [item for item in gold_canon if item]
        pred_canon = [canonical_smiles(item) for item in parsed]
        valid_pred_canon = [item for item in pred_canon if item]

        if not parsed:
            counters["empty_prediction"] += 1
        if pred.get("error"):
            counters["provider_error"] += 1
        valid_rate = len(valid_pred_canon) / max(len(parsed), 1)
        valid_rates.append(valid_rate)
        duplicate_flags.append(1.0 if len(set(valid_pred_canon)) < len(valid_pred_canon) else 0.0)

        matched = set(valid_pred_canon) & set(gold_canon)
        hit_at_1.append(1.0 if valid_pred_canon[:1] and valid_pred_canon[0] in set(gold_canon) else 0.0)
        hit_at_k.append(1.0 if matched else 0.0)
        recall_at_k.append(len(matched) / max(len(gold_canon), 1))
        all_gold_at_k.append(1.0 if len(matched) == len(gold_canon) and gold_canon else 0.0)

        pair_sims = [
            sim
            for pred_smiles in valid_pred_canon
            for gold_smiles in gold_canon
            for sim in [tanimoto(pred_smiles, gold_smiles)]
            if sim is not None
        ]
        if pair_sims:
            best_pair_tanimoto.append(max(pair_sims))
        gold_sims = []
        for gold_smiles in gold_canon:
            sims = [
                sim
                for pred_smiles in valid_pred_canon
                for sim in [tanimoto(pred_smiles, gold_smiles)]
                if sim is not None
            ]
            if sims:
                gold_sims.append(max(sims))
        if gold_sims:
            gold_coverage_tanimoto.append(sum(gold_sims) / len(gold_sims))

        local_flags = [
            flag
            for smiles in valid_pred_canon
            for flag in [local_constraint_pass(record.input_smiles, smiles, args)]
            if flag is not None
        ]
        if local_flags:
            local_pass_rates.append(sum(1 for flag in local_flags if flag) / len(local_flags))

        per_query.append(
            {
                "query_id": query_id,
                "gold_smiles": record.gold_smiles,
                "parsed_smiles": parsed,
                "valid_smiles": valid_pred_canon,
                "hit_at_k": bool(matched),
                "recall_at_k": recall_at_k[-1],
                "best_pair_tanimoto": max(pair_sims) if pair_sims else None,
            }
        )

    evaluated = len(per_query)
    metrics = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "data_dir": str(args.data_dir),
        "split": args.split,
        "predictions": str(args.predictions),
        "top_k": args.top_k,
        "rdkit_available": Chem is not None,
        "num_queries_in_split": len(records),
        "num_predictions": len(predictions),
        "num_evaluated": evaluated,
        "coverage": round(evaluated / len(records), 6) if records else 0.0,
        "valid_smiles_rate": mean(valid_rates),
        "duplicate_prediction_rate": mean(duplicate_flags),
        "exact_hit_at_1": mean(hit_at_1),
        f"exact_hit_at_{args.top_k}": mean(hit_at_k),
        f"exact_recall_at_{args.top_k}": mean(recall_at_k),
        f"all_gold_exact_at_{args.top_k}": mean(all_gold_at_k),
        "mean_best_pair_tanimoto": mean(best_pair_tanimoto),
        "mean_gold_coverage_tanimoto": mean(gold_coverage_tanimoto),
        "local_constraint_pass_rate": mean(local_pass_rates),
        "counters": dict(sorted(counters.items())),
        "per_query": per_query,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(metrics, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"out": str(args.out), "summary": {k: v for k, v in metrics.items() if k != "per_query"}}, indent=2))
    return 0


def mean(values: list[float]) -> float:
    return round(float(sum(values) / len(values)), 6) if values else 0.0


def main() -> int:
    init_rdkit()
    args = parse_args()
    if args.command == "export-prompts":
        return run_export_prompts(args)
    if args.command == "run":
        return run_model(args)
    if args.command == "evaluate":
        return evaluate_predictions(args)
    raise SystemExit(f"unknown command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())

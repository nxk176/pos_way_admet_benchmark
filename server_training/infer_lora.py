from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


JSON_RE = re.compile(r"\{.*\}", re.DOTALL)
FIELD_RE = re.compile(r'"(?P<key>stronger_smiles|hard_negative_smiles|smiles)"\s*:\s*"(?P<value>[^"]*)"')
ARRAY_FIELD_RE = re.compile(r'"(?P<key>edited_smiles|molecules|smiles)"\s*:\s*\[(?P<value>.*?)\]', re.DOTALL)
ARRAY_ITEM_RE = re.compile(r'"([^"]*)"')


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run inference with a base HF model plus a LoRA adapter.")
    parser.add_argument("--base-model", required=True)
    parser.add_argument("--adapter-dir", type=Path, required=True)
    parser.add_argument("--input-jsonl", type=Path, required=True)
    parser.add_argument("--output-jsonl", type=Path, required=True)
    parser.add_argument("--max-rows", type=int, default=0)
    parser.add_argument("--max-seq-length", type=int, default=2048)
    parser.add_argument("--max-new-tokens", type=int, default=192)
    parser.add_argument("--temperature", type=float, default=0.0)
    parser.add_argument("--top-p", type=float, default=1.0)
    parser.add_argument("--load-in-4bit", action="store_true")
    parser.add_argument("--bf16", action="store_true")
    parser.add_argument("--fp16", action="store_true")
    parser.add_argument("--trust-remote-code", action="store_true")
    parser.add_argument(
        "--device-map",
        choices=["single-gpu", "auto"],
        default="single-gpu",
        help="Use single-gpu for one selected CUDA device. Use auto only when intentional offload is configured.",
    )
    return parser.parse_args()


def read_jsonl(path: Path, max_rows: int = 0) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                rows.append(json.loads(line))
                if max_rows and len(rows) >= max_rows:
                    break
    return rows


def write_jsonl(path: Path, rows: list[dict[str, Any]], append: bool = True) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    mode = "a" if append else "w"
    with path.open(mode, encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False, separators=(",", ":")) + "\n")


def prompt_messages(row: dict[str, Any]) -> list[dict[str, str]]:
    messages = row["messages"]
    if messages and messages[-1].get("role") == "assistant":
        return messages[:-1]
    return messages


def parse_json_payload(text: str) -> Any:
    match = JSON_RE.search(text or "")
    if not match:
        return None
    candidate = match.group(0)
    try:
        return json.loads(candidate)
    except json.JSONDecodeError:
        pass

    repaired: dict[str, Any] = {}
    for field_match in FIELD_RE.finditer(candidate):
        repaired[field_match.group("key")] = field_match.group("value").strip()
    for field_match in ARRAY_FIELD_RE.finditer(candidate):
        values = [item.strip() for item in ARRAY_ITEM_RE.findall(field_match.group("value")) if item.strip()]
        if values:
            repaired[field_match.group("key")] = values
    return repaired or None


def main() -> int:
    args = parse_args()
    if args.output_jsonl.exists():
        args.output_jsonl.unlink()

    import torch
    from peft import PeftModel
    from transformers import AutoModelForCausalLM, AutoTokenizer, BitsAndBytesConfig

    dtype = torch.bfloat16 if args.bf16 else torch.float16
    quantization_config = None
    if args.load_in_4bit:
        quantization_config = BitsAndBytesConfig(
            load_in_4bit=True,
            bnb_4bit_quant_type="nf4",
            bnb_4bit_use_double_quant=True,
            bnb_4bit_compute_dtype=dtype,
        )

    tokenizer = AutoTokenizer.from_pretrained(args.base_model, trust_remote_code=args.trust_remote_code)
    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token

    if args.device_map == "single-gpu":
        if not torch.cuda.is_available():
            raise SystemExit("ERROR: --device-map single-gpu requires a CUDA GPU.")
        device_map: str | dict[str, int] = {"": 0}
    else:
        device_map = "auto"

    model = AutoModelForCausalLM.from_pretrained(
        args.base_model,
        device_map=device_map,
        quantization_config=quantization_config,
        torch_dtype=dtype if (args.bf16 or args.fp16) else "auto",
        trust_remote_code=args.trust_remote_code,
    )
    model = PeftModel.from_pretrained(model, str(args.adapter_dir))
    model.eval()

    rows = read_jsonl(args.input_jsonl, args.max_rows)
    for idx, row in enumerate(rows, start=1):
        messages = prompt_messages(row)
        if hasattr(tokenizer, "apply_chat_template"):
            inputs = tokenizer.apply_chat_template(
                messages,
                add_generation_prompt=True,
                tokenize=True,
                return_dict=True,
                return_tensors="pt",
            )
        else:
            prompt = "\n\n".join(f"{msg['role'].upper()}: {msg['content']}" for msg in messages) + "\n\nASSISTANT:"
            inputs = tokenizer(prompt, return_tensors="pt")
        inputs = {key: value.to(model.device) for key, value in inputs.items()}
        with torch.no_grad():
            output = model.generate(
                **inputs,
                max_new_tokens=args.max_new_tokens,
                do_sample=args.temperature > 0,
                temperature=max(args.temperature, 1e-6),
                top_p=args.top_p,
                pad_token_id=tokenizer.eos_token_id,
            )
        generated = output[0][inputs["input_ids"].shape[-1] :]
        text = tokenizer.decode(generated, skip_special_tokens=True).strip()
        out = {
            "id": row.get("id"),
            "task": row.get("task"),
            "source_dataset": row.get("source_dataset"),
            "input_smiles": row.get("input_smiles"),
            "gold_smiles": row.get("gold_smiles"),
            "hard_negative_smiles": row.get("hard_negative_smiles"),
            "raw_output": text,
            "parsed_json": parse_json_payload(text),
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
        }
        write_jsonl(args.output_jsonl, [out], append=True)
        if idx % 25 == 0:
            print(json.dumps({"event": "progress", "rows": idx, "out": str(args.output_jsonl)}))
    print(json.dumps({"event": "done", "rows": len(rows), "out": str(args.output_jsonl)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

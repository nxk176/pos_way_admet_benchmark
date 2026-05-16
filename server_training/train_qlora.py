from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train QLoRA adapter on POS-WAY chat JSONL data.")
    parser.add_argument("--base-model", required=True, help="HF model ID or local safetensors model directory.")
    parser.add_argument("--train-jsonl", type=Path, required=True)
    parser.add_argument("--val-jsonl", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-name", default="")
    parser.add_argument("--max-seq-length", type=int, default=2048)
    parser.add_argument("--max-train-samples", type=int, default=0)
    parser.add_argument("--max-val-samples", type=int, default=0)
    parser.add_argument("--num-train-epochs", type=float, default=1.0)
    parser.add_argument("--learning-rate", type=float, default=2e-4)
    parser.add_argument("--warmup-ratio", type=float, default=0.03)
    parser.add_argument("--weight-decay", type=float, default=0.0)
    parser.add_argument("--per-device-train-batch-size", type=int, default=1)
    parser.add_argument("--per-device-eval-batch-size", type=int, default=1)
    parser.add_argument("--gradient-accumulation-steps", type=int, default=16)
    parser.add_argument("--logging-steps", type=int, default=10)
    parser.add_argument("--eval-steps", type=int, default=250)
    parser.add_argument("--save-steps", type=int, default=250)
    parser.add_argument("--save-total-limit", type=int, default=2)
    parser.add_argument("--lora-r", type=int, default=16)
    parser.add_argument("--lora-alpha", type=int, default=32)
    parser.add_argument("--lora-dropout", type=float, default=0.05)
    parser.add_argument(
        "--target-modules",
        default="q_proj,k_proj,v_proj,o_proj,gate_proj,up_proj,down_proj",
        help="Comma-separated module names. Default fits Qwen/Qwen-derived models.",
    )
    parser.add_argument("--optim", default="paged_adamw_8bit")
    parser.add_argument("--seed", type=int, default=29)
    parser.add_argument("--bf16", action="store_true")
    parser.add_argument("--fp16", action="store_true")
    parser.add_argument("--load-in-4bit", action="store_true")
    parser.add_argument("--gradient-checkpointing", action="store_true")
    parser.add_argument("--trust-remote-code", action="store_true")
    parser.add_argument("--train-on-assistant-only", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--resume-from-checkpoint", default=None)
    return parser.parse_args()


def fail_if_missing(path: Path, label: str) -> None:
    if not path.is_file():
        raise SystemExit(f"ERROR: missing {label}: {path}")


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def chat_to_text(tokenizer: Any, messages: list[dict[str, str]], add_generation_prompt: bool) -> str:
    if hasattr(tokenizer, "apply_chat_template"):
        return tokenizer.apply_chat_template(
            messages,
            tokenize=False,
            add_generation_prompt=add_generation_prompt,
        )
    text = "\n\n".join(f"{msg['role'].upper()}: {msg['content']}" for msg in messages)
    if add_generation_prompt:
        text += "\n\nASSISTANT:"
    return text


def training_texts(tokenizer: Any, messages: list[dict[str, str]]) -> tuple[str, str]:
    if not messages or messages[-1].get("role") != "assistant":
        raise ValueError("Every SFT row must end with an assistant message.")
    full_text = chat_to_text(tokenizer, messages, add_generation_prompt=False)
    prefix_text = chat_to_text(tokenizer, messages[:-1], add_generation_prompt=True)
    if tokenizer.eos_token and not full_text.endswith(tokenizer.eos_token):
        full_text += tokenizer.eos_token
    return full_text, prefix_text


@dataclass
class AssistantOnlyCollator:
    tokenizer: Any

    def __call__(self, features: list[dict[str, Any]]) -> dict[str, Any]:
        labels = [item.pop("labels") for item in features]
        batch = self.tokenizer.pad(features, padding=True, return_tensors="pt")
        max_len = batch["input_ids"].shape[1]
        import torch

        padded_labels = []
        for label in labels:
            pad_len = max_len - len(label)
            padded_labels.append(label + [-100] * pad_len)
        batch["labels"] = torch.tensor(padded_labels, dtype=torch.long)
        return batch


def main() -> int:
    args = parse_args()
    fail_if_missing(args.train_jsonl, "train JSONL")
    fail_if_missing(args.val_jsonl, "validation JSONL")

    import torch
    from datasets import load_dataset
    from peft import LoraConfig, get_peft_model, prepare_model_for_kbit_training
    from transformers import AutoModelForCausalLM, AutoTokenizer, BitsAndBytesConfig, Trainer, TrainingArguments

    if args.bf16 and args.fp16:
        raise SystemExit("ERROR: choose only one of --bf16 or --fp16.")

    tokenizer = AutoTokenizer.from_pretrained(args.base_model, trust_remote_code=args.trust_remote_code)
    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token
    tokenizer.padding_side = "right"

    compute_dtype = torch.bfloat16 if args.bf16 else torch.float16
    quantization_config = None
    if args.load_in_4bit:
        quantization_config = BitsAndBytesConfig(
            load_in_4bit=True,
            bnb_4bit_quant_type="nf4",
            bnb_4bit_use_double_quant=True,
            bnb_4bit_compute_dtype=compute_dtype,
        )

    model = AutoModelForCausalLM.from_pretrained(
        args.base_model,
        device_map="auto",
        quantization_config=quantization_config,
        torch_dtype=compute_dtype if (args.bf16 or args.fp16) else "auto",
        trust_remote_code=args.trust_remote_code,
    )
    model.config.use_cache = False
    if args.gradient_checkpointing:
        model.gradient_checkpointing_enable(gradient_checkpointing_kwargs={"use_reentrant": False})
    if args.load_in_4bit:
        model = prepare_model_for_kbit_training(model)

    lora_config = LoraConfig(
        r=args.lora_r,
        lora_alpha=args.lora_alpha,
        lora_dropout=args.lora_dropout,
        bias="none",
        task_type="CAUSAL_LM",
        target_modules=[item.strip() for item in args.target_modules.split(",") if item.strip()],
    )
    model = get_peft_model(model, lora_config)
    model.print_trainable_parameters()

    data_files = {"train": str(args.train_jsonl), "validation": str(args.val_jsonl)}
    dataset = load_dataset("json", data_files=data_files)
    if args.max_train_samples:
        dataset["train"] = dataset["train"].select(range(min(args.max_train_samples, len(dataset["train"]))))
    if args.max_val_samples:
        dataset["validation"] = dataset["validation"].select(range(min(args.max_val_samples, len(dataset["validation"]))))

    def tokenize_row(example: dict[str, Any]) -> dict[str, Any]:
        messages = example["messages"]
        full_text, prefix_text = training_texts(tokenizer, messages)
        full = tokenizer(full_text, truncation=True, max_length=args.max_seq_length, padding=False)
        if args.train_on_assistant_only:
            prefix = tokenizer(prefix_text, truncation=True, max_length=args.max_seq_length, padding=False)
            prefix_len = min(len(prefix["input_ids"]), len(full["input_ids"]))
            labels = [-100] * prefix_len + full["input_ids"][prefix_len:]
        else:
            labels = list(full["input_ids"])
        full["labels"] = labels
        return full

    remove_columns = dataset["train"].column_names
    tokenized = dataset.map(tokenize_row, remove_columns=remove_columns, desc="Tokenizing")

    training_kwargs = {
        "output_dir": str(args.output_dir),
        "run_name": args.run_name or args.output_dir.name,
        "num_train_epochs": args.num_train_epochs,
        "learning_rate": args.learning_rate,
        "warmup_ratio": args.warmup_ratio,
        "weight_decay": args.weight_decay,
        "per_device_train_batch_size": args.per_device_train_batch_size,
        "per_device_eval_batch_size": args.per_device_eval_batch_size,
        "gradient_accumulation_steps": args.gradient_accumulation_steps,
        "logging_steps": args.logging_steps,
        "eval_steps": args.eval_steps,
        "save_steps": args.save_steps,
        "save_total_limit": args.save_total_limit,
        "save_strategy": "steps",
        "optim": args.optim,
        "bf16": args.bf16,
        "fp16": args.fp16,
        "report_to": "none",
        "seed": args.seed,
        "remove_unused_columns": False,
    }
    try:
        training_args = TrainingArguments(evaluation_strategy="steps", **training_kwargs)
    except TypeError:
        training_args = TrainingArguments(eval_strategy="steps", **training_kwargs)

    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=tokenized["train"],
        eval_dataset=tokenized["validation"],
        data_collator=AssistantOnlyCollator(tokenizer),
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_json(
        args.output_dir / "run_config.json",
        {
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "base_model": args.base_model,
            "train_jsonl": str(args.train_jsonl),
            "val_jsonl": str(args.val_jsonl),
            "args": vars(args),
        },
    )
    trainer.train(resume_from_checkpoint=args.resume_from_checkpoint)
    trainer.save_model(str(args.output_dir))
    tokenizer.save_pretrained(str(args.output_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Train a LoRA/QLoRA adapter on POS-WAY chat-style JSONL files."
    )
    parser.add_argument("--base-model", required=True, help="HF model ID or local safetensors model path.")
    parser.add_argument("--train-jsonl", type=Path, required=True)
    parser.add_argument("--val-jsonl", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--max-seq-length", type=int, default=2048)
    parser.add_argument("--num-train-epochs", type=float, default=1.0)
    parser.add_argument("--learning-rate", type=float, default=2e-4)
    parser.add_argument("--per-device-train-batch-size", type=int, default=1)
    parser.add_argument("--per-device-eval-batch-size", type=int, default=1)
    parser.add_argument("--gradient-accumulation-steps", type=int, default=16)
    parser.add_argument("--warmup-ratio", type=float, default=0.03)
    parser.add_argument("--logging-steps", type=int, default=10)
    parser.add_argument("--eval-steps", type=int, default=200)
    parser.add_argument("--save-steps", type=int, default=200)
    parser.add_argument("--save-total-limit", type=int, default=2)
    parser.add_argument("--optim", default="paged_adamw_8bit")
    parser.add_argument("--eval-strategy", choices=["no", "steps", "epoch"], default="steps")
    parser.add_argument("--save-strategy", choices=["no", "steps", "epoch"], default="steps")
    parser.add_argument("--seed", type=int, default=29)
    parser.add_argument("--lora-r", type=int, default=16)
    parser.add_argument("--lora-alpha", type=int, default=32)
    parser.add_argument("--lora-dropout", type=float, default=0.05)
    parser.add_argument(
        "--target-modules",
        default="q_proj,k_proj,v_proj,o_proj,gate_proj,up_proj,down_proj",
        help="Comma-separated LoRA target module names for Qwen/Llama-style decoders.",
    )
    parser.add_argument("--load-in-4bit", action="store_true", help="Use QLoRA 4-bit loading via bitsandbytes.")
    parser.add_argument("--gradient-checkpointing", action="store_true")
    parser.add_argument("--trust-remote-code", action="store_true")
    parser.add_argument("--bf16", action="store_true")
    parser.add_argument("--fp16", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if not args.train_jsonl.is_file():
        raise SystemExit(f"ERROR: missing train JSONL: {args.train_jsonl}")
    if not args.val_jsonl.is_file():
        raise SystemExit(f"ERROR: missing val JSONL: {args.val_jsonl}")

    import torch
    from datasets import load_dataset
    from peft import LoraConfig, get_peft_model, prepare_model_for_kbit_training
    from transformers import (
        AutoModelForCausalLM,
        AutoTokenizer,
        BitsAndBytesConfig,
        DataCollatorForLanguageModeling,
        Trainer,
        TrainingArguments,
    )

    tokenizer = AutoTokenizer.from_pretrained(args.base_model, trust_remote_code=args.trust_remote_code)
    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token
    tokenizer.padding_side = "right"

    quantization_config = None
    if args.load_in_4bit:
        quantization_config = BitsAndBytesConfig(
            load_in_4bit=True,
            bnb_4bit_quant_type="nf4",
            bnb_4bit_use_double_quant=True,
            bnb_4bit_compute_dtype=torch.bfloat16 if args.bf16 else torch.float16,
        )

    model = AutoModelForCausalLM.from_pretrained(
        args.base_model,
        device_map="auto",
        quantization_config=quantization_config,
        torch_dtype=torch.bfloat16 if args.bf16 else (torch.float16 if args.fp16 else "auto"),
        trust_remote_code=args.trust_remote_code,
    )
    if args.gradient_checkpointing:
        model.gradient_checkpointing_enable()
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

    dataset = load_dataset(
        "json",
        data_files={"train": str(args.train_jsonl), "validation": str(args.val_jsonl)},
    )

    def render_chat(example: dict) -> dict:
        messages = example["messages"]
        if hasattr(tokenizer, "apply_chat_template"):
            text = tokenizer.apply_chat_template(messages, tokenize=False, add_generation_prompt=False)
        else:
            text = "\n\n".join(f"{msg['role'].upper()}: {msg['content']}" for msg in messages)
        if tokenizer.eos_token and not text.endswith(tokenizer.eos_token):
            text += tokenizer.eos_token
        return {"text": text}

    dataset = dataset.map(render_chat, remove_columns=dataset["train"].column_names)

    def tokenize(example: dict) -> dict:
        return tokenizer(
            example["text"],
            truncation=True,
            max_length=args.max_seq_length,
            padding=False,
        )

    tokenized = dataset.map(tokenize, batched=False, remove_columns=["text"])

    training_kwargs = {
        "output_dir": str(args.out_dir),
        "num_train_epochs": args.num_train_epochs,
        "learning_rate": args.learning_rate,
        "per_device_train_batch_size": args.per_device_train_batch_size,
        "per_device_eval_batch_size": args.per_device_eval_batch_size,
        "gradient_accumulation_steps": args.gradient_accumulation_steps,
        "warmup_ratio": args.warmup_ratio,
        "logging_steps": args.logging_steps,
        "eval_steps": args.eval_steps,
        "save_steps": args.save_steps,
        "save_total_limit": args.save_total_limit,
        "save_strategy": args.save_strategy,
        "optim": args.optim,
        "report_to": "none",
        "seed": args.seed,
        "bf16": args.bf16,
        "fp16": args.fp16,
        "remove_unused_columns": False,
    }
    try:
        training_args = TrainingArguments(evaluation_strategy=args.eval_strategy, **training_kwargs)
    except TypeError:
        training_args = TrainingArguments(eval_strategy=args.eval_strategy, **training_kwargs)

    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=tokenized["train"],
        eval_dataset=tokenized["validation"],
        data_collator=DataCollatorForLanguageModeling(tokenizer=tokenizer, mlm=False),
    )
    trainer.train()
    trainer.save_model(str(args.out_dir))
    tokenizer.save_pretrained(str(args.out_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

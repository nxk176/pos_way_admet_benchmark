from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from huggingface_hub import snapshot_download


ALLOW_PATTERNS = [
    "*.json",
    "*.safetensors",
    "*.model",
    "*.txt",
    "*.py",
    "tokenizer*",
    "merges.txt",
    "vocab.json",
    "generation_config.json",
    "chat_template.jinja",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download HF safetensors models needed for POS-WAY LoRA runs.")
    parser.add_argument("--manifest", type=Path, default=Path("server_training/model_manifest.json"))
    parser.add_argument("--model", default="all", help="Model short name from manifest, or 'all'.")
    parser.add_argument("--revision", default=None)
    parser.add_argument("--token", default=None, help="Optional HF token. If omitted, huggingface_hub uses normal env/cache auth.")
    return parser.parse_args()


def load_manifest(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        raise SystemExit(f"ERROR: manifest not found: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    models = payload.get("models", [])
    if not isinstance(models, list) or not models:
        raise SystemExit(f"ERROR: manifest has no models: {path}")
    return models


def main() -> int:
    args = parse_args()
    models = load_manifest(args.manifest)
    if args.model != "all":
        models = [item for item in models if item.get("name") == args.model]
        if not models:
            raise SystemExit(f"ERROR: model {args.model!r} not found in {args.manifest}")

    for item in models:
        name = item["name"]
        hf_id = item["hf_id"]
        local_dir = Path(item["local_dir"])
        local_dir.mkdir(parents=True, exist_ok=True)
        print(json.dumps({"event": "download_start", "name": name, "hf_id": hf_id, "local_dir": str(local_dir)}))
        path = snapshot_download(
            repo_id=hf_id,
            revision=args.revision,
            local_dir=str(local_dir),
            allow_patterns=ALLOW_PATTERNS,
            token=args.token,
        )
        print(json.dumps({"event": "download_done", "name": name, "path": path}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

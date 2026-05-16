from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import zipfile
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

try:
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit.Contrib.SA_Score import sascorer
except Exception as exc:  # noqa: BLE001
    raise SystemExit(
        "ERROR: extract_pubchem_final_supported_observations.py requires RDKit. "
        "Run it with .\\myenv311\\Scripts\\python.exe"
    ) from exc


FIELDS = [
    "observation_uid",
    "source",
    "source_version",
    "property_family",
    "endpoint_name",
    "label_type",
    "experimental_only_flag",
    "molecule_chembl_id",
    "smiles_canon",
    "inchikey",
    "connectivity_key",
    "murcko_scaffold",
    "condition_bucket",
    "value_raw",
    "unit_raw",
    "relation_raw",
    "type_raw",
    "value_canonical",
    "unit_canonical",
    "confidence",
    "target_chembl_id",
    "target_label",
    "assay_chembl_id",
    "assay_description",
    "assay_type",
    "assay_organism",
    "assay_cell_type",
    "document_year",
    "doi",
    "pubmed_id",
    "mw",
    "logp",
    "tpsa",
    "hba",
    "hbd",
    "rotatable_bonds",
    "heavy_atoms",
    "qed",
    "sa_score",
]

CONCENTRATION_TYPES = {
    "AC50",
    "CC50",
    "EC50",
    "GI50",
    "IC50",
    "KD",
    "KI",
    "LC50",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract exact p-scale PubChem BioAssay observations needed for source-specific multi-answer rows."
    )
    parser.add_argument(
        "--zip-dir",
        type=Path,
        default=Path("pos_way_admet_benchmark/raw/public/pubchem/bioassay_csv"),
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/pubchem_bioassay/pubchem_final_supported_observations.csv"),
    )
    parser.add_argument(
        "--stats",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/pubchem_bioassay/pubchem_final_supported_stats.json"),
    )
    parser.add_argument("--source-version", default="pubchem_bioassay_csv_data_full")
    parser.add_argument("--progress-every", type=int, default=100)
    parser.add_argument("--start-aid", type=int, default=0, help="Resume from the zip shard containing this AID.")
    parser.add_argument("--append", action="store_true", help="Append to an existing CSV instead of overwriting it.")
    return parser.parse_args()


def clean_text(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value or "").strip())


def as_float(value: Any) -> float | None:
    try:
        parsed = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def aid_from_name(name: str) -> str:
    match = re.search(r"(\d+)\.csv\.gz$", name)
    return match.group(1) if match else ""


def zip_range(path: Path) -> tuple[int, int] | None:
    match = re.search(r"(\d+)_(\d+)\.zip$", path.name)
    if not match:
        return None
    return int(match.group(1)), int(match.group(2))


def canonical_type(raw: str) -> str:
    return clean_text(raw).upper().replace("-", "").replace("_", "").replace(" ", "")


def canonical_relation(raw: str) -> str:
    value = clean_text(raw)
    return value if value in {"=", "<", ">", "<=", ">="} else value


def pscale_value(std_type: str, value: float, unit: str) -> tuple[float, str] | None:
    type_norm = canonical_type(std_type)
    unit_norm = clean_text(unit).lower()
    if type_norm not in CONCENTRATION_TYPES or value <= 0:
        return None
    if unit_norm in {"nm", "nanomolar"}:
        return round(9.0 - math.log10(value), 6), f"p{type_norm}"
    if unit_norm in {"um", "Âµm", "micromolar"}:
        return round(6.0 - math.log10(value), 6), f"p{type_norm}"
    if unit_norm in {"mm", "millimolar"}:
        return round(3.0 - math.log10(value), 6), f"p{type_norm}"
    return None


def descriptor_record(smiles: str, cache: dict[str, dict[str, Any] | None]) -> dict[str, Any] | None:
    cached = cache.get(smiles)
    if smiles in cache:
        return cached
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        cache[smiles] = None
        return None
    smiles_canon = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    inchikey = Chem.MolToInchiKey(mol)
    scaffold = ""
    scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold_mol is not None and scaffold_mol.GetNumHeavyAtoms() > 0:
        scaffold = Chem.MolToSmiles(scaffold_mol, canonical=True, isomericSmiles=False)
    record = {
        "smiles_canon": smiles_canon,
        "inchikey": inchikey,
        "connectivity_key": inchikey.split("-")[0] if inchikey else "",
        "murcko_scaffold": scaffold,
        "mw": round(float(Descriptors.MolWt(mol)), 4),
        "logp": round(float(Crippen.MolLogP(mol)), 4),
        "tpsa": round(float(rdMolDescriptors.CalcTPSA(mol)), 4),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "rotatable_bonds": int(Lipinski.NumRotatableBonds(mol)),
        "heavy_atoms": int(mol.GetNumHeavyAtoms()),
        "qed": round(float(QED.qed(mol)), 6),
        "sa_score": round(float(sascorer.calculateScore(mol)), 6),
    }
    cache[smiles] = record
    return record


def data_rows_from_csv_gz(payload: bytes) -> Iterable[dict[str, str]]:
    text = gzip.decompress(payload).decode("utf-8", errors="replace").splitlines()
    if not text:
        return
    reader = csv.DictReader(text)
    for row in reader:
        tag = clean_text(row.get("PUBCHEM_RESULT_TAG"))
        if tag and not tag.startswith("RESULT_"):
            yield row


def iter_zip_rows(zip_path: Path, desc_cache: dict[str, dict[str, Any] | None], counters: Counter[str]) -> Iterable[dict[str, Any]]:
    with zipfile.ZipFile(zip_path) as zf:
        for info in zf.infolist():
            if not info.filename.endswith(".csv.gz"):
                continue
            counters["csv_gz_entries"] += 1
            aid = aid_from_name(info.filename)
            if not aid:
                counters["skipped_missing_aid"] += 1
                continue
            for row in data_rows_from_csv_gz(zf.read(info)):
                relation = canonical_relation(row.get("Standard Relation"))
                if relation != "=":
                    counters["skipped_non_exact_relation"] += 1
                    continue
                std_type_raw = clean_text(row.get("Standard Type"))
                std_type = canonical_type(std_type_raw)
                value = as_float(row.get("Standard Value"))
                unit = clean_text(row.get("Standard Units"))
                if value is None:
                    counters["skipped_missing_standard_value"] += 1
                    continue
                canonical = pscale_value(std_type, value, unit)
                if canonical is None:
                    counters["skipped_non_supported_endpoint"] += 1
                    continue
                smiles = clean_text(row.get("PUBCHEM_EXT_DATASOURCE_SMILES"))
                cid = clean_text(row.get("PUBCHEM_CID"))
                sid = clean_text(row.get("PUBCHEM_SID"))
                if not smiles or not cid:
                    counters["skipped_missing_cid_or_smiles"] += 1
                    continue
                desc = descriptor_record(smiles, desc_cache)
                if desc is None:
                    counters["skipped_invalid_smiles"] += 1
                    continue
                value_canonical, unit_canonical = canonical
                endpoint = f"PubChem_AID{aid}_{unit_canonical}"
                observation_uid = f"pubchem_aid{aid}_sid{sid or 'unknown'}_cid{cid}_{std_type}"
                counters["kept"] += 1
                counters[f"type:{std_type}"] += 1
                yield {
                    "observation_uid": observation_uid,
                    "source": "PubChem BioAssay",
                    "source_version": "pubchem_bioassay_csv_data_full",
                    "property_family": "pubchem_bioassay_activity",
                    "endpoint_name": endpoint,
                    "label_type": "experimental",
                    "experimental_only_flag": "true",
                    "molecule_chembl_id": "",
                    "smiles_canon": desc["smiles_canon"],
                    "inchikey": desc["inchikey"],
                    "connectivity_key": desc["connectivity_key"],
                    "murcko_scaffold": desc["murcko_scaffold"],
                    "condition_bucket": f"{endpoint}|{std_type}|AID:{aid}",
                    "value_raw": value,
                    "unit_raw": unit,
                    "relation_raw": relation,
                    "type_raw": std_type_raw or std_type,
                    "value_canonical": value_canonical,
                    "unit_canonical": unit_canonical,
                    "confidence": 0.8,
                    "target_chembl_id": "",
                    "target_label": "",
                    "assay_chembl_id": f"PubChem_AID:{aid}",
                    "assay_description": "",
                    "assay_type": "",
                    "assay_organism": "",
                    "assay_cell_type": "",
                    "document_year": "",
                    "doi": "",
                    "pubmed_id": "",
                    "mw": desc["mw"],
                    "logp": desc["logp"],
                    "tpsa": desc["tpsa"],
                    "hba": desc["hba"],
                    "hbd": desc["hbd"],
                    "rotatable_bonds": desc["rotatable_bonds"],
                    "heavy_atoms": desc["heavy_atoms"],
                    "qed": desc["qed"],
                    "sa_score": desc["sa_score"],
                }


def main() -> int:
    args = parse_args()
    zip_paths = sorted(args.zip_dir.glob("*.zip"))
    if args.start_aid:
        filtered = []
        for path in zip_paths:
            span = zip_range(path)
            if span is None or span[1] >= args.start_aid:
                filtered.append(path)
        zip_paths = filtered
    if not zip_paths:
        raise SystemExit(f"ERROR: no zip files found in {args.zip_dir}")
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    counters: Counter[str] = Counter()
    desc_cache: dict[str, dict[str, Any] | None] = {}
    unique_molecules: set[str] = set()
    unique_endpoints: set[str] = set()
    mode = "a" if args.append and args.out_csv.exists() else "w"
    with args.out_csv.open(mode, encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, extrasaction="ignore")
        if mode == "w":
            writer.writeheader()
        for idx, zip_path in enumerate(zip_paths, start=1):
            counters["zip_files"] += 1
            for out_row in iter_zip_rows(zip_path, desc_cache, counters):
                writer.writerow(out_row)
                unique_molecules.add(out_row["connectivity_key"])
                unique_endpoints.add(out_row["endpoint_name"])
            if args.progress_every and idx % args.progress_every == 0:
                print(
                    json.dumps(
                        {
                            "zip_files_done": idx,
                            "zip_files_total": len(zip_paths),
                            "kept": counters["kept"],
                            "unique_molecules": len(unique_molecules),
                            "unique_endpoints": len(unique_endpoints),
                        },
                        ensure_ascii=False,
                    ),
                    flush=True,
                )
    stats = {
        "source": "PubChem BioAssay CSV/Data full local crawl",
        "zip_dir": str(args.zip_dir),
        "out_csv": str(args.out_csv),
        "zip_files": len(zip_paths),
        "kept_rows": counters["kept"],
        "unique_molecules": len(unique_molecules),
        "unique_endpoints": len(unique_endpoints),
        "descriptor_cache_entries": len(desc_cache),
        "counters": dict(sorted(counters.items())),
        "row_policy": "Exact relation only; concentration-like endpoints only; p-scale canonical values only.",
    }
    args.stats.parent.mkdir(parents=True, exist_ok=True)
    args.stats.write_text(json.dumps(stats, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(stats, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


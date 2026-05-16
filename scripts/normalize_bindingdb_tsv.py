from __future__ import annotations

import argparse
import csv
import json
import math
import re
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

try:
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit.Contrib.SA_Score import sascorer
except Exception as exc:  # noqa: BLE001
    raise SystemExit(
        "ERROR: normalize_bindingdb_tsv.py requires RDKit. "
        "Run it with the project environment, for example .\\myenv311\\Scripts\\python.exe"
    ) from exc


FIELDNAMES = [
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

MEASUREMENT_COLUMNS = {
    "Ki (nM)": "pKi",
    "IC50 (nM)": "pIC50",
    "Kd (nM)": "pKd",
    "EC50 (nM)": "pEC50",
}
VALUE_RE = re.compile(r"([<>]=?|=)?\s*([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Normalize BindingDB TSV/ZIP into property_observations CSV schema.")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("pos_way_admet_benchmark/raw/public/bindingdb/BindingDB_BindingDB_Articles_202605_tsv.zip"),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/bindingdb/bindingdb_curated_observations.csv"),
    )
    parser.add_argument(
        "--stats",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/bindingdb/bindingdb_curated_observation_stats.json"),
    )
    parser.add_argument("--source-version", default="202605_curated_articles")
    parser.add_argument("--max-rows", type=int, default=0, help="0 means no row cap.")
    parser.add_argument("--max-rows-per-endpoint", type=int, default=0, help="0 means no endpoint cap.")
    parser.add_argument("--skip-rdkit-descriptors", action="store_true")
    parser.add_argument("--compute-sa-score", action="store_true")
    return parser.parse_args()


def clean_text(value: Any) -> str:
    return str(value or "").replace("\r", " ").replace("\n", " ").strip()


def slug(value: str, max_len: int = 64) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")
    return (cleaned or "target_unknown")[:max_len]


def connectivity_key(inchikey: str) -> str:
    return str(inchikey or "").split("-", 1)[0]


def read_bindingdb_rows(path: Path) -> Iterable[dict[str, str]]:
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as zf:
            entries = [entry for entry in zf.infolist() if not entry.is_dir()]
            if not entries:
                raise SystemExit(f"ERROR: zip has no files: {path}")
            with zf.open(entries[0], "r") as raw:
                text = (line.decode("utf-8", errors="replace") for line in raw)
                yield from csv.DictReader(text, delimiter="\t")
    else:
        with path.open("r", encoding="utf-8", newline="") as handle:
            yield from csv.DictReader(handle, delimiter="\t")


def descriptor_record(
    smiles: str,
    inchikey: str,
    compute_sa_score: bool,
    skip_rdkit_descriptors: bool,
    cache: dict[str, dict[str, Any]],
) -> dict[str, Any] | None:
    cache_key = inchikey or smiles
    if cache_key in cache:
        return cache[cache_key]
    if skip_rdkit_descriptors:
        record = {
            "smiles_canon": smiles,
            "murcko_scaffold": "",
            "mw": "",
            "logp": "",
            "tpsa": "",
            "hba": "",
            "hbd": "",
            "rotatable_bonds": "",
            "heavy_atoms": "",
            "qed": "",
            "sa_score": "",
        }
        cache[cache_key] = record
        return record

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    scaffold = ""
    try:
        scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
        if scaffold_mol is not None and scaffold_mol.GetNumHeavyAtoms() > 0:
            scaffold = Chem.MolToSmiles(scaffold_mol, canonical=True, isomericSmiles=False)
    except Exception:  # noqa: BLE001
        scaffold = ""
    record = {
        "smiles_canon": Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
        "murcko_scaffold": scaffold,
        "mw": round(float(Descriptors.MolWt(mol)), 4),
        "logp": round(float(Crippen.MolLogP(mol)), 4),
        "tpsa": round(float(rdMolDescriptors.CalcTPSA(mol)), 4),
        "hba": int(Lipinski.NumHAcceptors(mol)),
        "hbd": int(Lipinski.NumHDonors(mol)),
        "rotatable_bonds": int(Lipinski.NumRotatableBonds(mol)),
        "heavy_atoms": int(mol.GetNumHeavyAtoms()),
        "qed": round(float(QED.qed(mol)), 4),
        "sa_score": round(float(sascorer.calculateScore(mol)), 4) if compute_sa_score else "",
    }
    cache[cache_key] = record
    return record


def parse_nm_value(value: str) -> tuple[str, float, float] | None:
    text = clean_text(value).replace(",", "")
    if not text:
        return None
    match = VALUE_RE.search(text)
    if not match:
        return None
    relation = match.group(1) or "="
    nm_value = float(match.group(2))
    if not math.isfinite(nm_value) or nm_value <= 0.0:
        return None
    canonical = 9.0 - math.log10(nm_value)
    return relation, nm_value, canonical


def endpoint_from_row(row: dict[str, str], unit: str) -> tuple[str, str]:
    uniprot = clean_text(row.get("UniProt (SwissProt) Primary ID of Target Chain 1"))
    if uniprot:
        target_id = f"BDB_{slug(uniprot, 32)}"
    else:
        target_id = f"BDB_{slug(row.get('Target Name', ''), 64)}"
    return f"{target_id}_{unit}", target_id


def normalize(args: argparse.Namespace) -> dict[str, Any]:
    if not args.input.is_file():
        raise SystemExit(f"ERROR: BindingDB input not found: {args.input}")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.stats.parent.mkdir(parents=True, exist_ok=True)

    counters: Counter[str] = Counter()
    endpoint_counts: Counter[str] = Counter()
    descriptor_cache: dict[str, dict[str, Any]] = {}
    written = 0
    input_rows = 0

    with args.out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        for row in read_bindingdb_rows(args.input):
            if args.max_rows and input_rows >= args.max_rows:
                break
            input_rows += 1
            smiles = clean_text(row.get("Ligand SMILES"))
            inchikey = clean_text(row.get("Ligand InChI Key"))
            if not smiles:
                counters["skipped_missing_smiles"] += 1
                continue
            desc = descriptor_record(
                smiles,
                inchikey,
                args.compute_sa_score,
                args.skip_rdkit_descriptors,
                descriptor_cache,
            )
            if desc is None:
                counters["skipped_invalid_smiles"] += 1
                continue

            reactant_set_id = clean_text(row.get("BindingDB Reactant_set_id"))
            monomer_id = clean_text(row.get("BindingDB MonomerID"))
            molecule_id = f"BDBM{monomer_id}" if monomer_id else f"BDBR{reactant_set_id}"
            pubmed_id = clean_text(row.get("PMID"))
            doi = clean_text(row.get("Article DOI")) or clean_text(row.get("BindingDB Entry DOI"))
            pub_year = clean_text(row.get("Date of publication"))[-4:]
            document_year = pub_year if pub_year.isdigit() else ""
            target_label = clean_text(row.get("Target Name"))
            target_organism = clean_text(row.get("Target Source Organism According to Curator or DataSource"))
            pubchem_aid = clean_text(row.get("PubChem AID"))
            institution = clean_text(row.get("Institution"))

            for raw_column, unit in MEASUREMENT_COLUMNS.items():
                parsed = parse_nm_value(row.get(raw_column, ""))
                if parsed is None:
                    continue
                endpoint, target_id = endpoint_from_row(row, unit)
                if args.max_rows_per_endpoint and endpoint_counts[endpoint] >= args.max_rows_per_endpoint:
                    continue

                relation, nm_value, canonical = parsed
                condition_bucket = "|".join(
                    [
                        endpoint,
                        raw_column.removesuffix(" (nM)"),
                        f"pmid:{pubmed_id or 'unknown'}",
                        f"doi:{doi or 'unknown'}",
                    ]
                )
                observation_uid = f"bindingdb_{args.source_version}_{reactant_set_id}_{unit}"
                out = {
                    "observation_uid": observation_uid,
                    "source": "BindingDB",
                    "source_version": args.source_version,
                    "property_family": "bindingdb_target_bioactivity",
                    "endpoint_name": endpoint,
                    "label_type": "experimental",
                    "experimental_only_flag": "true",
                    "molecule_chembl_id": molecule_id,
                    "smiles_canon": desc["smiles_canon"],
                    "inchikey": inchikey,
                    "connectivity_key": connectivity_key(inchikey),
                    "murcko_scaffold": desc["murcko_scaffold"],
                    "condition_bucket": condition_bucket,
                    "value_raw": nm_value,
                    "unit_raw": "nM",
                    "relation_raw": relation,
                    "type_raw": raw_column.removesuffix(" (nM)"),
                    "value_canonical": round(float(canonical), 6),
                    "unit_canonical": unit,
                    "confidence": 0.9 if relation == "=" else 0.75,
                    "target_chembl_id": target_id,
                    "target_label": target_label,
                    "assay_chembl_id": pubchem_aid,
                    "assay_description": f"BindingDB curated literature measurement: {target_label}",
                    "assay_type": "binding",
                    "assay_organism": target_organism,
                    "assay_cell_type": "",
                    "document_year": document_year,
                    "doi": doi,
                    "pubmed_id": pubmed_id,
                    **{key: desc.get(key, "") for key in FIELDNAMES if key in desc},
                }
                writer.writerow(out)
                written += 1
                endpoint_counts[endpoint] += 1
                counters[f"{unit}:observations_written"] += 1

    stats = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": {
            "name": "BindingDB",
            "source_version": args.source_version,
            "input": str(args.input),
        },
        "output": str(args.out),
        "num_input_rows": input_rows,
        "num_observations": written,
        "num_endpoints": len(endpoint_counts),
        "parameters": {
            "max_rows": args.max_rows,
            "max_rows_per_endpoint": args.max_rows_per_endpoint,
            "skip_rdkit_descriptors": args.skip_rdkit_descriptors,
            "compute_sa_score": args.compute_sa_score,
        },
        "measurement_counts": {key: counters[key] for key in sorted(counters) if key.endswith(":observations_written")},
        "endpoint_counts_top50": dict(endpoint_counts.most_common(50)),
        "counters": dict(sorted(counters.items())),
    }
    args.stats.write_text(json.dumps(stats, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return stats


def main() -> int:
    args = parse_args()
    stats = normalize(args)
    print(json.dumps({k: stats[k] for k in ["num_input_rows", "num_observations", "num_endpoints"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

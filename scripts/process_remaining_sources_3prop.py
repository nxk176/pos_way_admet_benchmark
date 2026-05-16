from __future__ import annotations

import csv
import html
import json
import re
import zipfile
from pathlib import Path
from typing import Any
from xml.etree import ElementTree as ET


ROOT = Path("pos_way_admet_benchmark")
RAW = ROOT / "raw" / "public"
OUT = ROOT / "data" / "remaining_sources_3prop_2pos"


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def parse_apache_index(path: Path, source: str) -> list[dict[str, str]]:
    text = path.read_text(encoding="utf-8", errors="ignore")
    rows = []
    pattern = re.compile(
        r'<a href="(?P<href>[^"]+)">(?P<name>[^<]+)</a>\s+'
        r'(?P<modified>\d{4}-\d{2}-\d{2}\s+\d{2}:\d{2})\s+'
        r'(?P<size>[^\s<]+)',
        re.IGNORECASE,
    )
    for match in pattern.finditer(text):
        rows.append(
            {
                "source": source,
                "name": html.unescape(match.group("name")),
                "href": html.unescape(match.group("href")),
                "last_modified": match.group("modified"),
                "size": match.group("size"),
                "local_file": str(path),
            }
        )
    return rows


def parse_links(path: Path, source: str) -> list[dict[str, str]]:
    text = path.read_text(encoding="utf-8", errors="ignore")
    rows = []
    for href, label in re.findall(r'<a[^>]+href=["\']([^"\']+)["\'][^>]*>(.*?)</a>', text, flags=re.I | re.S):
        clean_label = re.sub(r"<[^>]+>", " ", label)
        clean_label = re.sub(r"\s+", " ", html.unescape(clean_label)).strip()
        rows.append({"source": source, "href": html.unescape(href), "label": clean_label, "local_file": str(path)})
    return rows


def col_index(cell_ref: str) -> int:
    letters = re.sub(r"[^A-Z]", "", cell_ref.upper())
    value = 0
    for char in letters:
        value = value * 26 + (ord(char) - ord("A") + 1)
    return max(value - 1, 0)


def ns(tag: str) -> str:
    return f"{{http://schemas.openxmlformats.org/spreadsheetml/2006/main}}{tag}"


def rel_ns(tag: str) -> str:
    return f"{{http://schemas.openxmlformats.org/package/2006/relationships}}{tag}"


def read_shared_strings(zf: zipfile.ZipFile) -> list[str]:
    try:
        data = zf.read("xl/sharedStrings.xml")
    except KeyError:
        return []
    root = ET.fromstring(data)
    strings = []
    for si in root.findall(ns("si")):
        text_parts = [node.text or "" for node in si.iter() if node.tag == ns("t")]
        strings.append("".join(text_parts))
    return strings


def workbook_sheets(zf: zipfile.ZipFile) -> list[tuple[str, str]]:
    workbook = ET.fromstring(zf.read("xl/workbook.xml"))
    rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
    rel_targets = {
        rel.attrib["Id"]: rel.attrib["Target"]
        for rel in rels.findall(rel_ns("Relationship"))
        if "Id" in rel.attrib and "Target" in rel.attrib
    }
    sheets = []
    for sheet in workbook.findall(f"{ns('sheets')}/{ns('sheet')}"):
        name = sheet.attrib.get("name", "sheet")
        rid = sheet.attrib.get("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id", "")
        target = rel_targets.get(rid, "")
        if target:
            if not target.startswith("xl/"):
                target = "xl/" + target.lstrip("/")
            sheets.append((name, target))
    return sheets


def cell_value(cell: ET.Element, shared_strings: list[str]) -> str:
    cell_type = cell.attrib.get("t", "")
    if cell_type == "inlineStr":
        return "".join(node.text or "" for node in cell.iter() if node.tag == ns("t"))
    value_node = cell.find(ns("v"))
    if value_node is None or value_node.text is None:
        return ""
    raw = value_node.text
    if cell_type == "s":
        try:
            return shared_strings[int(raw)]
        except (ValueError, IndexError):
            return raw
    return raw


def sheet_rows(zf: zipfile.ZipFile, sheet_path: str, shared_strings: list[str]) -> list[list[str]]:
    root = ET.fromstring(zf.read(sheet_path))
    rows = []
    for row in root.findall(f"{ns('sheetData')}/{ns('row')}"):
        values: list[str] = []
        for cell in row.findall(ns("c")):
            index = col_index(cell.attrib.get("r", "A1"))
            while len(values) <= index:
                values.append("")
            values[index] = cell_value(cell, shared_strings)
        rows.append(values)
    return rows


def convert_xlsx(path: Path, out_dir: Path) -> list[dict[str, Any]]:
    summary = []
    with zipfile.ZipFile(path) as zf:
        shared_strings = read_shared_strings(zf)
        for sheet_name, sheet_path in workbook_sheets(zf):
            rows = sheet_rows(zf, sheet_path, shared_strings)
            safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", sheet_name).strip("_") or "sheet"
            out_path = out_dir / f"{path.stem}__{safe_name}.csv"
            out_path.parent.mkdir(parents=True, exist_ok=True)
            max_cols = max((len(row) for row in rows), default=0)
            with out_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle)
                for row in rows:
                    writer.writerow(row + [""] * (max_cols - len(row)))
            header = rows[0] if rows else []
            summary.append(
                {
                    "source": "toxcast",
                    "workbook": str(path),
                    "sheet": sheet_name,
                    "rows_including_header": len(rows),
                    "data_rows": max(len(rows) - 1, 0),
                    "columns": max_cols,
                    "header_json": json.dumps(header, ensure_ascii=False, separators=(",", ":")),
                    "output_csv": str(out_path),
                }
            )
    return summary


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)

    source_status = []

    pubchem_bioassay = parse_apache_index(RAW / "pubchem" / "bioassay_csv_index.html", "pubchem_bioassay_csv_index")
    write_csv(OUT / "pubchem_bioassay_csv_index.csv", pubchem_bioassay, ["source", "name", "href", "last_modified", "size", "local_file"])
    source_status.append(
        {
            "source": "PubChem BioAssay",
            "processed_output": "pubchem_bioassay_csv_index.csv",
            "processable_into_final_rows": "no",
            "reason": "Only archive index HTML is available locally; assay CSV archives and CID-to-structure/activity payloads are not downloaded.",
        }
    )

    pubchem_compound = parse_apache_index(RAW / "pubchem" / "compound_sdf_index.html", "pubchem_compound_sdf_index")
    write_csv(OUT / "pubchem_compound_sdf_index.csv", pubchem_compound, ["source", "name", "href", "last_modified", "size", "local_file"])
    source_status.append(
        {
            "source": "PubChem Compound",
            "processed_output": "pubchem_compound_sdf_index.csv",
            "processable_into_final_rows": "no",
            "reason": "Only compound archive index and an incomplete partial SDF transfer are available locally; compounds alone are not property ground truth.",
        }
    )

    tox21_links = parse_links(RAW / "tox21" / "tox21_public_data_page.html", "tox21_public_data_page")
    write_csv(OUT / "tox21_public_page_links.csv", tox21_links, ["source", "href", "label", "local_file"])
    tox21_assay_path = RAW / "tox21" / "tox21_assays.json"
    if tox21_assay_path.exists():
        assays = json.loads(tox21_assay_path.read_text(encoding="utf-8"))
        assay_rows = [
            {
                "source": "tox21_assays_api",
                "protocol_name": row.get("PROTOCOL_NAME", ""),
                "local_file": str(tox21_assay_path),
            }
            for row in assays
            if isinstance(row, dict)
        ]
        write_csv(OUT / "tox21_assays.csv", assay_rows, ["source", "protocol_name", "local_file"])
        tox21_output = "tox21_public_page_links.csv and tox21_assays.csv"
        tox21_reason = (
            "Assay list was downloaded, but replicate/aggregated result tables were not downloaded. "
            "The assay list alone is not molecule-level activity data."
        )
    else:
        tox21_output = "tox21_public_page_links.csv"
        tox21_reason = "Only the public page HTML is available locally; replicate/aggregated assay data files are not downloaded."
    source_status.append(
        {
            "source": "Tox21",
            "processed_output": tox21_output,
            "processable_into_final_rows": "no",
            "reason": tox21_reason,
        }
    )

    toxcast_summary = []
    toxcast_out = OUT / "toxcast_metadata_csv"
    for workbook in sorted((RAW / "toxcast").glob("*.xlsx")):
        toxcast_summary.extend(convert_xlsx(workbook, toxcast_out))
    write_csv(
        OUT / "toxcast_workbook_summary.csv",
        toxcast_summary,
        ["source", "workbook", "sheet", "rows_including_header", "data_rows", "columns", "header_json", "output_csv"],
    )
    source_status.append(
        {
            "source": "ToxCast/invitroDB v4.3",
            "processed_output": "toxcast_workbook_summary.csv and toxcast_metadata_csv/*.csv",
            "processable_into_final_rows": "no",
            "reason": "Available files are annotations, target mappings, cytotoxicity annotations, and analytical QC metadata. They do not include a complete molecule-level activity matrix with structures needed to create 2-positive/1-negative molecular edit rows.",
        }
    )

    source_status.extend(
        [
            {
                "source": "DrugBank",
                "processed_output": "",
                "processable_into_final_rows": "no",
                "reason": "Not downloaded because redistribution requires authenticated license/access.",
            },
            {
                "source": "eTOX",
                "processed_output": "",
                "processable_into_final_rows": "no",
                "reason": "Not downloaded because access is controlled/proprietary.",
            },
        ]
    )
    write_csv(
        OUT / "source_processing_status.csv",
        source_status,
        ["source", "processed_output", "processable_into_final_rows", "reason"],
    )
    (OUT / "README.md").write_text(
        "# Remaining Sources 3-Property/2-Positive Processing\n\n"
        "This folder contains what can be processed locally from the remaining non-ZINC sources.\n"
        "None of these sources can currently produce final 3-property, 2-positive molecular edit rows from the locally downloaded files alone.\n"
        "See `source_processing_status.csv` for per-source reasons.\n",
        encoding="utf-8",
    )
    print(json.dumps({"out_dir": str(OUT), "sources": source_status}, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

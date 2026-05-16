from __future__ import annotations

import argparse
import csv
import html
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export two separate HTML tables from CSV and Markdown sources.")
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/normalized_csv/property_summary.csv"),
    )
    parser.add_argument(
        "--selection-md",
        type=Path,
        default=Path("pos_way_admet_benchmark/PROPERTY_SELECTION_TABLE.md"),
    )
    parser.add_argument(
        "--summary-html",
        type=Path,
        default=Path("pos_way_admet_benchmark/PROPERTY_SUMMARY_TABLE.html"),
    )
    parser.add_argument(
        "--selection-html",
        type=Path,
        default=Path("pos_way_admet_benchmark/PROPERTY_SELECTION_TABLE_VIEW.html"),
    )
    return parser.parse_args()


CSS = """
:root {
  --bg: #f6f7f9;
  --panel: #ffffff;
  --text: #17202a;
  --muted: #5f6b7a;
  --line: #d9dee7;
  --head: #eef2f7;
  --blue: #1f5eff;
  --green: #147a3f;
  --amber: #9a5b00;
  --red: #a83232;
  --slate: #435161;
}
* { box-sizing: border-box; }
body {
  margin: 0;
  background: var(--bg);
  color: var(--text);
  font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  line-height: 1.45;
}
main {
  max-width: 1320px;
  margin: 0 auto;
  padding: 32px 24px 48px;
}
h1 { margin: 0 0 8px; font-size: 30px; letter-spacing: 0; }
h2 { margin: 28px 0 10px; font-size: 20px; letter-spacing: 0; }
p, li { color: var(--muted); }
code {
  background: #eef1f4;
  border-radius: 5px;
  padding: 1px 5px;
  color: #263445;
}
.subtitle { margin: 0 0 20px; max-width: 900px; }
.table-wrap {
  overflow-x: auto;
  background: var(--panel);
  border: 1px solid var(--line);
  border-radius: 8px;
  margin: 14px 0 22px;
}
table {
  border-collapse: collapse;
  width: 100%;
  min-width: 980px;
  font-size: 13px;
}
thead th {
  background: var(--head);
  color: #263445;
  text-align: left;
  padding: 11px 12px;
  border-bottom: 1px solid var(--line);
  white-space: nowrap;
}
tbody td {
  padding: 10px 12px;
  border-bottom: 1px solid #edf0f4;
  vertical-align: middle;
}
tbody tr:hover { background: #f9fbff; }
.num {
  text-align: right;
  font-variant-numeric: tabular-nums;
  white-space: nowrap;
}
.endpoint { font-weight: 700; white-space: nowrap; }
.badge {
  display: inline-flex;
  align-items: center;
  min-height: 24px;
  border-radius: 999px;
  padding: 3px 9px;
  font-weight: 650;
  white-space: nowrap;
}
.strong { color: var(--green); background: #e7f5ec; }
.usable { color: var(--amber); background: #fff4df; }
.low { color: var(--red); background: #fdecec; }
.proxy { color: var(--slate); background: #eef1f4; }
.source-note { margin-top: 18px; font-size: 13px; color: var(--muted); }
@media (max-width: 820px) {
  main { padding: 24px 14px 36px; }
  h1 { font-size: 24px; }
}
"""


NUMERIC_COLUMNS = {
    "Samples",
    "Molecules",
    "Scaffolds",
    "Buckets",
    "Median",
    "Min",
    "Max",
    "Value",
    "sample_count",
    "unique_molecules",
    "unique_connectivity_keys",
    "unique_scaffolds",
    "unique_condition_buckets",
    "median_value",
    "min_value",
    "max_value",
}


def format_number(value: str) -> str:
    try:
        parsed = float(value)
    except ValueError:
        return value
    if parsed.is_integer():
        return f"{int(parsed):,}"
    return f"{parsed:,.4g}"


def badge_class(text: str) -> str:
    lower = text.lower()
    if "strong" in lower:
        return "strong"
    if "usable" in lower:
        return "usable"
    if "lower" in lower:
        return "low"
    if "proxy" in lower or "rdkit" in lower:
        return "proxy"
    return "proxy"


def render_page(title: str, subtitle: str, body: str, source: str) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(title)}</title>
  <style>{CSS}</style>
</head>
<body>
  <main>
    <h1>{html.escape(title)}</h1>
    <p class="subtitle">{html.escape(subtitle)}</p>
    {body}
    <p class="source-note">Source: <code>{html.escape(source)}</code></p>
  </main>
</body>
</html>
"""


def render_html_table(headers: list[str], rows: list[list[str]]) -> str:
    head = "".join(f"<th>{html.escape(header)}</th>" for header in headers)
    body_rows = []
    for row in rows:
        cells = []
        for idx, value in enumerate(row):
            header = headers[idx] if idx < len(headers) else ""
            display = format_number(value) if header in NUMERIC_COLUMNS else value
            css = []
            if header in NUMERIC_COLUMNS:
                css.append("num")
            if header.lower() in {"endpoint", "endpoint_name"}:
                css.append("endpoint")
            if header.lower() == "recommendation" or header.lower() == "priority":
                display_html = f'<span class="badge {badge_class(display)}">{html.escape(display)}</span>'
            else:
                display_html = html.escape(display)
            class_attr = f' class="{" ".join(css)}"' if css else ""
            cells.append(f"<td{class_attr}>{display_html}</td>")
        body_rows.append(f"<tr>{''.join(cells)}</tr>")
    return f"""
<div class="table-wrap">
  <table>
    <thead><tr>{head}</tr></thead>
    <tbody>{''.join(body_rows)}</tbody>
  </table>
</div>
""".strip()


def export_summary_csv(csv_path: Path, out_path: Path) -> None:
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        headers = list(reader.fieldnames or [])
        rows = [[row.get(header, "") for header in headers] for row in reader]

    body = render_html_table(headers, rows)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        render_page(
            title="Property Summary CSV Table",
            subtitle="Direct HTML rendering of the normalized property_summary.csv file.",
            body=body,
            source=str(csv_path),
        ),
        encoding="utf-8",
    )


def parse_inline_markdown(text: str) -> str:
    escaped = html.escape(text)
    parts = escaped.split("`")
    if len(parts) == 1:
        return escaped
    rendered = []
    for idx, part in enumerate(parts):
        rendered.append(f"<code>{part}</code>" if idx % 2 else part)
    return "".join(rendered)


def parse_md_table(lines: list[str], start: int) -> tuple[str, int]:
    headers = [cell.strip().strip("`") for cell in lines[start].strip().strip("|").split("|")]
    rows: list[list[str]] = []
    idx = start + 2
    while idx < len(lines) and lines[idx].lstrip().startswith("|"):
        rows.append([cell.strip().strip("`") for cell in lines[idx].strip().strip("|").split("|")])
        idx += 1
    return render_html_table(headers, rows), idx


def export_selection_markdown(md_path: Path, out_path: Path) -> None:
    lines = md_path.read_text(encoding="utf-8").splitlines()
    body_parts: list[str] = []
    title = "Property Selection Markdown Table"
    idx = 0
    while idx < len(lines):
        line = lines[idx].strip()
        if not line:
            idx += 1
            continue
        if line.startswith("# "):
            title = line[2:].strip()
            idx += 1
            continue
        if line.startswith("## "):
            body_parts.append(f"<h2>{html.escape(line[3:].strip())}</h2>")
            idx += 1
            continue
        if line.startswith("|") and idx + 1 < len(lines) and lines[idx + 1].lstrip().startswith("|---"):
            table_html, idx = parse_md_table(lines, idx)
            body_parts.append(table_html)
            continue
        if line.startswith("- "):
            items = []
            while idx < len(lines) and lines[idx].strip().startswith("- "):
                items.append(f"<li>{parse_inline_markdown(lines[idx].strip()[2:])}</li>")
                idx += 1
            body_parts.append(f"<ul>{''.join(items)}</ul>")
            continue
        body_parts.append(f"<p>{parse_inline_markdown(line)}</p>")
        idx += 1

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        render_page(
            title=title,
            subtitle="HTML rendering of PROPERTY_SELECTION_TABLE.md, including summary and recommendation tables.",
            body="\n".join(body_parts),
            source=str(md_path),
        ),
        encoding="utf-8",
    )


def main() -> int:
    args = parse_args()
    export_summary_csv(args.summary_csv, args.summary_html)
    export_selection_markdown(args.selection_md, args.selection_html)
    print(f"Wrote {args.summary_html}")
    print(f"Wrote {args.selection_html}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

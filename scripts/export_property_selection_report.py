from __future__ import annotations

import argparse
import csv
import html
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export a polished property-selection HTML report from property_summary.csv.")
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path("pos_way_admet_benchmark/data/normalized_csv/property_summary.csv"),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("pos_way_admet_benchmark/PROPERTY_SELECTION_REPORT.html"),
    )
    parser.add_argument(
        "--out-md",
        type=Path,
        default=Path("pos_way_admet_benchmark/PROPERTY_SELECTION_TABLE.md"),
    )
    return parser.parse_args()


def fmt_int(value: str) -> str:
    try:
        return f"{int(float(value)):,}"
    except ValueError:
        return value


def fmt_float(value: str) -> str:
    try:
        parsed = float(value)
    except ValueError:
        return value
    return f"{parsed:,.4g}"


def badge_class(recommendation: str) -> str:
    text = recommendation.lower()
    if "strong" in text:
        return "strong"
    if "usable" in text:
        return "usable"
    if "lower" in text:
        return "low"
    return "proxy"


def tier_class(tier: str) -> str:
    return "experimental" if tier == "experimental" else "proxy"


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return sorted(
        rows,
        key=lambda row: (
            0 if row["source_tier"] == "experimental" else 1,
            -int(float(row["sample_count"])),
            row["endpoint_name"],
        ),
    )


def render(rows: list[dict[str, str]]) -> str:
    experimental = [row for row in rows if row["source_tier"] == "experimental"]
    proxy = [row for row in rows if row["source_tier"] == "proxy"]
    total_samples = sum(int(float(row["sample_count"])) for row in rows)
    total_molecules = sum(int(float(row["unique_molecules"])) for row in experimental)
    strong = sum(1 for row in experimental if "Strong" in row["recommendation"])

    table_rows = []
    for row in rows:
        recommendation = html.escape(row["recommendation"])
        table_rows.append(
            f"""
            <tr class="{tier_class(row['source_tier'])}">
              <td class="endpoint">{html.escape(row['endpoint_name'])}</td>
              <td><span class="tier {tier_class(row['source_tier'])}">{html.escape(row['source_tier'])}</span></td>
              <td>{html.escape(row['property_family'])}</td>
              <td class="num">{fmt_int(row['sample_count'])}</td>
              <td class="num">{fmt_int(row['unique_molecules'])}</td>
              <td class="num">{fmt_int(row['unique_scaffolds'])}</td>
              <td class="num">{fmt_int(row['unique_condition_buckets'])}</td>
              <td>{html.escape(row['unit_canonical'])}</td>
              <td class="num">{fmt_float(row['median_value'])}</td>
              <td class="num">{fmt_float(row['min_value'])}</td>
              <td class="num">{fmt_float(row['max_value'])}</td>
              <td><span class="badge {badge_class(row['recommendation'])}">{recommendation}</span></td>
            </tr>
            """.strip()
        )

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Property Selection Report</title>
  <style>
    :root {{
      --bg: #f6f7f9;
      --panel: #ffffff;
      --text: #17202a;
      --muted: #5f6b7a;
      --line: #d9dee7;
      --blue: #1f5eff;
      --green: #147a3f;
      --amber: #9a5b00;
      --red: #a83232;
      --slate: #435161;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      background: var(--bg);
      color: var(--text);
      font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      line-height: 1.45;
    }}
    main {{
      max-width: 1320px;
      margin: 0 auto;
      padding: 32px 24px 48px;
    }}
    header {{
      margin-bottom: 20px;
    }}
    h1 {{
      margin: 0 0 8px;
      font-size: 30px;
      letter-spacing: 0;
    }}
    .subtitle {{
      margin: 0;
      color: var(--muted);
      max-width: 900px;
    }}
    .cards {{
      display: grid;
      grid-template-columns: repeat(4, minmax(0, 1fr));
      gap: 12px;
      margin: 22px 0;
    }}
    .card {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 14px 16px;
    }}
    .label {{
      color: var(--muted);
      font-size: 13px;
      margin-bottom: 6px;
    }}
    .value {{
      font-size: 24px;
      font-weight: 700;
    }}
    .table-wrap {{
      overflow-x: auto;
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
    }}
    table {{
      border-collapse: collapse;
      width: 100%;
      min-width: 1180px;
      font-size: 13px;
    }}
    thead th {{
      position: sticky;
      top: 0;
      background: #eef2f7;
      color: #263445;
      text-align: left;
      padding: 11px 12px;
      border-bottom: 1px solid var(--line);
      white-space: nowrap;
    }}
    tbody td {{
      padding: 10px 12px;
      border-bottom: 1px solid #edf0f4;
      vertical-align: middle;
    }}
    tbody tr:hover {{
      background: #f9fbff;
    }}
    .num {{
      text-align: right;
      font-variant-numeric: tabular-nums;
      white-space: nowrap;
    }}
    .endpoint {{
      font-weight: 700;
      white-space: nowrap;
    }}
    .tier,
    .badge {{
      display: inline-flex;
      align-items: center;
      min-height: 24px;
      border-radius: 999px;
      padding: 3px 9px;
      font-weight: 650;
      white-space: nowrap;
    }}
    .tier.experimental {{ color: var(--blue); background: #eaf0ff; }}
    .tier.proxy {{ color: var(--slate); background: #edf0f3; }}
    .badge.strong {{ color: var(--green); background: #e7f5ec; }}
    .badge.usable {{ color: var(--amber); background: #fff4df; }}
    .badge.low {{ color: var(--red); background: #fdecec; }}
    .badge.proxy {{ color: var(--slate); background: #eef1f4; }}
    .notes {{
      margin-top: 16px;
      color: var(--muted);
      font-size: 13px;
    }}
    @media (max-width: 820px) {{
      main {{ padding: 24px 14px 36px; }}
      .cards {{ grid-template-columns: repeat(2, minmax(0, 1fr)); }}
      h1 {{ font-size: 24px; }}
    }}
  </style>
</head>
<body>
  <main>
    <header>
      <h1>Property Selection Report</h1>
      <p class="subtitle">Normalized endpoint-level summary for deciding which properties should become paper-facing benchmark targets.</p>
    </header>

    <section class="cards" aria-label="Summary cards">
      <div class="card"><div class="label">Total endpoints</div><div class="value">{len(rows)}</div></div>
      <div class="card"><div class="label">Experimental endpoints</div><div class="value">{len(experimental)}</div></div>
      <div class="card"><div class="label">Strong experimental candidates</div><div class="value">{strong}</div></div>
      <div class="card"><div class="label">Experimental molecule rows</div><div class="value">{total_molecules:,}</div></div>
    </section>

    <section class="table-wrap" aria-label="Property selection table">
      <table>
        <thead>
          <tr>
            <th>Endpoint</th>
            <th>Tier</th>
            <th>Family</th>
            <th>Samples</th>
            <th>Molecules</th>
            <th>Scaffolds</th>
            <th>Buckets</th>
            <th>Unit</th>
            <th>Median</th>
            <th>Min</th>
            <th>Max</th>
            <th>Recommendation</th>
          </tr>
        </thead>
        <tbody>
          {''.join(table_rows)}
        </tbody>
      </table>
    </section>

    <p class="notes">
      Source: <code>data/normalized_csv/property_summary.csv</code>. Proxy RDKit descriptors are included as secondary constraints or silver signals;
      experimental endpoints are the primary candidates for gold benchmark targets.
    </p>
  </main>
</body>
</html>
"""


def render_markdown(rows: list[dict[str, str]]) -> str:
    experimental = [row for row in rows if row["source_tier"] == "experimental"]
    proxy = [row for row in rows if row["source_tier"] == "proxy"]
    strong = [row for row in experimental if "Strong" in row["recommendation"]]
    usable = [row for row in experimental if "Usable" in row["recommendation"]]
    lower = [row for row in experimental if "Lower" in row["recommendation"]]

    def table_line(row: dict[str, str]) -> str:
        return (
            f"| `{row['endpoint_name']}` | {row['source_tier']} | {row['property_family']} | "
            f"{fmt_int(row['sample_count'])} | {fmt_int(row['unique_molecules'])} | "
            f"{fmt_int(row['unique_scaffolds'])} | {fmt_int(row['unique_condition_buckets'])} | "
            f"`{row['unit_canonical']}` | {fmt_float(row['median_value'])} | "
            f"{fmt_float(row['min_value'])} | {fmt_float(row['max_value'])} | {row['recommendation']} |"
        )

    def endpoint_list(items: list[dict[str, str]]) -> str:
        return ", ".join(f"`{item['endpoint_name']}`" for item in items) if items else "None"

    lines = [
        "# Property Selection Table",
        "",
        "This table is generated from `data/normalized_csv/property_summary.csv` and is intended for deciding which properties should become benchmark targets.",
        "",
        "## Summary",
        "",
        f"| Metric | Value |",
        "|---|---:|",
        f"| Total endpoints | {len(rows)} |",
        f"| Experimental endpoints | {len(experimental)} |",
        f"| Proxy/RDKit endpoints | {len(proxy)} |",
        f"| Strong experimental candidates | {len(strong)} |",
        "",
        "## Recommended Experimental Targets",
        "",
        "| Priority | Endpoints | Notes |",
        "|---|---|---|",
        f"| Strong | {endpoint_list(strong)} | Best first-pass candidates based on coverage. |",
        f"| Usable | {endpoint_list(usable)} | Keep if the final task needs these endpoints; inspect assay buckets. |",
        f"| Lower coverage | {endpoint_list(lower)} | Use selectively or expand data before making it a main endpoint. |",
        "",
        "## Full Endpoint Table",
        "",
        "| Endpoint | Tier | Family | Samples | Molecules | Scaffolds | Buckets | Unit | Median | Min | Max | Recommendation |",
        "|---|---|---|---:|---:|---:|---:|---|---:|---:|---:|---|",
    ]
    lines.extend(table_line(row) for row in rows)
    lines.extend(
        [
            "",
            "## Reading Notes",
            "",
            "- `Samples`: number of normalized property observations.",
            "- `Molecules`: unique molecule count for that endpoint.",
            "- `Scaffolds`: unique Murcko scaffold count; this is a coverage statistic, not the fragment-matching rule.",
            "- `Buckets`: unique condition/assay buckets. More buckets usually means more assay heterogeneity to inspect.",
            "- Proxy RDKit endpoints are useful as secondary constraints or silver training signals, not as gold experimental endpoints.",
            "- The fragment multi-property dataset uses shared BRICS fragments rather than exact scaffold equality.",
            "",
            "HTML view: `PROPERTY_SELECTION_REPORT.html`.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    args = parse_args()
    rows = load_rows(args.summary)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(render(rows), encoding="utf-8")
    args.out_md.parent.mkdir(parents=True, exist_ok=True)
    args.out_md.write_text(render_markdown(rows), encoding="utf-8")
    print(f"Wrote {args.out}")
    print(f"Wrote {args.out_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""HTML, Markdown, and PDF report rendering utilities for BRD."""

from __future__ import annotations

import html
from pathlib import Path
from typing import Sequence
import pandas as pd


def df_to_html_table(df: pd.DataFrame, max_rows: int = 100) -> str:
    sub = df.head(max_rows)
    return sub.to_html(index=False, classes="brd-table", border=0, escape=True)


def wrap_html_report(title: str, subtitle: str, body_html: str) -> str:
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<title>{html.escape(title)}</title>
<style>
  body {{ font-family: 'Segoe UI', Arial, sans-serif; margin: 2rem auto; max-width: 1200px; color: #1e293b; line-height: 1.55; padding: 0 1.2rem; }}
  h1 {{ color: #0f172a; border-bottom: 3px solid #1d4ed8; padding-bottom: 0.4rem; }}
  h2 {{ color: #1e3a8a; margin-top: 2rem; border-bottom: 1px solid #cbd5e1; padding-bottom: 0.25rem; }}
  h3 {{ color: #334155; margin-top: 1.3rem; }}
  .subtitle {{ color: #475569; font-size: 1.05rem; margin-bottom: 1.5rem; }}
  .callout {{ background: #eff6ff; border-left: 4px solid #2563eb; padding: 0.85rem 1rem; margin: 1rem 0; border-radius: 4px; }}
  .warn {{ background: #fffbeb; border-left: 4px solid #d97706; padding: 0.85rem 1rem; margin: 1rem 0; border-radius: 4px; }}
  table.brd-table {{ border-collapse: collapse; width: 100%; margin: 1rem 0; font-size: 0.88rem; }}
  table.brd-table th {{ background: #1e3a8a; color: #ffffff; text-align: left; padding: 0.5rem 0.65rem; }}
  table.brd-table td {{ border-bottom: 1px solid #e2e8f0; padding: 0.45rem 0.65rem; vertical-align: top; }}
  table.brd-table tr:nth-child(even) {{ background: #f8fafc; }}
  code, pre {{ background: #f1f5f9; padding: 0.15rem 0.35rem; border-radius: 3px; font-family: Consolas, monospace; font-size: 0.88em; }}
  figure {{ margin: 1.5rem 0; text-align: center; }}
  figure img {{ max-width: 100%; border: 1px solid #cbd5e1; border-radius: 4px; }}
  figcaption {{ font-size: 0.88rem; color: #475569; margin-top: 0.4rem; }}
</style>
</head>
<body>
<h1>{html.escape(title)}</h1>
<div class="subtitle">{html.escape(subtitle)}</div>
{body_html}
</body>
</html>
"""


def write_text(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path

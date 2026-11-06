"""Supporting-Information fitting tables from result.json records.

One row per scattering path. Output formats: CSV, LaTeX (booktabs), Word
(.docx, needs the optional ``python-docx``).
"""

from __future__ import annotations

import csv
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path

from .uncertainty import format_value

COLUMNS = [
    ("sample", "Sample"),
    ("path", "Path"),
    ("cn", "CN"),
    ("r", "R (Å)"),
    ("sigma2", "σ² (10⁻³ Å²)"),
    ("e0", "ΔE₀ (eV)"),
    ("s02", "S₀²"),
    ("r_factor", "R-factor"),
    ("chi2_reduced", "Reduced χ²"),
]

LATEX_HEAD = {
    "sample": "Sample",
    "path": "Path",
    "cn": "CN",
    "r": r"$R$ (\AA)",
    "sigma2": r"$\sigma^2$ ($10^{-3}$ \AA$^2$)",
    "e0": r"$\Delta E_0$ (eV)",
    "s02": r"$S_0^2$",
    "r_factor": r"$R$-factor",
    "chi2_reduced": r"$\chi^2_\nu$",
}


@dataclass
class TableRow:
    cells: dict[str, str]
    notes: list[str] = field(default_factory=list)


def _fmt(rec: dict | None, style: str, scale: float = 1.0, fixed_mark: str = "*") -> tuple[str, bool]:
    """Format a {value, stderr, expr} record. Returns (text, is_fixed_or_constrained)."""
    if rec is None or rec.get("value") is None:
        return "–", False
    v = rec["value"] * scale
    s = rec.get("stderr")
    if s is None or s <= 0:
        return f"{v:.3g}{fixed_mark}", True
    return format_value(v, s * scale, style=style), False


def rows_from_results(results: Iterable[dict], style: str = "paren") -> list[TableRow]:
    """Build table rows from schema-v1 result dicts."""
    rows: list[TableRow] = []
    for res in results:
        fit = res["fit"]
        first = True
        s02_ref = res.get("s02_ref")
        for ds in res["datasets"]:
            for p in ds["paths"]:
                pp = p.get("params", {})
                notes = []
                cn, cn_fixed = _fmt(p.get("cn"), style)
                r, _ = _fmt(p.get("r"), style)
                ss, _ = _fmt(pp.get("sigma2"), style, scale=1e3)
                e0, _ = _fmt(pp.get("e0"), style)
                if s02_ref:
                    s02 = f"{s02_ref['value']:.2f}†"
                    notes.append("† S₀² fixed from a reference compound measured at the same edge and beamline.")
                else:
                    s02, _ = _fmt(pp.get("s02"), style)
                if cn_fixed:
                    notes.append("* fixed or constrained, not refined independently.")
                cells = {
                    "sample": res.get("sample", "") if first else "",
                    "path": p.get("label", ""),
                    "cn": cn, "r": r, "sigma2": ss, "e0": e0, "s02": s02,
                    "r_factor": f"{fit['r_factor']:.4f}" if first and fit.get("r_factor") is not None else "",
                    "chi2_reduced": f"{fit['chi2_reduced']:.1f}" if first and fit.get("chi2_reduced") is not None else "",
                }
                rows.append(TableRow(cells, notes))
                first = False
    return rows


def _footnotes(rows: list[TableRow]) -> list[str]:
    seen: list[str] = []
    for r in rows:
        for n in r.notes:
            if n not in seen:
                seen.append(n)
    seen.append("Uncertainties are 1σ from the fit covariance (Larch), given in parentheses on the last digit(s).")
    return seen


def to_csv(rows: list[TableRow], path: str | Path) -> Path:
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow([h for _, h in COLUMNS])
        for r in rows:
            w.writerow([r.cells.get(k, "") for k, _ in COLUMNS])
    return path


def _tex_escape(s: str) -> str:
    rep = {"&": r"\&", "%": r"\%", "_": r"\_", "#": r"\#", "σ": r"$\sigma$", "†": r"$^\dagger$",
           "±": r"$\pm$", "–": "--", "*": r"$^*$", "χ": r"$\chi$", "Å": r"\AA{}", "₀": r"$_0$", "²": r"$^2$",
           "⁻³": r"$^{-3}$", "Δ": r"$\Delta$"}
    for a, b in rep.items():
        s = s.replace(a, b)
    return s


def to_latex(rows: list[TableRow], caption: str = "EXAFS fitting results.", label: str = "tab:exafs") -> str:
    keys = [k for k, _ in COLUMNS]
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        rf"\caption{{{_tex_escape(caption)}}}",
        rf"\label{{{label}}}",
        r"\begin{tabular}{" + "l" * 2 + "c" * (len(keys) - 2) + "}",
        r"\toprule",
        " & ".join(LATEX_HEAD[k] for k in keys) + r" \\",
        r"\midrule",
    ]
    for r in rows:
        lines.append(" & ".join(_tex_escape(r.cells.get(k, "")) for k in keys) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    for n in _footnotes(rows):
        lines.append(r"\par\footnotesize " + _tex_escape(n))
    lines.append(r"\end{table}")
    return "\n".join(lines) + "\n"


def to_docx(rows: list[TableRow], path: str | Path, caption: str = "EXAFS fitting results.") -> Path:
    try:
        import docx
        from docx.shared import Pt
    except ImportError as exc:  # pragma: no cover
        raise ImportError("Word export needs python-docx (pip install 'xastk[docx]')") from exc
    document = docx.Document()
    document.add_paragraph(caption)
    table = document.add_table(rows=1, cols=len(COLUMNS))
    table.style = "Table Grid"
    for cell, (_, head) in zip(table.rows[0].cells, COLUMNS, strict=True):
        cell.text = head
    for r in rows:
        cells = table.add_row().cells
        for cell, (k, _) in zip(cells, COLUMNS, strict=True):
            cell.text = r.cells.get(k, "")
    for n in _footnotes(rows):
        para = document.add_paragraph(n)
        for run in para.runs:
            run.font.size = Pt(8)
    path = Path(path)
    document.save(path)
    return path

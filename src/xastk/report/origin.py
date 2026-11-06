"""OriginPro bridge, split into an OS-independent bundle and a Windows runner.

1. :func:`export_bundle` (any OS, safe in CI) writes one CSV per worksheet and a
   ``figure.json`` describing layers, columns, axis labels and templates.
2. :func:`run_bundle` (Windows + Origin 2021 or later + ``pip install originpro``)
   rebuilds the workbook and graph inside Origin and saves an ``.opju`` project
   that co-authors can fine-tune.

Status: ``export_bundle`` is tested. ``run_bundle`` follows the documented
``originpro`` API (new_sheet, from_df, new_graph, add_plot, save) but has not
yet been exercised on a Windows/Origin machine -- see docs/technical-report.md.
"""

from __future__ import annotations

import csv
import json
from collections.abc import Mapping, Sequence
from pathlib import Path

BUNDLE_VERSION = 1


def export_bundle(
    outdir: str | Path,
    sheets: Mapping[str, Mapping[str, Sequence[float]]],
    layers: Sequence[Mapping],
    title: str = "",
    template: str | None = None,
) -> Path:
    """Write a figure bundle.

    Parameters
    ----------
    sheets : ``{sheet_name: {column_name: values}}``; columns of one sheet must
        have equal length. The first column is the X column by convention.
    layers : list of layer specs, e.g.::

        {"sheet": "chik", "x": "k", "y": ["data", "fit"], "xlabel": "k (1/A)",
         "ylabel": "k^2 chi(k)", "offset": 0.0}

    template : name of an Origin graph template (.otpu) to use, if any
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for name, cols in sheets.items():
        names = list(cols)
        n = {len(cols[c]) for c in names}
        if len(n) != 1:
            raise ValueError(f"sheet {name!r}: columns have different lengths {sorted(n)}")
        with open(outdir / f"{name}.csv", "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(names)
            for row in zip(*(cols[c] for c in names), strict=True):
                w.writerow([f"{float(v):.10g}" for v in row])
    for layer in layers:
        if layer["sheet"] not in sheets:
            raise KeyError(f"layer refers to unknown sheet {layer['sheet']!r}")
        for col in [layer["x"], *layer["y"]]:
            if col not in sheets[layer["sheet"]]:
                raise KeyError(f"layer refers to unknown column {col!r} in sheet {layer['sheet']!r}")
    spec = {"bundle_version": BUNDLE_VERSION, "title": title, "template": template,
            "sheets": list(sheets), "layers": list(layers)}
    (outdir / "figure.json").write_text(json.dumps(spec, indent=2), encoding="utf-8")
    return outdir


def run_bundle(bundle_dir: str | Path, out_opju: str | Path | None = None,
               template_dir: str | Path | None = None, visible: bool = False):  # pragma: no cover
    """Build an Origin project from a bundle (Windows only, untested in CI)."""
    try:
        import originpro as op
        import pandas as pd
    except ImportError as exc:
        raise ImportError("run_bundle needs Windows, Origin 2021+, originpro and pandas") from exc

    bundle_dir = Path(bundle_dir)
    spec = json.loads((bundle_dir / "figure.json").read_text(encoding="utf-8"))
    op.set_show(visible)
    books = {}
    for name in spec["sheets"]:
        df = pd.read_csv(bundle_dir / f"{name}.csv")
        wks = op.new_sheet("w", lname=name)
        wks.from_df(df)
        books[name] = (wks, list(df.columns))
    template = spec.get("template")
    if template and template_dir:
        template = str(Path(template_dir) / template)
    graph = op.new_graph(template=template) if template else op.new_graph()
    for i, layer in enumerate(spec["layers"]):
        gl = graph[i] if i < len(graph) else graph.add_layer()
        wks, cols = books[layer["sheet"]]
        xi = cols.index(layer["x"])
        for ycol in layer["y"]:
            gl.add_plot(wks, cols.index(ycol), xi)
        gl.rescale()
        if layer.get("xlabel"):
            gl.axis("x").title = layer["xlabel"]
        if layer.get("ylabel"):
            gl.axis("y").title = layer["ylabel"]
    if out_opju is None:
        out_opju = bundle_dir / "figure.opju"
    op.save(str(Path(out_opju).resolve()))
    return Path(out_opju)

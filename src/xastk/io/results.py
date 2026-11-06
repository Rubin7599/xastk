"""The ``result.json`` record: one versioned file per fit, the single source
for SI tables, figures, sanity checks and the team audit log.

Schema: ``xastk/schemas/result-v1.schema.json`` (JSON Schema 2020-12).
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import platform
from collections.abc import Iterable
from pathlib import Path
from typing import Any

SCHEMA_ID = "xastk-result/1"
PATH_PARAMS = ("s02", "e0", "deltar", "sigma2", "third", "fourth", "ei")


def _vs(par) -> dict:
    if par is None:
        return {"value": None, "stderr": None}
    std = getattr(par, "stderr", None)
    out = {"value": float(par.value), "stderr": None if std is None else float(std)}
    expr = getattr(par, "expr", None)
    if expr:
        out["expr"] = str(expr)
    return out


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def software_versions() -> dict:
    from .. import __version__

    out = {"xastk": __version__, "python": platform.python_version()}
    for mod in ("larch", "numpy", "scipy", "lmfit"):
        try:
            m = __import__(mod)
            out[mod] = getattr(m, "__version__", "unknown")
        except ImportError:
            pass
    return out


def result_from_feffit(
    out,
    sample: str = "",
    s02_ref: tuple[float, float] | None = None,
    inputs: Iterable[str | Path] = (),
    notes: str = "",
    min_correl: float = 0.1,
) -> dict[str, Any]:
    """Convert a Larch ``feffit`` result into a schema-v1 dictionary.

    Parameters
    ----------
    out : return value of ``larch.xafs.feffit``
    sample : sample label for tables
    s02_ref : optional (value, stderr) of S0^2 from a reference compound. When
        given, each path's coordination number is reported as
        ``CN = degen * s02_path / S0^2_ref`` (the usual convention when the path
        amplitude is set to ``N * S0^2``); otherwise CN = degen and the fitted
        amplitude is reported separately.
    inputs : data/FEFF files whose SHA-256 is recorded for provenance
    """
    params = out.params
    variables = {
        n: {"value": float(params[n].value),
            "stderr": None if params[n].stderr is None else float(params[n].stderr),
            "init": None if getattr(params[n], "init_value", None) is None else float(params[n].init_value)}
        for n in out.var_names
    }
    correl: dict[str, float] = {}
    for n in out.var_names:
        for m, c in (getattr(params[n], "correl", None) or {}).items():
            key = "|".join(sorted((n, m)))
            if abs(c) >= min_correl and key not in correl:
                correl[key] = float(c)

    datasets = []
    for ds in out.datasets:
        tr = ds.transform
        paths = []
        for label, path in ds.paths.items():
            pp = {}
            for pn in PATH_PARAMS:
                par = path.params.get(path.pathpar_name(pn)) if path.params is not None else None
                if par is not None:
                    pp[pn] = _vs(par)
            dr = pp.get("deltar", {"value": 0.0, "stderr": None})
            r = {"value": float(path.reff + (dr["value"] or 0.0)), "stderr": dr.get("stderr")}
            amp = pp.get("s02", {"value": 1.0, "stderr": None})
            if s02_ref is not None:
                ref, ref_err = float(s02_ref[0]), float(s02_ref[1])
                cn = float(path.degen) * amp["value"] / ref
                rel = ((amp["stderr"] or 0.0) / amp["value"]) ** 2 + (ref_err / ref) ** 2 if amp["value"] else 0.0
                cn_rec = {"value": cn, "stderr": abs(cn) * rel ** 0.5, "convention": "degen*s02/s02_ref"}
            else:
                cn_rec = {"value": float(path.degen), "stderr": None, "convention": "degen (fixed)"}
            geom = getattr(path._feffdat, "geom", None) if hasattr(path, "_feffdat") else None
            paths.append({
                "label": str(label),
                "filename": str(getattr(path, "filename", "")),
                "degen": float(path.degen),
                "reff": float(path.reff),
                "nleg": int(path.nleg),
                "geometry": "-".join(str(g[0]) for g in geom) if geom else "",
                "params": pp,
                "r": r,
                "cn": cn_rec,
            })
        datasets.append({
            "label": str(getattr(ds, "hashkey", "")),
            "transform": {
                "kmin": float(tr.kmin), "kmax": float(tr.kmax),
                "rmin": float(tr.rmin), "rmax": float(tr.rmax),
                "kweight": [float(x) for x in (tr.kweight if hasattr(tr.kweight, "__iter__") else [tr.kweight])],
                "dk": float(tr.dk), "window": str(tr.window), "fitspace": str(getattr(tr, "fitspace", "r")),
            },
            "n_idp": float(ds.n_idp),
            "paths": paths,
        })

    def _f(name):
        v = getattr(out, name, None)
        return None if v is None else float(v)

    return {
        "schema": SCHEMA_ID,
        "sample": sample,
        "created": _dt.datetime.now(_dt.UTC).isoformat(timespec="seconds"),
        "software": software_versions(),
        "fit": {
            "n_idp": _f("n_independent"), "n_var": int(out.nvarys), "n_data": int(getattr(out, "ndata", 0)),
            "chi_square": _f("chi_square"), "chi2_reduced": _f("chi2_reduced"),
            "r_factor": _f("rfactor"), "aic": _f("aic"), "bic": _f("bic"),
        },
        "s02_ref": None if s02_ref is None else {"value": float(s02_ref[0]), "stderr": float(s02_ref[1])},
        "variables": variables,
        "correlations": correl,
        "datasets": datasets,
        "provenance": {"inputs": [{"path": str(p), "sha256": sha256_file(p)} for p in inputs]},
        "notes": notes,
    }


def save_result(result: dict, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return path


def load_result(path: str | Path) -> dict:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if data.get("schema") != SCHEMA_ID:
        raise ValueError(f"{path}: not an {SCHEMA_ID} record (schema={data.get('schema')!r})")
    return data


def schema() -> dict:
    """Return the JSON Schema for result records."""
    from importlib.resources import files

    return json.loads(files("xastk.schemas").joinpath("result-v1.schema.json").read_text(encoding="utf-8"))


def validate_result(result: dict) -> None:
    """Validate against the JSON Schema (needs the optional ``jsonschema`` package)."""
    import jsonschema

    jsonschema.validate(result, schema())

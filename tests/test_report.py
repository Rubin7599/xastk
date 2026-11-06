import json

import pytest

from xastk.io import validate_result
from xastk.reduce import Rules, check_result
from xastk.report import export_bundle, rows_from_results, to_csv, to_docx, to_latex


def test_schema_accepts_synthetic(synthetic_result):
    pytest.importorskip("jsonschema")
    validate_result(synthetic_result)


def test_rows_and_formats(synthetic_result, tmp_path):
    rows = rows_from_results([synthetic_result])
    assert len(rows) == 2
    first, second = rows[0].cells, rows[1].cells
    assert first["sample"] == "Fe-N-C" and second["sample"] == ""
    assert first["cn"] == "4.0(5)"
    assert first["r"] == "1.970(12)"
    assert first["sigma2"] == "7.1(12)"
    assert first["s02"] == "0.85†"
    assert second["cn"] == "8*"              # fixed CN marked
    assert first["r_factor"] == "0.0123"

    csv_path = to_csv(rows, tmp_path / "t.csv")
    assert "Fe-N" in csv_path.read_text(encoding="utf-8")

    tex = to_latex(rows)
    assert r"\toprule" in tex and r"\bottomrule" in tex
    assert "4.0(5)" in tex and r"$^\dagger$" in tex

    pytest.importorskip("docx")
    docx_path = to_docx(rows, tmp_path / "t.docx")
    assert docx_path.stat().st_size > 1000


def test_sanity_flags(synthetic_result):
    codes = {f.code for f in check_result(synthetic_result)}
    assert "correlation" in codes           # 0.93 > 0.90
    assert "sigma2_range" in codes          # Fe-C sigma2 = 0.0015 < 0.002
    assert "nidp_budget" not in codes       # 5 <= 2/3 * 12.4
    flags = check_result(synthetic_result, Rules(r_reference={"Fe-N": 2.10}))
    assert any(f.code == "r_reference" for f in flags)


def test_origin_bundle(tmp_path):
    out = export_bundle(
        tmp_path / "bundle",
        sheets={"chik": {"k": [3.0, 3.05, 3.1], "data": [0.1, 0.2, 0.1], "fit": [0.11, 0.19, 0.12]}},
        layers=[{"sheet": "chik", "x": "k", "y": ["data", "fit"], "xlabel": "k", "ylabel": "k2chi"}],
        title="Cu foil",
    )
    spec = json.loads((out / "figure.json").read_text())
    assert spec["sheets"] == ["chik"] and spec["layers"][0]["y"] == ["data", "fit"]
    assert (out / "chik.csv").read_text().splitlines()[0] == "k,data,fit"
    with pytest.raises(KeyError):
        export_bundle(tmp_path / "bad", sheets={"a": {"x": [1], "y": [2]}},
                      layers=[{"sheet": "a", "x": "x", "y": ["z"]}])
    with pytest.raises(ValueError):
        export_bundle(tmp_path / "bad2", sheets={"a": {"x": [1, 2], "y": [2]}}, layers=[])

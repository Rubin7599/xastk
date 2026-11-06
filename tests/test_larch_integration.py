"""Golden-value tests against the Larch Cu-foil example.

These pin both our own numbers and the parts of Larch's API we depend on
(feffit result attributes, path parameter names, the internal residual).
"""

import numpy as np
import pytest
from conftest import DATA, needs_larch

pytestmark = [needs_larch, pytest.mark.larch]


def test_cu_golden_values(cu_result):
    p = cu_result.params
    assert cu_result.n_independent == pytest.approx(15.26, abs=0.01)
    assert p["amp"].value == pytest.approx(0.931, abs=0.005)
    assert p["sig2"].value == pytest.approx(0.00867, abs=0.0001)
    assert p["del_r"].value == pytest.approx(-0.0060, abs=0.001)
    assert p["del_e0"].value == pytest.approx(4.36, abs=0.1)


def test_result_record_validates(cu_result, tmp_path):
    from xastk.io import load_result, result_from_feffit, save_result, validate_result

    rec = result_from_feffit(cu_result, sample="Cu foil", inputs=[DATA / "cu_metal_rt.xdi"])
    pytest.importorskip("jsonschema")
    validate_result(rec)
    path = rec["datasets"][0]["paths"][0]
    assert path["label"] == "Cu-Cu1"
    assert path["degen"] == 12.0
    assert path["r"]["value"] == pytest.approx(2.5478 - 0.0060, abs=0.002)
    assert path["params"]["sigma2"]["stderr"] == pytest.approx(0.00031, abs=0.00005)
    assert rec["correlations"]["amp|sig2"] == pytest.approx(0.93, abs=0.02)
    assert len(rec["provenance"]["inputs"][0]["sha256"]) == 64
    out = save_result(rec, tmp_path / "result.json")
    assert load_result(out)["sample"] == "Cu foil"


def test_cn_with_reference_s02(cu_result):
    from xastk.io import result_from_feffit

    rec = result_from_feffit(cu_result, s02_ref=(0.93, 0.04))
    cn = rec["datasets"][0]["paths"][0]["cn"]
    assert cn["value"] == pytest.approx(12.0, abs=0.1)
    assert cn["stderr"] == pytest.approx(12 * (2 ** 0.5) * 0.04 / 0.93, rel=0.1)


def test_hamilton_third_cumulant(cu_result, cu_result_c3):
    from xastk.stats import hamilton_from_results

    res = hamilton_from_results(cu_result, cu_result_c3)
    assert res.dof_num == 1
    assert res.confidence > 0.95


@pytest.mark.slow
def test_mcmc_matches_least_squares(cu_result):
    pytest.importorskip("emcee")
    from xastk.stats import sample_feffit

    best = {n: cu_result.params[n].value for n in cu_result.var_names}
    mc = sample_feffit(cu_result, nwalkers=32, nsteps=1500, burn=500, thin=5, seed=1)
    # parameters, and the propagated path-parameter uncertainties, are restored after sampling
    assert {n: cu_result.params[n].value for n in cu_result.var_names} == best
    from xastk.io import result_from_feffit

    path = result_from_feffit(cu_result)["datasets"][0]["paths"][0]
    assert path["params"]["sigma2"]["stderr"] == pytest.approx(0.00031, abs=0.00005)
    assert path["r"]["stderr"] == pytest.approx(0.0026, abs=0.0005)
    for n in mc.names:
        s = mc.summary[n]
        assert s["sigma"] / mc.stderr[n] == pytest.approx(1.0, abs=0.25)
        assert abs(s["median"] - mc.best[n]) < 0.5 * mc.stderr[n]
    assert any({a, b} == {"amp", "sig2"} for a, b, _ in mc.strong_correlations(0.85))


def test_rbkg_scan():
    from larch.io import read_ascii

    from xastk.reduce import rbkg_stability

    g = read_ascii(str(DATA / "cu_metal_rt.xdi"))
    scan = rbkg_stability(g.energy, g.mutrans, kmin=3, kmax=16, shell=(1.6, 2.8))
    assert scan.rbkg.size == 9
    assert scan.recommended is not None
    assert scan.recommended <= 0.5 * scan.peak_position.max()
    i = int(np.argmin(np.abs(scan.rbkg - scan.recommended)))
    assert scan.stable[i]
    allowed = scan.stable & (scan.rbkg <= 0.5 * scan.peak_position)
    assert scan.low_r_residue[i] == pytest.approx(scan.low_r_residue[allowed].min())
    assert "recommended" in scan.table()

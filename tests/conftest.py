from pathlib import Path

import numpy as np
import pytest

DATA = Path(__file__).parent / "data"


def _have_larch():
    try:
        import larch  # noqa: F401
        return True
    except ImportError:
        return False


HAVE_LARCH = _have_larch()
needs_larch = pytest.mark.skipif(not HAVE_LARCH, reason="xraylarch not installed")


def cu_fit(c3: bool = False):
    """First-shell fit of the Larch Cu-foil example (examples/feffit/doc_feffit1.py)."""
    from larch.fitting import guess, param, param_group
    from larch.io import read_ascii
    from larch.xafs import autobk, feffit, feffit_dataset, feffit_transform, feffpath

    cu = read_ascii(str(DATA / "cu_metal_rt.xdi"))
    autobk(cu.energy, cu.mutrans, group=cu, rbkg=1.0, kw=2)
    pars = param_group(amp=param(1.0, vary=True), del_e0=param(0.0, vary=True),
                       sig2=param(0.0, vary=True), del_r=guess(0.0, vary=True),
                       c3=param(0.0, vary=c3))
    path = feffpath(str(DATA / "feffcu01.dat"), label="Cu-Cu1", s02="amp", e0="del_e0",
                    sigma2="sig2", deltar="del_r", third="c3")
    trans = feffit_transform(kmin=3, kmax=17, kw=2, dk=4, window="kaiser", rmin=1.4, rmax=3.0)
    ds = feffit_dataset(data=cu, pathlist=[path], transform=trans)
    return feffit(pars, ds)


@pytest.fixture(scope="session")
def cu_result():
    if not HAVE_LARCH:
        pytest.skip("xraylarch not installed")
    return cu_fit(False)


@pytest.fixture(scope="session")
def cu_result_c3():
    if not HAVE_LARCH:
        pytest.skip("xraylarch not installed")
    return cu_fit(True)


@pytest.fixture
def synthetic_result():
    """A minimal schema-v1 record that does not need Larch."""
    return {
        "schema": "xastk-result/1",
        "sample": "Fe-N-C",
        "created": "2026-11-07T00:00:00+00:00",
        "software": {"xastk": "0.1.0.dev0", "python": "3.13"},
        "fit": {"n_idp": 12.4, "n_var": 5, "n_data": 80, "chi_square": 120.0,
                "chi2_reduced": 16.2, "r_factor": 0.0123, "aic": 30.1, "bic": 33.0},
        "s02_ref": {"value": 0.85, "stderr": 0.03},
        "variables": {"amp_n": {"value": 3.4, "stderr": 0.4}},
        "correlations": {"amp_n|sig2_n": 0.93},
        "datasets": [{
            "label": "ds1",
            "transform": {"kmin": 3.0, "kmax": 11.5, "rmin": 1.0, "rmax": 2.2, "kweight": [2.0],
                          "dk": 1.0, "window": "hanning", "fitspace": "r"},
            "n_idp": 7.5,
            "paths": [
                {"label": "Fe-N", "degen": 1.0, "reff": 2.0, "nleg": 2, "geometry": "Fe-N",
                 "params": {"s02": {"value": 3.4, "stderr": 0.4, "expr": "amp_n"},
                            "e0": {"value": 1.2, "stderr": 1.1}, "sigma2": {"value": 0.0071, "stderr": 0.0012},
                            "deltar": {"value": -0.03, "stderr": 0.012}},
                 "r": {"value": 1.97, "stderr": 0.012},
                 "cn": {"value": 4.0, "stderr": 0.5, "convention": "degen*s02/s02_ref"}},
                {"label": "Fe-C", "degen": 8.0, "reff": 3.0, "nleg": 2, "geometry": "Fe-C",
                 "params": {"s02": {"value": 0.85, "stderr": None, "expr": "s02_ref"},
                            "e0": {"value": 1.2, "stderr": 1.1}, "sigma2": {"value": 0.0015, "stderr": 0.003},
                            "deltar": {"value": 0.04, "stderr": 0.03}},
                 "r": {"value": 3.04, "stderr": 0.03},
                 "cn": {"value": 8.0, "stderr": None, "convention": "degen (fixed)"}},
            ],
        }],
        "provenance": {"inputs": []},
        "notes": "",
    }


@pytest.fixture(scope="session")
def operando_series():
    """Synthetic A -> B -> C XANES series (k1 = 0.12, k2 = 0.05 per minute, 0.4 % noise)."""
    from xastk.mcr import consecutive

    rng = np.random.default_rng(0)
    E = np.linspace(-20, 60, 400)
    t = np.linspace(0, 60, 120)

    def edge(e0, wl, w):
        return 0.5 + np.arctan((E - e0) / 1.5) / np.pi + wl * np.exp(-((E - e0 - 4) / w) ** 2)

    S = np.vstack([edge(0, 0.6, 3), edge(2.5, 1.1, 4), edge(5, 0.4, 6)])
    C = consecutive(t, 0.12, 0.05)
    D = C @ S + 0.004 * rng.standard_normal((t.size, E.size))
    return {"E": E, "t": t, "S": S, "C": C, "D": D}

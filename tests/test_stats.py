import math

import pytest
from scipy import stats

from xastk.stats import hamilton_test, n_idp, nvar_budget


def test_nidp_conventions():
    base = 2 * (17 - 3) * (3.0 - 1.4) / math.pi
    assert n_idp(3, 17, 1.4, 3.0) == pytest.approx(base + 1)
    assert n_idp(3, 17, 1.4, 3.0, "stern") == pytest.approx(base + 2)
    assert n_idp(3, 17, 1.4, 3.0, "nyquist") == pytest.approx(base)
    # Larch reports 15.26 for this window in its Cu example
    assert n_idp(3, 17, 1.4, 3.0) == pytest.approx(15.26, abs=0.01)


def test_nidp_rejects_bad_ranges():
    with pytest.raises(ValueError):
        n_idp(10, 3, 1, 2)
    with pytest.raises(ValueError):
        n_idp(3, 10, 1, 2, "nope")


def test_nvar_budget():
    b = nvar_budget(15.0, 11)
    assert not b["ok"] and b["dof"] == 4
    assert nvar_budget(15.0, 6)["ok"]


def test_hamilton_matches_f_distribution():
    res = hamilton_test(112.7, 4, 68.3, 5, 15.26)
    F = (112.7 - 68.3) / (68.3 / (15.26 - 5))
    assert res.F == pytest.approx(F)
    assert res.confidence == pytest.approx(stats.f.cdf(F, 1, 15.26 - 5))
    assert res.significant
    assert res.ratio == pytest.approx(math.sqrt(112.7 / 68.3))
    assert "significant" in str(res)


def test_hamilton_no_improvement():
    res = hamilton_test(100.0, 4, 99.9, 6, 20.0)
    assert not res.significant
    assert res.confidence < 0.5


@pytest.mark.parametrize("args", [(100, 5, 90, 5, 20), (100, 4, 90, 6, 6), (0, 4, 90, 5, 20)])
def test_hamilton_invalid(args):
    with pytest.raises(ValueError):
        hamilton_test(*args)

import numpy as np
import pytest

from xastk.mcr import CONSECUTIVE, consecutive, first_order, mcr_als, pca_rank


def test_kinetic_models_close_to_one():
    t = np.linspace(0, 50, 101)
    for C in (consecutive(t, 0.2, 0.05), consecutive(t, 0.1, 0.1), first_order(t, 0.3)):
        assert np.allclose(C.sum(axis=1), 1.0)
        assert C.min() > -1e-12


def test_equal_rates_limit_is_continuous():
    t = np.linspace(0, 50, 101)
    assert np.allclose(consecutive(t, 0.1, 0.1), consecutive(t, 0.1, 0.1 + 1e-7), atol=1e-6)


def test_pca_rank_finds_three_components(operando_series):
    est = pca_rank(operando_series["D"])
    assert est.n_components == 3
    assert "IND minimum at n = 3" in est.table()


def test_soft_mcr_fits_well_but_is_ambiguous(operando_series):
    D, C_true = operando_series["D"], operando_series["C"]
    res = mcr_als(D, S0=D[[0, 40, -1]], max_iter=300)
    assert res.lack_of_fit < 0.6            # fits the data to the noise level...
    assert np.abs(res.C - C_true).max() > 0.1   # ...yet the concentrations are wrong


def test_kinetic_mcr_recovers_truth(operando_series):
    D, t = operando_series["D"], operando_series["t"]
    res = mcr_als(D, S0=D[[0, 40, -1]], kinetic=CONSECUTIVE, t=t, max_iter=300)
    assert res.lack_of_fit < 0.6
    assert np.abs(res.C - operando_series["C"]).max() < 0.03
    assert np.abs(res.S - operando_series["S"]).max() < 0.05
    k1, k2 = res.kinetic_params
    assert k1 == pytest.approx(0.12, rel=0.05)
    assert k2 == pytest.approx(0.05, rel=0.05)


def test_mcr_argument_checks(operando_series):
    D = operando_series["D"]
    with pytest.raises(ValueError):
        mcr_als(D)
    with pytest.raises(ValueError):
        mcr_als(D, S0=D[[0, -1]], kinetic=CONSECUTIVE, t=operando_series["t"])
    with pytest.raises(ValueError):
        mcr_als(D, S0=D[[0, 40, -1]], kinetic=CONSECUTIVE)

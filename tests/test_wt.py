import numpy as np
import pytest

from xastk.wt import find_lobes, morlet_wt, plot_wt_composite

K = np.arange(0, 16.001, 0.05)


def single_shell(R=2.5):
    return np.sin(2 * K * R) * np.exp(-2 * 0.005 * K ** 2) / np.clip(K, 0.5, None)


@pytest.mark.parametrize("eta, sigma", [(5.0, 1.0), (10.0, 1.0)])
def test_ridge_at_shell_distance(eta, sigma):
    R = 2.5
    wt = morlet_wt(K, single_shell(R), kweight=2, eta=eta, sigma=sigma, rmin=0.5, rmax=5, nr=226)
    lobes = find_lobes(wt)
    assert len(lobes) == 1
    # the modulus maximum is biased low by about R / (2 (eta sigma)^2); it vanishes for eta*sigma >> 1
    bias = R / (2 * (eta * sigma) ** 2)
    assert lobes[0].r == pytest.approx(R - bias, abs=0.03)


def test_light_and_heavy_scatterers_separate_in_k():
    light = np.exp(-((K - 4) / 2) ** 2)
    heavy = np.exp(-((K - 10) / 2.5) ** 2)
    chi = (light + heavy) * np.sin(2 * K * 2.2) / np.clip(K, 0.5, None) ** 2
    wt = morlet_wt(K, chi, kweight=2, eta=7.5, sigma=0.5, rmin=0.5, rmax=5, nr=226)
    lobes = find_lobes(wt, threshold=0.3)
    ks = sorted(lb.k for lb in lobes[:2])
    assert ks[0] == pytest.approx(4.0, abs=0.5)
    assert ks[1] == pytest.approx(10.0, abs=0.5)
    assert all(lb.r == pytest.approx(2.15, abs=0.1) for lb in lobes[:2])


def test_nonuniform_input_is_resampled():
    k = np.sort(np.random.default_rng(1).uniform(0, 16, 500))
    chi = np.sin(2 * k * 2.5)
    wt = morlet_wt(k, chi, rmin=1, rmax=4, nr=31)
    assert np.allclose(np.diff(wt.k), 0.05)


def test_bad_parameters():
    with pytest.raises(ValueError):
        morlet_wt(K, single_shell(), rmin=0)
    with pytest.raises(ValueError):
        morlet_wt(K, single_shell(), eta=-1)


def test_composite_figure(tmp_path):
    import matplotlib

    matplotlib.use("Agg")
    wt = morlet_wt(K, single_shell(), rmin=0.5, rmax=5, nr=91)
    fig = plot_wt_composite(wt, lobes=find_lobes(wt), title="test")
    out = tmp_path / "wt.png"
    fig.savefig(out, dpi=80)
    assert out.stat().st_size > 1000
    dk, dr = wt.resolution(2.5)
    assert dk == pytest.approx(1.0) and dr == pytest.approx(0.5)

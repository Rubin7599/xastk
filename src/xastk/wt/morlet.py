"""Morlet wavelet transform of EXAFS chi(k).

Distance-calibrated form (Funke, Scheinost & Chukalina, Phys. Rev. B 71,
094110, 2005; same convention as the wtEXAFS tool), with x = 2 r (k - k') / eta::

    psi(x)   = 1 / (sqrt(2 pi) sigma) * exp(i eta x) * exp(-x^2 / (2 sigma^2))
    W(k', r) = sqrt(2 r / eta) * integral chi(k) k^n conj(psi(2 r (k - k') / eta)) dk

The carrier exp(i eta x) = exp(2 i r (k - k')) oscillates like sin(2 k R), so
the modulus |W| peaks at r = R (before any phase-shift correction), just as
|chi(R)| does. The Gaussian envelope has a k-width of about eta * sigma / (2 r):
a larger product eta * sigma gives better R resolution and worse k resolution.
Funke et al. call eta*sigma >~ 10 an "overview" wavelet and eta*sigma ~ 5 a
"detail" wavelet that separates light from heavy backscatterers in k.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class WaveletResult:
    k: np.ndarray            # (n_k,) wavelet position k'
    r: np.ndarray            # (n_r,) distance axis (Angstrom, not phase corrected)
    W: np.ndarray            # (n_r, n_k) complex coefficients
    eta: float
    sigma: float
    kweight: float
    meta: dict = field(default_factory=dict)

    @property
    def modulus(self) -> np.ndarray:
        return np.abs(self.W)

    def k_profile(self) -> np.ndarray:
        """Marginal over r (max of |W| along r) for the side panel of a composite plot."""
        return self.modulus.max(axis=0)

    def r_profile(self) -> np.ndarray:
        """Marginal over k (max of |W| along k)."""
        return self.modulus.max(axis=1)

    def resolution(self, r: float) -> tuple[float, float]:
        """Approximate (delta_k, delta_r) half-widths of the wavelet at distance r."""
        dk = self.eta * self.sigma / (2.0 * r)
        dr = 1.0 / (2.0 * dk)
        return dk, dr


def _uniform(k, chi, kstep):
    k = np.asarray(k, float)
    chi = np.asarray(chi, float)
    if kstep is None:
        steps = np.diff(k)
        if np.allclose(steps, steps[0], rtol=1e-3, atol=1e-6):
            return k, chi
        kstep = 0.05
    kk = np.arange(k[0], k[-1] + 0.5 * kstep, kstep)
    return kk, np.interp(kk, k, chi)


def morlet_wt(
    k: np.ndarray,
    chi: np.ndarray,
    kweight: float = 2.0,
    eta: float = 5.0,
    sigma: float = 1.0,
    rmin: float = 0.2,
    rmax: float = 6.0,
    nr: int = 291,
    kmin: float | None = None,
    kmax: float | None = None,
    kstep: float | None = None,
) -> WaveletResult:
    """Continuous Morlet wavelet transform of k^n chi(k).

    Parameters
    ----------
    k, chi : EXAFS on a k grid (Larch groups provide ``group.k`` and ``group.chi``)
    kweight : n in k^n chi(k)
    eta, sigma : Morlet parameters (see module docstring)
    rmin, rmax, nr : distance axis; rmin must be > 0
    kmin, kmax : optional k range of the input to use (outside -> zero)
    kstep : resample to this uniform step if k is not uniform (default 0.05)
    """
    if rmin <= 0:
        raise ValueError("rmin must be > 0")
    if eta <= 0 or sigma <= 0:
        raise ValueError("eta and sigma must be positive")
    kk, cc = _uniform(k, chi, kstep)
    signal = cc * kk ** kweight
    if kmin is not None:
        signal = np.where(kk >= kmin, signal, 0.0)
    if kmax is not None:
        signal = np.where(kk <= kmax, signal, 0.0)
    dk = kk[1] - kk[0]
    r = np.linspace(rmin, rmax, nr)
    norm = 1.0 / (np.sqrt(2.0 * np.pi) * sigma)
    W = np.empty((nr, kk.size), complex)
    diff = kk[None, :] - kk[:, None]          # [k', k] -> k - k'
    for i, ri in enumerate(r):
        x = 2.0 * ri * diff / eta
        kern = norm * np.exp(-1j * eta * x) * np.exp(-(x * x) / (2.0 * sigma * sigma))  # conj(psi)
        W[i] = np.sqrt(2.0 * ri / eta) * (kern @ signal) * dk
    return WaveletResult(k=kk, r=r, W=W, eta=eta, sigma=sigma, kweight=kweight,
                         meta={"kmin": kmin, "kmax": kmax, "kstep": float(dk)})

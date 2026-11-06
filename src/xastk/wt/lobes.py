"""Detect and locate lobe maxima in a wavelet modulus map."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage

from .morlet import WaveletResult


@dataclass(frozen=True)
class Lobe:
    k: float          # 1/Angstrom, sub-pixel
    r: float          # Angstrom, sub-pixel (not phase corrected)
    amplitude: float  # |W| at the grid maximum
    relative: float   # amplitude / global maximum


def _parabolic(y_m, y_0, y_p):
    den = y_m - 2.0 * y_0 + y_p
    if den == 0:
        return 0.0
    return 0.5 * (y_m - y_p) / den


def find_lobes(
    wt: WaveletResult,
    threshold: float = 0.2,
    size: int = 9,
    rmin: float | None = None,
    rmax: float | None = None,
) -> list[Lobe]:
    """Local maxima of |W| above ``threshold`` x global max, sorted by amplitude.

    ``size`` is the neighbourhood (pixels) for the maximum filter. Positions are
    refined with a 3-point parabola in k and in r independently.
    """
    M = wt.modulus.copy()
    if rmin is not None:
        M[wt.r < rmin] = 0.0
    if rmax is not None:
        M[wt.r > rmax] = 0.0
    gmax = M.max()
    if gmax <= 0:
        return []
    peaks = (M == ndimage.maximum_filter(M, size=size, mode="nearest")) & (M >= threshold * gmax)
    lobes = []
    dr = wt.r[1] - wt.r[0]
    dk = wt.k[1] - wt.k[0]
    for ir, ik in zip(*np.nonzero(peaks), strict=True):
        r = wt.r[ir]
        k = wt.k[ik]
        if 0 < ir < M.shape[0] - 1:
            r += _parabolic(M[ir - 1, ik], M[ir, ik], M[ir + 1, ik]) * dr
        if 0 < ik < M.shape[1] - 1:
            k += _parabolic(M[ir, ik - 1], M[ir, ik], M[ir, ik + 1]) * dk
        lobes.append(Lobe(k=float(k), r=float(r), amplitude=float(M[ir, ik]), relative=float(M[ir, ik] / gmax)))
    return sorted(lobes, key=lambda lb: -lb.amplitude)

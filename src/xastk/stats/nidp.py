"""Number of independent points in an EXAFS fit."""

from __future__ import annotations

import math


def n_idp(kmin: float, kmax: float, rmin: float, rmax: float, convention: str = "larch") -> float:
    """Number of independent points of a Fourier-filtered EXAFS fit.

    Parameters
    ----------
    kmin, kmax : fit range in k (1/Angstrom)
    rmin, rmax : fit range in R (Angstrom)
    convention : ``"larch"`` -> 1 + 2 dk dR / pi   (what Larch ``feffit`` uses)
                 ``"stern"`` -> 2 + 2 dk dR / pi   (Stern, Phys. Rev. B 48, 9825, 1993)
                 ``"nyquist"`` -> 2 dk dR / pi

    Returns
    -------
    float
    """
    if kmax <= kmin or rmax <= rmin:
        raise ValueError("require kmax > kmin and rmax > rmin")
    base = 2.0 * (kmax - kmin) * (rmax - rmin) / math.pi
    offsets = {"larch": 1.0, "stern": 2.0, "nyquist": 0.0}
    try:
        return base + offsets[convention]
    except KeyError as exc:
        raise ValueError(f"unknown convention {convention!r}; use one of {sorted(offsets)}") from exc


def nvar_budget(nidp: float, nvar: int, max_fraction: float = 2.0 / 3.0) -> dict:
    """Compare the number of variables with the information budget.

    ``max_fraction`` is a lab rule of thumb (default 2/3 of N_idp), not a law.
    Returns a dict with the ratio, the degrees of freedom and an ``ok`` flag.
    """
    return {
        "n_idp": nidp,
        "n_var": nvar,
        "dof": nidp - nvar,
        "ratio": nvar / nidp,
        "ok": nvar <= max_fraction * nidp,
        "max_fraction": max_fraction,
    }

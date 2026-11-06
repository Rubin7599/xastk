"""AUTOBK R_bkg stability scan.

Runs Larch ``autobk`` over a grid of R_bkg values and records, for each:
  - the height and position of the first-shell peak of |chi(R)| inside ``shell``
  - the low-R residue: mean |chi(R)| in a FIXED window [0, ``low_rmax``], so
    that values are comparable across R_bkg

Recommendation rule
-------------------
1. *stable*: first-shell peak height within ``tol`` (relative) of the median
   over the whole grid -- the background choice is not changing the structural
   signal;
2. *allowed*: R_bkg <= ``max_fraction`` x first-shell peak position (the
   spline must not be able to follow the first shell);
3. among stable and allowed values, pick the one with the smallest low-R residue.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class RbkgScan:
    rbkg: np.ndarray
    peak_height: np.ndarray
    peak_position: np.ndarray
    low_r_residue: np.ndarray
    stable: np.ndarray
    recommended: float | None
    reason: str

    def table(self) -> str:
        lines = [f"{'R_bkg':>6s} {'peak |chi(R)|':>14s} {'peak R':>7s} {'low-R residue':>14s} {'stable':>7s}"]
        for i in range(self.rbkg.size):
            mark = "  <- recommended" if self.recommended is not None and np.isclose(self.rbkg[i], self.recommended) else ""
            lines.append(f"{self.rbkg[i]:6.2f} {self.peak_height[i]:14.4f} {self.peak_position[i]:7.3f} "
                         f"{self.low_r_residue[i]:14.4f} {'yes' if self.stable[i] else 'no':>7s}{mark}")
        lines.append(self.reason)
        return "\n".join(lines)


def rbkg_stability(
    energy: np.ndarray,
    mu: np.ndarray,
    rbkg_values=None,
    kmin: float = 2.0,
    kmax: float = 14.0,
    kweight: float = 2.0,
    dk: float = 1.0,
    window: str = "hanning",
    shell: tuple[float, float] = (1.0, 3.2),
    low_rmax: float = 0.9,
    tol: float = 0.05,
    max_fraction: float = 0.5,
    e0: float | None = None,
) -> RbkgScan:
    """Scan R_bkg with Larch ``autobk`` + ``xftf`` and recommend a value.

    ``shell`` is the R window (Angstrom, not phase corrected) of the first-shell
    peak; ``low_rmax`` the upper edge of the fixed low-R window.
    """
    try:
        from larch import Group
        from larch.xafs import autobk, pre_edge, xftf
    except ImportError as exc:  # pragma: no cover
        raise ImportError("rbkg_stability needs xraylarch") from exc

    rbkg_values = np.round(np.arange(0.6, 1.41, 0.1), 3) if rbkg_values is None else np.asarray(rbkg_values, float)
    heights, positions, residues = [], [], []
    for rb in rbkg_values:
        g = Group(energy=np.asarray(energy, float), mu=np.asarray(mu, float))
        pre_edge(g.energy, g.mu, group=g, e0=e0)
        autobk(g.energy, g.mu, group=g, rbkg=float(rb), kweight=kweight, e0=e0)
        xftf(g.k, g.chi, group=g, kmin=kmin, kmax=kmax, kweight=kweight, dk=dk, window=window)
        r, mag = g.r, g.chir_mag
        sel = (r >= shell[0]) & (r <= shell[1])
        i = np.argmax(mag[sel])
        heights.append(float(mag[sel][i]))
        positions.append(float(r[sel][i]))
        residues.append(float(np.mean(mag[r <= low_rmax])))
    heights = np.array(heights)
    positions = np.array(positions)
    residues = np.array(residues)

    med = np.median(heights)
    stable = np.abs(heights - med) / med < tol
    allowed = stable & (rbkg_values <= max_fraction * positions)
    if allowed.any():
        idx = np.flatnonzero(allowed)
        best = idx[np.argmin(residues[idx])]
        rec = float(rbkg_values[best])
        reason = (f"recommended: lowest mean |chi(R)| below {low_rmax:g} A among R_bkg values whose first-shell peak "
                  f"is within {tol:.0%} of the grid median and R_bkg <= {max_fraction:g} x peak R")
    else:
        rec = None
        reason = ("no stable R_bkg <= max_fraction x peak R: inspect mu(E) for glitches, "
                  "normalization problems or an unsuitable shell window")
    return RbkgScan(rbkg_values, heights, positions, residues, stable, rec, reason)

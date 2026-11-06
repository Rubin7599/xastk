"""How many components? PCA eigenvalues and Malinowski's IND function."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class RankEstimate:
    singular_values: np.ndarray
    explained: np.ndarray        # fraction of variance per component
    ind: np.ndarray              # Malinowski indicator function, index n-1 -> n components
    n_components: int            # argmin of IND

    def table(self, nmax: int = 8) -> str:
        lines = [f"{'n':>3s} {'sing. value':>12s} {'explained':>10s} {'IND':>12s}"]
        for i in range(min(nmax, len(self.ind))):
            lines.append(f"{i + 1:3d} {self.singular_values[i]:12.5g} {self.explained[i]:10.5f} {self.ind[i]:12.5g}")
        lines.append(f"-> IND minimum at n = {self.n_components}")
        return "\n".join(lines)


def pca_rank(D: np.ndarray, center: bool = False) -> RankEstimate:
    """Estimate the number of spectroscopically distinct components in D.

    Uses Malinowski's indicator function (Anal. Chem. 49, 612, 1977)::

        RE(n)  = sqrt( sum_{j>n} lambda_j / (r (c - n)) )
        IND(n) = RE(n) / (c - n)^2

    with r >= c the matrix dimensions and lambda_j = s_j^2. XAS data are
    normally not mean-centred because the components are real spectra.
    """
    D = np.asarray(D, float)
    if center:
        D = D - D.mean(axis=0, keepdims=True)
    s = np.linalg.svd(D, compute_uv=False)
    lam = s ** 2
    r, c = max(D.shape), min(D.shape)
    ind = np.empty(c - 1)
    for n in range(1, c):
        re = np.sqrt(lam[n:].sum() / (r * (c - n)))
        ind[n - 1] = re / (c - n) ** 2
    return RankEstimate(
        singular_values=s,
        explained=lam / lam.sum(),
        ind=ind,
        n_components=int(np.argmin(ind)) + 1,
    )

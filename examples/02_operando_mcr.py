"""Synthetic operando XANES series A -> B -> C: why MCR-ALS needs physical constraints.

    python examples/02_operando_mcr.py

Prints the PCA rank estimate and compares soft-constrained MCR-ALS with the
kinetics-constrained (hard-soft) variant; saves examples/output/mcr_concentrations.png.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from xastk.mcr import CONSECUTIVE, consecutive, mcr_als, pca_rank

OUT = Path(__file__).resolve().parents[1] / "examples" / "output"
OUT.mkdir(exist_ok=True)

rng = np.random.default_rng(0)
E = np.linspace(-20, 60, 400)
t = np.linspace(0, 60, 120)


def edge(e0, wl, w):
    return 0.5 + np.arctan((E - e0) / 1.5) / np.pi + wl * np.exp(-((E - e0 - 4) / w) ** 2)


S_true = np.vstack([edge(0, 0.6, 3), edge(2.5, 1.1, 4), edge(5, 0.4, 6)])
C_true = consecutive(t, 0.12, 0.05)
D = C_true @ S_true + 0.004 * rng.standard_normal((t.size, E.size))

print(pca_rank(D).table(nmax=5))

S0 = D[[0, 40, -1]]   # common practice: first, middle, last spectrum
soft = mcr_als(D, S0=S0, max_iter=300)
hard = mcr_als(D, S0=S0, kinetic=CONSECUTIVE, t=t, max_iter=300)
for name, res in (("soft (non-negativity + closure)", soft), ("+ kinetic A->B->C", hard)):
    print(f"{name:34s} lack of fit = {res.lack_of_fit:.3f} %   "
          f"max |C error| = {np.abs(res.C - C_true).max():.3f}   "
          f"max |S error| = {np.abs(res.S - S_true).max():.3f}")
print(f"recovered k1, k2 = {hard.kinetic_params[0]:.4f}, {hard.kinetic_params[1]:.4f} (true 0.12, 0.05)")

fig, axes = plt.subplots(1, 2, figsize=(8, 3.2), sharey=True, constrained_layout=True)
for ax, res, title in ((axes[0], soft, "Soft constraints only"), (axes[1], hard, "With kinetic constraint")):
    for j, lab in enumerate("ABC"):
        ax.plot(t, C_true[:, j], color=f"C{j}", lw=2.5, alpha=0.35)
        ax.plot(t, res.C[:, j], color=f"C{j}", lw=1.2, ls="--", label=f"{lab}")
    ax.set_title(f"{title}\nlack of fit {res.lack_of_fit:.3f} %", fontsize=9)
    ax.set_xlabel("time (min)")
axes[0].set_ylabel("fraction")
axes[1].legend(title="thick: true, dashed: MCR", fontsize=8, title_fontsize=8)
fig.savefig(OUT / "mcr_concentrations.png", dpi=200)
print(OUT / "mcr_concentrations.png")

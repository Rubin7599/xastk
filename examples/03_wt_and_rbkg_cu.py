"""Cu foil: R_bkg stability scan and a Morlet WT-EXAFS composite figure.

    python examples/03_wt_and_rbkg_cu.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
from larch.io import read_ascii
from larch.xafs import autobk

from xastk.reduce import rbkg_stability
from xastk.wt import find_lobes, morlet_wt, plot_wt_composite

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "examples" / "output"
OUT.mkdir(exist_ok=True)

cu = read_ascii(str(ROOT / "tests" / "data" / "cu_metal_rt.xdi"))

print("== R_bkg stability scan")
scan = rbkg_stability(cu.energy, cu.mutrans, kmin=3, kmax=16, shell=(1.6, 2.8))
print(scan.table())

autobk(cu.energy, cu.mutrans, group=cu, rbkg=scan.recommended or 1.0, kweight=2)
print("\n== Morlet WT (eta = 10, sigma = 1, k^2)")
wt = morlet_wt(cu.k, cu.chi, kweight=2, eta=10, sigma=1, rmin=0.5, rmax=6, kmin=3, kmax=16)
lobes = find_lobes(wt, threshold=0.25)
for lb in lobes:
    print(f"lobe at k = {lb.k:.2f} 1/A, R = {lb.r:.3f} A (relative {lb.relative:.2f})")
fig = plot_wt_composite(wt, lobes=lobes, title="Cu foil, Cu K-edge", kmax_plot=16, rmax_plot=6)
fig.savefig(OUT / "wt_cu.png", dpi=200)
print(OUT / "wt_cu.png")

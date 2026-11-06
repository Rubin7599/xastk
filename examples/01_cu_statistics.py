"""Cu foil, first shell: fit with Larch, then validate with xastk.

Run from the repository root:
    python examples/01_cu_statistics.py

Writes examples/output/cu_result.json, cu_table.{csv,tex} and prints the
Hamilton test, the MCMC comparison and the sanity flags.
"""

from pathlib import Path

from larch.fitting import guess, param, param_group
from larch.io import read_ascii
from larch.xafs import autobk, feffit, feffit_dataset, feffit_transform, feffpath

from xastk.io import result_from_feffit, save_result
from xastk.reduce import check_result
from xastk.report import rows_from_results, to_csv, to_latex
from xastk.stats import hamilton_from_results, sample_feffit

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data"
OUT = ROOT / "examples" / "output"
OUT.mkdir(exist_ok=True)

cu = read_ascii(str(DATA / "cu_metal_rt.xdi"))
autobk(cu.energy, cu.mutrans, group=cu, rbkg=1.0, kw=2)
trans = feffit_transform(kmin=3, kmax=17, kw=2, dk=4, window="kaiser", rmin=1.4, rmax=3.0)


def fit(with_c3: bool):
    pars = param_group(amp=param(1.0, vary=True), del_e0=param(0.0, vary=True),
                       sig2=param(0.0, vary=True), del_r=guess(0.0, vary=True),
                       c3=param(0.0, vary=with_c3))
    path = feffpath(str(DATA / "feffcu01.dat"), label="Cu-Cu1", s02="amp", e0="del_e0",
                    sigma2="sig2", deltar="del_r", third="c3")
    return feffit(pars, feffit_dataset(data=cu, pathlist=[path], transform=trans))


out_a = fit(False)
out_b = fit(True)

print("== Hamilton test: does a third cumulant (C3) improve the fit?")
print(hamilton_from_results(out_a, out_b))

print("\n== MCMC posterior vs Larch least squares (model A)")
mc = sample_feffit(out_a, nwalkers=32, nsteps=3000, burn=1000, thin=10, seed=1)
print(mc.table())
for a, b, c in mc.strong_correlations(0.85):
    print(f"strong correlation: {a} - {b}: {c:+.2f}")

rec = result_from_feffit(out_a, sample="Cu foil (RT)", inputs=[DATA / "cu_metal_rt.xdi", DATA / "feffcu01.dat"])
save_result(rec, OUT / "cu_result.json")

print("\n== Sanity checks")
for flag in check_result(rec) or ["no flags"]:
    print(flag)

rows = rows_from_results([rec])
to_csv(rows, OUT / "cu_table.csv")
(OUT / "cu_table.tex").write_text(to_latex(rows, caption="Cu foil first-shell fit."), encoding="utf-8")
print("\n== SI table row")
print(rows[0].cells)

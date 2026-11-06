# Getting started

## 1. Install

```bash
git clone https://github.com/Rubin7599/xastk && cd xastk
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -e ".[all]"
pytest
```

| Extra | Adds | Needed for |
| --- | --- | --- |
| `larch` | xraylarch | `result_from_feffit`, `rbkg_stability`, the Larch tests |
| `stats` | xraylarch, emcee | `sample_feffit` |
| `docx` | python-docx | `to_docx` |
| `schema` | jsonschema | `validate_result` |
| `origin` | originpro, pandas (Windows only) | `run_bundle` |
| `all` | larch, stats, docx, schema | |
| `dev` | all + pytest, ruff | development |

## 2. Statistics on a Larch fit

`examples/01_cu_statistics.py` runs end to end on the Cu foil example in `tests/data/`:

```bash
python examples/01_cu_statistics.py
```

It fits model A (S0², ΔE₀, σ², ΔR) and model B (A + third cumulant), then:

1. `hamilton_from_results(out_a, out_b)` tests whether B is a real improvement;
2. `sample_feffit(out_a)` samples the posterior and prints it next to Larch's least-squares values;
3. `result_from_feffit(out_a, ...)` builds the `result.json` record (with SHA-256 of the inputs);
4. `check_result(rec)` lists physical sanity warnings;
5. `rows_from_results` + `to_csv` / `to_latex` write the SI table.

**Tips**

- Use the same transform (k, R ranges, window, k-weights) for the fits you compare with the
  Hamilton test. `hamilton_from_results` refuses fits whose N_idp differs.
- Pass physical `bounds` to `sample_feffit`, e.g. `bounds={"sig2": (0, 0.05)}`. Variables with
  `min`/`max` set in Larch are bounded automatically.
- Call `result_from_feffit` after sampling if you like: `sample_feffit` restores the best-fit
  values and the propagated path-parameter uncertainties.

## 3. Operando series with MCR-ALS

```python
from xastk.mcr import pca_rank, mcr_als, CONSECUTIVE

print(pca_rank(D).table())                      # D: rows = spectra in time order
res = mcr_als(D, S0=D[[0, n_mid, -1]],          # initial spectra
              kinetic=CONSECUTIVE, t=t_minutes)  # A -> B -> C hard model
res.C, res.S, res.kinetic_params, res.lack_of_fit
```

Without a physical constraint, treat MCR results as one member of a family of equally good
solutions. See `examples/02_operando_mcr.py` and technical report §5.4. You can write your own
`KineticModel(func, param_names, p0, n_components)` for other mechanisms; `func(t, *params)` must
return an (n_t × n) matrix whose rows sum to one.

## 4. Wavelet transform

```python
from xastk.wt import morlet_wt, find_lobes, plot_wt_composite

wt = morlet_wt(group.k, group.chi, kweight=2, eta=10, sigma=1, rmax=6, kmin=3, kmax=14)
lobes = find_lobes(wt, threshold=0.25)
fig = plot_wt_composite(wt, lobes=lobes, title="sample", kmax_plot=14, rmax_plot=6)
fig.savefig("wt.pdf")
```

Choosing η and σ: ησ ≳ 10 gives sharp R and broad k ("overview"). ησ ≈ 5 separates light (low k)
from heavy (high k) scatterers at the same distance ("detail"). Lobe R is biased low by about
R/(2(ησ)²); the parameters are printed on the figure. See report §5.5.

## 5. R_bkg stability scan

```bash
xastk rbkg sample.xdi --mu mutrans --kmin 3 --kmax 14 --shell-rmin 1.2 --shell-rmax 2.4
```

Pick the shell window around the first-shell peak of |χ(R)| (not phase corrected). The table shows
the peak height, its stability and the low-R residue for each R_bkg, and marks the recommendation.

## 6. Tables from many fits

```bash
xastk table runs/*/result.json --to csv latex docx --out si_table --caption "EXAFS fitting results."
xastk check runs/*/result.json --strict      # non-zero exit status if any warning: use in CI
```

## 7. OriginPro

On any OS, write a figure bundle:

```python
from xastk.report import export_bundle
export_bundle("fig_chik", sheets={"chik": {"k": k, "data": kchi, "fit": kchi_fit}},
              layers=[{"sheet": "chik", "x": "k", "y": ["data", "fit"],
                       "xlabel": "k (1/Å)", "ylabel": "k²χ(k)"}],
              template="exafs_k.otpu")
```

On Windows with Origin 2021+ and `pip install "xastk[origin]"`:

```python
from xastk.report import run_bundle
run_bundle("fig_chik", template_dir="templates/origin")   # writes fig_chik/figure.opju
```

`run_bundle` has not yet been verified on an Origin installation; please report problems.

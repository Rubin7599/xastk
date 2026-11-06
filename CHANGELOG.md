# Changelog

## 0.1.0.dev0 — 2026-11-07

First public scaffold.

### Added
- `xastk.stats`: `n_idp` (Larch, Stern, Nyquist conventions), `nvar_budget`, Hamilton F-test
  (`hamilton_test`, `hamilton_from_results`), MCMC posterior for Larch `feffit` results
  (`sample_feffit`, likelihood weighted by N_idp/N_pts, noise scale sampled).
- `xastk.mcr`: PCA rank with Malinowski IND, MCR-ALS with non-negativity, closure and kinetic
  hard-model constraints (`FIRST_ORDER`, `CONSECUTIVE`, custom `KineticModel`).
- `xastk.wt`: distance-calibrated Morlet wavelet transform, sub-pixel lobe detection, composite
  figure with marginals.
- `xastk.reduce`: AUTOBK R_bkg stability scan, post-fit sanity rules.
- `xastk.io`: `result.json` schema v1 and `result_from_feffit` with SHA-256 provenance.
- `xastk.report`: PDG-rule uncertainty notation with decimal half-up rounding, SI tables
  (CSV, LaTeX, Word), OriginPro figure bundles (`run_bundle` not yet verified on Origin).
- CLI: `nidp`, `hamilton`, `table`, `check`, `wt`, `rbkg`.
- Technical report, getting-started guide, examples, 49 tests, CI.

### Notes
- Drafted with AI assistance (Claude, Anthropic); under maintainer review.

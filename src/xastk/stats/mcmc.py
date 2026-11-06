"""MCMC posterior sampling for a finished Larch ``feffit`` result.

Likelihood
----------
Larch evaluates the residual on an oversampled grid (k or R points), so the
N_pts residual values carry only about N_idp independent pieces of
information. The log-likelihood is therefore down-weighted by
w = N_idp / N_pts -- the same factor Larch uses to rescale chi-square::

    ln L(theta, s) = -1/2 * w * [ sum_i r_i(theta)^2 / s^2 + N_pts * ln s^2 ]

The noise scale ``s`` is sampled jointly (Jeffreys-like flat prior on ln s),
so the posterior width does not depend on Larch's epsilon estimate being
exact. Validated on the Larch Cu-foil example: posterior medians agree with
the least-squares values to within 0.1 sigma, and posterior 1-sigma widths are
7-10 % larger than Larch's least-squares standard errors.

Priors are uniform within each variable's ``min``/``max`` (from the lmfit
Parameter) intersected with any user ``bounds``.

Notes
-----
This module calls Larch's internal residual function
(``larch.xafs.feffit._feffit_resid``). It is pinned by the integration test
so that a Larch API change fails CI instead of silently changing results.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


def _residual_function():
    try:
        from larch.xafs.feffit import _feffit_resid
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise ImportError("xastk.stats.sample_feffit needs xraylarch (pip install xraylarch)") from exc
    return _feffit_resid


@dataclass
class MCMCResult:
    names: list[str]
    chain: np.ndarray                 # (n_samples, n_var + 1); last column is ln(s)
    best: dict[str, float]            # least-squares values from Larch
    stderr: dict[str, float]          # least-squares standard errors from Larch
    acceptance: float
    weight: float                     # w = N_idp / N_pts
    n_idp: float
    n_pts: int
    summary: dict[str, dict] = field(default_factory=dict)
    correlation: np.ndarray | None = None

    def __post_init__(self):
        nv = len(self.names)
        for i, name in enumerate(self.names):
            q16, q50, q84 = np.percentile(self.chain[:, i], [16, 50, 84])
            self.summary[name] = {
                "median": float(q50),
                "lo68": float(q16),
                "hi68": float(q84),
                "sigma": float((q84 - q16) / 2),
                "mean": float(np.mean(self.chain[:, i])),
            }
        if nv > 1:
            self.correlation = np.corrcoef(self.chain[:, :nv].T)

    def table(self) -> str:
        head = f"{'param':10s} {'best (LSQ)':>12s} {'LSQ 1sig':>10s} {'MCMC median':>12s} {'MCMC 1sig':>10s} {'ratio':>6s}"
        lines = [head, "-" * len(head)]
        for name in self.names:
            s = self.summary[name]
            se = self.stderr.get(name) or float("nan")
            lines.append(
                f"{name:10s} {self.best[name]:12.5g} {se:10.3g} {s['median']:12.5g} "
                f"{s['sigma']:10.3g} {s['sigma'] / se:6.2f}"
            )
        lines.append(f"acceptance = {self.acceptance:.2f}, w = N_idp/N_pts = {self.n_idp:.2f}/{self.n_pts}")
        return "\n".join(lines)

    def strong_correlations(self, threshold: float = 0.9) -> list[tuple[str, str, float]]:
        out = []
        if self.correlation is None:
            return out
        for i in range(len(self.names)):
            for j in range(i + 1, len(self.names)):
                c = float(self.correlation[i, j])
                if abs(c) >= threshold:
                    out.append((self.names[i], self.names[j], c))
        return out


def sample_feffit(
    result,
    datasets=None,
    nwalkers: int = 32,
    nsteps: int = 3000,
    burn: int = 1000,
    thin: int = 10,
    bounds: dict[str, tuple[float, float]] | None = None,
    init_scale: float = 0.05,
    seed: int | None = 0,
    progress: bool = False,
) -> MCMCResult:
    """Sample the posterior of the variables of a Larch ``feffit`` result.

    Parameters
    ----------
    result : output of ``larch.xafs.feffit``
    datasets : list of ``FeffitDataSet``; defaults to ``result.datasets``
    nwalkers, nsteps, burn, thin : emcee settings
    bounds : optional ``{name: (lo, hi)}`` uniform priors, e.g. ``{"sig2": (0, 0.05)}``
    init_scale : walkers start at best + init_scale * stderr * N(0, 1)
    seed : random seed for reproducible chains

    The Larch parameters are restored to their best-fit values on return.
    """
    try:
        import emcee
    except ImportError as exc:  # pragma: no cover
        raise ImportError("sample_feffit needs emcee (pip install 'xastk[stats]')") from exc

    resid = _residual_function()
    datasets = list(datasets if datasets is not None else result.datasets)
    params = result.params
    names = [n for n in result.var_names if params[n].vary]
    best = {n: float(params[n].value) for n in names}
    stderr = {n: float(params[n].stderr) if params[n].stderr is not None else float("nan") for n in names}

    lo = np.array([params[n].min if np.isfinite(params[n].min) else -np.inf for n in names], float)
    hi = np.array([params[n].max if np.isfinite(params[n].max) else np.inf for n in names], float)
    for name, (blo, bhi) in (bounds or {}).items():
        if name not in names:
            raise KeyError(f"bound given for unknown or fixed variable {name!r}")
        i = names.index(name)
        lo[i], hi[i] = max(lo[i], blo), min(hi[i], bhi)

    n_idp = float(result.n_independent)
    # Every Larch residual call re-creates the path parameters (create_path_params),
    # which drops the propagated stderr of path parameters. Snapshot them here.
    saved_path_stderr = {}
    for ds in datasets:
        for path in ds.paths.values():
            if path.params is None:
                continue
            for pname in ("s02", "e0", "deltar", "sigma2", "third", "fourth", "ei"):
                key = path.pathpar_name(pname)
                par = path.params.get(key)
                if par is not None and par.stderr is not None:
                    saved_path_stderr[key] = par.stderr
    try:
        r0 = resid(params, datasets=datasets)
        npts = int(r0.size)
        w = n_idp / npts

        def log_prob(theta):
            x, lns = theta[:-1], theta[-1]
            if np.any(x < lo) or np.any(x > hi) or not (-50 < lns < 50):
                return -np.inf
            for n, v in zip(names, x, strict=True):
                params[n].value = v
            r = resid(params, datasets=datasets)
            if not np.all(np.isfinite(r)):
                return -np.inf
            s2 = np.exp(2.0 * lns)
            return -0.5 * w * (np.sum(r * r) / s2 + npts * np.log(s2))

        rng = np.random.default_rng(seed)
        x0 = np.array([best[n] for n in names])
        scale = np.array([stderr[n] if np.isfinite(stderr[n]) and stderr[n] > 0 else 1e-3 * max(abs(best[n]), 1e-3)
                          for n in names])
        lns0 = 0.5 * np.log(np.mean(r0 * r0))
        p0 = np.column_stack([
            x0 + init_scale * scale * rng.standard_normal((nwalkers, len(names))),
            lns0 + 0.01 * rng.standard_normal(nwalkers),
        ])
        p0[:, :-1] = np.clip(p0[:, :-1], lo + 1e-12, hi - 1e-12)

        sampler = emcee.EnsembleSampler(nwalkers, len(names) + 1, log_prob)
        if seed is not None:  # emcee keeps its own RandomState; seed it for reproducible chains
            sampler.random_state = np.random.RandomState(seed).get_state()
        sampler.run_mcmc(p0, nsteps, progress=progress)
        chain = sampler.get_chain(discard=burn, thin=thin, flat=True)
        acceptance = float(np.mean(sampler.acceptance_fraction))
    finally:
        for n in names:  # restore best fit: the sampler mutates the shared Parameters
            params[n].value = best[n]
        resid(params, datasets=datasets)   # model arrays back at the best fit
        for key, std in saved_path_stderr.items():
            if key in params:
                params[key].stderr = std

    return MCMCResult(
        names=names, chain=chain, best=best, stderr=stderr, acceptance=acceptance,
        weight=w, n_idp=n_idp, n_pts=npts,
    )

"""Post-fit physical sanity checks on a result.json record.

Every rule produces a warning, not a hard failure: the analyst decides. The
defaults are lab conventions and can be overridden per fit template.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Flag:
    level: str      # "warning" or "info"
    code: str
    message: str

    def __str__(self) -> str:
        return f"[{self.level}] {self.code}: {self.message}"


@dataclass
class Rules:
    max_var_fraction: float = 2.0 / 3.0          # N_var <= fraction * N_idp
    sigma2_range: tuple[float, float] = (0.002, 0.03)   # Angstrom^2
    e0_abs_max: float = 10.0                     # eV
    s02_range: tuple[float, float] = (0.7, 1.05)
    max_abs_correlation: float = 0.90
    max_rel_cn_error: float = 0.25
    r_reference: dict[str, float] = field(default_factory=dict)  # path label -> expected R
    r_tolerance: float = 0.10                    # Angstrom


def check_result(result: dict, rules: Rules | None = None) -> list[Flag]:
    """Run all checks; returns a (possibly empty) list of flags."""
    rules = rules or Rules()
    flags: list[Flag] = []
    fit = result["fit"]
    nidp, nvar = fit["n_idp"], fit["n_var"]
    if nvar > rules.max_var_fraction * nidp:
        flags.append(Flag("warning", "nidp_budget",
                          f"N_var = {nvar} exceeds {rules.max_var_fraction:.2f} x N_idp = {rules.max_var_fraction * nidp:.1f}"))
    for key, c in (result.get("correlations") or {}).items():
        if abs(c) > rules.max_abs_correlation:
            a, b = key.split("|")
            flags.append(Flag("warning", "correlation", f"|corr({a}, {b})| = {abs(c):.2f} > {rules.max_abs_correlation}"))
    s02_ref = result.get("s02_ref")
    for ds in result["datasets"]:
        for p in ds["paths"]:
            lab = p["label"]
            pp = p.get("params", {})
            s2 = (pp.get("sigma2") or {}).get("value")
            if s2 is not None and not (rules.sigma2_range[0] <= s2 <= rules.sigma2_range[1]):
                lvl = "warning" if s2 <= 0 else "info"
                flags.append(Flag(lvl, "sigma2_range", f"{lab}: sigma2 = {s2:.4g} outside {rules.sigma2_range}"))
            e0 = (pp.get("e0") or {}).get("value")
            if e0 is not None and abs(e0) > rules.e0_abs_max:
                flags.append(Flag("warning", "e0_range", f"{lab}: |dE0| = {abs(e0):.2f} eV > {rules.e0_abs_max}"))
            if not s02_ref:
                amp = (pp.get("s02") or {}).get("value")
                if amp is not None and not (rules.s02_range[0] <= amp <= rules.s02_range[1]):
                    flags.append(Flag("info", "s02_range",
                                      f"{lab}: amplitude = {amp:.3f} outside {rules.s02_range} "
                                      "(expected if the amplitude encodes N x S0^2)"))
            cn = p.get("cn") or {}
            if cn.get("value") and cn.get("stderr"):
                rel = cn["stderr"] / abs(cn["value"])
                if rel > rules.max_rel_cn_error:
                    flags.append(Flag("warning", "cn_uncertainty", f"{lab}: CN relative error {rel:.0%} > {rules.max_rel_cn_error:.0%}"))
            if lab in rules.r_reference:
                dr = p["r"]["value"] - rules.r_reference[lab]
                if abs(dr) > rules.r_tolerance:
                    flags.append(Flag("warning", "r_reference",
                                      f"{lab}: R = {p['r']['value']:.3f} differs from expected {rules.r_reference[lab]:.3f} by {dr:+.3f} A"))
    return flags

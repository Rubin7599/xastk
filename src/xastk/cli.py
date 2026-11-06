"""Command-line interface: ``xastk <command> ...``."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__


def _cmd_nidp(a):
    from .stats import n_idp

    print(f"{n_idp(a.kmin, a.kmax, a.rmin, a.rmax, a.convention):.3f}")


def _cmd_hamilton(a):
    from .stats import hamilton_test

    print(hamilton_test(a.chi2_a, a.nvar_a, a.chi2_b, a.nvar_b, a.nidp, a.level))


def _cmd_table(a):
    from .io import load_result
    from .report import rows_from_results, to_csv, to_docx, to_latex

    rows = rows_from_results([load_result(p) for p in a.results], style=a.style)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    for fmt in a.to:
        if fmt == "csv":
            print(to_csv(rows, out.with_suffix(".csv")))
        elif fmt == "latex":
            p = out.with_suffix(".tex")
            p.write_text(to_latex(rows, caption=a.caption), encoding="utf-8")
            print(p)
        elif fmt == "docx":
            print(to_docx(rows, out.with_suffix(".docx"), caption=a.caption))


def _cmd_check(a):
    from .io import load_result
    from .reduce import check_result

    worst = 0
    for p in a.results:
        flags = check_result(load_result(p))
        print(f"{p}: {len(flags)} flag(s)")
        for f in flags:
            print(f"  {f}")
            if f.level == "warning":
                worst = 1
    return worst if a.strict else 0


def _cmd_wt(a):
    import numpy as np

    from .wt import find_lobes, morlet_wt, plot_wt_composite

    data = np.loadtxt(a.chi_file, comments="#")
    k, chi = data[:, a.kcol], data[:, a.chicol]
    wt = morlet_wt(k, chi, kweight=a.kweight, eta=a.eta, sigma=a.sigma, rmin=a.rmin, rmax=a.rmax,
                   kmin=a.kmin, kmax=a.kmax)
    lobes = find_lobes(wt, threshold=a.threshold)
    for lb in lobes:
        print(f"lobe: k = {lb.k:.2f} 1/A, R = {lb.r:.3f} A, relative amplitude = {lb.relative:.2f}")
    if a.out:
        fig = plot_wt_composite(wt, lobes=lobes, title=a.title or Path(a.chi_file).stem, rmax_plot=a.rmax)
        fig.savefig(a.out, dpi=300)
        print(a.out)


def _cmd_rbkg(a):
    from larch.io import read_ascii

    from .reduce import rbkg_stability

    g = read_ascii(a.data)
    energy = getattr(g, a.energy)
    if a.mu:
        mu = getattr(g, a.mu)
    else:
        import numpy as np

        mu = np.log(getattr(g, a.i0) / getattr(g, a.it))
    scan = rbkg_stability(energy, mu, kmin=a.kmin, kmax=a.kmax, kweight=a.kweight,
                          shell=(a.shell_rmin, a.shell_rmax))
    print(scan.table())


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="xastk", description="Reproducible, statistically validated XAS analysis.")
    p.add_argument("--version", action="version", version=f"xastk {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("nidp", help="number of independent points")
    for name in ("kmin", "kmax", "rmin", "rmax"):
        s.add_argument(name, type=float)
    s.add_argument("--convention", default="larch", choices=["larch", "stern", "nyquist"])
    s.set_defaults(func=_cmd_nidp)

    s = sub.add_parser("hamilton", help="Hamilton F-test for nested fits A (simpler) and B")
    s.add_argument("--chi2-a", type=float, required=True)
    s.add_argument("--nvar-a", type=int, required=True)
    s.add_argument("--chi2-b", type=float, required=True)
    s.add_argument("--nvar-b", type=int, required=True)
    s.add_argument("--nidp", type=float, required=True)
    s.add_argument("--level", type=float, default=0.95)
    s.set_defaults(func=_cmd_hamilton)

    s = sub.add_parser("table", help="SI fitting table from result.json files")
    s.add_argument("results", nargs="+")
    s.add_argument("--to", nargs="+", default=["csv", "latex"], choices=["csv", "latex", "docx"])
    s.add_argument("--out", default="exafs_table")
    s.add_argument("--style", default="paren", choices=["paren", "pm"])
    s.add_argument("--caption", default="EXAFS fitting results.")
    s.set_defaults(func=_cmd_table)

    s = sub.add_parser("check", help="physical sanity checks on result.json files")
    s.add_argument("results", nargs="+")
    s.add_argument("--strict", action="store_true", help="exit with status 1 if any warning is raised")
    s.set_defaults(func=_cmd_check)

    s = sub.add_parser("wt", help="Morlet wavelet transform of a two-column k, chi(k) file")
    s.add_argument("chi_file")
    s.add_argument("--kcol", type=int, default=0)
    s.add_argument("--chicol", type=int, default=1)
    s.add_argument("--kweight", type=float, default=2.0)
    s.add_argument("--eta", type=float, default=5.0)
    s.add_argument("--sigma", type=float, default=1.0)
    s.add_argument("--rmin", type=float, default=0.5)
    s.add_argument("--rmax", type=float, default=5.0)
    s.add_argument("--kmin", type=float, default=None)
    s.add_argument("--kmax", type=float, default=None)
    s.add_argument("--threshold", type=float, default=0.2)
    s.add_argument("--title", default=None)
    s.add_argument("--out", default=None, help="figure file (png, pdf, svg)")
    s.set_defaults(func=_cmd_wt)

    s = sub.add_parser("rbkg", help="AUTOBK R_bkg stability scan (needs xraylarch)")
    s.add_argument("data")
    s.add_argument("--energy", default="energy")
    s.add_argument("--mu", default=None, help="column holding mu(E); otherwise ln(i0/it)")
    s.add_argument("--i0", default="i0")
    s.add_argument("--it", default="itrans")
    s.add_argument("--kmin", type=float, default=2.0)
    s.add_argument("--kmax", type=float, default=14.0)
    s.add_argument("--kweight", type=float, default=2.0)
    s.add_argument("--shell-rmin", type=float, default=1.0)
    s.add_argument("--shell-rmax", type=float, default=3.2)
    s.set_defaults(func=_cmd_rbkg)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    rc = args.func(args)
    return int(rc or 0)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())

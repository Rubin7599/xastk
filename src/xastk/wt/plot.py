"""Publication-style WT-EXAFS composite: 2D map with marginal k and R projections."""

from __future__ import annotations

from collections.abc import Sequence

from .lobes import Lobe
from .morlet import WaveletResult


def plot_wt_composite(
    wt: WaveletResult,
    lobes: Sequence[Lobe] | None = None,
    cmap: str = "viridis",
    title: str | None = None,
    kmax_plot: float | None = None,
    rmax_plot: float | None = None,
    figsize: tuple[float, float] = (6.0, 5.2),
    annotate_parameters: bool = True,
):
    """Return a matplotlib Figure with the |W(k, R)| map and its marginals.

    Uses a perceptually uniform colormap by default. Lobe maxima, if given, are
    marked with crosses and labelled with (k, R).
    """
    import matplotlib.pyplot as plt

    M = wt.modulus
    fig = plt.figure(figsize=figsize, constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=(4, 1), height_ratios=(1, 4))
    ax = fig.add_subplot(gs[1, 0])
    axk = fig.add_subplot(gs[0, 0], sharex=ax)
    axr = fig.add_subplot(gs[1, 1], sharey=ax)

    mesh = ax.pcolormesh(wt.k, wt.r, M, shading="auto", cmap=cmap, rasterized=True)
    ax.set_xlabel(r"$k$ ($\mathrm{\AA^{-1}}$)")
    ax.set_ylabel(r"$R$ ($\mathrm{\AA}$, not phase corrected)")
    if kmax_plot:
        ax.set_xlim(wt.k[0], kmax_plot)
    if rmax_plot:
        ax.set_ylim(wt.r[0], rmax_plot)

    axk.plot(wt.k, wt.k_profile(), color="0.2", lw=1.0)
    axk.set_ylabel(r"max$_R$ |W|")
    axk.tick_params(labelbottom=False)
    axr.plot(wt.r_profile(), wt.r, color="0.2", lw=1.0)
    axr.set_xlabel(r"max$_k$ |W|")
    axr.tick_params(labelleft=False)

    for lb in lobes or []:
        ax.plot(lb.k, lb.r, marker="+", color="white", ms=10, mew=1.5)
        ax.annotate(f"({lb.k:.1f}, {lb.r:.2f})", (lb.k, lb.r), xytext=(6, 6),
                    textcoords="offset points", color="white", fontsize=8)
    if annotate_parameters:
        ax.text(0.02, 0.97, rf"Morlet $\eta$={wt.eta:g}, $\sigma$={wt.sigma:g}, $k^{{{wt.kweight:g}}}$",
                transform=ax.transAxes, va="top", ha="left", color="white", fontsize=8)
    if title:
        axk.set_title(title)
    fig.colorbar(mesh, ax=axr, label="|W|", shrink=0.8)
    return fig

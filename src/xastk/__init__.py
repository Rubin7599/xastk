"""xastk: reproducible, statistically validated XAS analysis on top of xraylarch.

Subpackages
-----------
stats   N_idp, Hamilton (F) test, MCMC posteriors for Larch feffit results
mcr     PCA rank estimation and MCR-ALS with closure and kinetic constraints
wt      Morlet wavelet transform of EXAFS, lobe detection, composite figures
report  uncertainty notation, SI fitting tables, OriginPro figure bundles
reduce  AUTOBK R_bkg stability scan and post-fit physical sanity checks
io      versioned result.json schema and conversion from Larch feffit output

Heavy dependencies (xraylarch, emcee, python-docx, originpro) are imported
lazily inside the functions that need them, so ``import xastk`` stays light.
"""

from ._version import __version__

__all__ = ["__version__"]

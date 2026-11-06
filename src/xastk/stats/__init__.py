"""Statistical validation of EXAFS fits.

- :func:`n_idp` -- number of independent points (Larch / Stern conventions)
- :func:`hamilton_test`, :func:`hamilton_from_results` -- nested-model F-test
- :func:`sample_feffit` -- MCMC posterior for a finished Larch ``feffit`` result
"""

from .hamilton import HamiltonResult, hamilton_from_results, hamilton_test
from .mcmc import MCMCResult, sample_feffit
from .nidp import n_idp, nvar_budget

__all__ = [
    "n_idp",
    "nvar_budget",
    "HamiltonResult",
    "hamilton_test",
    "hamilton_from_results",
    "MCMCResult",
    "sample_feffit",
]

"""Uncertainty notation for tables: 1.95(4) and 1.95 +/- 0.04.

Rounding uses :mod:`decimal` with ROUND_HALF_UP on the *decimal* string of
each number, so 0.975 rounds to 0.98 (binary floats would give 0.97).

Rules for the number of significant digits of the uncertainty
-------------------------------------------------------------
``"pdg"`` (default, Particle Data Group): look at the three leading digits h
of the uncertainty. 100-354 -> two digits; 355-949 -> one digit; 950-999 ->
round up to 1000 and keep two digits.
``"one"`` / ``"two"``: always one / two significant digits.
"""

from __future__ import annotations

import math
from decimal import ROUND_FLOOR, ROUND_HALF_UP, Decimal


def _dec(x) -> Decimal:
    return x if isinstance(x, Decimal) else Decimal(repr(float(x)))


def round_uncertainty(u, rule: str = "pdg") -> tuple[Decimal, int]:
    """Round an uncertainty; return (rounded value, exponent of its last digit)."""
    ud = _dec(u)
    if not (ud > 0) or not math.isfinite(float(ud)):
        raise ValueError("uncertainty must be a positive finite number")
    e = ud.adjusted()                                    # exponent of leading digit
    h = int((ud.scaleb(2 - e)).to_integral_value(rounding=ROUND_FLOOR))  # three leading digits
    if rule == "pdg":
        if h >= 950:
            last = e
            return Decimal(1).scaleb(e + 1).quantize(Decimal(1).scaleb(last)), last
        nsig = 2 if h <= 354 else 1
    elif rule == "one":
        nsig = 1
    elif rule == "two":
        nsig = 2
    else:
        raise ValueError(f"unknown rule {rule!r}")
    last = e - nsig + 1
    ur = ud.quantize(Decimal(1).scaleb(last), rounding=ROUND_HALF_UP)
    if ur.adjusted() > e and rule != "pdg":            # e.g. 0.096 -> 0.10 with rule "one"
        last += 1
        ur = ur.quantize(Decimal(1).scaleb(last), rounding=ROUND_HALF_UP)
    return ur, last


def _fmt_fixed(d: Decimal, last: int) -> str:
    q = d.quantize(Decimal(1).scaleb(last), rounding=ROUND_HALF_UP)
    s = format(q, "f")
    if s.startswith("-") and float(q) == 0:
        s = s[1:]
    return s


def _paren(ur: Decimal, last: int) -> str:
    if last < 0:
        return str(int(ur.scaleb(-last)))
    return format(ur, "f")


def format_value(value, unc, style: str = "paren", rule: str = "pdg") -> str:
    """Format value with one uncertainty.

    >>> format_value(1.9512, 0.0412)
    '1.95(4)'
    >>> format_value(1.9512, 0.0412, style="pm")
    '1.95 ± 0.04'
    """
    if unc is None or not (float(unc) > 0) or not math.isfinite(float(unc)):
        return f"{float(value):g}"
    ur, last = round_uncertainty(unc, rule)
    vs = _fmt_fixed(_dec(value), last)
    if style == "paren":
        return f"{vs}({_paren(ur, last)})"
    if style == "pm":
        return f"{vs} ± {_fmt_fixed(ur, last)}"
    if style == "latex":
        return rf"{vs} \pm {_fmt_fixed(ur, last)}"
    raise ValueError(f"unknown style {style!r}")


def format_two(value, unc_stat, unc_choice, rule: str = "pdg") -> str:
    """Value with statistical and choice (analyst-sensitivity) uncertainties: 4.1(3)(5).

    Precision is set by the smaller of the two uncertainties.
    """
    us, ls = round_uncertainty(unc_stat, rule)
    uc, lc = round_uncertainty(unc_choice, rule)
    last = min(ls, lc)
    vs = _fmt_fixed(_dec(value), last)
    us2 = _dec(unc_stat).quantize(Decimal(1).scaleb(last), rounding=ROUND_HALF_UP)
    uc2 = _dec(unc_choice).quantize(Decimal(1).scaleb(last), rounding=ROUND_HALF_UP)
    return f"{vs}({_paren(us2, last)})({_paren(uc2, last)})"

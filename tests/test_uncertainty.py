from decimal import Decimal

import pytest

from xastk.report import format_two, format_value, round_uncertainty


@pytest.mark.parametrize(
    "value, unc, paren, pm",
    [
        (1.9512, 0.0412, "1.95(4)", "1.95 ± 0.04"),
        (8.6, 1.0, "8.6(10)", "8.6 ± 1.0"),
        (4.2112, 0.51135, "4.2(5)", "4.2 ± 0.5"),
        (0.00910, 0.00031, "0.00910(31)", "0.00910 ± 0.00031"),   # PDG: 310 -> two digits
        (12.3, 0.096, "12.30(10)", "12.30 ± 0.10"),               # PDG: 960 -> round up, two digits
        (0.975, 0.0391, "0.98(4)", "0.98 ± 0.04"),                # half-up on the decimal value
        (-0.00546, 0.00261, "-0.0055(26)", "-0.0055 ± 0.0026"),
        (1234.0, 56.0, "1230(60)", "1230 ± 60"),                  # PDG: 560 -> one digit
        (1234.0, 16.0, "1234(16)", "1234 ± 16"),
        (12345.0, 640.0, "12300(600)", "12300 ± 600"),
    ],
)
def test_format_value(value, unc, paren, pm):
    assert format_value(value, unc) == paren
    assert format_value(value, unc, style="pm") == pm


def test_round_uncertainty_rules():
    assert round_uncertainty(0.0412, "two") == (Decimal("0.041"), -3)
    assert round_uncertainty(0.096, "one")[0] == Decimal("0.1")
    with pytest.raises(ValueError):
        round_uncertainty(0)


def test_missing_uncertainty_falls_back():
    assert format_value(0.85, None) == "0.85"
    assert format_value(0.85, 0.0) == "0.85"


def test_format_two():
    assert format_two(4.13, 0.3, 0.5, rule="one") == "4.1(3)(5)"
    assert format_two(4.13, 0.3, 0.5) == "4.13(30)(50)"      # PDG keeps two digits for 0.30
    assert format_two(2.0134, 0.012, 0.004) == "2.013(12)(4)"


def test_latex_style():
    assert format_value(1.9512, 0.0412, style="latex") == r"1.95 \pm 0.04"

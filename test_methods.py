import math
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from integration import trapezoidal, simpson
from optimization import golden_section_search
from airmass import airmass_kasten_young, observation_quality
from astronomy import altitude_deg, max_altitude_deg
from objective import total_quality, make_objective


# ---------- integration ----------
def test_trapezoidal_linear_exact():
    assert trapezoidal(lambda x: 2 * x + 1, 0, 3, 4) == pytest.approx(12.0)


def test_trapezoidal_sin():
    assert trapezoidal(math.sin, 0, math.pi, 1000) == pytest.approx(2.0, abs=1e-5)


def test_simpson_cubic_exact():
    # Simpson is exact for polynomials up to degree 3: int_0^2 x^3 = 4
    assert simpson(lambda x: x ** 3, 0, 2, 2) == pytest.approx(4.0)


def test_simpson_sin():
    assert simpson(math.sin, 0, math.pi, 100) == pytest.approx(2.0, abs=1e-7)


def test_simpson_more_accurate_than_trapezoid():
    exact = math.e - 1
    et = abs(trapezoidal(math.exp, 0, 1, 20) - exact)
    es = abs(simpson(math.exp, 0, 1, 20) - exact)
    assert es < et


def test_zero_width_interval():
    assert trapezoidal(math.sin, 1, 1, 10) == 0.0
    assert simpson(math.sin, 1, 1, 10) == 0.0


@pytest.mark.parametrize("rule", [trapezoidal, simpson])
def test_invalid_inputs(rule):
    with pytest.raises(ValueError):
        rule(math.sin, 0, 1, 0)
    with pytest.raises(ValueError):
        rule(math.sin, 0, 1, -4)
    with pytest.raises(TypeError):
        rule(math.sin, 0, 1, 2.5)
    with pytest.raises(TypeError):
        rule("not a function", 0, 1, 4)
    with pytest.raises(ValueError):
        rule(math.sin, 2, 1, 4)
    with pytest.raises(ValueError):
        rule(math.sin, 0, float("inf"), 4)


def test_simpson_odd_n():
    with pytest.raises(ValueError):
        simpson(math.sin, 0, 1, 3)


# ---------- optimization ----------
def test_golden_section_quadratic():
    x, fx = golden_section_search(lambda x: -(x - 2) ** 2 + 3, 0, 5, tol=1e-8)
    assert x == pytest.approx(2.0, abs=1e-6)
    assert fx == pytest.approx(3.0, abs=1e-10)


def test_golden_section_sin():
    x, fx = golden_section_search(math.sin, 0, math.pi, tol=1e-8)
    assert x == pytest.approx(math.pi / 2, abs=1e-6)
    assert fx == pytest.approx(1.0)


def test_golden_section_invalid():
    with pytest.raises(ValueError):
        golden_section_search(math.sin, 3, 1)
    with pytest.raises(ValueError):
        golden_section_search(math.sin, 0, 1, tol=0)


# ---------- airmass / astronomy ----------
def test_airmass_values():
    assert airmass_kasten_young(90) == pytest.approx(1.0, abs=1e-3)
    assert airmass_kasten_young(30) == pytest.approx(2.0, abs=0.01)
    assert airmass_kasten_young(0) == pytest.approx(38.0, abs=0.5)


def test_airmass_invalid():
    with pytest.raises(ValueError):
        airmass_kasten_young(-5)


def test_quality_properties():
    assert observation_quality(90) == pytest.approx(1.0, abs=1e-3)
    assert observation_quality(5) == 0.0                   # below min altitude
    assert observation_quality(60) > observation_quality(20)


def test_altitude_at_transit():
    # h_max = 90 - |phi - delta|
    assert altitude_deg(0, 40, 20) == pytest.approx(70.0)
    assert max_altitude_deg(40, 20) == pytest.approx(70.0)


# ---------- objective ----------
def test_objective_methods_agree():
    a = total_quality(-1, 2, 40, 20, method="trapezoidal", n=400)
    b = total_quality(-1, 2, 40, 20, method="simpson", n=400)
    assert a == pytest.approx(b, rel=1e-5)


def test_objective_invalid_method():
    with pytest.raises(ValueError):
        total_quality(0, 1, 40, 20, method="bad")


def test_golden_section_on_objective_centres_on_transit():
    # Quality is symmetric about transit, so the best exposure is centred on it:
    # t_start* = t_transit - duration / 2
    T = 2.0
    J = make_objective(T, 40.0, 20.0, t_transit=0.0, method="simpson", n=100)
    t_opt, J_opt = golden_section_search(J, -6.0, 6.0 - T, tol=1e-6)
    assert t_opt == pytest.approx(-T / 2, abs=1e-3)
    assert J_opt > 0

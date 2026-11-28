import math, os, sys
import pytest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import experiment as ex
from optimization import golden_section_search


def test_simpson_order_four_and_trapezoid_order_two_on_smooth_case():
    ref, rows = ex.convergence(20.0, -1.0, 2.0, [2 ** k for k in range(1, 9)])
    assert ex.observed_order(rows, 2) == pytest.approx(2.0, abs=0.15)
    assert ex.observed_order(rows, 3) == pytest.approx(4.0, abs=0.3)


def test_order_drops_to_about_one_when_cutoff_inside_interval():
    ref, rows = ex.convergence(-20.0, 2.5, 2.0, [2 ** k for k in range(1, 9)])
    assert 0.6 < ex.observed_order(rows, 2) < 1.5
    assert 0.6 < ex.observed_order(rows, 3) < 1.6


def test_golden_section_agrees_with_brute_force():
    t_gs, J_gs, _ = ex.optimize(20.0, 2.0, "simpson")
    t_bf, J_bf = ex.brute_force(20.0, 2.0)
    assert t_gs == pytest.approx(t_bf, abs=0.011)
    assert J_gs >= J_bf - 1e-9


@pytest.mark.parametrize("dec", [-20, 0, 20, 40, 60])
@pytest.mark.parametrize("T", [1.0, 2.0, 3.0])
def test_optimal_start_is_centred_on_transit(dec, T):
    for m in ("trapezoidal", "simpson"):
        t, _, _ = ex.optimize(dec, T, m)
        assert t == pytest.approx(-T / 2, abs=1e-4)

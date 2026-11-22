import math
import pytest

from root_finding import bisection, newton, secant


def test_bisection_sqrt2():
    root, it = bisection(lambda x: x**2 - 2, 0, 2)
    assert root == pytest.approx(math.sqrt(2), abs=1e-7)
    assert it > 0


def test_bisection_cos_x_minus_x():
    root, _ = bisection(lambda x: math.cos(x) - x, 0, 1)
    assert root == pytest.approx(0.7390851332, abs=1e-7)


def test_bisection_no_sign_change():
    with pytest.raises(ValueError):
        bisection(lambda x: x**2 + 1, -1, 1)


def test_bisection_invalid_bracket():
    with pytest.raises(ValueError):
        bisection(lambda x: x, 2, 1)


def test_bisection_invalid_tolerance():
    with pytest.raises(ValueError):
        bisection(lambda x: x**2 - 2, 0, 2, tol=-1)


def test_bisection_max_iter_exceeded():
    with pytest.raises(RuntimeError):
        bisection(lambda x: x**2 - 2, 0, 2, tol=1e-15, max_iter=3)


def test_newton_sqrt2():
    root, _ = newton(lambda x: x**2 - 2, lambda x: 2 * x, 1.0)
    assert root == pytest.approx(math.sqrt(2), abs=1e-9)


def test_newton_cubic():
    root, _ = newton(lambda x: x**3 - x - 2, lambda x: 3 * x**2 - 1, 1.5)
    assert root**3 - root - 2 == pytest.approx(0, abs=1e-8)


def test_newton_zero_derivative():
    with pytest.raises(ValueError):
        newton(lambda x: x**2 - 2, lambda x: 2 * x, 0.0)


def test_newton_invalid_max_iter():
    with pytest.raises(ValueError):
        newton(lambda x: x, lambda x: 1, 1.0, max_iter=0)


def test_secant_sqrt2():
    root, _ = secant(lambda x: x**2 - 2, 1.0, 2.0)
    assert root == pytest.approx(math.sqrt(2), abs=1e-9)


def test_secant_cos_x_minus_x():
    root, _ = secant(lambda x: math.cos(x) - x, 0.0, 1.0)
    assert root == pytest.approx(0.7390851332, abs=1e-9)


def test_secant_same_starting_points():
    with pytest.raises(ValueError):
        secant(lambda x: x**2 - 2, 1.0, 1.0)


def test_secant_equal_function_values():
    with pytest.raises(ValueError):
        secant(lambda x: 5, 1.0, 2.0)


def test_methods_agree():
    f = lambda x: x**3 - x - 2
    df = lambda x: 3 * x**2 - 1
    r1, _ = bisection(f, 1, 2)
    r2, _ = newton(f, df, 1.5)
    r3, _ = secant(f, 1.0, 2.0)
    assert r1 == pytest.approx(r2, abs=1e-6)
    assert r2 == pytest.approx(r3, abs=1e-6)

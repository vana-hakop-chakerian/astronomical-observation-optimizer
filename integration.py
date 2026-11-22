"""Numerical integration implemented from scratch."""
import math


def _validate(f, a, b, n):
    if not callable(f):
        raise TypeError("f must be callable")
    if isinstance(n, bool) or not isinstance(n, int):
        raise TypeError("n must be an integer")
    if n <= 0:
        raise ValueError("n must be a positive integer")
    if not (math.isfinite(a) and math.isfinite(b)):
        raise ValueError("integration limits must be finite")
    if a > b:
        raise ValueError("require a <= b")


def trapezoidal(f, a, b, n):
    """Composite trapezoidal rule with n subintervals.

    I ~ h * [ f(a)/2 + f(x_1) + ... + f(x_{n-1}) + f(b)/2 ],  h = (b-a)/n
    Error is O(h^2).
    """
    _validate(f, a, b, n)
    if a == b:
        return 0.0
    h = (b - a) / n
    s = 0.5 * (f(a) + f(b))
    for i in range(1, n):
        s += f(a + i * h)
    return h * s


def simpson(f, a, b, n):
    """Composite Simpson's 1/3 rule with n subintervals (n must be even).

    I ~ h/3 * [ f(x_0) + 4 f(x_1) + 2 f(x_2) + ... + 4 f(x_{n-1}) + f(x_n) ]
    Error is O(h^4).
    """
    _validate(f, a, b, n)
    if n % 2 != 0:
        raise ValueError("Simpson's rule requires an even number of subintervals")
    if a == b:
        return 0.0
    h = (b - a) / n
    s = f(a) + f(b)
    for i in range(1, n):
        s += (4 if i % 2 == 1 else 2) * f(a + i * h)
    return h * s / 3.0

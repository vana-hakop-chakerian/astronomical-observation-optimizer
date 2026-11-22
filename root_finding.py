"""Root-finding methods: Bisection, Newton's method and Secant method.

Each function returns a tuple (root, iterations) so that convergence
can be compared later (Week 3).
"""

import math


def _validate_common(tol, max_iter):
    if not isinstance(tol, (int, float)) or tol <= 0:
        raise ValueError("tol must be a positive number")
    if not isinstance(max_iter, int) or max_iter <= 0:
        raise ValueError("max_iter must be a positive integer")


def bisection(f, a, b, tol=1e-8, max_iter=200):
    """Find a root of f in [a, b] using the Bisection method.

    Requires f(a) and f(b) to have opposite signs.
    Returns (root, iterations).
    """
    _validate_common(tol, max_iter)
    if not a < b:
        raise ValueError("a must be smaller than b")

    fa, fb = f(a), f(b)
    if fa == 0:
        return a, 0
    if fb == 0:
        return b, 0
    if fa * fb > 0:
        raise ValueError("f(a) and f(b) must have opposite signs")

    for i in range(1, max_iter + 1):
        mid = (a + b) / 2
        fm = f(mid)
        if fm == 0 or (b - a) / 2 < tol:
            return mid, i
        if fa * fm < 0:
            b, fb = mid, fm
        else:
            a, fa = mid, fm

    raise RuntimeError("Bisection did not converge within max_iter iterations")


def newton(f, df, x0, tol=1e-10, max_iter=100):
    """Find a root of f using Newton's method, starting at x0.

    df is the derivative of f. Returns (root, iterations).
    """
    _validate_common(tol, max_iter)

    x = x0
    for i in range(1, max_iter + 1):
        fx = f(x)
        dfx = df(x)
        if dfx == 0:
            raise ValueError("Derivative is zero; Newton's method cannot continue")
        x_new = x - fx / dfx
        if abs(x_new - x) < tol:
            return x_new, i
        x = x_new

    raise RuntimeError("Newton's method did not converge within max_iter iterations")


def secant(f, x0, x1, tol=1e-10, max_iter=100):
    """Find a root of f using the Secant method with starting points x0, x1.

    Returns (root, iterations).
    """
    _validate_common(tol, max_iter)
    if x0 == x1:
        raise ValueError("x0 and x1 must be different")

    f0, f1 = f(x0), f(x1)
    for i in range(1, max_iter + 1):
        if f1 == f0:
            raise ValueError("f(x0) equals f(x1); Secant method cannot continue")
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        if abs(x2 - x1) < tol:
            return x2, i
        x0, f0 = x1, f1
        x1, f1 = x2, f(x2)

    raise RuntimeError("Secant method did not converge within max_iter iterations")

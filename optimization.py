"""Golden-section search (maximization of a unimodal function)."""
import math

INV_PHI = (math.sqrt(5.0) - 1.0) / 2.0   # 1/phi ~ 0.6180339887


def golden_section_search(f, a, b, tol=1e-6, max_iter=500):
    """Find the maximum of a unimodal function f on [a, b].

    Parameters
    ----------
    f : callable      objective function to MAXIMIZE
    a, b : float      search interval, a < b
    tol : float       stop when the interval length b - a < tol
    max_iter : int    safety limit on iterations

    Returns
    -------
    (x_opt, f_opt) : estimated maximizer and the function value there.
    """
    if not callable(f):
        raise TypeError("f must be callable")
    if not (math.isfinite(a) and math.isfinite(b)):
        raise ValueError("interval limits must be finite")
    if a >= b:
        raise ValueError("require a < b")
    if tol <= 0:
        raise ValueError("tol must be positive")

    c = b - INV_PHI * (b - a)
    d = a + INV_PHI * (b - a)
    fc, fd = f(c), f(d)
    it = 0
    while (b - a) > tol and it < max_iter:
        if fc > fd:            # maximum lies in [a, d]
            b, d, fd = d, c, fc
            c = b - INV_PHI * (b - a)
            fc = f(c)
        else:                  # maximum lies in [c, b]
            a, c, fc = c, d, fd
            d = a + INV_PHI * (b - a)
            fd = f(d)
        it += 1
    x = 0.5 * (a + b)
    return x, f(x)

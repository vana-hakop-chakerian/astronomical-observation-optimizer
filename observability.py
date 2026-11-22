"""Observability constraints: minimum altitude, astronomical night and
observable time intervals.

This module does NOT recompute the altitude formula. Pass in the altitude
function already defined in astronomy.py, as a function of time t
(for example hours relative to transit) returning altitude in degrees.
"""

from root_finding import bisection

ASTRONOMICAL_NIGHT_SUN_ALTITUDE = -18.0  # degrees


def is_astronomical_night(sun_altitude_deg):
    """True if the Sun is more than 18 degrees below the horizon."""
    return sun_altitude_deg < ASTRONOMICAL_NIGHT_SUN_ALTITUDE


def altitude_constraint(altitude_fn, h_min):
    """Return g(t) = altitude(t) - h_min. Its roots are the times when the
    object crosses the minimum altitude."""
    return lambda t: altitude_fn(t) - h_min


def find_boundary_crossings(altitude_fn, h_min, t_start, t_end,
                            samples=200, tol=1e-8):
    """Find the times in [t_start, t_end] where altitude crosses h_min.

    The interval is scanned on a grid to locate sign changes, then each
    one is refined with the Bisection method. Returns a sorted list.
    """
    if not t_start < t_end:
        raise ValueError("t_start must be smaller than t_end")
    if samples < 2:
        raise ValueError("samples must be at least 2")

    g = altitude_constraint(altitude_fn, h_min)
    step = (t_end - t_start) / samples
    crossings = []
    for k in range(samples):
        a = t_start + k * step
        b = a + step
        ga, gb = g(a), g(b)
        if ga == 0:
            crossings.append(a)
        elif ga * gb < 0:
            root, _ = bisection(g, a, b, tol=tol)
            crossings.append(root)
    if g(t_end) == 0:
        crossings.append(t_end)
    return sorted(set(crossings))


def observable_intervals(altitude_fn, h_min, t_start, t_end, samples=200):
    """Return a list of (start, end) intervals where altitude >= h_min."""
    crossings = find_boundary_crossings(altitude_fn, h_min, t_start, t_end, samples)
    points = [t_start] + crossings + [t_end]
    intervals = []
    for a, b in zip(points[:-1], points[1:]):
        if b > a and altitude_fn((a + b) / 2) >= h_min:
            intervals.append((a, b))
    return intervals

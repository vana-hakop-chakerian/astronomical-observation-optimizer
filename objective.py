"""Objective function: total observation quality over a fixed exposure.

    J(t_s) = integral_{t_s}^{t_s + T}  Q( h(t) ) dt

where h(t) is the object's altitude (astronomy.py) and Q the quality based on
airmass (airmass.py). J is evaluated by numerical integration and will be
maximized over t_s with golden-section search (optimization.py).
"""
from astronomy import altitude_deg
from airmass import observation_quality, DEFAULT_K, DEFAULT_MIN_ALT
from integration import trapezoidal, simpson

METHODS = {"trapezoidal": trapezoidal, "simpson": simpson}


def quality_at_time(t, latitude_deg, declination_deg, t_transit=0.0,
                    k=DEFAULT_K, min_altitude_deg=DEFAULT_MIN_ALT):
    """Q(h(t)) at time t (hours)."""
    h = altitude_deg(t, latitude_deg, declination_deg, t_transit)
    return observation_quality(h, k, min_altitude_deg)


def total_quality(t_start, duration, latitude_deg, declination_deg,
                  t_transit=0.0, k=DEFAULT_K, min_altitude_deg=DEFAULT_MIN_ALT,
                  method="simpson", n=100):
    """Total observation quality J(t_start) for an exposure of `duration` hours."""
    if method not in METHODS:
        raise ValueError(f"method must be one of {list(METHODS)}")
    if duration <= 0:
        raise ValueError("duration must be positive")

    def q(t):
        return quality_at_time(t, latitude_deg, declination_deg, t_transit,
                               k, min_altitude_deg)

    return METHODS[method](q, t_start, t_start + duration, n)


def make_objective(duration, latitude_deg, declination_deg, t_transit=0.0,
                   k=DEFAULT_K, min_altitude_deg=DEFAULT_MIN_ALT,
                   method="simpson", n=100):
    """Return J(t_start) as a one-argument function, ready for golden_section_search."""
    def J(t_start):
        return total_quality(t_start, duration, latitude_deg, declination_deg,
                             t_transit, k, min_altitude_deg, method, n)
    return J

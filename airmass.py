"""Airmass and observation-quality model.

Airmass X(h): relative path length of light through the atmosphere
(X = 1 at the zenith). We use the Kasten & Young (1989) formula, valid down to
the horizon (h in degrees):

    X(h) = 1 / ( sin(h) + 0.50572 * (h + 6.07995)^(-1.6364) )

Extinction: the transmitted flux is  T(X) = 10^(-0.4 * k * X)  where k is the
extinction coefficient (mag / airmass). Relative to the zenith, the observation
quality is

    Q(h) = 10^(-0.4 * k * (X(h) - 1))     for h >= h_min
    Q(h) = 0                              for h <  h_min

so Q = 1 at the zenith (best) and decreases as the object gets lower.
"""
import math

DEFAULT_K = 0.20          # extinction coefficient, mag/airmass (typical V band)
DEFAULT_MIN_ALT = 10.0    # degrees, below this the object is not observable


def airmass_kasten_young(altitude_deg):
    """Airmass from altitude in degrees (0 <= h <= 90)."""
    if not (0.0 <= altitude_deg <= 90.0):
        raise ValueError("altitude must be in [0, 90] degrees")
    h = altitude_deg
    return 1.0 / (math.sin(math.radians(h)) + 0.50572 * (h + 6.07995) ** -1.6364)


def airmass_plane_parallel(altitude_deg):
    """Simple approximation X = sec(z) = 1/sin(h); invalid near the horizon."""
    if not (0.0 < altitude_deg <= 90.0):
        raise ValueError("altitude must be in (0, 90] degrees")
    return 1.0 / math.sin(math.radians(altitude_deg))


def observation_quality(altitude_deg, k=DEFAULT_K, min_altitude_deg=DEFAULT_MIN_ALT):
    """Relative observation quality Q in [0, 1] for an object at given altitude."""
    if k < 0:
        raise ValueError("extinction coefficient k must be non-negative")
    if altitude_deg < min_altitude_deg:
        return 0.0
    X = airmass_kasten_young(min(altitude_deg, 90.0))
    return 10.0 ** (-0.4 * k * (X - 1.0))

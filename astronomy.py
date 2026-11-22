"""Basic astronomical geometry: altitude of an object as a function of time.

Conventions
-----------
* Time ``t`` is in hours; ``t_transit`` is the time (hours) at which the object
  crosses the local meridian (hour angle H = 0).
* Latitude, declination and altitude are in degrees.
* The hour angle grows by 15 degrees per hour (Earth's rotation, sidereal
  vs. solar difference neglected).
"""
import math

DEG_PER_HOUR = 15.0


def hour_angle_deg(t, t_transit=0.0):
    """Hour angle H (degrees) at time t (hours)."""
    return DEG_PER_HOUR * (t - t_transit)


def altitude_deg(t, latitude_deg, declination_deg, t_transit=0.0):
    """Altitude h (degrees above the horizon) of the object at time t.

    sin(h) = sin(phi) sin(delta) + cos(phi) cos(delta) cos(H)
    """
    phi = math.radians(latitude_deg)
    delta = math.radians(declination_deg)
    H = math.radians(hour_angle_deg(t, t_transit))
    sin_h = math.sin(phi) * math.sin(delta) + math.cos(phi) * math.cos(delta) * math.cos(H)
    sin_h = max(-1.0, min(1.0, sin_h))  # guard against round-off
    return math.degrees(math.asin(sin_h))


def max_altitude_deg(latitude_deg, declination_deg):
    """Altitude at transit (the highest the object ever gets)."""
    return altitude_deg(0.0, latitude_deg, declination_deg, 0.0)

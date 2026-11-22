import pytest

from observability import (
    is_astronomical_night,
    altitude_constraint,
    find_boundary_crossings,
    observable_intervals,
)


def fake_altitude(t):
    """Simple test altitude: peaks at 60 deg at t = 0, falls 5 deg per hour^2."""
    return 60 - 5 * t**2


def test_astronomical_night_true():
    assert is_astronomical_night(-20.0)


def test_astronomical_night_boundary_is_false():
    assert not is_astronomical_night(-18.0)


def test_astronomical_night_false():
    assert not is_astronomical_night(-5.0)


def test_altitude_constraint():
    g = altitude_constraint(fake_altitude, 30)
    assert g(0) == pytest.approx(30)


def test_boundary_crossings():
    # 60 - 5 t^2 = 30  ->  t = +/- sqrt(6)
    crossings = find_boundary_crossings(fake_altitude, 30, -5, 5)
    assert len(crossings) == 2
    assert crossings[0] == pytest.approx(-(6 ** 0.5), abs=1e-6)
    assert crossings[1] == pytest.approx(6 ** 0.5, abs=1e-6)


def test_observable_interval():
    intervals = observable_intervals(fake_altitude, 30, -5, 5)
    assert len(intervals) == 1
    start, end = intervals[0]
    assert start == pytest.approx(-(6 ** 0.5), abs=1e-6)
    assert end == pytest.approx(6 ** 0.5, abs=1e-6)


def test_never_observable():
    assert observable_intervals(fake_altitude, 70, -5, 5) == []


def test_invalid_time_range():
    with pytest.raises(ValueError):
        find_boundary_crossings(fake_altitude, 30, 5, -5)

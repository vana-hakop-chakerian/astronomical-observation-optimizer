
# Astronomical Observation Optimizer

Optimizing the observation time of an astronomical object using numerical methods.

Given an observer's location and an object's coordinates, this project finds **when the object can be observed** (minimum altitude and astronomical night) and **the best time to start an exposure** (maximum integrated observation quality).

## Team

| Member | Part | Methods |
|---|---|---|
| Vana Hakop Chakerian | Optimal observation time | Numerical integration (Trapezoidal Rule, Simpson's Rule), optimization (Golden-Section Search) |
| Angelina Nersesyan | Observable time interval | Root finding (Bisection, Newton, Secant), observability constraints |

## How the Two Parts Fit Together

1. **Observable window (Angelina):** Root-finding methods are used to solve the altitude constraint `altitude(t) - h_min = 0` and determine when the object crosses the minimum required altitude. Combined with the astronomical night condition (the Sun is below -18 degrees), this determines the time interval during which the object is observable.

2. **Optimal observation time (Vana):** Within the observable interval, numerical integration calculates the total observation quality over an exposure duration `T`. Golden-Section Search then finds the exposure start time `t_s` that maximizes the integrated quality `J(t_s)`.

## Project Structure

```text
.
├── astronomy.py          # Altitude of an object (latitude, declination, hour angle)
├── airmass.py            # Kasten-Young airmass formula
├── integration.py        # Composite Trapezoidal Rule and Composite Simpson's Rule
├── optimization.py       # Golden-Section Search
├── objective.py          # Observation quality Q(h) and objective function J(t_s)
├── root_finding.py       # Bisection, Newton's method, Secant method
├── observability.py      # Minimum altitude, astronomical night, observable intervals
├── tests/
│   ├── test_observability.py
│   └── test_root_finding.py
└── README.md
```

## Modules

### Optimal Observation Time (Vana)

- **`astronomy.py`**: Calculates the object's altitude using the observer's latitude, the object's declination, and its hour angle.
- **`airmass.py`**: Implements the Kasten-Young airmass formula `X(h)`.
- **`objective.py`**: Defines the observation quality function `Q(h) = 10^(-0.4 * k * (X(h) - 1))`, where `k` is the extinction coefficient, and the objective function `J(t_s)`, which integrates `Q(h(t))` over the exposure interval from `t_s` to `t_s + T`.
- **`integration.py`**: Implements the Composite Trapezoidal Rule and Composite Simpson's Rule, including input validation.
- **`optimization.py`**: Implements Golden-Section Search with a stopping tolerance to estimate the optimal exposure start time and its function value.

### Observable Time Interval (Angelina)

- **`root_finding.py`**:
  - `bisection(f, a, b, tol, max_iter)`
  - `newton(f, df, x0, tol, max_iter)`
  - `secant(f, x0, x1, tol, max_iter)`

  Each method returns `(root, iterations)` and validates its inputs.

- **`observability.py`**:
  - `is_astronomical_night(sun_altitude_deg)`: Checks whether the Sun is below -18 degrees.
  - `altitude_constraint(altitude_fn, h_min)`: Defines the altitude constraint function `altitude(t) - h_min`.
  - `find_boundary_crossings(altitude_fn, h_min, t_start, t_end)`: Finds the times when the object crosses the minimum altitude.
  - `observable_intervals(altitude_fn, h_min, t_start, t_end)`: Determines intervals during which the object's altitude is at least `h_min`.

## Getting Started

Clone the repository:

```bash
git clone https://github.com/vana-hakop-chakerian/astronomical-observation-optimizer.git
cd astronomical-observation-optimizer
```

Install the required dependencies:

```bash
pip install -r requirements.txt
```

## Running the Tests

Run the test suite using:

```bash
python -m pytest -v
```

## Current Status

| Component | Status |
|---|---|
| Altitude model, airmass, quality function, objective | Done, tested |
| Trapezoidal and Simpson integration | Done, tested |
| Golden-Section Search | Done, tested |
| Bisection, Newton, Secant | Done, tested |
| Observability constraints and intervals | Done, tested |
| Final experiment with a real target and observatory | Planned |
| Error comparison, convergence comparison, plots and tables | Planned |
| Final report | Planned |

## Planned Work

- Run the final experiment with a real astronomical target and observatory.
- Compare the errors of the Trapezoidal and Simpson's Rules for different numbers of subintervals.
- Compare the convergence and accuracy of Bisection, Newton, and Secant methods.
- Calculate boundary crossing times and analyze the observable intervals.
- Apply Golden-Section Search to find optimal exposure start times under different conditions.
- Generate plots and tables and prepare the final report.

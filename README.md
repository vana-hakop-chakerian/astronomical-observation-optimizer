# Optimizing the Observation Time of an Astronomical Object

Numerical Analysis course project. This repository contains the **Optimal Observation Time**
component (optimization + numerical integration), status: **Second weekly progress**.

## Problem
Given an object (declination δ), an observatory (latitude φ) and a fixed exposure duration T,
choose the exposure **start time** t_s that maximizes the total observation quality.

## Mathematical model
| Quantity | Formula |
|---|---|
| Hour angle | H(t) = 15° · (t − t_transit), t in hours |
| Altitude | sin h = sin φ sin δ + cos φ cos δ cos H |
| Airmass (Kasten–Young 1989) | X(h) = 1 / (sin h + 0.50572 (h + 6.07995)^−1.6364) |
| Quality | Q(h) = 10^(−0.4 k (X(h) − 1)) if h ≥ h_min, else 0 |
| **Objective** | **J(t_s) = ∫ from t_s to t_s+T of Q(h(t)) dt** |
| Optimum | t_s* = argmax J(t_s) |

Assumptions: Q = 1 at zenith; extinction coefficient k = 0.20 mag/airmass; minimum usable
altitude h_min = 10°; solar/sidereal time difference, refraction and moon/weather neglected;
J is unimodal over the search interval.

## Modules
| File | Content |
|---|---|
| `astronomy.py` | hour angle and altitude of an object |
| `airmass.py` | airmass models and observation-quality function |
| `integration.py` | composite trapezoidal and Simpson rules (from scratch, with input validation) |
| `optimization.py` | golden-section search (maximization, from scratch) |
| `objective.py` | J(t_s), integration method selectable (`"trapezoidal"` / `"simpson"`) |
| `tests/test_methods.py` | tests with known integrals / known maximum |
| `main_demo.py` | small end-to-end demo |

## Usage
```bash
pip install -r requirements.txt
python -m pytest -v
python main_demo.py
```

```python
from objective import make_objective
from optimization import golden_section_search

J = make_objective(duration=2.0, latitude_deg=40, declination_deg=20, method="simpson", n=100)
t_opt, J_opt = golden_section_search(J, -6, 4, tol=1e-6)
```

## Current status
- [x] Mathematical model and objective function
- [x] Composite trapezoidal and Simpson rules + tests (sin, x³, exp)
- [x] Golden-section search + tests
- [x] Objective function linking integration and optimization
- [ ] Final astronomical experiment, error comparison, plots, tables (third reporting period)

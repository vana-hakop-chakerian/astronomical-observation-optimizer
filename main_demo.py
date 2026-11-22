"""Quick demonstration of the full pipeline (not the final experiment)."""
from objective import make_objective
from optimization import golden_section_search

LAT, DEC, T_TRANSIT, DURATION = 40.0, 20.0, 0.0, 2.0   # deg, deg, h, h

for method in ("trapezoidal", "simpson"):
    J = make_objective(DURATION, LAT, DEC, T_TRANSIT, method=method, n=100)
    t_opt, J_opt = golden_section_search(J, T_TRANSIT - 6, T_TRANSIT + 6 - DURATION, tol=1e-6)
    print(f"{method:12s}: optimal start = {t_opt:+.5f} h, total quality = {J_opt:.6f}")

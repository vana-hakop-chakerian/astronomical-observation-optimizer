"""Week 3 experiment: error comparison, optimization results, plots and tables.

Run:  python experiment.py      (needs numpy + matplotlib only for the plots)
Outputs go to the results/ folder (PNG figures and CSV tables).
"""
import csv
import math
import os

from astronomy import altitude_deg, max_altitude_deg
from airmass import airmass_kasten_young, observation_quality
from integration import trapezoidal, simpson
from objective import make_objective, total_quality
from optimization import golden_section_search, INV_PHI

OUT = "results"
LAT = 40.335            # Byurakan Observatory, Armenia (deg N)
T_TRANSIT = 0.0
K, H_MIN = 0.20, 10.0
PHI = 1 / INV_PHI


def q_of_t(dec):
    return lambda t: observation_quality(altitude_deg(t, LAT, dec, T_TRANSIT), K, H_MIN)


def write_csv(name, header, rows):
    with open(os.path.join(OUT, name), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


# ---------------------------------------------------------------- 1. convergence
def convergence(dec, t_s, T, ns, ref_n=2 ** 17):
    q = q_of_t(dec)
    ref = simpson(q, t_s, t_s + T, ref_n)
    rows = []
    for n in ns:
        et = abs(trapezoidal(q, t_s, t_s + T, n) - ref)
        es = abs(simpson(q, t_s, t_s + T, n) - ref)
        rows.append((n, (T) / n, et, es))
    return ref, rows


def observed_order(rows, col):
    """Slope of log(error) vs log(h) over the points with error above round-off."""
    pts = [(math.log(r[1]), math.log(r[col])) for r in rows if r[col] > 1e-13]
    xs, ys = zip(*pts)
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    return sum((x - mx) * (y - my) for x, y in pts) / sum((x - mx) ** 2 for x in xs)


# ---------------------------------------------------------------- 2. optimization
def optimize(dec, T, method, tol=1e-6, n=200):
    calls = [0]
    J0 = make_objective(T, LAT, dec, T_TRANSIT, K, H_MIN, method, n)

    def J(t):
        calls[0] += 1
        return J0(t)

    t_opt, J_opt = golden_section_search(J, -6.0, 6.0 - T, tol=tol)
    return t_opt, J_opt, calls[0]


def brute_force(dec, T, method="simpson", n=200, step=0.01):
    J = make_objective(T, LAT, dec, T_TRANSIT, K, H_MIN, method, n)
    pts = [(-6 + i * step, J(-6 + i * step)) for i in range(int((12 - T) / step) + 1)
           if -6 + i * step <= 6 - T + 1e-12]
    return max(pts, key=lambda p: p[1])


def main():
    os.makedirs(OUT, exist_ok=True)
    ns = [2 ** k for k in range(1, 11)]

    # ---- convergence: smooth case and case with the altitude cutoff inside the interval
    cases = {
        "smooth (dec=+20, t_s=-1, T=2)": (20.0, -1.0, 2.0),
        "cutoff inside (dec=-20, t_s=2.5, T=2)": (-20.0, 2.5, 2.0),
    }
    conv = {}
    for label, (dec, t_s, T) in cases.items():
        ref, rows = convergence(dec, t_s, T, ns)
        conv[label] = (ref, rows)
        write_csv("convergence_" + ("smooth" if "smooth" in label else "cutoff") + ".csv",
                  ["n", "h", "error_trapezoidal", "error_simpson"], rows)
        print(f"\n{label}: reference J = {ref:.12f}")
        print(f"  observed order: trapezoidal {observed_order(rows, 2):.2f}, "
              f"simpson {observed_order(rows, 3):.2f}")

    # ---- optimization table
    decs = [-20, 0, 20, 40, 60]
    Ts = [1.0, 2.0, 3.0, 4.0]
    opt_rows = []
    for dec in decs:
        for T in Ts:
            ts_t, J_t, c_t = optimize(dec, T, "trapezoidal")
            ts_s, J_s, c_s = optimize(dec, T, "simpson")
            opt_rows.append((dec, T, ts_t, ts_s, -T / 2, J_t, J_s, c_s))
    write_csv("optimization_results.csv",
              ["dec_deg", "T_h", "t_opt_trap", "t_opt_simpson", "t_opt_theory",
               "J_opt_trap", "J_opt_simpson", "evaluations"], opt_rows)
    print("\n dec    T   t_opt(trap) t_opt(simp) theory     J_opt(simp)  evals")
    for r in opt_rows:
        print(f"{r[0]:4d} {r[1]:4.1f} {r[2]:+11.6f} {r[3]:+11.6f} {r[4]:+7.2f} {r[6]:12.6f} {r[7]:5d}")

    # ---- tolerance study (iterations vs tol) and brute-force validation
    tol_rows = []
    for tol in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8]:
        t, J, calls = optimize(20.0, 2.0, "simpson", tol=tol)
        L = 12 - 2.0
        theory = math.ceil(math.log(L / tol) / math.log(PHI))
        tol_rows.append((tol, abs(t + 1.0), calls, theory))
    write_csv("tolerance_study.csv", ["tol", "abs_error_t_opt", "J_evaluations", "theory_iterations"], tol_rows)
    print("\n tol      |t_opt - theory|   evals  theory-iter")
    for r in tol_rows:
        print(f"{r[0]:8.0e} {r[1]:14.2e} {r[2]:8d} {r[3]:8d}")

    bf_t, bf_J = brute_force(20.0, 2.0)
    gs_t, gs_J, gs_c = optimize(20.0, 2.0, "simpson")
    print(f"\nBrute force (step 0.01 h): t={bf_t:+.3f}, J={bf_J:.6f}   "
          f"(~1000 evaluations)\nGolden section          : t={gs_t:+.6f}, J={gs_J:.6f} ({gs_c} evaluations)")

    make_plots(conv, decs)


# ---------------------------------------------------------------- plots
def make_plots(conv, decs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    dec = 20.0
    ts = [-6 + 0.05 * i for i in range(241)]
    alt = [altitude_deg(t, LAT, dec, T_TRANSIT) for t in ts]
    X = [airmass_kasten_young(max(a, 0)) if a > 0 else float("nan") for a in alt]
    Q = [observation_quality(a, K, H_MIN) for a in alt]

    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
    ax[0].plot(ts, alt); ax[0].axhline(H_MIN, ls="--", c="gray"); ax[0].set(title="Altitude h(t)", xlabel="time from transit (h)", ylabel="degrees")
    ax[1].plot(ts, X); ax[1].set(title="Airmass X(t)", xlabel="time from transit (h)", ylim=(0.9, 8))
    ax[2].plot(ts, Q); ax[2].set(title="Quality Q(t)", xlabel="time from transit (h)", ylim=(0, 1.05))
    fig.suptitle(f"Byurakan (lat {LAT}°), object at dec = +{dec:.0f}°"); fig.tight_layout()
    fig.savefig(f"{OUT}/fig1_altitude_airmass_quality.png", dpi=160); plt.close(fig)

    # J(t_s) with optimum
    T = 2.0
    J = make_objective(T, LAT, dec, T_TRANSIT, K, H_MIN, "simpson", 200)
    grid = [-6 + 0.05 * i for i in range(int((12 - T) / 0.05) + 1)]
    t_opt, J_opt = golden_section_search(J, -6, 6 - T, 1e-6)
    fig, a = plt.subplots(figsize=(6.5, 4))
    a.plot(grid, [J(t) for t in grid]); a.plot([t_opt], [J_opt], "ro", label=f"optimum t_s = {t_opt:.3f} h")
    a.set(title="Objective J(t_s), T = 2 h", xlabel="exposure start t_s (h from transit)", ylabel="total quality (h)")
    a.legend(); fig.tight_layout(); fig.savefig(f"{OUT}/fig2_objective.png", dpi=160); plt.close(fig)

    # convergence
    fig, ax = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for a, (label, (ref, rows)) in zip(ax, conv.items()):
        h = [r[1] for r in rows]
        a.loglog(h, [max(r[2], 1e-17) for r in rows], "o-", label="Trapezoidal")
        a.loglog(h, [max(r[3], 1e-17) for r in rows], "s-", label="Simpson")
        h0 = h[3]
        a.loglog(h, [rows[3][2] * (x / h0) ** 2 for x in h], "k:", lw=1, label="slope 2")
        a.loglog(h, [rows[3][3] * (x / h0) ** 4 for x in h], "k--", lw=1, label="slope 4")
        a.loglog(h, [rows[3][2] * (x / h0) for x in h], "gray", lw=1, ls="-.", label="slope 1")
        a.set(title=label, xlabel="step size h (h)"); a.grid(True, which="both", alpha=.3)
    ax[0].set_ylabel("absolute error"); ax[0].legend(); ax[0].set_ylim(1e-16, 1)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig3_convergence.png", dpi=160); plt.close(fig)

    # optimum quality vs declination
    fig, a = plt.subplots(figsize=(6.5, 4))
    for T in (1.0, 2.0, 4.0):
        vals = [optimize(d, T, "simpson")[1] for d in decs]
        a.plot(decs, vals, "o-", label=f"T = {T:g} h")
    a.set(title="Best achievable total quality vs declination", xlabel="declination (deg)", ylabel="J at optimum")
    a.legend(); a.grid(alpha=.3); fig.tight_layout(); fig.savefig(f"{OUT}/fig4_quality_vs_declination.png", dpi=160); plt.close(fig)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Figure 2 -- Bound-mode spectrum (v1.1 beautification pass).

Source data (unchanged from v1.0, per figures_manifest.json):
  data/spectral/spectral_validation.csv  (61 rows: profile,W0,a,sqrtW0_a,n,
                                           omega2_n,omega2_matching,kappa_n)
  data/spectral/mode_counts.csv          (19 rows: profile,W0,a,sqrtW0_a,
                                           N_numeric,N_analytic)
All 61 + 19 rows are read verbatim from the archived CSVs below; nothing
is simulated, re-fit, or extrapolated beyond what those rows already
contain.

Two numbers ARE derived here, both direct arithmetic on archived values,
not new scientific calculations:
  - relative eigenvalue deviation = |omega2_n - omega2_matching| / omega2_matching,
    i.e. exactly the quantity the main text already reports in prose
    ("agree to better than 2e-6 relative").
  - the analytic curves in the middle panel are PIECEWISE-LINEAR through the
    archived omega2_matching points (sorted by sqrtW0_a, per mode branch) --
    deliberately NOT a fresh root-find of the transcendental matching
    condition, to guarantee the plotted analytic curve passes through
    exactly the archived, already-verified values with no new calculation
    introduced. At the sweep's 0.05-in-W0 spacing (9-19 points per branch
    depending on mode) this is visually smooth for modes 1-4; mode 5 has
    only 2 archived points and is shown as the short segment the data
    supports, not extrapolated further.

v1.1 changes (readability only, no scientific content change):
  - Mode legend (modes 1-5) moved outside the plotting area (single shared
    legend above all three panels) instead of sitting inside panel B.
  - Figure promoted from single-column (\\columnwidth) to a two-column
    spanning figure* in main.tex, giving each of the three panels roughly
    double its previous width -- this is the fix for "three-panel layout
    is visually dense" (documented in CHANGELOG.md / FIGURE_QA_REPORT.md).
  - Colorblind-safe, grayscale-distinguishable Okabe-Ito palette with a
    distinct marker shape per mode (belt-and-suspenders for print).
  - The near-threshold mode-count discrepancy at W0=0.10 is kept fully
    visible (explicit marker + annotation in panel A); not concealed.
"""
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import common_style as cs


@cs.isolated_style
def make_figure(data_dir, out_dir, synthetic=False):
    mpl.rcParams.update({
        "font.size": 10,
        "axes.labelsize": 11,
        "legend.fontsize": 8,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "font.family": "serif",
        "mathtext.fontset": "cm",
    })

    DATA = f"{data_dir}/spectral"

    # ---- load archived data verbatim -------------------------------------------
    rows = []
    with open(f"{DATA}/spectral_validation.csv") as f:
        for r in csv.DictReader(f):
            rows.append({
                "W0": float(r["W0"]), "sqrtW0_a": float(r["sqrtW0_a"]),
                "n": int(r["n"]), "omega2_n": float(r["omega2_n"]),
                "omega2_matching": float(r["omega2_matching"]),
            })

    counts = []
    with open(f"{DATA}/mode_counts.csv") as f:
        for r in csv.DictReader(f):
            counts.append({
                "W0": float(r["W0"]), "sqrtW0_a": float(r["sqrtW0_a"]),
                "N_numeric": int(r["N_numeric"]), "N_analytic": int(r["N_analytic"]),
            })

    if len(rows) != 61 or len(counts) != 19:
        raise ValueError("Expected original spectral data: 61 eigenvalue rows and 19 mode counts")

    # Okabe-Ito colorblind-safe palette + distinct marker per mode
    MODE_STYLE = {
        1: dict(color="#0072B2", marker="o"),
        2: dict(color="#D55E00", marker="s"),
        3: dict(color="#009E73", marker="^"),
        4: dict(color="#E69F00", marker="D"),
        5: dict(color="#CC79A7", marker="*"),
    }

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(7.0, 2.75))

    # ---- Panel A: bound-mode count staircase -----------------------------------
    x_analytic = np.linspace(0.3, 15.0, 2000)
    thresholds = [(2 * n - 1) * np.pi / 2 for n in range(1, 7)]
    N_analytic_curve = np.searchsorted(thresholds, x_analytic, side="right")
    axA.step(x_analytic, N_analytic_curve, where="post", color="0.45", lw=1.3,
              label="analytic (Eq.~A4)", zorder=2)
    xs = [c["sqrtW0_a"] for c in counts]
    ys = [c["N_numeric"] for c in counts]
    axA.plot(xs, ys, "o", ms=4.2, color="black", label="numerical", zorder=3)
    # explicitly flag the one mismatch rather than concealing it
    mism = [c for c in counts if c["N_numeric"] != c["N_analytic"]]
    for c in mism:
        axA.annotate("mismatch\n" + r"($W_0=0.10$)",
                     xy=(c["sqrtW0_a"], c["N_numeric"]),
                     xytext=(c["sqrtW0_a"] + 1.3, c["N_numeric"] - 0.9),
                     fontsize=7, color="#D55E00", ha="left",
                     arrowprops=dict(arrowstyle="->", color="#D55E00", lw=0.9))
    axA.set_xlabel(r"$\sqrt{W_0}\,a$")
    axA.set_ylabel("bound-mode count $N$")
    axA.set_xlim(0, 15)
    axA.set_ylim(-0.3, 5.8)
    axA.legend(loc="upper left", frameon=False, fontsize=7.5, handlelength=1.4)

    # ---- Panel B: eigenvalue trajectories ---------------------------------------
    for n in range(1, 6):
        pts = sorted([r for r in rows if r["n"] == n], key=lambda r: r["sqrtW0_a"])
        if not pts:
            continue
        x = [p["sqrtW0_a"] for p in pts]
        y_num = [p["omega2_n"] for p in pts]
        y_ana = [p["omega2_matching"] for p in pts]
        st = MODE_STYLE[n]
        axB.plot(x, y_ana, "-", color=st["color"], lw=1.1, alpha=0.55, zorder=1)
        axB.plot(x, y_num, st["marker"], color=st["color"], ms=4.0,
                  mfc=st["color"], mec="none", zorder=3)
    axB.set_xlabel(r"$\sqrt{W_0}\,a$")
    axB.set_ylabel(r"$\omega_n^2/m_0^2$")
    axB.set_xlim(0, 15)
    axB.set_ylim(-0.02, 1.05)

    # ---- Panel C: relative eigenvalue deviation ---------------------------------
    for n in range(1, 6):
        pts = sorted([r for r in rows if r["n"] == n], key=lambda r: r["sqrtW0_a"])
        if not pts:
            continue
        x = [p["sqrtW0_a"] for p in pts]
        dev = [abs(p["omega2_n"] - p["omega2_matching"]) / p["omega2_matching"] for p in pts]
        st = MODE_STYLE[n]
        axC.plot(x, dev, st["marker"], color=st["color"], ms=4.0,
                  mfc=st["color"], mec="none", zorder=3)
    axC.set_yscale("log")
    axC.set_xlabel(r"$\sqrt{W_0}\,a$")
    axC.set_ylabel("rel. eigenvalue deviation")
    axC.set_xlim(0, 15)

    # single shared mode legend, fully outside all three plotting areas
    handles = [plt.Line2D([0], [0], marker=MODE_STYLE[n]["marker"], color=MODE_STYLE[n]["color"],
                            mfc=MODE_STYLE[n]["color"], mec="none", lw=0, ms=6,
                            label=f"mode {n}") for n in range(1, 6)]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.13),
               ncol=5, frameon=True, fontsize=8.5, handletextpad=0.4, columnspacing=1.2)

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    return cs.save_figure(fig, "figure2_spectrum", out_dir,
                          source_data=[f"{DATA}/spectral_validation.csv", f"{DATA}/mode_counts.csv"],
                          run_ids=[], synthetic=synthetic)


if __name__ == "__main__":
    cs.standalone(make_figure)

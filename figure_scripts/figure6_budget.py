#!/usr/bin/env python3
"""
Figure 6 -- Energy budget and closure diagnostics (v1.1 beautification pass).

Source data (unchanged from v1.0, per figures_manifest.json):
  data/baseline/time_series.csv   (800 data rows: t,E_box,E_esc,E_mode_1,E_mode_2,
                                    f_cap,f_esc,f_res,closure,cons_err)
  data/baseline/run_metadata.json (status=REJECTED, t_clear=178.7234...)
All curves are the archived columns, read verbatim; nothing is
recomputed except the sign-stripped |closure| = |f_esc+f_cap-1| already
used as the bottom-panel diagnostic in v1.0 (direct absolute value of an
archived column, not a new calculation).

v1.1 changes (readability only, no scientific content change):
  - "run status: REJECTED", which sat directly on the top plot border in
    v1.0, is removed from the figure body; the status is carried in the
    caption instead (one of the two options the brief allows).
  - "radiation cleared" is now an arrow-annotated label pointing at the
    actual t_clear crossing, placed in empty space above the curves,
    instead of sitting on top of the diagnostic curve.
  - Both legends (top panel: f_esc/f_res/f_cap/sum; bottom panel:
    |f_esc+f_cap-1|/f_res/epsilon_budget) are pulled out of the data
    region: the top legend sits above the top panel; the bottom legend
    sits in the one genuinely empty region of the log-scale panel
    (late time, high value -- nothing is ever there because both
    diagnostics have already dropped by then).
  - Does NOT imply the internal closure diagnostic independently
    resolves the T6 flux-reintegration limitation (unchanged from v1.0;
    reinforced in the updated caption).
"""
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import common_style as cs
import common_io as cio


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

    DATA = f"{data_dir}/baseline"

    t, f_cap, f_esc, f_res, closure = [], [], [], [], []
    with open(f"{DATA}/time_series.csv") as f:
        for r in csv.DictReader(f):
            t.append(float(r["t"])); f_cap.append(float(r["f_cap"]))
            f_esc.append(float(r["f_esc"])); f_res.append(float(r["f_res"]))
            closure.append(float(r["closure"]))
    t = np.array(t); f_cap = np.array(f_cap); f_esc = np.array(f_esc)
    f_res = np.array(f_res); closure = np.array(closure)
    meta = cio.load_run_metadata(f"{DATA}/run_metadata.json")
    t_clear = float(meta["outcomes"]["t_clear"])
    eps_budget = float(meta["gates"]["eps_budget"])
    if len(t) != 800:
        raise ValueError("Expected original baseline: 800 data rows")

    fig, (axT, axB) = plt.subplots(2, 1, figsize=(3.4, 4.6),
                                     gridspec_kw={"hspace": 0.32})

    # ---- top panel: stacked energy-fraction budget -----------------------------
    axT.fill_between(t, 0, f_esc, facecolor="#D55E00", alpha=0.35, hatch="//",
                       edgecolor="#D55E00", lw=0, label=r"$f_{\rm esc}$", zorder=2)
    axT.fill_between(t, f_esc, f_esc + f_res, facecolor="#CC79A7", lw=0,
                       label=r"$f_{\rm res}$", zorder=2)
    axT.fill_between(t, f_esc + f_res, f_esc + f_res + f_cap, facecolor="#0072B2",
                       lw=0, label=r"$f_{\rm cap}$", zorder=2)
    axT.plot(t, f_esc + f_res + f_cap, "-", color="black", lw=1.0,
              label=r"$f_{\rm esc}+f_{\rm res}+f_{\rm cap}$", zorder=3)
    axT.set_xlim(0, 400)
    axT.set_ylim(0, 1.03)
    axT.set_ylabel(r"fraction of $E_{\rm in}$")
    axT.set_xlabel(r"$t$ [model units]")
    axT.legend(loc="upper center", bbox_to_anchor=(0.5, 1.30), ncol=2,
                frameon=True, fontsize=7.5, handlelength=1.3, columnspacing=1.0,
                handletextpad=0.4, borderaxespad=0.3)

    # ---- bottom panel: closure diagnostics (log scale) -------------------------
    abs_closure = np.abs(closure)
    (hB1,) = axB.plot(t, abs_closure, "-", color="black", lw=1.1,
                        label=r"$|f_{\rm esc}+f_{\rm cap}-1|$", zorder=3)
    (hB2,) = axB.plot(t, f_res, "-.", color="#CC79A7", lw=1.1,
                        label=r"$f_{\rm res}$", zorder=3)
    hB3 = axB.axhline(eps_budget, color="#E69F00", lw=1.1, ls="--",
                        label=r"$\varepsilon_{\rm budget}$", zorder=2)
    axB.set_yscale("log")
    axB.set_xlim(0, 400)
    axB.set_ylim(1e-4, 2.0)
    axB.set_xlabel(r"$t$ [model units]")
    axB.set_ylabel("closure diagnostics")

    # radiation-cleared: arrow to the actual crossing, label in empty space above
    axB.annotate("radiation\ncleared", xy=(t_clear, eps_budget),
                 xytext=(t_clear + 70, 0.5), fontsize=7.5, color="0.2", ha="left",
                 arrowprops=dict(arrowstyle="->", color="0.3", lw=0.9))

    # bottom-panel legend placed in the one genuinely empty region: late time,
    # high value -- both diagnostics have already dropped well below there
    axB.legend(handles=[hB1, hB2, hB3], loc="lower left", frameon=True,
                fontsize=7, handlelength=1.6, borderaxespad=0.4)

    fig.align_ylabels([axT, axB])
    return cs.save_figure(fig, "figure6_budget", out_dir,
                          source_data=[f"{DATA}/time_series.csv", f"{DATA}/run_metadata.json"],
                          run_ids=[meta["run_id"]], synthetic=synthetic)


if __name__ == "__main__":
    cs.standalone(make_figure)

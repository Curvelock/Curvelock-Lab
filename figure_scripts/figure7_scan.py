#!/usr/bin/env python3
"""
Figure 7 -- Robustness scan + per-run status matrix (v1.1 beautification
pass, including the Figure-7 addendum redesign).

Source data (unchanged from v1.0, per figures_manifest.json):
  data/scan_summary.csv  (7 rows: baseline + the 6-row original robustness
                           subset -- profile x depth, a=5 throughout)
This file contains ONLY the original six-configuration robustness subset
(plus one duplicate "baseline" row -- see note below). It does NOT contain
control_k0_zero or the four survey_k* runs: those are excluded from this
figure by construction of the source data, not by a plotting choice, and
their explanations live in the main text (Sec. 6.6) and Supplemental
Table S6, per the addendum's instruction not to imply they belong to this
six-run matrix.

Baseline-cell investigation (addendum Sec. 1 & 4): scan_summary.csv
contains two rows describing the SAME physical configuration (profile=
square, W0=0.9, a=5): one tagged run_id="baseline" (in_robustness_set=
False) and one tagged run_id="rob_square_W0p9" (in_robustness_set=True).
Both rows agree exactly on net_f_cap (0.0023572250382...) and on every
gate/status field (status=REJECTED, seed_clean=True, numerics_valid=True,
radiation_cleared=True, budget_closed=True). There is therefore NO
ambiguity or inconsistency in the underlying status classification for
this cell -- only a presentational one in the v1.0 figure (an unexplained
lighter border on this one cell). Per the addendum's own decision rule,
this is treated as a plotting-style inconsistency to correct, not as
evidence of a different scientific classification: the cell is now drawn
with the SAME styling as every other rejected cell, with an explicit,
subtle "(baseline)" text label added -- rather than retaining an
unexplained border difference that could be mistaken for a different
status.

v1.1 / addendum changes (readability only, no scientific content change):
  - Figure promoted from single-column (\\columnwidth) to a two-column
    spanning figure* in main.tex (same fix as Figure 2, for the same
    "visually dense" / legend-crowding reason).
  - Left panel: legends (profile shapes, status) moved out of the
    plotting area; a small horizontal jitter (+/-0.003 in W0) separates
    the three profiles at each depth for visual clarity -- the true W0
    values are exactly 0.85 / 0.90 and are unaffected; this is disclosed
    here and in FIGURE_QA_REPORT.md.
  - Left panel: the H2 threshold (f_min=0.02=2%) is stated in a text
    annotation placed in the empty lower part of the axis, since 0.02 is
    roughly 8-10x above the top of the zoomed axis range and cannot be
    drawn as an on-scale line without either compressing all the real
    data to a sliver or drawing a line at a misleading position.
  - Right panel: completely redesigned per the addendum. Each of the six
    cells now shows an explicit "REJECTED" label and its net-capture
    percentage; a salmon/red fill (never green, per the addendum) is
    used consistently for all six (all REJECTED) cells; the previously
    redundant five-category hatch-pattern legend is removed from this
    panel (every cell here is the same status), and the square/W0=0.90
    cell is explicitly labeled "(baseline)" rather than carrying an
    unexplained different border.
"""
import csv
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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

    DATA = f"{data_dir}"

    rows = []
    with open(f"{DATA}/scan_summary.csv") as f:
        for r in csv.DictReader(f):
            rows.append(r)
    print(f"rows loaded from scan_summary.csv: {len(rows)} (expect 7: baseline + 6 robustness)")

    # the six-configuration robustness subset proper (drop the duplicate
    # "baseline" row for the scatter -- it is the same point as rob_square_W0p9;
    # plotting both would silently double-count one configuration)
    robustness = [r for r in rows if r["in_robustness_set"] == "True"]
    if len(rows) != 7 or len(robustness) != 6:
        raise ValueError("Expected baseline plus six original robustness rows")
    for r in robustness:
        if r["status"] != "REJECTED":
            raise ValueError("Original robustness status mismatch")

    PROFILE_STYLE = {
        "square": dict(color="#0072B2", marker="s", jitter=-0.003),
        "tanh":   dict(color="#D55E00", marker="o", jitter=0.0),
        "sech2":  dict(color="#009E73", marker="^", jitter=0.003),
    }
    PROFILE_LABEL = {"square": "square", "tanh": "tanh", "sech2": r"sech$^2$"}

    fig = plt.figure(figsize=(7.0, 3.3))
    axL = fig.add_axes([0.085, 0.16, 0.40, 0.72])
    axR = fig.add_axes([0.58, 0.08, 0.40, 0.84])

    # ---- left panel: robustness scatter -----------------------------------
    for r in robustness:
        st = PROFILE_STYLE[r["profile"]]
        x = float(r["W0"]) + st["jitter"]
        y = float(r["net_f_cap"])
        axL.plot(x, y, st["marker"], mfc="none", mec=st["color"], mew=1.6,
                  ms=8, zorder=3)

    handles = [plt.Line2D([0], [0], marker=PROFILE_STYLE[p]["marker"], color=PROFILE_STYLE[p]["color"],
                            mfc="none", mec=PROFILE_STYLE[p]["color"], mew=1.6, lw=0, ms=8,
                            label=PROFILE_LABEL[p]) for p in ["square", "tanh", "sech2"]]
    axL.legend(handles=handles, title="profile (all REJECTED)", loc="upper center",
                bbox_to_anchor=(0.5, 1.30), ncol=3, frameon=True, fontsize=7.5,
                title_fontsize=7.5, handletextpad=0.3, columnspacing=0.9, borderaxespad=0.3)

    axL.set_xlim(0.84, 0.91)
    axL.set_ylim(0.0, 0.0028)
    axL.set_xticks([0.85, 0.90])
    axL.set_xlabel(r"$W_0$  [model units$^2$]")
    axL.set_ylabel(r"net $f_{\rm cap}$")
    axL.text(0.845, 0.0004,
          "$H2$ threshold $f_{\\min}=0.02$ (2%)\nis $\\sim$8--10$\\times$ the observed\n"
          "capture and is off-scale.",
              fontsize=7, color="0.25", va="bottom", ha="left")

    # ---- right panel: redesigned per-run status matrix --------------------
    row_order = ["sech2", "tanh", "square"]  # top to bottom
    col_order = ["0.85", "0.9"]              # left to right
    cell = {(r["profile"], r["W0"]): r for r in robustness}

    axR.set_xlim(0, 2)
    axR.set_ylim(0, 3)
    axR.set_xticks([0.5, 1.5])
    axR.set_xticklabels([r"$W_0=0.85$", r"$W_0=0.90$"])
    axR.set_yticks([2.5, 1.5, 0.5])
    axR.set_yticklabels([PROFILE_LABEL[p] for p in row_order])
    axR.tick_params(length=0)
    for spine in axR.spines.values():
        spine.set_visible(False)
    axR.set_title("Per-run status (original robustness subset)", fontsize=8.5, pad=10)

    REJECTED_FILL = "#F5CFC9"     # salmon -- never green, per the addendum
    REJECTED_EDGE = "#8B1A1A"     # dark red border

    for i, prof in enumerate(row_order):
        for j, w0s in enumerate(col_order):
            r = cell[(prof, w0s)]
            pct = float(r["net_f_cap"]) * 100
            x0, y0 = j, 2 - i
            axR.add_patch(mpatches.Rectangle((x0 + 0.05, y0 + 0.05), 0.9, 0.9,
                                               facecolor=REJECTED_FILL, edgecolor=REJECTED_EDGE,
                                               lw=1.4, zorder=2))
            label = f"REJECTED\n{pct:.3f}%"
            if prof == "square" and w0s == "0.9":
                label += "\n(baseline)"
            axR.text(x0 + 0.5, y0 + 0.5, label, ha="center", va="center",
                      fontsize=7.3, color="#5A0F0F", zorder=3, linespacing=1.4)

    return cs.save_figure(fig, "figure7_scan", out_dir,
                          source_data=[f"{DATA}/scan_summary.csv"],
                          run_ids=[r["run_id"] for r in robustness], synthetic=synthetic)


if __name__ == "__main__":
    cs.standalone(make_figure)

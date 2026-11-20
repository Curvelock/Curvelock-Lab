#!/usr/bin/env python3
"""
Figure 1 -- Prescribed potential well schematic (v1.1 beautification pass).

Source data: NONE (per figures_manifest.json, figure1_schematic has
source_data=[]). This figure is, and was in the original v1.0 package, a
purely analytic schematic of the three prescribed well profiles defined
in main.tex Eqs. (squarewell), (tanhwell), (sech2well), and the effective
mass relation Eq. (meff) = m0^2 + W(r). No simulation output is plotted;
nothing here is "invented data" -- these are the closed-form profile
definitions already published in the manuscript, evaluated at the
baseline parameters (W0=0.9, a=5, Delta=0.5, m0=1, same as every other
baseline figure in this package).

v1.1 changes (readability only, no scientific content change):
  - The "-W0 (positive depth W0)" floating annotation, which overlapped
    the legend and curves in v1.0, is removed; the depth value is now
    carried in the legend text instead (no information lost).
  - "exterior threshold" and "bound-mode window" annotations are moved
    into the genuinely empty upper-left quadrant of the bottom panel
    (r small, m_eff^2/m0^2 large), where no curve, legend, or border
    is nearby.
  - The legend is pulled out of the plotting area entirely, as a single
    shared horizontal legend above the top panel.
  - Both panels, both profiles, and the W(r)-vs-m_eff^2(r) distinction
    are all preserved unchanged.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import common_style as cs


@cs.isolated_style
def make_figure(data_dir, out_dir, synthetic=False):
    mpl.rcParams.update({
        "font.size": 10,
        "axes.labelsize": 11,
        "legend.fontsize": 9,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "font.family": "serif",
        "mathtext.fontset": "cm",
    })

    # ---- baseline parameters (same as every other baseline-run figure) --------
    W0 = 0.9
    a = 5.0
    Delta = 0.5
    m0_sq = 1.0

    r = np.linspace(0.0, 15.0, 2000)

    def W_square(r):
        return np.where(r < a, -W0, 0.0)

    def W_tanh(r):
        return -(W0 / 2.0) * (1.0 - np.tanh((r - a) / Delta))

    def W_sech2(r):
        return -W0 / np.cosh(r / a) ** 2

    profiles = [
        ("square", W_square(r), "-", 1.8),
        ("tanh", W_tanh(r), "--", 1.8),
        (r"sech$^2$", W_sech2(r), "-.", 1.8),
    ]
    colors = {"square": "black", "tanh": "#1f77b4", r"sech$^2$": "#d2691e"}

    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(3.4, 4.6), sharex=True,
                                          gridspec_kw={"hspace": 0.12})

    # interior shading (same gray band both panels)
    for ax in (ax_top, ax_bot):
        ax.axvspan(0, a, color="0.88", zorder=0)
        ax.axvline(a, color="0.4", lw=0.8, ls=":", zorder=1)

    handles = []
    for name, y, ls, lw in profiles:
        c = colors[name]
        (h,) = ax_top.plot(r, y, ls, color=c, lw=lw, label=name, zorder=3)
        handles.append(h)
        ax_bot.plot(r, m0_sq + y, ls, color=c, lw=lw, zorder=3)

    ax_top.set_ylabel(r"$W(r)$  [model units$^2$]")
    ax_top.set_ylim(-1.05, 0.12)
    ax_top.tick_params(labelbottom=False)

    ax_bot.set_ylabel(r"$m_{\rm eff}^2(r)/m_0^2$")
    ax_bot.set_xlabel(r"$r$  [model units]")
    ax_bot.set_ylim(-0.03, 1.12)
    ax_bot.set_xlim(0, 15)
    ax_bot.axhline(1.0, color="0.3", lw=0.8, ls=":", zorder=1)
    ax_bot.axhline(0.0, color="0.3", lw=0.6, ls=":", zorder=1)

    # -- annotations relocated into genuinely empty whitespace (upper-left of
    #    the bottom panel: small r, large m_eff^2/m0^2 -- no curve passes here) --
    ax_bot.text(0.4, 0.97, "exterior threshold\n$m_0^2$",
                ha="left", va="top", fontsize=8, color="0.25")
    ax_bot.text(7.0, 0.45, "bound-mode window\n" r"$0<\omega^2<m_0^2$",
                ha="left", va="top", fontsize=8, color="0.25")

    # shared legend, fully outside both axes, with the common (W0, a) values
    # carried in the legend title so no separate floating depth label is needed
    fig.legend(handles, [h.get_label() for h in handles],
               title=r"Profiles ($W_0=0.9$, $a=5$):",
               loc="upper center", bbox_to_anchor=(0.54, 1.045),
               ncol=3, frameon=True, fontsize=8, title_fontsize=8,
               columnspacing=1.0, handletextpad=0.5, borderaxespad=0.3)

    fig.align_ylabels([ax_top, ax_bot])
    return cs.save_figure(fig, "figure1_schematic", out_dir,
                          source_data=[],
                          run_ids=[], synthetic=synthetic)


if __name__ == "__main__":
    cs.standalone(make_figure)

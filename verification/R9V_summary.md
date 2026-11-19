# R9-V verification summary — Paper 1

**Frozen solver:** `curvature_well_capture_v0.3.py` (v0.3.0), SHA-256 `036dde53b227729fd20eb41e6ab82b6020828bf26a184e3f012cc28e45282491`.

This is a compact transcription of existing archived results and the November 17, 2026 C1/C2 closure review. No new simulation, data fitting, gate change, or independent full-field reintegration was performed to prepare this repository record.

## Scientific verdict

| Question | Disposition |
| --- | --- |
| H1: existence of localized bound modes in a prescribed well | Supported |
| Dynamic capture: nonlinear population of localized modes by the incoming packet | Supported at finite time in the tested regime |
| H2: operational capture at or above the original 2% threshold | Rejected in the tested regime |
| C1: fundamental-mode dominance of captured energy | Accepted with disclosure |
| C2: original/native basis comparison | Accepted with disclosure |

The clean fundamental eigenvalue is approximately 0.364043 (omega squared); the shallow second mode, approximately 0.999236, is box-affected. All six original robustness members were rejected by the 2% H2 gate, despite satisfying numerical, clearance and budget criteria. The original six-member set spans 0.1985–0.2491% net capture. The solver's archived gates are eps_lin = 0.002, f_min = 0.02, eps_res = 0.01, eps_budget = 0.01, and rtol = 0.3 (fractions). A solver status of REJECTED is the operational H2 outcome, not a denial of H1 or weak finite-time modal occupation.

## Verification tests

- **R0 — cross-platform baseline:** original Linux net capture = 0.002357225038249329; Windows rerun = 0.0023344811261043884. The spectral eigenvalues are identical in these two archived outputs; the capture differs by about 0.97% relative to the Linux value.
- **T1 — seed identity:** passed, with both array equality and identical SHA-256 seed hashes for nonlinear and linear-null initial conditions. Seed hash: `5ed5b04b9e2ba6cfffe4ffd3a18c802c759279119820462f28c3442f60384175`.
- **T2 — independent basis comparison:** nonlinear maximum absolute difference in the full trajectory = 0.10054564827258936 (10.0546 percentage points); final-window maximum absolute difference = 0.00022152314530342419 (0.02215 percentage points). The latter is about 8.7% when normalized by the original-basis endpoint 0.002539149086897734, not a maximum pointwise relative error. The native-basis endpoint is 0.0024072741169051. The diagnostic used 600 snapshots rather than the baseline's 81; the native estimator reconstructed velocities with central finite differences. At the diagnostic cadence, the expected amplitude loss is approximately 73% for the carrier and 2.7% for the fundamental. This is quantitatively consistent with the transient discrepancy, but does not prove its root cause. The primary reported occupation uses the original basis.
- **T3 — grid refinement:** net capture at dr = 0.05, 0.025 and 0.0125 is 0.2334%, 0.2765% and 0.2881%. The finite-time occupation survives refinement; this is not proof of a fully converged continuum capture fraction.
- **T4 — domain/absorber sensitivity:** r_max = 150 gives 0.2334%; sponge start fraction = 0.85 gives 0.2339%. The baseline r_max = 100, sponge start = 0.75, strength = 3 is reused in the archive as the reference for each corresponding sweep. Strength = 3 is a reference, not evidence of an independent strength variation.
- **T5 — nonlinear-coupling sweep:** at lambda_chi = 0 the net fraction is 2.2695707389754556e-16. At lambda_chi = 0.05, 0.1, 0.25, 0.5 and 1 the respective capture percentages are 0.0818%, 0.1539%, 0.4563%, 0.5285% and 0.2334%. The dependence is nonmonotonic; no point passes 2%. At 0.05 and 0.1 the closure residual is comparable to or larger than the captured fraction. The mode-share diagnostic is withheld by the solver below its 0.2% capture guard.
- **T6 — independent energy-flux audit:** NOT_IDENTIFIABLE on R0, the finest T3 grid, T4 r_max = 150, and T5 lambda_chi = 0. R0 snapshots are decimated and float32; the other audited runs lack full field snapshots. This is a missing-evidence limitation, not a detected conservation failure. The structural checks report relative self-adjointness error at most 4.85e-16, modal orthogonality error at most 1.3e-17, and origin-stencil error decreasing from 3.70e-8 to 2.32e-9 under refinement. These checks do not substitute for independent flux reintegration.

## C1: mode-resolved capture

The nominal final analysis window is t = 375–400 (last five characteristic times); radiation clearance precedes this window. Exact sampled window endpoints are in the CSV.

| Configuration | dr | r_max | Sponge start | Net capture / incident energy | Mode 1 / captured energy | Mode 2 / captured energy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| baseline | 0.05 | 100 | 0.75 | 0.2334% | 97.6529% | 2.5443% |
| t3_dr_0.0125 | 0.0125 | 100 | 0.75 | 0.2881% | 96.2958% | 3.6867% |
| t3_dr_0.025 | 0.025 | 100 | 0.75 | 0.2765% | 96.7747% | 3.0354% |
| t4_rmax_150 | 0.05 | 150 | 0.75 | 0.2334% | 97.7054% | 2.4041% |
| t4_spongefrac_0.85 | 0.05 | 100 | 0.85 | 0.2339% | 97.6176% | 2.4058% |
| Original Linux baseline | 0.05 | 100 | 0.75 | 0.2357% | 97.5242% | 2.5138% |

The fundamental accounts for 96.30–97.71% of the captured fraction across the six independent C1 configurations, exceeding the >90% criterion. The closely spaced percentages in this table represent different denominators: net capture is relative to incident energy, whereas mode shares are relative to captured energy. Mode shares use independently computed medians and need not sum to exactly 100%. The second mode's share increases modestly with grid refinement; it is not established to vanish in the continuum limit.

## Data provenance

`R9V_results.csv` contains 16 archival metadata records (15 R9-V records plus the original Linux baseline). Five R9-V reference rows reuse R0; they are not independent runs. `c1_included` identifies the six independent C1 configurations, not all runs with a nonblank mode share.

Each numeric field was taken from archived `run_metadata.json`. Fraction-valued columns are not percentages: multiply by 100 to express percent. Blank entries are missing or solver-withheld values, never imputed zero. Archive paths use normalized forward slashes for readability. Source archives remain separate from this minimal repository:

| Archive | SHA-256 of supplied archive |
| --- | --- |
| data_r9v.zip | `e8c8facf601599d20d29c4ae0af9d60b599d2462e66ae1cf87df29ee35ee1373` |
| paper1_v0.3_data_and_figures(3).zip | `eb3b51ae920918e77916c1f4944278944203ce9f54fe9e5c9e19fd70dd16d4f4` |

Other sources: `T1/t1_t2_report.json`, `T6/t6_audit/` diagnostic tables and JSON, original `data/scan_summary.csv`, and the R9-V C1/C2 closure review dated November 17, 2026. The data tables are not a replacement for the raw archives.

## Scope and limitations

- Single deterministic realization per configuration; no statistical confidence intervals or claim of stationarity. Residual final-window mode variation is recorded in the CSV.
- The basis discrepancy and non-identifiability of independent flux closure are explicitly disclosed.
- H2 rejection is specific to the tested prescribed-well models, parameter ranges, and packet family. It is not a universal no-capture theorem.
- Finite-time modal occupation up to t = 400 does not establish indefinite persistence, asymptotic stability, or radiative lifetime.
- The study uses a prescribed potential/well in model units, not a self-consistent dynamical spacetime geometry or a particle-physics model.

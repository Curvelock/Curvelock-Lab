# Original publication figure inputs

These files are byte-for-byte copies of the corresponding files in the original July 23, 2026 Paper 1 data archive (provided as `data.zip`). They are a deliberately small subset of the original archive, not R9-V output and not the complete raw simulation record. The frozen solver is v0.3.0. The archived baseline metadata records the original Linux run, its configuration and seed hashes, numerical environment, parameters and adjudication.

| Figure | Source files (relative to this directory) |
| --- | --- |
| 1 | None (analytic profile equations) |
| 2 | `spectral/spectral_validation.csv`, `spectral/mode_counts.csv` |
| 3 | `baseline/mode_profiles.csv`, `baseline/run_metadata.json` |
| 4 | `baseline/field_evolution.npz`, `baseline/run_metadata.json` |
| 5 | `baseline/time_series.csv`, `null_linear/time_series.csv`, `baseline/run_metadata.json` |
| 6 | `baseline/time_series.csv`, `baseline/run_metadata.json` |
| 7 | `scan_summary.csv` |

`null_linear/run_metadata.json` is included for matched-control provenance. The field snapshot archive contains 80 saved time snapshots at 501 radial positions in float32, not full-resolution fields. The baseline and null time series each contain 800 data rows. The scan summary has six original robustness rows plus a duplicated baseline row; the latter is not counted as an extra robustness configuration. The original baseline net capture fraction is 0.002357225038249329; the R9-V verification results are separately maintained in `verification/`.

Run `python figure_scripts/generate_all_figures.py` from the repository root to regenerate the figures and their source SHA-256 manifest. No numerical run is repeated by that command.

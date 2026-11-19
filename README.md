# Curvelock Lab — Paper 1

**Localized Standing Waves in Prescribed Curvature Wells**

This repository provides the frozen numerical solver and a compact verification record for Paper 1 of the Curvelock Research Program (Danai Leangtragul, independent researcher).

The study asks whether a prescribed inward-curvature well can support localized standing-wave modes and whether an initially orthogonalized incoming wave packet can populate those modes through nonlinear evolution. The prescribed well is a model input, not a self-consistent solution of gravitational field equations. The results do not establish matter formation, particle masses, or a mechanism of quantum chromodynamics.

## Repository contents

| File | Purpose |
| --- | --- |
| [curvature_well_capture_v0.3.py](curvature_well_capture_v0.3.py) | Frozen standalone numerical solver (v0.3.0) |
| [requirements.txt](requirements.txt) | Python dependencies |
| [LICENSE](LICENSE) | MIT license |
| [verification/R9V_summary.md](verification/R9V_summary.md) | Verification findings, acceptance gates and limitations |
| [verification/R9V_results.csv](verification/R9V_results.csv) | Extracted run-level numerical results, including reused reference rows |

The manuscript and full raw simulation archive are maintained separately. This repository is the minimal code and verification reference, not a complete archive of every analysis script or figure-generation workflow.

## Installation

```bash
git clone https://github.com/Curvelock/Curvelock-Lab.git
cd Curvelock-Lab
python -m venv .venv
```

Activate the virtual environment (Windows PowerShell: `.venv\Scripts\Activate.ps1`; Linux/macOS: `source .venv/bin/activate`), then:

```bash
python -m pip install -r requirements.txt
python curvature_well_capture_v0.3.py --help
```

The archived environments were Linux (Python 3.12.3, NumPy 2.4.4, SciPy 1.17.1) and Windows (Python 3.14.6, NumPy 2.5.0, SciPy 1.18.0). These are recorded versions, not minimum supported versions. Dependencies are unpinned to permit compatible installations; use the archived versions when matching the numerical environment closely. Last-digit differences across platforms are possible.

## Verify solver integrity

Frozen solver SHA-256:

```text
036dde53b227729fd20eb41e6ab82b6020828bf26a184e3f012cc28e45282491
```

Windows PowerShell:

```powershell
Get-FileHash .\curvature_well_capture_v0.3.py -Algorithm SHA256
```

Linux:

```bash
sha256sum curvature_well_capture_v0.3.py
```

The original frozen solver is also available at [this immutable Git commit](https://github.com/Curvelock/Curvelock-Lab/blob/f5063d0788805d564e1b98f00265a8ce6d886eb9/curvature_well_capture_v0.3.py).

## Reproduce the baseline and linear null

**Important:** The solver's built-in default configurations are placeholders, not the publication configuration. Save the JSON below as `baseline.json` in the repository root. It reconstructs the frozen baseline parameters from the archived T1/T2 report and includes the original six robustness configurations and selected R9-V diagnostic overrides. It is not a byte-identical copy of the original configuration file.

```json
{
  "baseline": {
    "profile": "square",
    "W0": 0.9,
    "a": 5,
    "Delta": 0.5,
    "lambda_chi": 1,
    "amp": 0.2,
    "r0_packet": 55,
    "packet_width": 2,
    "packet_kind": "chi_carrier",
    "k0": 3.5,
    "orthogonalize": true,
    "r_max": 100,
    "dr": 0.05,
    "cfl": 0.2,
    "t_end": 400,
    "n_out": 600,
    "n_snapshots": 81,
    "sponge_start_frac": 0.75,
    "sponge_strength": 3,
    "field_r_stride": 4,
    "field_n_snapshots": 80,
    "run_id": "baseline",
    "in_robustness_set": false
  },
  "scan": {
    "runs": [
      {
        "run_id": "rob_square_W0p85",
        "profile": "square",
        "W0": 0.85,
        "in_robustness_set": true
      },
      {
        "run_id": "rob_square_W0p9",
        "profile": "square",
        "W0": 0.9,
        "in_robustness_set": true
      },
      {
        "run_id": "rob_tanh_W0p85",
        "profile": "tanh",
        "W0": 0.85,
        "in_robustness_set": true
      },
      {
        "run_id": "rob_tanh_W0p9",
        "profile": "tanh",
        "W0": 0.9,
        "in_robustness_set": true
      },
      {
        "run_id": "rob_sech2_W0p85",
        "profile": "sech2",
        "W0": 0.85,
        "in_robustness_set": true
      },
      {
        "run_id": "rob_sech2_W0p9",
        "profile": "sech2",
        "W0": 0.9,
        "in_robustness_set": true
      }
    ]
  },
  "controls": {
    "runs": [
      {
        "run_id": "t3_dr_0.0125",
        "dr": 0.0125
      },
      {
        "run_id": "t3_dr_0.025",
        "dr": 0.025
      },
      {
        "run_id": "t4_rmax_150",
        "r_max": 150
      },
      {
        "run_id": "t4_spongefrac_0.85",
        "sponge_start_frac": 0.85
      },
      {
        "run_id": "t5_lambda_0",
        "lambda_chi": 0
      },
      {
        "run_id": "t5_lambda_0.05",
        "lambda_chi": 0.05
      },
      {
        "run_id": "t5_lambda_0.1",
        "lambda_chi": 0.1
      },
      {
        "run_id": "t5_lambda_0.25",
        "lambda_chi": 0.25
      },
      {
        "run_id": "t5_lambda_0.5",
        "lambda_chi": 0.5
      }
    ]
  }
}
```

Run:

```bash
python curvature_well_capture_v0.3.py baseline --config baseline.json --out reproduce
python curvature_well_capture_v0.3.py null --config baseline.json --out reproduce
python curvature_well_capture_v0.3.py scan --config baseline.json --out reproduce
python curvature_well_capture_v0.3.py controls --config baseline.json --out reproduce
python curvature_well_capture_v0.3.py adjudicate --out reproduce
```

The first two runs establish the baseline and identical-seed linear control; `scan` runs the original robustness set, and `controls` reproduces the selected T3/T4/T5 parameter variations using the core solver. The baseline is not counted as an additional robustness member. R9-V-specific T1/T2 and T6 post-processing routines are not included; the results of those checks are summarized in the verification record. The archived R9-V control settings are reconstructed here from recorded parameters and do not reproduce the diagnostic wrapper's provenance fields or reference-row reuse. These commands have been checked against the solver's CLI and configuration schema, not newly executed for this repository update.

An optional built-in smoke test is available:

```bash
python curvature_well_capture_v0.3.py selftest --scratch selftest_scratch
```

### Outputs

- `reproduce/baseline/run_metadata.json` and `time_series.csv`: baseline spectral and dynamic results.
- `reproduce/null_linear/run_metadata.json`: linear control, including the seed identity and null-floor diagnostics.
- `reproduce/scan/` and `reproduce/controls/`: robustness and diagnostic configurations.
- `reproduce/capture_budget.csv` and `reproduce/scan_summary.csv`: adjudication outputs.

## Expected findings

| Measurement | Archived result |
| --- | --- |
| Baseline fundamental eigenvalue, omega_1 squared | 0.3640429192 |
| Baseline net capture, Linux | 0.235723% |
| Baseline net capture, Windows R0 | 0.233448% |
| Linear-null endpoint, original Linux | 2.27e-16 (fraction) |
| Original six-member robustness set, net capture | 0.1985–0.2491% |
| R9-V C1 configurations, net capture | 0.2334–0.2881% |
| R9-V C1 fundamental share of captured energy | 96.30–97.71% |
| Operational H2 capture threshold | 2% of incident energy |

**H1 supported:** localized bound modes exist within the prescribed-well model. **Dynamic capture supported at finite time:** weak nonlinear modal occupation occurs in the tested regime; most of the small captured fraction resides in the fundamental mode. **H2 rejected:** robust capture at the original 2% threshold is not supported by the tested configurations. A large share *of captured energy* is not a large fraction *of incident energy*. No claim of indefinite stability or experimental verification is made.

See [R9-V verification summary](verification/R9V_summary.md) for the limitations, acceptance decisions and definitions.

## Citation and license

If using the solver, cite Danai Leangtragul, *Localized Standing Waves in Prescribed Curvature Wells*, and reference this repository and the frozen solver version/commit. A publication DOI can be added once assigned.

The software and repository documentation are distributed under the [MIT License](LICENSE). No license assertion is made here about external manuscripts or separately archived data.

---
name: insar-gnss-validation
description: Validate InSAR LOS and vertical-land-motion velocity products against GNSS ENU rates using explicit look-vector conventions, strict pixel sampling, independent-station quality control, spatially grouped cross-validation, PyGMT diagnostics, and evidence-tiered reporting. Use for geodetic accuracy assessment; not for generic InSAR processing, VLM inversion, or GNSS time-series estimation.
---

# InSAR–GNSS validation

Use this skill to produce an auditable comparison, not merely a correlation plot. Preserve the distinction between internal consistency, out-of-sample fitting checks, and independent external accuracy.

## Route the task

1. For any new station table or raster bundle, read [data contracts](references/data-contracts.md) and run `scripts/preflight_inputs.py` before computing metrics.
2. For ENU-to-LOS projection, VLM-versus-Up comparison, exclusions, or cross-validation, read [methods](references/methods.md) and use `scripts/validate.py` with the relevant subcommand.
3. For maps, histograms, scatter plots, or RMS distributions, read [plotting](references/plotting.md) and use `scripts/plot_validation.py`. Plotting is optional; the numerical core must remain usable without PyGMT.
4. Before drafting a technical report or reviewer response, read [report guide](references/report-guide.md). Use `scripts/build_report_stub.py` only as a metrics-grounded starting point, then apply scientific judgment.

The scripts are portable reference implementations. Inspect `--help`, use explicit input and output paths, and adapt column flags instead of editing constants.

## Non-negotiable checks

- Require declared velocity units, vertical-positive convention, reference frame, observation window, LOS sign, incidence definition, and heading convention. Never infer them from filenames.
- For heading clockwise from north and incidence from vertical, use the documented look-vector equation exactly. A different convention requires an explicit conversion.
- Audit LOS, incidence, and heading grids for identical shape, CRS, affine transform, and pixel alignment before sampling.
- Sample the station-containing pixel. Do not silently substitute a nearby finite pixel. A display-only fallback must record its distance and remain outside strict statistics.
- Treat raster masks and NoData metadata as invalid; do not assume that a numerical zero is NoData.
- Define every residual as `InSAR − GNSS`. Retain the mean residual in headline bias and raw RMS. Report mean-removed RMS only as `cRMS`, a spatial-scatter diagnostic.
- Derive exclusions from declared rules. Never hard-code a station name, expected sample count, study region, product path, or validation result.
- Audit any supplied `used-in-fitting` status and fitting-site distance against a frozen fitting-station list or documented crosswalk; those columns are assertions, not self-validating evidence.
- Keep sites used to fit an InSAR correction separate from independent sites. Comparisons with fitting stations are internal consistency, not independent accuracy.
- Preserve the full ENU vector when projecting to LOS. Do not assume north motion is zero unless the user explicitly requests and labels that sensitivity test.
- Record all selection thresholds, exclusions, random seeds, software versions, hashes, sample counts, and product metadata in provenance output.

## Standard workflow

### 1. Preflight

Run the read-only audit first:

```bash
python scripts/preflight_inputs.py los \
  --stations stations.csv \
  --los final_los.tif \
  --incidence incidence.tif \
  --heading heading.tif
```

Stop for clarification when units, frame, sign, angle meaning, or coordinate CRS are missing. Do not auto-download data or install dependencies.

### 2. Calculate comparisons

Use one of:

```bash
python scripts/validate.py los ...
python scripts/validate.py vlm ...
python scripts/validate.py cross-validate ...
```

Each comparison writes deterministic tables plus `provenance.json`. By default, output directories must not already exist. An explicit `--overwrite` authorizes replacement of generated outputs only, never source data.

For LOS and VLM accuracy, report at least `N`, bias, raw RMS, cRMS, MAE, Pearson correlation, and the fraction within the declared tolerance when calculable. Keep excluded and unavailable stations in a disposition table with reasons.

For quadratic-surface validation:

- group all observations sharing a raster cell or physical site before splitting;
- use group-aware leave-one-group-out validation as the primary overfitting check;
- pool all held-out observation residuals once for the headline LOOCV RMS;
- use repeated group holdout only as a complementary stability test;
- call only held-out results out of sample;
- reject rank-deficient or ill-conditioned fits rather than returning an unstable surface.

### 3. Plot

Plots consume saved station-level and metric tables; they do not recalculate or redefine the formal sample. Keep raw residual statistics in the figure or caption even when a centered residual map is used for spatial diagnosis.

### 4. Interpret and report

Build an evidence ladder:

1. grouped out-of-sample validation of the GNSS-fitted correction;
2. internal comparison with fitting-period GNSS;
3. independent long-record GNSS LOS comparison;
4. independent long-record GNSS Up versus InSAR VLM;
5. optional local analyses with explicit spatial and temporal limitations.

State the strongest supported result and its limitation in the same paragraph. Never call GNSS ground truth, use cRMS to hide a bias, or claim validation beyond the sampled spatial and temporal domain.

## Verification

Run before delivery:

```bash
python scripts/run_self_test.py
python /path/to/skill-creator/scripts/quick_validate.py .
```

Also inspect every generated figure and verify that the numerical annotations equal the saved raw metrics. A passing schema test does not resolve ambiguous geodetic conventions.

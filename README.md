# InSAR–GNSS validation

A reusable agent skill and command-line toolkit for validating InSAR
line-of-sight (LOS) and vertical-land-motion (VLM) velocity products against
GNSS rates. Numerical analysis is independent of the optional PyGMT figures.

这是一个可复用的 InSAR–GNSS 精度验证 skill：包含 LOS 投影、严格像元采样、
独立站点筛选、分组交叉验证、PyGMT 绘图及技术报告规范。
本仓库仅分享通用方法、代码和合成测试，不包含论文原始数据、审稿材料或报告。

## Scope

- Project the complete GNSS East/North/Up vector into LOS using explicitly
  declared incidence, heading and sign conventions.
- Compare the station-containing InSAR pixel with GNSS; retain missing-pixel
  and quality-control exclusions rather than silently moving stations.
- Report bias, raw RMS, centered RMS, MAE and correlation with archived
  station-level residuals and provenance.
- Test quadratic spatial corrections with grouped leave-one-out validation
  and repeated group holdout (Monte Carlo cross-validation).
- Produce optional validation maps, histograms, scatter plots with GNSS
  uncertainty bars, and a metrics-grounded Markdown report starting point.

This is **not** an InSAR processing pipeline, GNSS time-series estimator, or
ready-made archive of the results of a particular study. Fitting-network
consistency, held-out fitting checks and independent accuracy remain separate
evidence tiers. Users must supply their own data and verify its provenance.

## Get the repository

```bash
git clone https://github.com/wangtwtusm/InSAR-GNSS-validation.git insar-gnss-validation
cd insar-gnss-validation
```

For agent use, [SKILL.md](SKILL.md) is the entry point. The repository root is
the skill directory; keep `scripts/`, `references/`, and `agents/` with it.
All numerical commands can also be used directly without an agent.

## Python environment and checks

Use Python 3.10 or later. An isolated environment is recommended:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/run_self_test.py
```

The deterministic tests generate synthetic inputs in temporary directories;
they do not download or require research data. Raster tests require rasterio
and otherwise report a skip, so check the test summary. The complete 17-test
suite was checked with Python 3.14.5, NumPy 2.4.6, pandas 3.0.2 and rasterio
1.5.0 for this public-sharing update. This is a tested environment, not a
claim that every dependency version has been tested.

Optional figures additionally need PyGMT, its compatible GMT installation,
xarray and rasterio. See [plotting](references/plotting.md) and run:

```bash
python scripts/plot_validation.py preflight
```

No command installs dependencies automatically. `build_report_stub.py` creates
Markdown, not a finished Word document; document assembly and visual review
remain separate steps.

## Input preparation

Read [data contracts](references/data-contracts.md) and
[methods](references/methods.md) before running an analysis. In particular:

- Use geographic longitude/latitude in degrees and a matching geographic
  raster CRS (for example EPSG:4326); explicitly prepare aligned inputs
  beforehand. The CLI checks CRS consistency and does not reproject coordinates.
- Declare units, vertical sign, reference realization, observation windows and
  whether stations were used in fitting. Audit fitting membership separately;
  the software cannot establish independence from a supplied flag alone.
- Supply the required ENU columns even for the VLM command. Optional uncertainty
  fields must contain **rate** uncertainties, not position scatter.
- Define selection thresholds before examining residuals. Example commands
  below do not by themselves establish an independent sample.

Metadata is supplied as a JSON file. The common required keys are
`velocity_units` (`mm/year` or `mm/yr`), `coordinate_crs`,
`vertical_positive` (`up`), `reference_frame`, `observation_window`, and
`comparison_role` (`internal` or `independent`). Use documented source values,
not placeholders, in an actual analysis. When using ENU sigma fields, also
declare `uncertainty_semantics` as
`one_standard_deviation_rate_uncertainty`.

LOS metadata additionally requires `los_sign`, `incidence_definition`,
`heading_convention`, `angle_units`, and `insar_product_stage`. The implemented
geometry uses `from_vertical`,
`satellite_flight_heading_clockwise_from_north`, and `degrees`. For LOS
uncertainty propagation, the current CLI requires the explicit diagonal
approximation `los_uncertainty_model` =
`diagonal_enu_rate_covariance_zero_cross_terms`; it does not use full covariance
columns automatically.

## Command examples

Paths below are examples referring to your own prepared inputs. Run commands
from the repository root and inspect each subcommand's `--help` for selection
options.

### VLM versus GNSS Up

```bash
python scripts/preflight_inputs.py vlm \
  --stations inputs/stations.csv \
  --vlm inputs/vlm.tif \
  --metadata inputs/metadata.json

python scripts/validate.py vlm \
  --stations inputs/stations.csv \
  --vlm inputs/vlm.tif \
  --metadata inputs/metadata.json \
  --output-dir outputs/vlm
```

### LOS versus projected GNSS

```bash
python scripts/preflight_inputs.py los \
  --stations inputs/stations.csv \
  --los inputs/los.tif \
  --incidence inputs/incidence.tif \
  --heading inputs/heading.tif \
  --metadata inputs/los_metadata.json

python scripts/validate.py los \
  --stations inputs/stations.csv \
  --los inputs/los.tif \
  --incidence inputs/incidence.tif \
  --heading inputs/heading.tif \
  --metadata inputs/los_metadata.json \
  --output-dir outputs/los
```

### Group-aware correction validation

```bash
python scripts/validate.py cross-validate \
  --observations inputs/correction_observations.csv \
  --metadata inputs/correction_metadata.json \
  --trials 100 --holdout-fraction 0.20 --seed 0 \
  --output-dir outputs/cross_validation
```

The observations table must explicitly define co-pixel/physical-site groups,
fixed normalized coordinates, and pre-fit residuals. See
[methods](references/methods.md) for grouping and residual-sign requirements.

### Optional VLM figure

```bash
python scripts/plot_validation.py vlm \
  --raster inputs/vlm.tif \
  --station-values outputs/vlm/station_values.csv \
  --metrics outputs/vlm/metrics.csv \
  --output outputs/vlm/validation.pdf
```

## Repository guide

| Location | Purpose |
| --- | --- |
| `SKILL.md` | Agent workflow and scientific safeguards |
| `scripts/preflight_inputs.py` | Input schema and raster alignment audit |
| `scripts/validate.py` | LOS/VLM metrics and grouped cross-validation |
| `scripts/validation_core.py` | Numerical functions |
| `scripts/plot_validation.py` | Optional PyGMT diagnostics |
| `scripts/build_report_stub.py` | Metrics-grounded Markdown report starter |
| `references/` | Data contracts, methods, plotting and report conventions |
| `tests/` | Offline synthetic unit and CLI tests |

## Sharing and limitations

The public repository contains generic code and documentation, not private
station tables, rasters, manuscripts or reviewer correspondence. When sharing
your own outputs, inspect metadata and provenance for sensitive information.
Record the commit used in a scientific analysis (`git rev-parse HEAD`).

The scripts support reproducibility; passing tests does not establish the
scientific suitability of a new dataset, justify extrapolation beyond the
validated domain, or resolve reference-frame and observation-window mismatch.

License: not yet specified by the repository owner.

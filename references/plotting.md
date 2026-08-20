# Optional PyGMT validation figures

Use this reference only when the user asks for a map or validation figure.
Numerical validation should already be complete: this plotting command verifies
standardized outputs, but it does not refit a correction, change station
selection, sample a nearest valid pixel, or remove an offset.

## Environment and safety

Run the command with the Python interpreter that contains NumPy, PyGMT, GMT 6,
and (for map commands) rasterio plus xarray. The command imports these
dependencies only after parsing arguments, so `--help` remains available in an
incomplete environment. It never installs packages.

GMT requires writable state directories. The command creates them under the
system temporary directory by default. Select another writable location by
placing the global option before the subcommand:

```bash
python scripts/plot_validation.py --gmt-workdir /path/to/work preflight
```

All input and output paths are explicit. An existing output is never replaced
unless `--overwrite` is supplied. Tests and examples should use a temporary
workspace, never a source-data directory.

## Standardized inputs

Column names and rate units are strict. Longitudes and latitudes are decimal
degrees; velocity and uncertainty values are millimetres per year. The raster
must be a single-band geographic GeoTIFF with valid CRS and nodata metadata.
Numeric zero is accepted as a valid velocity.

### VLM `station_values.csv`

Required columns:

- `station_id`
- `longitude_deg`, `latitude_deg`
- `insar_vlm_mm_per_year`
- `gnss_up_mm_per_year`
- `residual_mm_per_year`
- `sigma_up_mm_per_year`
- `include_in_metrics`

### LOS `station_values.csv`

Required columns:

- `station_id`
- `longitude_deg`, `latitude_deg`
- `insar_los_mm_per_year`
- `gnss_los_mm_per_year`
- `residual_mm_per_year`
- `sigma_los_mm_per_year`
- `include_in_metrics`

For both modes, `residual_mm_per_year` must equal the raw value
`InSAR - GNSS`. `include_in_metrics` accepts `true/false` or `1/0`. Excluded
rows may remain visible for diagnosis but never enter the histogram, scatter,
or metric check. An excluded diagnostic row may have a missing coordinate,
InSAR sample, GNSS value, residual, or rate uncertainty; the renderer skips any
symbol it cannot place. Every included row must be complete and finite. A
negative rate uncertainty is invalid.

The map preflight samples the station-containing raster pixel and compares it
with the InSAR value in the CSV. An included station on a masked pixel fails
preflight. There is deliberately no nearest-pixel fallback.

### `metrics.csv`

Required columns:

- `n`
- `bias_mm_per_year`
- `rms_mm_per_year`
- `centered_rms_mm_per_year`
- `mae_mm_per_year`
- `correlation`

The file may contain additional selection or product columns. If it contains
more than one row, select exactly one with repeated
`--metric-filter COLUMN=VALUE` options. The command recomputes all metrics from
the included station rows and also verifies
`RMS^2 = bias^2 + centered_RMS^2`. The annotation calls RMS **raw RMS** and
labels cRMS as a centered diagnostic; no offset is removed.

### `trials.csv`

Required columns:

- `trial`
- `orbit` (the default; another shared column may be selected with
  `--group-column`)
- `all_rms_mm_per_year`
- `holdout_rms_mm_per_year`

Each `(group, trial)` pair must be unique. The optional `loocv_metrics.csv` has
one row per group and requires the shared group column plus
`rms_mm_per_year`. Only the held-out distribution is fully out of sample;
all-station RMS combines fitting and held-out observations.

## Preflight

Dependency-only check:

```bash
python scripts/plot_validation.py preflight
```

Full VLM input check without plotting:

```bash
python scripts/plot_validation.py preflight \
  --mode vlm \
  --raster /path/to/vlm.tif \
  --station-values /path/to/station_values.csv \
  --metrics /path/to/metrics.csv
```

Use `--mode los` for LOS inputs. Use `--mode cv-rms --trials ...` for
cross-validation trials. A successful preflight reports the dependency
versions, exact-pixel check, station count, and raw metrics as JSON.

## VLM comparison: 2 x 2

```bash
python scripts/plot_validation.py vlm \
  --raster /path/to/vlm.tif \
  --station-values /path/to/vlm_station_values.csv \
  --metrics /path/to/vlm_metrics.csv \
  --output /path/to/vlm_validation.pdf
```

The panels are fixed by scientific role:

- A: the InSAR VLM raster with GNSS Up squares;
- B: the same VLM raster with raw InSAR-minus-GNSS residual squares;
- C: the raw residual histogram for included stations;
- D: InSAR on x, GNSS on y, with vertical GNSS-rate uncertainty bars.

Maps draw the velocity raster first, then light-blue water and coastlines, then
station symbols. The VLM and residual colour maps are separate, symmetric
about zero, and use `vik`. Supply `--vlm-limit` and `--residual-limit` for fixed
cross-figure comparability; otherwise the command uses the configured robust
absolute-value percentile. `--scatter-range MIN MAX` fixes the same range on
both scatter axes.

The default water fill deliberately masks deformation pixels over mapped
water. Use `--water-fill none` only for a documented diagnostic raster that is
intended to remain visible offshore; the coastline is still drawn.

Map circles encode GNSS rate uncertainty as radius. The option
`--uncertainty-radius-scale` is centimetres of circle radius per
millimetre-per-year uncertainty; the script doubles that value for GMT's
diameter-based circle size. Scatter error bars are vertical because GNSS is on
the y axis.

## LOS validation map

```bash
python scripts/plot_validation.py los \
  --raster /path/to/los.tif \
  --station-values /path/to/los_station_values.csv \
  --metrics /path/to/los_metrics.csv \
  --title "Ascending LOS validation" \
  --output /path/to/los_validation.png
```

The LOS raster is the bottom layer. Raw residual squares and GNSS LOS-rate
uncertainty circles are drawn above the coast. The in-map annotation reports
the retained bias and raw RMS, with cRMS explicitly marked as diagnostic.
Background and residual colour bars are distinct and both symmetric about
zero. No centered residual is plotted.

## Repeated-holdout RMS distributions

```bash
python scripts/plot_validation.py cv-rms \
  --trials /path/to/trials.csv \
  --metrics /path/to/loocv_metrics.csv \
  --group-column orbit \
  --output /path/to/cv_rms.pdf
```

Blue dots show all-station RMS, orange dots show held-out RMS, black bars show
medians, and optional magenta dashed segments show group-wise LOOCV RMS. The
point jitter is visual only and deterministic under `--seed`.

## Figure invariants

- Panel letters and descriptive subtitles are bold and outside the coordinate
  frame; frame titles are not used.
- Frames are plain, drawn once, with annotations on west and south edges.
- Displacement and residual colour scales use zero-centred symmetric `vik`.
- Colour bars are horizontal and below the map axes.
- Water is light blue; land is not painted over the InSAR raster.
- The raster is always drawn before coastlines and station symbols.
- All accuracy annotations retain the residual mean. Centering is permitted
  only as a clearly labelled diagnostic statistic already present in the
  standardized metrics.

After every production render, inspect the rasterized output at 100% zoom for
clipped external panel titles, overlapping legends or colour bars, double
frames, invisible uncertainty circles, and aliasing. Increase panel size or
DPI rather than accepting undersampled dense InSAR texture.

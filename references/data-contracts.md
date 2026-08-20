# Data Contracts

Use this reference when ingesting, normalizing, sampling, or archiving InSAR and GNSS data. The objective is to make every plotted point and reported statistic traceable to an immutable source and an explicit transformation.

## Contract Principles

- Preserve raw downloads and supplied rasters unchanged. Write normalized and derived products separately.
- Store machine-readable tables in UTF-8 CSV or Parquet with one header row and no implicit index.
- Give every raster, station solution, comparison sample, and validation run a stable identifier.
- Record units, sign convention, terrestrial-frame realization, estimator, and observation interval as data, not prose-only assumptions.
- Use `null` or IEEE `NaN` for unavailable values. Never encode missing measurements as zero.
- Keep excluded observations in derived tables with explicit flags and reasons; do not silently delete them.
- Make headline statistics reproducible from archived per-observation values.

## Coordinate Velocity and Angle Conventions

Unless the source requires another convention, normalize to the following and retain the original convention in provenance metadata.

| Quantity | Required convention |
|---|---|
| Longitude | Decimal degrees, east-positive, range `[-180, 180)` |
| Latitude | Decimal degrees, north-positive, range `[-90, 90]` |
| Geographic CRS | EPSG:4326 or a fully specified equivalent |
| East velocity | Positive east, `mm yr^-1` |
| North velocity | Positive north, `mm yr^-1` |
| Up velocity or VLM | Positive upward, `mm yr^-1` |
| Displacement or residual scatter | `mm`; never label it as a rate |
| Incidence angle | Degrees from vertical unless explicitly declared otherwise |
| Heading angle | Degrees clockwise from north unless explicitly declared otherwise |
| LOS velocity | `mm yr^-1`; declare whether positive is toward or away from the radar |
| Residual | InSAR minus GNSS by default; record the sign definition in every result table |
| Rate uncertainty | Published or estimated one-standard-deviation value in `mm yr^-1` |

Do not infer an angle convention from a filename. Verify it from metadata, code, or a known-vector projection test. Sample LOS, incidence, and heading from the same raster cell unless the product documentation explicitly defines a different geometry grid.

## Time and Frame Metadata

Every velocity solution must carry:

- `start_epoch` and `end_epoch`, preferably ISO dates; decimal years may be retained in additional fields;
- `duration_years` calculated by a documented convention;
- `n_epochs` and, when relevant, cadence or sampling density;
- `velocity_estimator` and modeled terms, offsets, discontinuities, and noise assumptions;
- `frame_realization`, such as an IGS or ITRF realization, and reference epoch when applicable;
- `processing_center`, product version, and retrieval snapshot time.

Never silently compare or merge rates from different realization, time window, estimator, or discontinuity models. The comparison may still be valid, but those differences must remain explicit covariates and caveats.

## Raster Inventory Schema

Maintain one row per raster or raster band.

| Field | Type | Requirement |
|---|---|---|
| `product_id` | string | Stable unique identifier |
| `product_role` | enum | `los_velocity`, `vlm`, `east_velocity`, `north_velocity`, `incidence`, `heading`, `mask`, or documented extension |
| `track_id` | string or null | Orbit or mosaic identifier |
| `orbit_direction` | enum or null | `ascending`, `descending`, or `not_applicable` |
| `data_locator` | string | Relative archive locator or stable URI; no credentials |
| `sha256` | string | Hash of archived bytes |
| `crs` | string | EPSG code or complete WKT identifier |
| `width`, `height` | integer | Pixel dimensions |
| `transform` | string | Full affine transform or equivalent georeferencing |
| `pixel_size_x`, `pixel_size_y` | number | CRS units per pixel |
| `nodata_value` | number or null | Explicit nodata encoding |
| `units` | string | Physical units |
| `positive_direction` | string | Required for LOS and vertical products |
| `start_epoch`, `end_epoch` | date or null | Observation interval |
| `frame_realization` | string or null | Terrestrial-frame realization |
| `correction_state` | string | Raw, long-wavelength corrected, decomposed, or other precise state |
| `parent_product_ids` | list-like string | Inputs used to derive the raster |
| `processing_record_id` | string | Link to provenance record |

For a matched raster stack, additionally verify and record whether CRS, dimensions, affine transform, pixel registration, and nodata mask are identical. Do not assume one-to-one pixels from similar bounds alone.

## GNSS Station Rate Schema

Maintain one row per station solution and velocity-estimator version.

| Field | Type | Units or values |
|---|---|---|
| `solution_id` | string | Stable unique identifier |
| `station_id` | string | Provider identifier |
| `physical_site_id` | string | Audited co-location or monument cluster |
| `longitude_deg`, `latitude_deg` | number | EPSG:4326 |
| `east_mm_per_year`, `north_mm_per_year`, `up_mm_per_year` | number or null | ENU rates |
| `sigma_east_mm_per_year`, `sigma_north_mm_per_year`, `sigma_up_mm_per_year` | number or null | One-standard-deviation rate uncertainty |
| `cov_en`, `cov_eu`, `cov_nu` | number or null | Velocity covariance in `(mm yr^-1)^2` |
| `start_epoch`, `end_epoch`, `duration_years` | date and number | Observation support |
| `n_epochs` | integer or null | Number of position epochs |
| `frame_realization` | string | Exact realization |
| `velocity_estimator` | string | Named estimator and version |
| `processing_center` | string | Provider or processing branch |
| `used_in_fitting` | boolean | Whether this solution constrained InSAR correction or decomposition |
| `source_record_id` | string | Link to provenance |

If an archive contains fields named generically such as `std_e`, `std_n`, or `std_u`, determine whether they are rate uncertainty, coordinate scatter, or post-fit displacement-residual standard deviation. Rename normalized fields to encode both quantity and units. Never promote a displacement scatter in `mm` to a velocity uncertainty in `mm yr^-1`.

## Station Raster Sample Schema

Maintain one row per station, product, and sampling rule.

| Field | Type | Meaning |
|---|---|---|
| `sample_id` | string | Stable unique identifier |
| `solution_id`, `product_id` | string | Parent records |
| `row_index`, `column_index` | integer or null | Sampled cell |
| `pixel_longitude_deg`, `pixel_latitude_deg` | number or null | Cell center |
| `sample_distance_m` | number or null | Station-to-cell-center distance |
| `sampling_rule` | enum | `station_cell`, `nearest_valid_display_only`, or documented extension |
| `station_cell_finite` | boolean | Exact station cell is valid |
| `inside_raster` | boolean | Coordinate is inside the grid bounds |
| `inside_auxiliary_hull` | boolean or null | Required when interpolation support matters |
| `insar_value_mm_per_year` | number or null | Sampled velocity |
| `incidence_deg`, `heading_deg` | number or null | Co-located geometry |
| `gnss_projected_los_mm_per_year` | number or null | Full-vector projection |
| `sigma_projected_los_mm_per_year` | number or null | Propagated one-standard-deviation uncertainty |
| `residual_insar_minus_gnss_mm_per_year` | number or null | Signed residual |
| `included_in_statistics` | boolean | Final statistical inclusion |
| `exclusion_reasons` | list-like string | All applicable reasons |

Nearest-valid sampling must not replace exact-pixel sampling in strict validation. A nearest-valid point may be retained for display only if its distance and exclusion from statistics are explicit.

## Vertical Comparison Schema

Use a separate view or table when the comparison is specifically VLM versus GNSS Up.

Required fields are:

```text
solution_id
product_id
station_id
physical_site_id
longitude_deg
latitude_deg
gnss_up_mm_per_year
sigma_up_mm_per_year
insar_vlm_mm_per_year
residual_insar_minus_gnss_mm_per_year
station_cell_finite
inside_auxiliary_hull
used_in_fitting
included_in_statistics
exclusion_reasons
```

## Cross Validation Records

Archive both assignments and predictions. A summary alone is insufficient.

### Fold or Trial Table

```text
validation_run_id
method
orbit_or_product_id
random_seed
repeat_id
fold_id
group_id
solution_id
assignment
n_fitting_groups
n_held_out_groups
model_specification_id
```

The `assignment` field is `fitting` or `held_out`.

### Prediction Table

```text
validation_run_id
repeat_id
fold_id
group_id
solution_id
observed_gnss_mm_per_year
uncorrected_insar_mm_per_year
predicted_correction_mm_per_year
corrected_insar_mm_per_year
held_out_residual_mm_per_year
```

Each leave-one-group-out observation must appear exactly once as held out. For repeated holdout, the assignment table must make the random split reproducible from the recorded seed and grouping universe.

## Metric Summary Schema

Use one row per product, orbit, sample selection, and metric scope.

```text
validation_run_id
product_id
orbit_direction
selection_name
evidence_tier
n_observations
n_unique_physical_sites
bias_mm_per_year
rms_mm_per_year
centered_rms_mm_per_year
mae_mm_per_year
nmad_mm_per_year
pearson_r
fraction_abs_residual_within_threshold
threshold_mm_per_year
offset_removed
notes
```

The primary result must retain the residual mean unless the method explicitly defines an offset removal. Centered metrics require a separate field and label.

## Provenance Record Schema

Every raw or derived artifact should link to a provenance record containing:

```text
record_id
artifact_id
source_provider
source_product
source_version
source_url_or_citation
retrieved_at_utc
license_or_terms
sha256
byte_size
parent_artifact_ids
processing_script
processing_version_or_commit
parameters
runtime_environment
created_at_utc
```

Store no passwords, tokens, cookies, private keys, or signed URLs in a manifest or report.

## Validation Gates

Before analysis:

- confirm required fields, uniqueness constraints, coordinate bounds, and units;
- verify raster CRS, transform, registration, and nodata behavior;
- test LOS projection with a known ENU vector and the declared angle and sign convention;
- audit station identifiers, coordinates, and physical-site duplicates;
- verify time interval, estimator, realization, and uncertainty semantics;
- confirm fitting membership independently of the validation selection;
- ensure each retained station has all data required by the stated selection rule.

Before reporting:

- recompute every table metric from archived per-observation residuals;
- verify that excluded observations contribute to no statistic designated as strict;
- confirm residual sign and offset policy in data and caption text;
- confirm that sample counts agree across tables, figures, and prose;
- hash the final inputs, derived tables, figures, and report.

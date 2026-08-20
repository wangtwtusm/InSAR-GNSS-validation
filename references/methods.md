# Validation Methods

Use this reference when designing or auditing InSAR LOS or VLM validation against GNSS. Read [Data Contracts](data-contracts.md) first when inputs have not yet been normalized.

## Evidence Tiers

Keep the following evidence classes distinct in methods, tables, figures, and claims.

1. **Out-of-sample validation of the fitted correction** tests whether a GNSS-derived spatial model generalizes to fitting-network sites withheld from the fit.
2. **Internal consistency** compares a final product with observations that contributed to its correction or construction. It is useful, but not independent accuracy.
3. **Independent external validation** uses quality-screened stations not used in fitting, with explicit spatial separation and de-duplication. It bounds the accuracy of the final product.
4. **Local descriptive analysis** evaluates a named subregion or gradient. Small local samples usually support pattern direction more strongly than absolute amplitude or cause.

Never merge tiers into one headline metric. Report the tier with every result.

## LOS Projection

For east-positive `E`, north-positive `N`, up-positive `U`, incidence `i` measured from vertical, heading `h` clockwise from north, and LOS positive toward the radar, use

```math
v_LOS = -cos(h) sin(i) E + sin(h) sin(i) N + cos(i) U.
```

Define the look coefficients

```math
l_E = -cos(h) sin(i),
l_N =  sin(h) sin(i),
l_U =  cos(i).
```

If the source defines LOS positive away from the radar, reverse the full sign consistently. Do not change a single coefficient in isolation.

Sample LOS, heading, and incidence from the same pixel. If geometry grids differ, resample once with a documented method and archive the aligned grids before station sampling.

## LOS Uncertainty Projection

With GNSS velocity covariance matrix `C_ENU` and look vector `l`, the projected variance is

```math
sigma_LOS^2 = l^T C_ENU l.
```

Expanded,

```math
sigma_LOS^2 = l_E^2 sigma_E^2 + l_N^2 sigma_N^2 + l_U^2 sigma_U^2
              + 2 l_E l_N cov_EN + 2 l_E l_U cov_EU + 2 l_N l_U cov_NU.
```

If covariance terms are unavailable, a diagonal approximation may be used only if the independence assumption is stated. Do not substitute coordinate scatter or post-fit displacement-residual standard deviation for formal velocity covariance.

## North Motion Correction and Two Look Decomposition

When ascending and descending LOS are decomposed into East and Up, setting North to zero can leak real North motion into the recovered components. Given an independently specified or interpolated North field, remove its predicted contribution first:

```math
b_k = v_LOS,k - l_N,k N,  for k in {ascending, descending}.
```

Then solve at each jointly valid pixel:

```math
[l_E,a  l_U,a] [E] = [b_a]
[l_E,d  l_U,d] [U]   [b_d].
```

Requirements:

- use geometry from the corresponding LOS pixel for each orbit;
- require both LOS values and both geometry pairs to be finite;
- document interpolation support and flag extrapolated pixels or sites;
- do not use independent validation rates or validation-derived offsets in the recalculation;
- quantify the correction increment separately from product uncertainty;
- assess determinant or condition number and mask unstable geometry.

After defining the recalculated product once, use a stable name such as `InSAR VLM` throughout the report.

## Spatial Correction Model

A common long-wavelength correction is a six-parameter quadratic in normalized coordinates:

```math
q(x, y) = beta_0 + beta_1 x + beta_2 y + beta_3 xy + beta_4 x^2 + beta_5 y^2.
```

Normalize coordinates using constants fixed for the whole analysis, not separately within each fold. Fit by least squares unless another estimator or weighting scheme is pre-specified. Include the intercept in every fold. Define and archive whether the correction is added to or subtracted from the InSAR field; never infer the sign from the fitted coefficients.

Stations mapped to the same raster cell are one validation group. This prevents nearly duplicate observations from leaking between fitting and validation samples.

## Group Aware Leave One Out Validation

For each unique raster-cell group `g`:

1. remove every observation in `g`;
2. fit the full model specification to all remaining groups;
3. predict the correction at the held-out locations;
4. construct the held-out corrected InSAR value using the pipeline's declared correction sign;
5. compute residual `e_j = corrected InSAR_j - GNSS_j` for each held-out observation.

After every group has been left out once, calculate

```math
bias_LOOCV = (1/N) sum_j e_j,
RMS_LOOCV  = sqrt((1/N) sum_j e_j^2).
```

Invariants:

- each observation appears exactly once in the combined held-out result;
- all co-pixel observations are withheld together;
- coordinate normalization, model terms, and correction sign are fixed across folds;
- no held-out observation influences the fit, intercept, hyperparameters, or offset;
- do not subtract fold means or a final mean residual;
- aggregate squared residuals over held-out observations; do not average fold-level RMS values;
- report both full-fit and LOOCV metrics, but label only LOOCV out of sample.

Leave-one-site-out validation assesses generalization within the spatial support of the fitting network. It does not validate regions with no nearby constraints or prove a localized amplitude correct.

## Repeated Site Wise Monte Carlo Cross Validation

Use repeated group holdout as a complementary stability test, not as a replacement for group-aware leave-one-out validation.

For each repeat:

- sample a pre-specified fraction of unique site groups without replacement for validation;
- hold out every observation belonging to those groups;
- fit the unchanged model to the remaining groups;
- compute metrics separately for held-out observations and, if useful, all observations;
- use a recorded random seed and archive every assignment.

The held-out distribution is out of sample. The all-observation distribution mixes fitting and held-out observations and must not be described as out of sample.

For repeat-level held-out RMS values `RMS_r`, report their mean, median, percentile interval, and sample standard deviation:

```math
s_RMS = sqrt((1/(R - 1)) sum_r (RMS_r - mean(RMS))^2).
```

This is split-to-split variability in validation performance. It is not the standard deviation of individual residuals, a formal InSAR uncertainty, or a confidence interval on the mean unless a justified inferential procedure defines one.

Because groups are sampled without replacement, call the method **repeated site-wise Monte Carlo cross-validation**, not bootstrap validation.

## Independent Station Selection

Define selection before inspecting residuals. Apply filters in an auditable order and retain every flag.

1. **Not used in fitting**: exclude every station solution that constrained correction, decomposition, or auxiliary fields relevant to the comparison.
2. **Coordinate audit**: reconcile station identity, monument changes, and suspicious coordinate discrepancies.
3. **Minimum duration**: require the declared record length using the source solution interval.
4. **Component-wise uncertainty**: require finite ENU uncertainties and apply the threshold to each component, not only Up or projected LOS.
5. **Spatial independence**: require the declared minimum distance from all fitting sites.
6. **Physical-site de-duplication**: within a declared radius, retain one representative by a pre-specified rule such as longest valid record.
7. **Exact coverage**: require a finite station-containing LOS or VLM pixel and all necessary geometry.
8. **Interpolation support**: identify sites outside an auxiliary-field convex hull and report a strict sensitivity sample without them.

Do not choose an uncertainty threshold, distance, clustering radius, or outlier list after examining which choice improves RMS. High-uncertainty sites may remain visible for diagnosis, but they must contribute to none of the strict metrics.

## Metric Definitions

For residuals `e_j = InSAR_j - GNSS_j`:

```math
bias = mean(e)
RMS  = sqrt(mean(e^2))
cRMS = sqrt(mean((e - bias)^2))
MAE  = mean(abs(e))
NMAD = 1.4826 median(abs(e - median(e))).
```

Also report Pearson correlation between InSAR and GNSS and the fraction satisfying a pre-specified absolute-residual threshold.

Use raw bias and RMS as the primary absolute-agreement measures. cRMS is a diagnostic for spatial scatter after removing the sample mean. A high correlation or low cRMS does not erase a non-zero bias.

Do not remove a constant offset in the principal accuracy calculation unless offset removal is part of the pre-declared measurement model. If a centered residual map is useful for spatial diagnosis, label it mean-removed and pair it with raw statistics.

## Internal Vertical Consistency

When GNSS stations contributed to the InSAR correction or VLM construction, compare GNSS Up and VLM only as internal consistency.

- require the exact finite station cell and jointly valid decomposition inputs;
- explain why other fitting-network stations are not comparable;
- report the same residual sign and raw metrics used in external validation;
- do not call archived displacement-residual scatter a rate uncertainty;
- do not infer pointwise VLM accuracy from agreement with reused constraints.

## Reference Realization and Time Window

A realization change alters the representation of terrestrial-frame origin, scale, orientation, rates, and associated processing models. It is not a physical station displacement. Differences between two velocity products may also contain observation-window, estimator, equipment, discontinuity, antenna-model, and epoch effects.

Therefore:

- identify both realizations and observation intervals;
- do not label a practical product difference as a pure frame transformation without performing that transformation;
- do not automatically interpret a coherent component difference as geophysical VLM;
- treat a long-record GNSS rate as an external benchmark with its own temporal support, not timeless ground truth.

## Local and Regional Analysis

For a focused subregion:

- define the geographic boundary before selecting stations;
- list finite exact-pixel stations and NoData sites;
- report small-sample metrics as descriptive, with `N` prominent;
- compare point GNSS support with the InSAR pixel or effective spatial support;
- distinguish the sign of a gradient from its absolute amplitude;
- consider time-window mismatch and localized deformation as alternative explanations;
- do not assert a physical mechanism without independent coherence, stability, or time-series evidence.

## Method Audit Checklist

- LOS sign and angle convention verified by a known-vector test.
- All ENU components projected; North is not silently set to zero.
- Same-cell LOS and geometry sampling verified.
- Co-pixel station grouping enforced in every cross-validation split.
- Random seed, assignments, model specification, and normalization archived.
- Independent sample not used in any fitting or validation-derived correction.
- Quality, duration, distance, cluster, exact-pixel, and hull rules applied consistently.
- Raw bias and RMS retained; centered results labeled diagnostic.
- Every statistic recomputable from archived per-observation residuals.

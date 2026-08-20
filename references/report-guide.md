# Validation Report Guide

Use this reference when preparing a technical report, manuscript supplement, or reviewer response from completed InSAR-GNSS validation. Read [Validation Methods](methods.md) for metric and selection rules and [Data Contracts](data-contracts.md) for auditable schemas.

## Writing Objective

Organize the report by **independence of evidence**, not by the order in which analyses were performed. A reader should be able to distinguish model generalization, internal consistency, independent external accuracy, and local descriptive evidence without reading code.

Every major result paragraph should contain:

1. the result and sample size;
2. the evidence tier;
3. the spatial and temporal scope;
4. the principal limitation.

Lead with the outcome, then give the method detail needed to interpret it.

## Recommended Report Structure

### Title Block

Include a report label, descriptive title, one-line purpose, analysis snapshot date, and the principal product versions, frame realizations, and observation windows.

### Methodological Update

Use this section only when the evaluated product was recalculated or corrected. State:

- the original assumption;
- the revised equation or processing step;
- which data entered the recalculation;
- which independent validation data and validation-derived offsets did not enter;
- the quantitative magnitude of the increment;
- whether the principal patterns or conclusions changed.

Define the final product name once. Thereafter use a stable term such as `InSAR VLM` rather than repeatedly advertising a correction label.

### Executive Summary

Frame the concern, identify the evidence ladder, state the strongest defensible conclusion, give the independent external accuracy bound, and name the main remaining limitation. Add a compact headline table with columns such as `Test`, `Sample`, `Result`, and `Interpretation`.

### Scope Products and Statistical Conventions

Define data roles, observation periods, frame realizations, LOS projection, sampling rule, residual sign, metric formulas, quality screen, and offset policy before showing results. Use a definitions table when several sections reuse the same terms.

### Fitting Network Validation

Present group-aware leave-one-site-out validation first because it directly addresses overfitting. Follow with repeated site-wise Monte Carlo cross-validation as a stability test. End with internal VLM consistency and clearly state its non-independence.

### Independent Long Record Validation

Introduce the independent product and estimator, cross-match and duration screen, realization and time-window differences, pre-specified uncertainty exclusions, independent LOS results, and independent VLM results. Include a strict sensitivity sample when interpolation or extrapolation support differs.

Verify non-use and fitting-site distances against an archived fitting-station
list or crosswalk. If the analysis only receives precomputed status or distance
columns, identify them as provider assertions rather than independently
recomputed facts. When station-specific start and end epochs are unavailable,
state explicitly that temporal overlap could not be audited station by station.

### Focused Local Analysis

Make this optional. Include it only when the user requests it or the available sample can answer a defined local question. Report sparse local evidence descriptively and distinguish gradient direction from amplitude and cause.

### Data Provenance and Reproducibility

End with an inventory of inputs, roles, archived versions or snapshot dates, derived tables, analysis scripts, and figure sources. State whether reported numbers were read from saved summaries or recomputed from saved per-observation data.

## Terminology

Prefer these terms:

- **fitting sites**, **fitting observations**, and **fitting constraints**; avoid `training`;
- **InSAR-period GNSS** for the network used in correction or decomposition;
- **independent GNSS** only after non-use, quality, spatial-separation, and de-duplication rules are satisfied;
- **group-aware station-wise leave-one-out validation**;
- **repeated site-wise Monte Carlo cross-validation (MCCV)** for without-replacement group holdout;
- **raw RMS** or **datum-retaining RMS** for the primary statistic;
- **mean-removed RMS (cRMS)** for diagnostic centered scatter;
- **internal consistency**, **external-station check**, **qualified support**, and **external accuracy bound**.

Distinguish `GNSS Up` from `InSAR VLM`. Spell out **standard deviation (SD)** at first use and name the quantity whose SD is reported.

## Numeric and Unit Style

- Use `mm yr^-1` in source data and code; render the final document as `mm yr⁻¹` or the equivalent OMML structure.
- Use signed values for bias and signed residual summaries.
- Use a true minus sign in typeset mathematics.
- Use en dashes for date spans and numeric ranges in prose.
- Choose precision from measurement support; rates commonly use two decimals, very small biases or correlations may justify more.
- Put shared units in the caption or note rather than repeating them in every cell.
- Always show `N` and define whether it counts observations, raster-cell groups, or unique physical sites.

## Table Conventions

Place the caption above the table as `Table N. Descriptive title.` Place a concise note below when needed to define units, residual sign, offset policy, sample selection, or exclusions.

Recommended validation columns are:

```text
Selection or orbit | N | Bias | RMS | cRMS | MAE | r | |e| within threshold
```

Additional conventions:

- use an `Interpretation` column in the headline table;
- report raw bias and RMS before centered metrics;
- never show cRMS as the only accuracy statistic;
- distinguish the broad eligible sample from the formal independent quality-controlled sample;
- identify uncertainty as rate uncertainty, temporal displacement-residual SD, spatial SD, or trial-to-trial SD;
- list excluded stations or observations with the pre-specified reason in a QC table when the list is short enough to audit;
- include a provenance table with `Item`, `Role`, and `Version or archived source`;
- use content-based column widths, repeat header rows, prevent row splitting, and keep numeric columns aligned consistently.

## Figure Conventions

Place each figure with its caption immediately after the relevant method or table, followed by the interpretation paragraph. A caption must define the panels, background, marker values, residual type, uncertainty encoding, quality-control outlines, time window, regional boundary, and any display-only exception.

### Maps

- Draw the raster background before coastlines and symbols.
- Use the same background and color scale for paired value and residual maps when direct comparison is intended.
- Use a separate, symmetric residual color scale centered at zero when appropriate.
- Use colored squares for station values or residuals and outlined circles for uncertainty; include a quantitative size legend.
- Keep excluded high-uncertainty sites visible only when diagnostically useful, distinguish them with a special outline, and state that they are excluded from statistics.
- Do not label excluded stations individually unless a specific diagnostic requires it.
- Mark NoData sites explicitly rather than moving them to a valid pixel.
- Panel letters and subtitles belong outside the coordinate frame.

### Scatter and Distribution Panels

- Prefer InSAR on the horizontal axis and GNSS on the vertical axis.
- Show vertical GNSS rate-uncertainty bars when formal uncertainties exist.
- Include a one-to-one reference line and report the same sample selection used in the table.
- Plot the raw residual histogram for absolute accuracy. If a mean-removed distribution is shown, label it diagnostic and retain raw metrics nearby.

### Time Series

- Show individual epochs in a neutral color.
- Add a robust local summary only where sampling supports it.
- Distinguish published velocity trends from rates refitted for the report.
- Annotate documented discontinuities, equipment events, or long gaps when they explain quality flags.
- Shade the InSAR analysis interval so temporal overlap is visible.

### Layout and Accessibility

- Keep image and caption together and use deliberate page breaks for large multipanel figures.
- Avoid redundant colorbars and legends.
- Provide meaningful alternative text describing the analytical purpose, not only the filename.
- Verify that symbols, legends, error bars, and panel titles remain readable at final page size.

## Word Math and OMML

In a Word deliverable, equations and mathematical symbols must be editable Office Math Markup Language.

- Use inline `m:oMath` for equations embedded in prose, captions, notes, and table cells.
- Use `m:oMathPara` only for a genuinely standalone display equation.
- Build fractions, radicals, sums, hats, subscripts, and superscripts structurally; do not simulate them with ASCII, Unicode superscripts, or images.
- Keep variables in mathematical italic and functions, labels, and units upright.
- Preserve the surrounding paragraph, run formatting, table geometry, revision wrappers, fields, comments, and unrelated text.
- Do not convert page-number or cross-reference fields into mathematics.
- After an OOXML edit, verify package integrity, formula counts, visible text equivalence, and a full Word-to-PDF or PNG render.

Examples of content that should be OMML in the final Word file include LOS projection, RMS definitions, cross-validation formulas, quality inequalities, `N`, `r`, sigma with component subscripts, plus-minus expressions, and rate units with an inverse-year exponent.

## Captions and Notes

A table note or figure caption should answer the questions a skeptical reader would otherwise ask:

- What is the residual sign?
- Was a mean or constant offset removed?
- Which sample and quality rule were used?
- Are values raw or centered?
- What does marker size or outline mean?
- Are shown exceptions included in statistics?
- What observation interval and frame apply?

Do not hide a material selection rule or caveat only in the methods section.

## Claims and Anti Claims

### Supported Claim Patterns

- “Held-out tests reveal no evidence of substantial overfitting within the spatial domain sampled by the fitting network.”
- “Independent long-record stations provide an external accuracy bound for the final GNSS-assisted product.”
- “The correction improves the physical fidelity of the decomposition without materially changing the principal spatial patterns.”
- “The local observations support the sign of the gradient more strongly than its absolute amplitude.”
- “The internal comparison quantifies consistency with the observations used in processing; it is not independent validation.”

### Claims to Avoid

Do not state or imply that:

- cross-validation proves overfitting impossible;
- validation within the station network tests areas without spatial support;
- every pixel or localized amplitude is independently validated;
- reused fitting stations provide independent absolute accuracy;
- an external-station check of a GNSS-assisted mosaic validates raw GNSS-independent InSAR;
- high correlation or low cRMS cancels a non-zero bias;
- coordinate scatter or displacement-residual SD is velocity uncertainty;
- a long-record GNSS product is ground truth or automatically superior;
- a realization change is physical land motion or an untransformed product difference is a pure frame transformation;
- nearest-valid display samples belong in exact-pixel statistics;
- a plausible local process is demonstrated without independent temporal or spatial evidence.

## Reviewer Response Pattern

Use calibrated language that answers the criticism directly:

> We tested the fitted correction with group-aware held-out predictions and obtained results consistent with repeated site-wise cross-validation. This provides no evidence of material overfitting within the sampled domain. Independent long-record GNSS then supplies an external accuracy bound for the final product. These tests do not establish that every localized amplitude is correct, so the revised interpretation preserves the regional pattern while qualifying local magnitude and long-term representativeness.

Replace generic assertions with the actual sample, metrics, and limitations. Do not copy this wording when the evidence tier or result differs.

## Report Quality Checklist

- Evidence tiers are explicit and not blended.
- Terminology uses `fitting`, not `training`.
- Every result includes sample definition and `N`.
- Residual signs, units, offset policy, and quality rule are consistent everywhere.
- Raw bias and RMS are visible wherever centered metrics appear.
- Internal comparisons are not described as independent accuracy.
- Observation intervals, estimators, and frame realizations are stated.
- Exclusions and display-only exceptions are auditable.
- Tables and prose reproduce archived per-observation calculations.
- Captions are self-contained and panel labels are outside axes.
- Word mathematics is OMML and the final document has passed structural and visual QA.

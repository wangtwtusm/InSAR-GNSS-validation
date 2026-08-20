#!/usr/bin/env python3
"""Numerical core for reproducible InSAR--GNSS validation.

The module deliberately contains no project paths, file I/O, plotting, or
import-time side effects.  Callers are responsible for loading and quality
controlling their data, recording units/reference frames, and serializing the
returned results.

The LOS convention implemented here is explicit: ``heading`` is the satellite
flight heading measured clockwise from north, incidence is measured from the
upward vertical, and a positive ``toward_satellite`` LOS has coefficients

``[-cos(heading) sin(incidence), sin(heading) sin(incidence), cos(incidence)]``

for East, North, and Up.  Selecting ``away_from_satellite`` reverses all three
coefficients.  Requiring the convention and positive direction as keyword
arguments prevents a silent sign/convention choice.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Hashable, Iterable, Mapping, Sequence

import numpy as np


HEADING_CONVENTION = "satellite_heading_clockwise_from_north"
POSITIVE_DIRECTIONS = frozenset({"toward_satellite", "away_from_satellite"})
QUADRATIC_TERM_NAMES = ("x", "y", "1", "xy", "x2", "y2")


class ValidationCoreError(ValueError):
    """Base class for invalid numerical inputs or unsupported conventions."""


class QuadraticFitError(ValidationCoreError):
    """Raised when a six-term quadratic fit is rank deficient or ill-conditioned."""


@dataclass(frozen=True)
class QuadraticModel:
    """An unweighted six-term quadratic least-squares model.

    Coefficients follow :data:`QUADRATIC_TERM_NAMES`, namely
    ``[x, y, 1, x*y, x**2, y**2]``.
    """

    coefficients: tuple[float, float, float, float, float, float]
    rank: int
    condition_number: float


@dataclass(frozen=True)
class SiteDeduplicationResult:
    """Connected site clusters and one deterministic representative per cluster.

    ``component_by_row[i]`` is the zero-based component number for source row
    ``i``.  Components are ordered by their smallest station identifier.
    """

    components: tuple[tuple[int, ...], ...]
    representative_indices: tuple[int, ...]
    component_by_row: tuple[int, ...]


@dataclass(frozen=True)
class GroupSplitResult:
    """One fitted group split and its station-level residual summaries."""

    split_index: int
    heldout_groups: tuple[Hashable, ...]
    fitting_indices: tuple[int, ...]
    heldout_indices: tuple[int, ...]
    model: QuadraticModel
    all_metrics: Mapping[str, float | int]
    fitting_metrics: Mapping[str, float | int]
    heldout_metrics: Mapping[str, float | int]


@dataclass(frozen=True)
class GroupLOOCVResult:
    """Aggregate station-wise predictions from leave-one-group-out validation."""

    predicted_correction: tuple[float, ...]
    corrected_residual: tuple[float, ...]
    folds: tuple[GroupSplitResult, ...]
    metrics: Mapping[str, float | int]


@dataclass(frozen=True)
class RepeatedGroupHoldoutResult:
    """Deterministic repeated group-holdout results."""

    seed: int
    holdout_fraction: float
    trials: tuple[GroupSplitResult, ...]


def _as_1d_float(name: str, values: Sequence[float] | np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1:
        raise ValidationCoreError(f"{name} must be one-dimensional, got {array.shape}")
    return array


def _require_same_length(**arrays: Sequence[object] | np.ndarray) -> int:
    lengths = {name: len(values) for name, values in arrays.items()}
    if len(set(lengths.values())) != 1:
        raise ValidationCoreError(f"Input lengths differ: {lengths}")
    return next(iter(lengths.values()), 0)


def _validate_angles(incidence_deg: np.ndarray, heading_deg: np.ndarray) -> None:
    if not np.isfinite(incidence_deg).all() or not np.isfinite(heading_deg).all():
        raise ValidationCoreError("Incidence and heading must be finite")
    if np.any((incidence_deg < 0.0) | (incidence_deg > 90.0)):
        raise ValidationCoreError("Incidence must lie in the closed interval [0, 90] degrees")


def los_unit_vector(
    incidence_deg: Sequence[float] | np.ndarray | float,
    heading_deg: Sequence[float] | np.ndarray | float,
    *,
    heading_convention: str,
    positive_direction: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return East/North/Up LOS coefficients under an explicit convention.

    Parameters
    ----------
    incidence_deg
        Radar incidence angle in degrees from the upward vertical.
    heading_deg
        Satellite flight heading in degrees, clockwise from geographic north.
    heading_convention
        Must equal ``"satellite_heading_clockwise_from_north"``.  The explicit
        argument is intentional so callers cannot silently confuse flight
        heading with look azimuth.
    positive_direction
        ``"toward_satellite"`` uses a positive Up coefficient;
        ``"away_from_satellite"`` reverses the complete vector.

    Returns
    -------
    tuple of numpy.ndarray
        Broadcast arrays ``(look_east, look_north, look_up)``.
    """

    if heading_convention != HEADING_CONVENTION:
        raise ValidationCoreError(
            f"Unsupported heading convention {heading_convention!r}; "
            f"expected {HEADING_CONVENTION!r}"
        )
    if positive_direction not in POSITIVE_DIRECTIONS:
        raise ValidationCoreError(
            f"Unsupported positive direction {positive_direction!r}; "
            f"expected one of {sorted(POSITIVE_DIRECTIONS)}"
        )

    incidence, heading = np.broadcast_arrays(
        np.asarray(incidence_deg, dtype=np.float64),
        np.asarray(heading_deg, dtype=np.float64),
    )
    _validate_angles(incidence, heading)
    incidence_rad = np.deg2rad(incidence)
    heading_rad = np.deg2rad(heading)
    look_east = -np.cos(heading_rad) * np.sin(incidence_rad)
    look_north = np.sin(heading_rad) * np.sin(incidence_rad)
    look_up = np.cos(incidence_rad)
    if positive_direction == "away_from_satellite":
        look_east, look_north, look_up = -look_east, -look_north, -look_up
    return look_east, look_north, look_up


def enu_to_los(
    east: Sequence[float] | np.ndarray | float,
    north: Sequence[float] | np.ndarray | float,
    up: Sequence[float] | np.ndarray | float,
    incidence_deg: Sequence[float] | np.ndarray | float,
    heading_deg: Sequence[float] | np.ndarray | float,
    *,
    heading_convention: str,
    positive_direction: str,
) -> np.ndarray:
    """Project ENU velocities onto LOS using broadcast NumPy arithmetic.

    The output has the common broadcast shape of ENU, incidence, and heading.
    Input and output velocity units are identical.
    """

    look_east, look_north, look_up = los_unit_vector(
        incidence_deg,
        heading_deg,
        heading_convention=heading_convention,
        positive_direction=positive_direction,
    )
    east_array, north_array, up_array, look_east, look_north, look_up = (
        np.broadcast_arrays(
            np.asarray(east, dtype=np.float64),
            np.asarray(north, dtype=np.float64),
            np.asarray(up, dtype=np.float64),
            look_east,
            look_north,
            look_up,
        )
    )
    return look_east * east_array + look_north * north_array + look_up * up_array


def project_diagonal_sigma(
    sigma_east: Sequence[float] | np.ndarray | float,
    sigma_north: Sequence[float] | np.ndarray | float,
    sigma_up: Sequence[float] | np.ndarray | float,
    incidence_deg: Sequence[float] | np.ndarray | float,
    heading_deg: Sequence[float] | np.ndarray | float,
    *,
    heading_convention: str,
    positive_direction: str,
) -> np.ndarray:
    """Propagate independent ENU one-sigma values into LOS one-sigma.

    This implements ``sqrt((lE*sE)^2 + (lN*sN)^2 + (lU*sU)^2)``.  It is valid
    only for a diagonal ENU covariance matrix.  Callers with covariance terms
    must use the full ``sqrt(l.T @ covariance @ l)`` expression instead.
    """

    look_east, look_north, look_up = los_unit_vector(
        incidence_deg,
        heading_deg,
        heading_convention=heading_convention,
        positive_direction=positive_direction,
    )
    sigma_east, sigma_north, sigma_up, look_east, look_north, look_up = (
        np.broadcast_arrays(
            np.asarray(sigma_east, dtype=np.float64),
            np.asarray(sigma_north, dtype=np.float64),
            np.asarray(sigma_up, dtype=np.float64),
            look_east,
            look_north,
            look_up,
        )
    )
    if not (
        np.isfinite(sigma_east).all()
        and np.isfinite(sigma_north).all()
        and np.isfinite(sigma_up).all()
    ):
        raise ValidationCoreError("ENU sigmas must be finite")
    if np.any((sigma_east < 0.0) | (sigma_north < 0.0) | (sigma_up < 0.0)):
        raise ValidationCoreError("ENU sigmas must be non-negative")
    return np.sqrt(
        (look_east * sigma_east) ** 2
        + (look_north * sigma_north) ** 2
        + (look_up * sigma_up) ** 2
    )


def summarize_residuals(
    residuals: Sequence[float] | np.ndarray,
    *,
    tolerance_mm_per_year: float = 2.0,
    drop_nonfinite: bool = True,
) -> dict[str, float | int]:
    """Return datum-retaining and centered summaries for supplied residuals.

    ``rms`` is the raw ``sqrt(mean(r**2))`` and therefore retains the residual
    mean.  ``centered_rms`` subtracts the residual mean only for diagnosis; it
    does not alter the raw values.  NMAD is ``1.4826*median(|r-median(r)|)``.
    """

    if not np.isfinite(tolerance_mm_per_year) or tolerance_mm_per_year < 0.0:
        raise ValidationCoreError(
            "tolerance_mm_per_year must be finite and non-negative"
        )
    residual = _as_1d_float("residuals", residuals)
    finite = np.isfinite(residual)
    if not finite.all():
        if not drop_nonfinite:
            raise ValidationCoreError("Residuals contain non-finite values")
        residual = residual[finite]
    if residual.size == 0:
        raise ValidationCoreError("Cannot summarize an empty residual set")

    bias = float(np.mean(residual))
    median = float(np.median(residual))
    return {
        "n": int(residual.size),
        "bias_mm_per_year": bias,
        "rms_mm_per_year": float(np.sqrt(np.mean(residual**2))),
        "centered_rms_mm_per_year": float(
            np.sqrt(np.mean((residual - bias) ** 2))
        ),
        "mae_mm_per_year": float(np.mean(np.abs(residual))),
        "median_residual_mm_per_year": median,
        "nmad_mm_per_year": float(
            1.4826 * np.median(np.abs(residual - median))
        ),
        "fraction_within_tolerance": float(
            np.mean(np.abs(residual) <= tolerance_mm_per_year)
        ),
        "tolerance_mm_per_year": float(tolerance_mm_per_year),
        "minimum_residual_mm_per_year": float(np.min(residual)),
        "maximum_residual_mm_per_year": float(np.max(residual)),
    }


def accuracy_metrics(
    insar: Sequence[float] | np.ndarray,
    gnss: Sequence[float] | np.ndarray,
    *,
    tolerance_mm_per_year: float = 2.0,
    drop_nonfinite: bool = True,
) -> dict[str, float | int]:
    """Compare paired values with residual ``InSAR - GNSS``.

    Pairwise non-finite observations are dropped by default.  Pearson's
    correlation is returned as NaN when fewer than two values remain or either
    series is constant, avoiding warnings and misleading infinities.
    """

    insar_array = _as_1d_float("insar", insar)
    gnss_array = _as_1d_float("gnss", gnss)
    _require_same_length(insar=insar_array, gnss=gnss_array)
    finite = np.isfinite(insar_array) & np.isfinite(gnss_array)
    if not finite.all():
        if not drop_nonfinite:
            raise ValidationCoreError("Comparison contains non-finite pairs")
        insar_array = insar_array[finite]
        gnss_array = gnss_array[finite]
    if insar_array.size == 0:
        raise ValidationCoreError("Cannot compare an empty paired sample")

    result = summarize_residuals(
        insar_array - gnss_array,
        tolerance_mm_per_year=tolerance_mm_per_year,
        drop_nonfinite=False,
    )
    if (
        insar_array.size < 2
        or np.ptp(insar_array) == 0.0
        or np.ptp(gnss_array) == 0.0
    ):
        pearson_r = float("nan")
    else:
        pearson_r = float(np.corrcoef(insar_array, gnss_array)[0, 1])
    result.update(
        {
            "mean_insar": float(np.mean(insar_array)),
            "mean_gnss": float(np.mean(gnss_array)),
            "correlation": pearson_r,
        }
    )
    return result


def comparison_metrics(
    insar: Sequence[float] | np.ndarray,
    gnss: Sequence[float] | np.ndarray,
    *,
    tolerance_mm_per_year: float = 2.0,
    drop_nonfinite: bool = True,
) -> dict[str, float | int]:
    """Descriptive alias for :func:`accuracy_metrics`."""

    return accuracy_metrics(
        insar,
        gnss,
        tolerance_mm_per_year=tolerance_mm_per_year,
        drop_nonfinite=drop_nonfinite,
    )


def haversine_m(
    longitude_1_deg: Sequence[float] | np.ndarray | float,
    latitude_1_deg: Sequence[float] | np.ndarray | float,
    longitude_2_deg: Sequence[float] | np.ndarray | float,
    latitude_2_deg: Sequence[float] | np.ndarray | float,
    *,
    earth_radius_m: float = 6_371_008.8,
) -> np.ndarray | float:
    """Return great-circle distance(s) in metres on a spherical Earth."""

    if not np.isfinite(earth_radius_m) or earth_radius_m <= 0.0:
        raise ValidationCoreError("earth_radius_m must be finite and positive")
    lon1, lat1, lon2, lat2 = np.broadcast_arrays(
        np.asarray(longitude_1_deg, dtype=np.float64),
        np.asarray(latitude_1_deg, dtype=np.float64),
        np.asarray(longitude_2_deg, dtype=np.float64),
        np.asarray(latitude_2_deg, dtype=np.float64),
    )
    if not (
        np.isfinite(lon1).all()
        and np.isfinite(lat1).all()
        and np.isfinite(lon2).all()
        and np.isfinite(lat2).all()
    ):
        raise ValidationCoreError("Coordinates must be finite")
    if np.any(np.abs(lat1) > 90.0) or np.any(np.abs(lat2) > 90.0):
        raise ValidationCoreError("Latitude must lie within [-90, 90] degrees")

    phi1 = np.deg2rad(lat1)
    phi2 = np.deg2rad(lat2)
    delta_phi = phi2 - phi1
    delta_lambda = np.deg2rad(lon2 - lon1)
    term = (
        np.sin(delta_phi / 2.0) ** 2
        + np.cos(phi1) * np.cos(phi2) * np.sin(delta_lambda / 2.0) ** 2
    )
    distance = 2.0 * earth_radius_m * np.arcsin(np.sqrt(np.clip(term, 0.0, 1.0)))
    return float(distance) if distance.ndim == 0 else distance


def deduplicate_connected_sites(
    longitude_deg: Sequence[float] | np.ndarray,
    latitude_deg: Sequence[float] | np.ndarray,
    station_ids: Sequence[str],
    duration_years: Sequence[float] | np.ndarray,
    max_component_sigma: Sequence[float] | np.ndarray,
    *,
    threshold_m: float = 100.0,
) -> SiteDeduplicationResult:
    """Cluster nearby station IDs and select one deterministic representative.

    Rows are joined by single-link connected components: an edge exists when
    the great-circle separation is strictly less than ``threshold_m``.  Thus a
    chain can join endpoints farther apart than the threshold.  Within each
    component the representative is selected by longest duration, then lowest
    maximum component sigma, then lexicographically smallest station ID, then
    source-row index.  This matches a transparent physical-site deduplication
    rule while retaining deterministic ties.
    """

    longitude = _as_1d_float("longitude_deg", longitude_deg)
    latitude = _as_1d_float("latitude_deg", latitude_deg)
    duration = _as_1d_float("duration_years", duration_years)
    sigma = _as_1d_float("max_component_sigma", max_component_sigma)
    station_id_tuple = tuple(str(value) for value in station_ids)
    n_rows = _require_same_length(
        longitude=longitude,
        latitude=latitude,
        station_ids=station_id_tuple,
        duration=duration,
        sigma=sigma,
    )
    if n_rows == 0:
        raise ValidationCoreError("Cannot deduplicate an empty station table")
    if len(set(station_id_tuple)) != n_rows:
        raise ValidationCoreError("station_ids must be unique")
    if not np.isfinite(longitude).all() or not np.isfinite(latitude).all():
        raise ValidationCoreError("Station coordinates must be finite")
    if not np.isfinite(threshold_m) or threshold_m <= 0.0:
        raise ValidationCoreError("threshold_m must be finite and positive")

    parent = list(range(n_rows))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        root_left, root_right = find(left), find(right)
        if root_left != root_right:
            parent[root_right] = root_left

    for left in range(n_rows):
        for right in range(left):
            if (
                haversine_m(
                    longitude[left], latitude[left], longitude[right], latitude[right]
                )
                < threshold_m
            ):
                union(left, right)

    raw_components: dict[int, list[int]] = {}
    for index in range(n_rows):
        raw_components.setdefault(find(index), []).append(index)
    components = sorted(
        (tuple(sorted(indices)) for indices in raw_components.values()),
        key=lambda indices: min(station_id_tuple[index] for index in indices),
    )

    representatives: list[int] = []
    component_by_row = [-1] * n_rows
    for component_number, indices in enumerate(components):
        for index in indices:
            component_by_row[index] = component_number
        representative = min(
            indices,
            key=lambda index: (
                -duration[index] if np.isfinite(duration[index]) else np.inf,
                sigma[index] if np.isfinite(sigma[index]) else np.inf,
                station_id_tuple[index],
                index,
            ),
        )
        representatives.append(representative)

    return SiteDeduplicationResult(
        components=tuple(components),
        representative_indices=tuple(representatives),
        component_by_row=tuple(component_by_row),
    )


def normalize_pixel_coordinates(
    pixel_rows_zero_based: Sequence[float] | np.ndarray,
    pixel_columns_zero_based: Sequence[float] | np.ndarray,
    *,
    raster_height: int,
    raster_width: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Normalize full-raster pixel coordinates to the historical quadratic space.

    The returned ``x`` and ``y`` are equivalent to converting zero-based Python
    indices to one-based pixel centers and applying
    ``((index+1)-(size+1)/2)/(size/2)``.  They span the same polynomial space
    as an unnormalized six-term pixel-coordinate quadratic.
    """

    if raster_height <= 0 or raster_width <= 0:
        raise ValidationCoreError("Raster height and width must be positive")
    rows = _as_1d_float("pixel_rows_zero_based", pixel_rows_zero_based)
    columns = _as_1d_float("pixel_columns_zero_based", pixel_columns_zero_based)
    _require_same_length(rows=rows, columns=columns)
    if not np.isfinite(rows).all() or not np.isfinite(columns).all():
        raise ValidationCoreError("Pixel coordinates must be finite")
    if np.any((rows < 0) | (rows >= raster_height)):
        raise ValidationCoreError("Pixel row lies outside the raster")
    if np.any((columns < 0) | (columns >= raster_width)):
        raise ValidationCoreError("Pixel column lies outside the raster")
    x = ((columns + 1.0) - (raster_width + 1.0) / 2.0) / (raster_width / 2.0)
    y = ((rows + 1.0) - (raster_height + 1.0) / 2.0) / (raster_height / 2.0)
    return x, y


def quadratic_basis(
    x: Sequence[float] | np.ndarray,
    y: Sequence[float] | np.ndarray,
) -> np.ndarray:
    """Return the fixed six-term design ``[x, y, 1, xy, x², y²]``."""

    x_array = _as_1d_float("x", x)
    y_array = _as_1d_float("y", y)
    _require_same_length(x=x_array, y=y_array)
    if not np.isfinite(x_array).all() or not np.isfinite(y_array).all():
        raise ValidationCoreError("Quadratic coordinates must be finite")
    return np.column_stack(
        (
            x_array,
            y_array,
            np.ones_like(x_array),
            x_array * y_array,
            x_array * x_array,
            y_array * y_array,
        )
    )


def quadratic_design_matrix(
    x: Sequence[float] | np.ndarray,
    y: Sequence[float] | np.ndarray,
) -> np.ndarray:
    """Backward-compatible descriptive alias for :func:`quadratic_basis`."""

    return quadratic_basis(x, y)


def fit_quadratic(
    x: Sequence[float] | np.ndarray,
    y: Sequence[float] | np.ndarray,
    response: Sequence[float] | np.ndarray,
    *,
    maximum_condition_number: float = 1.0e12,
    rcond: float | None = None,
) -> QuadraticModel:
    """Fit an unweighted six-term quadratic with rank and condition guards."""

    if not np.isfinite(maximum_condition_number) or maximum_condition_number <= 1.0:
        raise ValidationCoreError("maximum_condition_number must be finite and > 1")
    response_array = _as_1d_float("response", response)
    design = quadratic_basis(x, y)
    _require_same_length(design=design, response=response_array)
    if design.shape[0] < design.shape[1]:
        raise QuadraticFitError(
            f"Quadratic fit requires at least 6 observations, found {design.shape[0]}"
        )
    if not np.isfinite(response_array).all():
        raise ValidationCoreError("Quadratic response must be finite")

    coefficients, _, rank, singular_values = np.linalg.lstsq(
        design, response_array, rcond=rcond
    )
    if int(rank) != 6:
        raise QuadraticFitError(f"Quadratic design rank is {rank}, expected 6")
    if singular_values.size != 6 or singular_values[-1] <= 0.0:
        condition_number = float("inf")
    else:
        condition_number = float(singular_values[0] / singular_values[-1])
    if not np.isfinite(condition_number) or condition_number > maximum_condition_number:
        raise QuadraticFitError(
            "Quadratic design is ill-conditioned: "
            f"condition={condition_number:.6g}, limit={maximum_condition_number:.6g}"
        )
    coefficient_tuple = tuple(float(value) for value in coefficients)
    return QuadraticModel(
        coefficients=coefficient_tuple,  # type: ignore[arg-type]
        rank=int(rank),
        condition_number=condition_number,
    )


def predict_quadratic(
    model: QuadraticModel,
    x: Sequence[float] | np.ndarray,
    y: Sequence[float] | np.ndarray,
) -> np.ndarray:
    """Evaluate a :class:`QuadraticModel` at supplied coordinates."""

    return quadratic_basis(x, y) @ np.asarray(
        model.coefficients, dtype=np.float64
    )


def _normalized_groups(groups: Sequence[Hashable] | np.ndarray) -> tuple[Hashable, ...]:
    normalized: list[Hashable] = []
    for group in groups:
        value = group.item() if isinstance(group, np.generic) else group
        try:
            hash(value)
        except TypeError as error:
            raise ValidationCoreError(f"Group label is not hashable: {value!r}") from error
        normalized.append(value)
    return tuple(normalized)


def _ordered_unique(values: Iterable[Hashable]) -> tuple[Hashable, ...]:
    return tuple(sorted(set(values), key=lambda value: (type(value).__name__, repr(value))))


def _prepare_group_inputs(
    x: Sequence[float] | np.ndarray,
    y: Sequence[float] | np.ndarray,
    response: Sequence[float] | np.ndarray,
    groups: Sequence[Hashable] | np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[Hashable, ...], tuple[Hashable, ...]]:
    x_array = _as_1d_float("x", x)
    y_array = _as_1d_float("y", y)
    response_array = _as_1d_float("response", response)
    group_tuple = _normalized_groups(groups)
    n_rows = _require_same_length(
        x=x_array, y=y_array, response=response_array, groups=group_tuple
    )
    if n_rows == 0:
        raise ValidationCoreError("Grouped validation requires observations")
    if not (
        np.isfinite(x_array).all()
        and np.isfinite(y_array).all()
        and np.isfinite(response_array).all()
    ):
        raise ValidationCoreError("Grouped validation inputs must be finite")
    unique_groups = _ordered_unique(group_tuple)
    if len(unique_groups) < 2:
        raise ValidationCoreError("Grouped validation requires at least two groups")
    return x_array, y_array, response_array, group_tuple, unique_groups


def _split_result(
    split_index: int,
    heldout_groups: tuple[Hashable, ...],
    fitting_indices: np.ndarray,
    heldout_indices: np.ndarray,
    model: QuadraticModel,
    corrected_residual_all: np.ndarray,
) -> GroupSplitResult:
    return GroupSplitResult(
        split_index=split_index,
        heldout_groups=heldout_groups,
        fitting_indices=tuple(int(value) for value in fitting_indices),
        heldout_indices=tuple(int(value) for value in heldout_indices),
        model=model,
        all_metrics=summarize_residuals(corrected_residual_all),
        fitting_metrics=summarize_residuals(corrected_residual_all[fitting_indices]),
        heldout_metrics=summarize_residuals(corrected_residual_all[heldout_indices]),
    )


def grouped_loocv(
    x: Sequence[float] | np.ndarray,
    y: Sequence[float] | np.ndarray,
    response: Sequence[float] | np.ndarray,
    groups: Sequence[Hashable] | np.ndarray,
    *,
    maximum_condition_number: float = 1.0e12,
) -> GroupLOOCVResult:
    """Leave each group out once and aggregate all held-out residuals.

    ``response`` is normally the pre-correction ``InSAR - GNSS`` residual.
    Each fold fits the quadratic only on other groups, predicts a correction
    for the held-out observations, and stores ``response - prediction``.  Rows
    sharing a group are always withheld together.  The final RMS is calculated
    once over all station-level held-out residuals, not averaged across folds.
    """

    x_array, y_array, response_array, group_tuple, unique_groups = _prepare_group_inputs(
        x, y, response, groups
    )
    predicted = np.full(response_array.shape, np.nan, dtype=np.float64)
    corrected = np.full(response_array.shape, np.nan, dtype=np.float64)
    folds: list[GroupSplitResult] = []

    for split_index, heldout_group in enumerate(unique_groups):
        heldout_mask = np.fromiter(
            (value == heldout_group for value in group_tuple),
            dtype=bool,
            count=len(group_tuple),
        )
        fitting_indices = np.flatnonzero(~heldout_mask)
        heldout_indices = np.flatnonzero(heldout_mask)
        model = fit_quadratic(
            x_array[fitting_indices],
            y_array[fitting_indices],
            response_array[fitting_indices],
            maximum_condition_number=maximum_condition_number,
        )
        prediction_all = predict_quadratic(model, x_array, y_array)
        residual_all = response_array - prediction_all
        predicted[heldout_indices] = prediction_all[heldout_indices]
        corrected[heldout_indices] = residual_all[heldout_indices]
        folds.append(
            _split_result(
                split_index,
                (heldout_group,),
                fitting_indices,
                heldout_indices,
                model,
                residual_all,
            )
        )

    if not np.isfinite(predicted).all() or not np.isfinite(corrected).all():
        raise RuntimeError("Internal error: not every observation received one LOOCV prediction")
    return GroupLOOCVResult(
        predicted_correction=tuple(float(value) for value in predicted),
        corrected_residual=tuple(float(value) for value in corrected),
        folds=tuple(folds),
        metrics=summarize_residuals(corrected),
    )


def repeated_group_holdout(
    x: Sequence[float] | np.ndarray,
    y: Sequence[float] | np.ndarray,
    response: Sequence[float] | np.ndarray,
    groups: Sequence[Hashable] | np.ndarray,
    *,
    trials: int = 100,
    holdout_fraction: float = 0.20,
    seed: int = 0,
    maximum_condition_number: float = 1.0e12,
) -> RepeatedGroupHoldoutResult:
    """Run seeded Monte Carlo group holdout without station leakage.

    For each trial, ``round(holdout_fraction*n_groups)`` unique groups are
    sampled without replacement (with a minimum of one).  The quadratic is
    fitted to all observations in the remaining groups.  Metrics are returned
    for all, fitting, and held-out observations; only the held-out metrics are
    fully out of sample.
    """

    if trials < 1:
        raise ValidationCoreError("trials must be at least 1")
    if not np.isfinite(holdout_fraction) or not 0.0 < holdout_fraction < 1.0:
        raise ValidationCoreError("holdout_fraction must lie strictly between 0 and 1")
    x_array, y_array, response_array, group_tuple, unique_groups = _prepare_group_inputs(
        x, y, response, groups
    )
    n_holdout = max(1, int(round(holdout_fraction * len(unique_groups))))
    if n_holdout >= len(unique_groups):
        raise ValidationCoreError("Holdout selection leaves no fitting group")

    rng = np.random.default_rng(seed)
    trial_results: list[GroupSplitResult] = []
    for trial_index in range(trials):
        selected_positions = np.sort(
            rng.choice(len(unique_groups), size=n_holdout, replace=False)
        )
        heldout_groups = tuple(unique_groups[int(index)] for index in selected_positions)
        heldout_set = set(heldout_groups)
        heldout_mask = np.fromiter(
            (value in heldout_set for value in group_tuple),
            dtype=bool,
            count=len(group_tuple),
        )
        fitting_indices = np.flatnonzero(~heldout_mask)
        heldout_indices = np.flatnonzero(heldout_mask)
        try:
            model = fit_quadratic(
                x_array[fitting_indices],
                y_array[fitting_indices],
                response_array[fitting_indices],
                maximum_condition_number=maximum_condition_number,
            )
        except QuadraticFitError as error:
            raise QuadraticFitError(
                f"Trial {trial_index + 1} failed for held-out groups "
                f"{heldout_groups!r}: {error}"
            ) from error
        residual_all = response_array - predict_quadratic(model, x_array, y_array)
        trial_results.append(
            _split_result(
                trial_index,
                heldout_groups,
                fitting_indices,
                heldout_indices,
                model,
                residual_all,
            )
        )

    return RepeatedGroupHoldoutResult(
        seed=int(seed),
        holdout_fraction=float(holdout_fraction),
        trials=tuple(trial_results),
    )


__all__ = [
    "HEADING_CONVENTION",
    "POSITIVE_DIRECTIONS",
    "QUADRATIC_TERM_NAMES",
    "GroupLOOCVResult",
    "GroupSplitResult",
    "QuadraticFitError",
    "QuadraticModel",
    "RepeatedGroupHoldoutResult",
    "SiteDeduplicationResult",
    "ValidationCoreError",
    "accuracy_metrics",
    "comparison_metrics",
    "deduplicate_connected_sites",
    "enu_to_los",
    "fit_quadratic",
    "grouped_loocv",
    "haversine_m",
    "los_unit_vector",
    "normalize_pixel_coordinates",
    "predict_quadratic",
    "project_diagonal_sigma",
    "quadratic_basis",
    "quadratic_design_matrix",
    "repeated_group_holdout",
    "summarize_residuals",
]

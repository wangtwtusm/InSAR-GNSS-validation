#!/usr/bin/env python3
"""Portable command-line validation of InSAR LOS/VLM against GNSS rates."""

from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import sys
import tempfile
from pathlib import Path
from typing import Any, Iterable

import numpy as np

try:
    import pandas as pd
except ImportError as exc:  # pragma: no cover - dependency error path
    raise SystemExit(
        "pandas is required by validate.py; use a compatible scientific Python "
        "environment. The skill does not install packages automatically."
    ) from exc

from preflight_inputs import audit_coregistration, load_metadata, sha256_file
from validation_core import (
    HEADING_CONVENTION,
    accuracy_metrics,
    deduplicate_connected_sites,
    enu_to_los,
    fit_quadratic,
    grouped_loocv,
    predict_quadratic,
    project_diagonal_sigma,
    repeated_group_holdout,
    summarize_residuals,
)


BASE_COLUMNS = {
    "station_id",
    "longitude_deg",
    "latitude_deg",
    "east_mm_per_year",
    "north_mm_per_year",
    "up_mm_per_year",
}
SIGMA_COLUMNS = {
    "sigma_east_mm_per_year",
    "sigma_north_mm_per_year",
    "sigma_up_mm_per_year",
}
CV_COLUMNS = {
    "orbit",
    "station_id",
    "group_id",
    "normalized_pixel_x",
    "normalized_pixel_y",
    "residual_before_fit_mm_per_year",
}
GENERATED_NAMES = {
    "station_values.csv",
    "metrics.csv",
    "exclusions.csv",
    "provenance.json",
    "full_fit_metrics.csv",
    "loocv_station_values.csv",
    "loocv_metrics.csv",
    "loocv_folds.csv",
    "cv_trials.csv",
    "cv_assignments.csv",
    "cv_distribution.csv",
}


def require_columns(frame: pd.DataFrame, required: set[str], label: str) -> None:
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing columns: {', '.join(missing)}")


def require_unique_station_ids(frame: pd.DataFrame) -> None:
    if frame["station_id"].astype(str).duplicated().any():
        duplicates = sorted(
            frame.loc[
                frame["station_id"].astype(str).duplicated(keep=False), "station_id"
            ].astype(str).unique()
        )
        raise ValueError(f"station_id values must be unique: {duplicates[:10]}")


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as stream:
            stream.write(text)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def atomic_write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        frame.to_csv(temporary, index=False, lineterminator="\n")
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def prepare_output(directory: Path, overwrite: bool) -> None:
    if directory.exists() and not directory.is_dir():
        raise NotADirectoryError(directory)
    if directory.exists() and not overwrite:
        conflicts = sorted(name for name in GENERATED_NAMES if (directory / name).exists())
        if conflicts:
            raise FileExistsError(
                f"Refusing to replace existing generated files in {directory}: {conflicts}"
            )
    directory.mkdir(parents=True, exist_ok=True)


def parse_required_status(values: Iterable[str]) -> list[tuple[str, str]]:
    parsed: list[tuple[str, str]] = []
    for value in values:
        if "=" not in value:
            raise ValueError(f"Use COLUMN=VALUE for --require-status, got {value!r}")
        column, expected = value.split("=", 1)
        if not column or not expected:
            raise ValueError(f"Use COLUMN=VALUE for --require-status, got {value!r}")
        parsed.append((column, expected))
    return parsed


def add_reason(reasons: list[list[str]], mask: np.ndarray, reason: str) -> None:
    for index in np.flatnonzero(mask):
        reasons[int(index)].append(reason)


def station_selection(frame: pd.DataFrame, args: argparse.Namespace) -> tuple[np.ndarray, list[str]]:
    """Apply declared non-raster rules and return include mask plus reasons."""

    n_rows = len(frame)
    reasons: list[list[str]] = [[] for _ in range(n_rows)]
    coordinates = frame[["longitude_deg", "latitude_deg"]].apply(
        pd.to_numeric, errors="coerce"
    ).to_numpy(dtype=float)
    velocities = frame[
        ["east_mm_per_year", "north_mm_per_year", "up_mm_per_year"]
    ].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    add_reason(reasons, ~np.isfinite(coordinates).all(axis=1), "nonfinite_coordinate")
    add_reason(reasons, ~np.isfinite(velocities).all(axis=1), "nonfinite_enu_rate")

    station_ids = frame["station_id"].astype(str).to_numpy()
    excluded_ids = set(str(value) for value in args.exclude_station_id)
    if excluded_ids:
        add_reason(
            reasons,
            np.asarray([value in excluded_ids for value in station_ids]),
            "explicit_station_exclusion",
        )

    for column, expected in parse_required_status(args.require_status):
        if column not in frame.columns:
            raise ValueError(f"--require-status refers to missing column {column!r}")
        add_reason(
            reasons,
            frame[column].astype(str).to_numpy() != expected,
            f"status_{column}_not_{expected}",
        )

    if args.min_duration_years is not None:
        if "duration_years" not in frame.columns:
            raise ValueError("--min-duration-years requires duration_years")
        duration = pd.to_numeric(frame["duration_years"], errors="coerce").to_numpy()
        add_reason(
            reasons,
            ~np.isfinite(duration) | (duration < args.min_duration_years),
            "record_too_short",
        )

    if args.min_distance_from_fitting_m is not None:
        if "nearest_fitting_distance_m" not in frame.columns:
            raise ValueError(
                "--min-distance-from-fitting-m requires nearest_fitting_distance_m"
            )
        distance = pd.to_numeric(
            frame["nearest_fitting_distance_m"], errors="coerce"
        ).to_numpy()
        add_reason(
            reasons,
            ~np.isfinite(distance) | (distance <= args.min_distance_from_fitting_m),
            "too_close_to_fitting_site",
        )

    has_sigmas = SIGMA_COLUMNS.issubset(frame.columns)
    if args.max_component_sigma is not None and not has_sigmas:
        raise ValueError(
            "--max-component-sigma requires sigma_east/north/up_mm_per_year"
        )
    max_sigma = np.full(n_rows, np.nan)
    if has_sigmas:
        sigma_matrix = frame[
            [
                "sigma_east_mm_per_year",
                "sigma_north_mm_per_year",
                "sigma_up_mm_per_year",
            ]
        ].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        max_sigma = np.max(sigma_matrix, axis=1)
        add_reason(
            reasons,
            ~np.isfinite(sigma_matrix).all(axis=1) | (sigma_matrix < 0).any(axis=1),
            "invalid_rate_uncertainty",
        )

    # Physical-site deduplication precedes the formal sigma threshold. The
    # representative is chosen by record length, then sigma, then ID.
    if args.deduplicate_within_m is not None:
        if "duration_years" not in frame.columns or not has_sigmas:
            raise ValueError(
                "--deduplicate-within-m requires duration_years and all ENU sigmas"
            )
        eligible_before_dedup = np.asarray([not value for value in reasons], dtype=bool)
        candidate_indices = np.flatnonzero(eligible_before_dedup)
        if candidate_indices.size:
            result = deduplicate_connected_sites(
                coordinates[candidate_indices, 0],
                coordinates[candidate_indices, 1],
                station_ids[candidate_indices],
                pd.to_numeric(frame["duration_years"], errors="coerce")
                .to_numpy(dtype=float)[candidate_indices],
                max_sigma[candidate_indices],
                threshold_m=args.deduplicate_within_m,
            )
            representatives = {
                int(candidate_indices[index]) for index in result.representative_indices
            }
            duplicate_mask = np.zeros(n_rows, dtype=bool)
            for source_index in candidate_indices:
                if int(source_index) not in representatives:
                    duplicate_mask[int(source_index)] = True
            add_reason(reasons, duplicate_mask, "nonrepresentative_physical_site_record")

    if args.max_component_sigma is not None:
        add_reason(
            reasons,
            ~np.isfinite(max_sigma) | (max_sigma > args.max_component_sigma),
            "component_rate_uncertainty_above_threshold",
        )

    include = np.asarray([not value for value in reasons], dtype=bool)
    return include, [";".join(value) if value else "eligible" for value in reasons]


def import_rasterio() -> Any:
    try:
        import rasterio
    except ImportError as exc:  # pragma: no cover - dependency error path
        raise RuntimeError(
            "rasterio/GDAL is required for strict pixel sampling. Use a compatible "
            "scientific Python environment; the skill does not install packages."
        ) from exc
    return rasterio


def crs_matches(dataset_crs: Any, declared_crs: str) -> bool:
    rasterio = import_rasterio()
    try:
        return dataset_crs == rasterio.crs.CRS.from_user_input(declared_crs)
    except Exception as exc:
        raise ValueError(f"Invalid coordinate_crs {declared_crs!r}") from exc


def sample_dataset(
    dataset: Any, longitude: np.ndarray, latitude: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[str]]:
    """Sample containing pixels without nearest-neighbour fallback."""

    n_rows = len(longitude)
    values = np.full(n_rows, np.nan)
    rows = np.full(n_rows, -1, dtype=int)
    columns = np.full(n_rows, -1, dtype=int)
    center_x = np.full(n_rows, np.nan)
    center_y = np.full(n_rows, np.nan)
    status = ["out_of_bounds"] * n_rows
    inside_indices: list[int] = []
    coordinates: list[tuple[float, float]] = []
    for index, (x_value, y_value) in enumerate(zip(longitude, latitude, strict=True)):
        if not np.isfinite(x_value) or not np.isfinite(y_value):
            status[index] = "invalid_coordinate"
            continue
        if not (
            dataset.bounds.left <= x_value < dataset.bounds.right
            and dataset.bounds.bottom < y_value <= dataset.bounds.top
        ):
            continue
        row, column = dataset.index(float(x_value), float(y_value))
        if not (0 <= row < dataset.height and 0 <= column < dataset.width):
            continue
        rows[index], columns[index] = int(row), int(column)
        pixel_x, pixel_y = dataset.xy(row, column)
        center_x[index], center_y[index] = float(pixel_x), float(pixel_y)
        inside_indices.append(index)
        coordinates.append((float(x_value), float(y_value)))
        status[index] = "no_data"

    for index, sample in zip(
        inside_indices,
        dataset.sample(coordinates, indexes=1, masked=True),
        strict=True,
    ):
        value = sample[0]
        if not np.ma.is_masked(value) and np.isfinite(float(value)):
            values[index] = float(value)
            status[index] = "finite"
    return values, rows, columns, center_x, center_y, status


def metadata_and_sigma_guard(frame: pd.DataFrame, metadata: dict[str, Any]) -> bool:
    has_any_sigma = bool(SIGMA_COLUMNS & set(frame.columns))
    has_sigmas = SIGMA_COLUMNS.issubset(frame.columns)
    if has_any_sigma and not has_sigmas:
        raise ValueError(
            "Rate-uncertainty columns must provide all three ENU components"
        )
    if has_any_sigma and metadata.get("uncertainty_semantics") != (
        "one_standard_deviation_rate_uncertainty"
    ):
        raise ValueError(
            "ENU sigma columns may be used only when metadata declares "
            "uncertainty_semantics='one_standard_deviation_rate_uncertainty'"
        )
    return has_sigmas


def metrics_frame(
    metrics: dict[str, float | int], *, label: str, role: str
) -> pd.DataFrame:
    return pd.DataFrame([{"selection": label, "comparison_role": role, **metrics}])


def disposition_table(frame: pd.DataFrame) -> pd.DataFrame:
    excluded = frame.loc[~frame["include_in_metrics"].astype(bool)].copy()
    columns = [
        name
        for name in (
            "station_id",
            "longitude_deg",
            "latitude_deg",
            "selection_reason",
            "pixel_status",
        )
        if name in excluded.columns
    ]
    return excluded[columns]


def input_manifest(paths: dict[str, Path]) -> dict[str, Any]:
    return {
        name: {"file": path.name, "sha256": sha256_file(path)}
        for name, path in paths.items()
    }


def provenance(
    *,
    mode: str,
    metadata: dict[str, Any],
    inputs: dict[str, Path],
    args: argparse.Namespace,
    counts: dict[str, int],
) -> dict[str, Any]:
    rasterio = import_rasterio() if mode in {"los", "vlm"} else None
    return {
        "mode": mode,
        "metadata": metadata,
        "inputs": input_manifest(inputs),
        "selection": {
            "exclude_station_ids": sorted(set(args.exclude_station_id))
            if hasattr(args, "exclude_station_id")
            else [],
            "require_status": list(args.require_status)
            if hasattr(args, "require_status")
            else [],
            "min_duration_years": getattr(args, "min_duration_years", None),
            "min_distance_from_fitting_m": getattr(
                args, "min_distance_from_fitting_m", None
            ),
            "deduplicate_within_m": getattr(args, "deduplicate_within_m", None),
            "max_component_sigma": getattr(args, "max_component_sigma", None),
            "sampling": "station_containing_pixel",
            "nearest_valid_fallback": "none",
            "constant_offset_removed": False,
        },
        "counts": counts,
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "rasterio": rasterio.__version__ if rasterio is not None else None,
        },
    }


def save_comparison(
    frame: pd.DataFrame,
    metrics: dict[str, float | int],
    output_dir: Path,
    provenance_record: dict[str, Any],
    *,
    label: str,
    role: str,
) -> None:
    atomic_write_csv(frame, output_dir / "station_values.csv")
    atomic_write_csv(metrics_frame(metrics, label=label, role=role), output_dir / "metrics.csv")
    atomic_write_csv(disposition_table(frame), output_dir / "exclusions.csv")
    atomic_write_text(
        output_dir / "provenance.json",
        json.dumps(provenance_record, indent=2, sort_keys=True) + "\n",
    )


def validate_los(args: argparse.Namespace) -> None:
    metadata = load_metadata(args.metadata, "los")
    frame = pd.read_csv(args.stations)
    require_columns(frame, BASE_COLUMNS, args.stations.name)
    require_unique_station_ids(frame)
    has_sigmas = metadata_and_sigma_guard(frame, metadata)
    include, reasons = station_selection(frame, args)

    audit_coregistration([args.los, args.incidence, args.heading])
    rasterio = import_rasterio()
    longitude = pd.to_numeric(frame["longitude_deg"], errors="coerce").to_numpy()
    latitude = pd.to_numeric(frame["latitude_deg"], errors="coerce").to_numpy()
    with (
        rasterio.open(args.los) as los_dataset,
        rasterio.open(args.incidence) as incidence_dataset,
        rasterio.open(args.heading) as heading_dataset,
    ):
        if not crs_matches(los_dataset.crs, metadata["coordinate_crs"]):
            raise ValueError(
                f"Station coordinate CRS {metadata['coordinate_crs']} does not match "
                f"raster CRS {los_dataset.crs}; reproject explicitly before validation"
            )
        insar_los, rows, columns, center_x, center_y, status = sample_dataset(
            los_dataset, longitude, latitude
        )
        incidence, *_ = sample_dataset(incidence_dataset, longitude, latitude)
        heading, *_ = sample_dataset(heading_dataset, longitude, latitude)

    geometry_finite = np.isfinite(incidence) & np.isfinite(heading)
    for index in np.flatnonzero(~geometry_finite & (np.asarray(status) == "finite")):
        status[int(index)] = "geometry_no_data"
    exact_finite = np.isfinite(insar_los) & geometry_finite

    gnss_los = np.full(len(frame), np.nan)
    valid_enu = frame[
        ["east_mm_per_year", "north_mm_per_year", "up_mm_per_year"]
    ].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    if exact_finite.any():
        gnss_los[exact_finite] = enu_to_los(
            valid_enu[exact_finite, 0],
            valid_enu[exact_finite, 1],
            valid_enu[exact_finite, 2],
            incidence[exact_finite],
            heading[exact_finite],
            heading_convention=HEADING_CONVENTION,
            positive_direction=metadata["los_sign"],
        )
    sigma_los = np.full(len(frame), np.nan)
    if has_sigmas and exact_finite.any():
        sigma = frame[
            [
                "sigma_east_mm_per_year",
                "sigma_north_mm_per_year",
                "sigma_up_mm_per_year",
            ]
        ].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        finite_sigma = exact_finite & np.isfinite(sigma).all(axis=1) & (sigma >= 0).all(axis=1)
        sigma_los[finite_sigma] = project_diagonal_sigma(
            sigma[finite_sigma, 0],
            sigma[finite_sigma, 1],
            sigma[finite_sigma, 2],
            incidence[finite_sigma],
            heading[finite_sigma],
            heading_convention=HEADING_CONVENTION,
            positive_direction=metadata["los_sign"],
        )

    raster_missing = ~exact_finite
    final_reasons: list[str] = []
    for base_reason, missing, pixel_status in zip(reasons, raster_missing, status, strict=True):
        values = [] if base_reason == "eligible" else base_reason.split(";")
        if missing:
            values.append(f"raster_{pixel_status}")
        final_reasons.append(";".join(values) if values else "included")
    include = include & exact_finite & np.isfinite(gnss_los)

    result = frame.copy()
    result["pixel_row_zero_based"] = rows
    result["pixel_column_zero_based"] = columns
    result["pixel_center_x"] = center_x
    result["pixel_center_y"] = center_y
    result["pixel_status"] = status
    result["incidence_deg"] = incidence
    result["heading_deg"] = heading
    result["insar_los_mm_per_year"] = insar_los
    result["gnss_los_mm_per_year"] = gnss_los
    result["sigma_los_mm_per_year"] = sigma_los
    result["residual_mm_per_year"] = insar_los - gnss_los
    result["selection_reason"] = final_reasons
    result["include_in_metrics"] = include
    result["group_id"] = [
        f"r{row}_c{column}" if row >= 0 and column >= 0 else ""
        for row, column in zip(rows, columns, strict=True)
    ]

    formal = result.loc[result["include_in_metrics"]]
    metrics = accuracy_metrics(
        formal["insar_los_mm_per_year"],
        formal["gnss_los_mm_per_year"],
        tolerance_mm_per_year=args.tolerance,
        drop_nonfinite=False,
    )
    record = provenance(
        mode="los",
        metadata=metadata,
        inputs={
            "stations": args.stations,
            "metadata": args.metadata,
            "los": args.los,
            "incidence": args.incidence,
            "heading": args.heading,
        },
        args=args,
        counts={"input_rows": len(result), "formal_rows": len(formal)},
    )
    save_comparison(
        result,
        metrics,
        args.output_dir,
        record,
        label=args.label,
        role=metadata["comparison_role"],
    )


def validate_vlm(args: argparse.Namespace) -> None:
    metadata = load_metadata(args.metadata, "vlm")
    frame = pd.read_csv(args.stations)
    require_columns(frame, BASE_COLUMNS, args.stations.name)
    require_unique_station_ids(frame)
    metadata_and_sigma_guard(frame, metadata)
    include, reasons = station_selection(frame, args)

    rasterio = import_rasterio()
    longitude = pd.to_numeric(frame["longitude_deg"], errors="coerce").to_numpy()
    latitude = pd.to_numeric(frame["latitude_deg"], errors="coerce").to_numpy()
    with rasterio.open(args.vlm) as dataset:
        if dataset.count != 1:
            raise ValueError("VLM raster must contain one band")
        if not crs_matches(dataset.crs, metadata["coordinate_crs"]):
            raise ValueError(
                f"Station coordinate CRS {metadata['coordinate_crs']} does not match "
                f"raster CRS {dataset.crs}; reproject explicitly before validation"
            )
        insar_vlm, rows, columns, center_x, center_y, status = sample_dataset(
            dataset, longitude, latitude
        )
    gnss_up = pd.to_numeric(frame["up_mm_per_year"], errors="coerce").to_numpy()
    exact_finite = np.isfinite(insar_vlm) & np.isfinite(gnss_up)
    final_reasons: list[str] = []
    for base_reason, finite, pixel_status in zip(reasons, exact_finite, status, strict=True):
        values = [] if base_reason == "eligible" else base_reason.split(";")
        if not finite:
            values.append(f"raster_{pixel_status}")
        final_reasons.append(";".join(values) if values else "included")
    include = include & exact_finite

    result = frame.copy()
    result["pixel_row_zero_based"] = rows
    result["pixel_column_zero_based"] = columns
    result["pixel_center_x"] = center_x
    result["pixel_center_y"] = center_y
    result["pixel_status"] = status
    result["insar_vlm_mm_per_year"] = insar_vlm
    result["gnss_up_mm_per_year"] = gnss_up
    result["sigma_up_mm_per_year"] = (
        pd.to_numeric(frame["sigma_up_mm_per_year"], errors="coerce").to_numpy()
        if "sigma_up_mm_per_year" in frame.columns
        else np.nan
    )
    result["residual_mm_per_year"] = insar_vlm - gnss_up
    result["selection_reason"] = final_reasons
    result["include_in_metrics"] = include
    result["group_id"] = [
        f"r{row}_c{column}" if row >= 0 and column >= 0 else ""
        for row, column in zip(rows, columns, strict=True)
    ]

    formal = result.loc[result["include_in_metrics"]]
    metrics = accuracy_metrics(
        formal["insar_vlm_mm_per_year"],
        formal["gnss_up_mm_per_year"],
        tolerance_mm_per_year=args.tolerance,
        drop_nonfinite=False,
    )
    record = provenance(
        mode="vlm",
        metadata=metadata,
        inputs={
            "stations": args.stations,
            "metadata": args.metadata,
            "vlm": args.vlm,
        },
        args=args,
        counts={"input_rows": len(result), "formal_rows": len(formal)},
    )
    save_comparison(
        result,
        metrics,
        args.output_dir,
        record,
        label=args.label,
        role=metadata["comparison_role"],
    )


def prefixed(prefix: str, values: dict[str, float | int]) -> dict[str, float | int]:
    return {f"{prefix}_{key}": value for key, value in values.items()}


def distribution_rows(trials: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for orbit, subset in trials.groupby("orbit", sort=True):
        for scope in ("all", "fitting", "holdout"):
            values = subset[f"{scope}_rms_mm_per_year"].to_numpy(dtype=float)
            rows.append(
                {
                    "orbit": orbit,
                    "scope": scope,
                    "trials": len(values),
                    "mean_rms_mm_per_year": np.mean(values),
                    "standard_deviation_of_rms_mm_per_year": np.std(values, ddof=1)
                    if len(values) > 1
                    else np.nan,
                    "p2_5_rms_mm_per_year": np.percentile(values, 2.5),
                    "median_rms_mm_per_year": np.median(values),
                    "p97_5_rms_mm_per_year": np.percentile(values, 97.5),
                }
            )
    return pd.DataFrame(rows)


def cross_validate(args: argparse.Namespace) -> None:
    metadata = load_metadata(args.metadata, "cross-validate")
    observations = pd.read_csv(args.observations)
    require_columns(observations, CV_COLUMNS, args.observations.name)
    numeric_columns = [
        "normalized_pixel_x",
        "normalized_pixel_y",
        "residual_before_fit_mm_per_year",
    ]
    numeric = observations[numeric_columns].apply(pd.to_numeric, errors="coerce")
    if not np.isfinite(numeric.to_numpy()).all():
        raise ValueError("Cross-validation coordinates and responses must be finite")

    loocv_station_rows: list[dict[str, Any]] = []
    loocv_metric_rows: list[dict[str, Any]] = []
    fold_rows: list[dict[str, Any]] = []
    full_fit_rows: list[dict[str, Any]] = []
    trial_rows: list[dict[str, Any]] = []
    assignment_rows: list[dict[str, Any]] = []
    orbit_names = sorted(observations["orbit"].astype(str).unique())
    seed_sequences = np.random.SeedSequence(args.seed).spawn(len(orbit_names))

    for orbit, seed_sequence in zip(orbit_names, seed_sequences, strict=True):
        subset = observations.loc[observations["orbit"].astype(str) == orbit].copy()
        x = pd.to_numeric(subset["normalized_pixel_x"]).to_numpy(dtype=float)
        y = pd.to_numeric(subset["normalized_pixel_y"]).to_numpy(dtype=float)
        response = pd.to_numeric(
            subset["residual_before_fit_mm_per_year"]
        ).to_numpy(dtype=float)
        groups = subset["group_id"].astype(str).to_numpy()

        full_model = fit_quadratic(
            x,
            y,
            response,
            maximum_condition_number=args.maximum_condition_number,
        )
        full_residual = response - predict_quadratic(full_model, x, y)
        full_fit_rows.append(
            {
                "orbit": orbit,
                "n_groups": len(set(groups)),
                "condition_number": full_model.condition_number,
                **summarize_residuals(full_residual, tolerance_mm_per_year=args.tolerance),
            }
        )

        loocv = grouped_loocv(
            x,
            y,
            response,
            groups,
            maximum_condition_number=args.maximum_condition_number,
        )
        loocv_summary = summarize_residuals(
            loocv.corrected_residual,
            tolerance_mm_per_year=args.tolerance,
        )
        loocv_metric_rows.append(
            {
                "orbit": orbit,
                "n_groups": len(set(groups)),
                **loocv_summary,
            }
        )
        for source, prediction, residual in zip(
            subset.to_dict("records"),
            loocv.predicted_correction,
            loocv.corrected_residual,
            strict=True,
        ):
            loocv_station_rows.append(
                {
                    **source,
                    "predicted_correction_mm_per_year": prediction,
                    "residual_mm_per_year": residual,
                }
            )
        for fold in loocv.folds:
            fold_residual_all = response - predict_quadratic(fold.model, x, y)
            fold_fitting_metrics = summarize_residuals(
                fold_residual_all[list(fold.fitting_indices)],
                tolerance_mm_per_year=args.tolerance,
            )
            fold_heldout_metrics = summarize_residuals(
                fold_residual_all[list(fold.heldout_indices)],
                tolerance_mm_per_year=args.tolerance,
            )
            fold_rows.append(
                {
                    "orbit": orbit,
                    "fold": fold.split_index + 1,
                    "heldout_groups": ";".join(map(str, fold.heldout_groups)),
                    "n_fitting_observations": len(fold.fitting_indices),
                    "n_heldout_observations": len(fold.heldout_indices),
                    "condition_number": fold.model.condition_number,
                    **prefixed("fitting", fold_fitting_metrics),
                    **prefixed("heldout", fold_heldout_metrics),
                }
            )

        orbit_seed = int(seed_sequence.generate_state(1, dtype=np.uint64)[0])
        repeated = repeated_group_holdout(
            x,
            y,
            response,
            groups,
            trials=args.trials,
            holdout_fraction=args.holdout_fraction,
            seed=orbit_seed,
            maximum_condition_number=args.maximum_condition_number,
        )
        for trial in repeated.trials:
            trial_residual_all = response - predict_quadratic(trial.model, x, y)
            trial_all_metrics = summarize_residuals(
                trial_residual_all,
                tolerance_mm_per_year=args.tolerance,
            )
            trial_fitting_metrics = summarize_residuals(
                trial_residual_all[list(trial.fitting_indices)],
                tolerance_mm_per_year=args.tolerance,
            )
            trial_heldout_metrics = summarize_residuals(
                trial_residual_all[list(trial.heldout_indices)],
                tolerance_mm_per_year=args.tolerance,
            )
            trial_rows.append(
                {
                    "trial": trial.split_index + 1,
                    "orbit": orbit,
                    "group": orbit,
                    "heldout_groups": ";".join(map(str, trial.heldout_groups)),
                    "n_fitting_observations": len(trial.fitting_indices),
                    "n_heldout_observations": len(trial.heldout_indices),
                    "condition_number": trial.model.condition_number,
                    **prefixed("all", trial_all_metrics),
                    **prefixed("fitting", trial_fitting_metrics),
                    **prefixed("holdout", trial_heldout_metrics),
                }
            )
            for group in trial.heldout_groups:
                assignment_rows.append(
                    {
                        "trial": trial.split_index + 1,
                        "orbit": orbit,
                        "group_id": group,
                    }
                )

    trials = pd.DataFrame(trial_rows)
    atomic_write_csv(pd.DataFrame(full_fit_rows), args.output_dir / "full_fit_metrics.csv")
    atomic_write_csv(pd.DataFrame(loocv_station_rows), args.output_dir / "loocv_station_values.csv")
    atomic_write_csv(pd.DataFrame(loocv_metric_rows), args.output_dir / "loocv_metrics.csv")
    atomic_write_csv(pd.DataFrame(fold_rows), args.output_dir / "loocv_folds.csv")
    atomic_write_csv(trials, args.output_dir / "cv_trials.csv")
    atomic_write_csv(pd.DataFrame(assignment_rows), args.output_dir / "cv_assignments.csv")
    atomic_write_csv(distribution_rows(trials), args.output_dir / "cv_distribution.csv")
    record = {
        "mode": "cross-validate",
        "metadata": metadata,
        "inputs": input_manifest(
            {"observations": args.observations, "metadata": args.metadata}
        ),
        "model": {
            "basis_order": ["x", "y", "1", "xy", "x2", "y2"],
            "split_unit": "group_id",
            "fit": "unweighted_least_squares",
            "loocv_summary": "pooled_heldout_station_observations",
            "trials": args.trials,
            "holdout_fraction": args.holdout_fraction,
            "seed": args.seed,
            "maximum_condition_number": args.maximum_condition_number,
            "constant_offset_removed_after_prediction": False,
        },
        "counts": {
            "input_rows": len(observations),
            "orbits": len(orbit_names),
            "unique_groups": int(observations[["orbit", "group_id"]].drop_duplicates().shape[0]),
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
    }
    atomic_write_text(
        args.output_dir / "provenance.json",
        json.dumps(record, indent=2, sort_keys=True) + "\n",
    )


def add_selection_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--exclude-station-id", action="append", default=[])
    parser.add_argument(
        "--require-status",
        action="append",
        default=[],
        metavar="COLUMN=VALUE",
    )
    parser.add_argument("--min-duration-years", type=float)
    parser.add_argument("--min-distance-from-fitting-m", type=float)
    parser.add_argument("--deduplicate-within-m", type=float)
    parser.add_argument("--max-component-sigma", type=float)
    parser.add_argument("--tolerance", type=float, default=2.0)
    parser.add_argument("--label", default="formal")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate InSAR LOS/VLM products against GNSS rates."
    )
    subparsers = parser.add_subparsers(dest="mode", required=True)

    los = subparsers.add_parser("los", help="Project ENU to LOS and compare exact pixels")
    los.add_argument("--stations", type=Path, required=True)
    los.add_argument("--metadata", type=Path, required=True)
    los.add_argument("--los", type=Path, required=True)
    los.add_argument("--incidence", type=Path, required=True)
    los.add_argument("--heading", type=Path, required=True)
    add_selection_arguments(los)

    vlm = subparsers.add_parser("vlm", help="Compare VLM raster with GNSS Up")
    vlm.add_argument("--stations", type=Path, required=True)
    vlm.add_argument("--metadata", type=Path, required=True)
    vlm.add_argument("--vlm", type=Path, required=True)
    add_selection_arguments(vlm)

    cv = subparsers.add_parser(
        "cross-validate", help="Run grouped quadratic LOOCV and repeated holdout"
    )
    cv.add_argument("--observations", type=Path, required=True)
    cv.add_argument("--metadata", type=Path, required=True)
    cv.add_argument("--trials", type=int, default=100)
    cv.add_argument("--holdout-fraction", type=float, default=0.20)
    cv.add_argument("--seed", type=int, default=0)
    cv.add_argument("--maximum-condition-number", type=float, default=1.0e12)
    cv.add_argument("--tolerance", type=float, default=2.0)
    cv.add_argument("--output-dir", type=Path, required=True)
    cv.add_argument("--overwrite", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    for path_name in (
        "stations",
        "metadata",
        "los",
        "incidence",
        "heading",
        "vlm",
        "observations",
    ):
        path = getattr(args, path_name, None)
        if path is not None and not path.is_file():
            raise FileNotFoundError(path)
    if hasattr(args, "tolerance") and args.tolerance < 0:
        raise ValueError("--tolerance must be non-negative")
    for option_name in (
        "min_duration_years",
        "min_distance_from_fitting_m",
        "max_component_sigma",
    ):
        option_value = getattr(args, option_name, None)
        if option_value is not None and option_value < 0:
            option_flag = "--" + option_name.replace("_", "-")
            raise ValueError(f"{option_flag} must be non-negative")
    deduplication_distance = getattr(args, "deduplicate_within_m", None)
    if deduplication_distance is not None and deduplication_distance <= 0:
        raise ValueError("--deduplicate-within-m must be positive")
    if args.mode == "cross-validate":
        if args.trials < 1:
            raise ValueError("--trials must be at least 1")
        if not 0 < args.holdout_fraction < 1:
            raise ValueError("--holdout-fraction must lie between 0 and 1")
    prepare_output(args.output_dir, args.overwrite)
    if args.mode == "los":
        validate_los(args)
    elif args.mode == "vlm":
        validate_vlm(args)
    else:
        cross_validate(args)
    print(args.output_dir)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc

#!/usr/bin/env python3
"""Optional publication-quality PyGMT plots for InSAR–GNSS validation.

The command is deliberately presentation-only.  It reads standardized outputs
from the validation workflow, verifies their internal consistency, and never
fits or removes an offset.  All paths are supplied explicitly; no study area,
station, orbit, or product is built into this module.

Run ``plot_validation.py --help`` for the available subcommands.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import os
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence


VLM_COLUMNS = {
    "station_id",
    "longitude_deg",
    "latitude_deg",
    "insar_vlm_mm_per_year",
    "gnss_up_mm_per_year",
    "residual_mm_per_year",
    "sigma_up_mm_per_year",
    "include_in_metrics",
}
LOS_COLUMNS = {
    "station_id",
    "longitude_deg",
    "latitude_deg",
    "insar_los_mm_per_year",
    "gnss_los_mm_per_year",
    "residual_mm_per_year",
    "sigma_los_mm_per_year",
    "include_in_metrics",
}
METRIC_COLUMNS = {
    "n",
    "bias_mm_per_year",
    "rms_mm_per_year",
    "centered_rms_mm_per_year",
    "mae_mm_per_year",
    "correlation",
}
TRIAL_COLUMNS = {
    "trial",
    "all_rms_mm_per_year",
    "holdout_rms_mm_per_year",
}

OCEAN_COLOR = "#DCEFF7"
INCLUDED_PEN = "0.55p,black"
EXCLUDED_HALO_PEN = "1.0p,white"
ERROR_PEN = "0.55p,gray30"


class ValidationError(RuntimeError):
    """Raised when a standardized input violates the plotting contract."""


class DependencyError(RuntimeError):
    """Raised when the selected Python environment cannot run PyGMT."""


@dataclass(frozen=True)
class Dependencies:
    np: Any
    pygmt: Any
    rasterio: Any | None
    resampling: Any | None
    xarray: Any | None


@dataclass(frozen=True)
class StationValues:
    station_id: list[str]
    longitude: Any
    latitude: Any
    insar: Any
    gnss: Any
    residual: Any
    sigma: Any
    include: Any

    @property
    def count(self) -> int:
        return len(self.station_id)


@dataclass(frozen=True)
class RasterInfo:
    path: Path
    region: tuple[float, float, float, float]
    finite_values: Any
    crs: str
    width: int
    height: int


def default_gmt_workdir() -> Path:
    return Path(tempfile.gettempdir()) / "insar-gnss-validation-gmt"


def configure_gmt(workdir: Path) -> dict[str, str]:
    """Set writable GMT locations before importing PyGMT."""
    root = workdir.expanduser().resolve()
    userdir = root / "user"
    tmpdir = root / "tmp"
    userdir.mkdir(parents=True, exist_ok=True)
    tmpdir.mkdir(parents=True, exist_ok=True)
    os.environ["GMT_USERDIR"] = str(userdir)
    os.environ["GMT_TMPDIR"] = str(tmpdir)
    os.environ.setdefault("GMT_SESSION_NAME", f"insar_gnss_{os.getpid()}")
    return {"GMT_USERDIR": str(userdir), "GMT_TMPDIR": str(tmpdir)}


def load_dependencies(*, require_raster: bool, gmt_workdir: Path) -> Dependencies:
    locations = configure_gmt(gmt_workdir)
    try:
        import numpy as np
    except Exception as exc:  # pragma: no cover - depends on host environment
        raise DependencyError(
            "NumPy is required. Select a Python environment containing numpy."
        ) from exc
    try:
        import pygmt
    except Exception as exc:  # pragma: no cover - depends on host environment
        raise DependencyError(
            "PyGMT could not be imported. Use a Python environment with PyGMT "
            "and a compatible GMT 6 shared library. Writable GMT directories "
            f"were configured at {locations['GMT_USERDIR']!r} and "
            f"{locations['GMT_TMPDIR']!r}. Original error: {exc}"
        ) from exc

    rasterio = None
    resampling = None
    xarray = None
    if require_raster:
        try:
            import rasterio
            import xarray
            from rasterio.enums import Resampling
        except Exception as exc:  # pragma: no cover - depends on host environment
            raise DependencyError(
                "Raster plotting requires rasterio and xarray in the same Python "
                "environment as PyGMT. No package installation was attempted."
            ) from exc
        resampling = Resampling
    return Dependencies(
        np=np,
        pygmt=pygmt,
        rasterio=rasterio,
        resampling=resampling,
        xarray=xarray,
    )


def read_csv_rows(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    if not path.is_file():
        raise ValidationError(f"CSV file does not exist: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValidationError(f"CSV has no header: {path}")
        fields = [field.strip() for field in reader.fieldnames]
        if len(fields) != len(set(fields)):
            raise ValidationError(f"CSV has duplicate column names: {path}")
        rows = []
        for raw in reader:
            rows.append({str(key).strip(): (value or "").strip() for key, value in raw.items()})
    if not rows:
        raise ValidationError(f"CSV contains no data rows: {path}")
    return rows, fields


def require_columns(path: Path, fields: Iterable[str], required: set[str]) -> None:
    missing = required - set(fields)
    if missing:
        raise ValidationError(
            f"{path} is missing standardized column(s): {', '.join(sorted(missing))}"
        )


def parse_float(value: str, *, context: str, allow_nan: bool = False) -> float:
    if allow_nan and not value.strip():
        return math.nan
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"Expected a numeric value for {context}, found {value!r}") from exc
    if allow_nan and math.isnan(result):
        return result
    if not math.isfinite(result):
        raise ValidationError(f"Expected a finite value for {context}, found {value!r}")
    return result


def parse_bool(value: str, *, context: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"1", "true", "yes", "y", "include", "included"}:
        return True
    if normalized in {"0", "false", "no", "n", "exclude", "excluded"}:
        return False
    raise ValidationError(
        f"Expected a Boolean for {context}; use true/false or 1/0, found {value!r}"
    )


def load_station_values(path: Path, *, mode: str, np: Any) -> StationValues:
    rows, fields = read_csv_rows(path)
    if mode == "vlm":
        required = VLM_COLUMNS
        insar_column = "insar_vlm_mm_per_year"
        gnss_column = "gnss_up_mm_per_year"
        sigma_column = "sigma_up_mm_per_year"
    elif mode == "los":
        required = LOS_COLUMNS
        insar_column = "insar_los_mm_per_year"
        gnss_column = "gnss_los_mm_per_year"
        sigma_column = "sigma_los_mm_per_year"
    else:  # defensive; callers constrain this
        raise ValidationError(f"Unsupported station-value mode: {mode}")
    require_columns(path, fields, required)

    station_ids = [row["station_id"] for row in rows]
    if any(not station_id for station_id in station_ids):
        raise ValidationError(f"{path}: station_id cannot be empty")
    duplicates = sorted(value for value, count in Counter(station_ids).items() if count > 1)
    if duplicates:
        raise ValidationError(
            f"{path}: station_id must be unique; duplicate(s): {', '.join(duplicates[:8])}"
        )

    def values(column: str, *, allow_nan: bool = False) -> Any:
        return np.asarray(
            [
                parse_float(
                    row[column],
                    context=f"{path.name} row {index + 2} column {column}",
                    allow_nan=allow_nan,
                )
                for index, row in enumerate(rows)
            ],
            dtype=float,
        )

    # The validation workflow deliberately retains excluded diagnostic rows.
    # Such rows may have no exact raster sample; only included rows must be
    # complete and finite.
    longitude = values("longitude_deg", allow_nan=True)
    latitude = values("latitude_deg", allow_nan=True)
    insar = values(insar_column, allow_nan=True)
    gnss = values(gnss_column, allow_nan=True)
    residual = values("residual_mm_per_year", allow_nan=True)
    sigma = values(sigma_column, allow_nan=True)
    include = np.asarray(
        [
            parse_bool(
                row["include_in_metrics"],
                context=f"{path.name} row {index + 2} include_in_metrics",
            )
            for index, row in enumerate(rows)
        ],
        dtype=bool,
    )
    if np.any(np.isfinite(sigma) & (sigma < 0)):
        raise ValidationError(f"{path}: uncertainty values cannot be negative")
    if not np.any(include):
        raise ValidationError(f"{path}: no station is included in metrics")

    complete = (
        np.isfinite(longitude)
        & np.isfinite(latitude)
        & np.isfinite(insar)
        & np.isfinite(gnss)
        & np.isfinite(residual)
    )
    if np.any(include & ~complete):
        index = int(np.flatnonzero(include & ~complete)[0])
        raise ValidationError(
            f"{path}: included station {station_ids[index]} has a non-finite "
            "coordinate or comparison value"
        )

    expected_residual = insar - gnss
    comparable = np.isfinite(residual) & np.isfinite(expected_residual)
    if not np.allclose(
        residual[comparable],
        expected_residual[comparable],
        rtol=1e-7,
        atol=1e-7,
    ):
        errors = np.where(comparable, np.abs(residual - expected_residual), -np.inf)
        index = int(np.argmax(errors))
        raise ValidationError(
            f"{path}: residual must be raw InSAR minus GNSS. Largest mismatch is "
            f"{errors[index]:.6g} mm/yr at station {station_ids[index]}; no "
            "offset-removed residual may be supplied as residual_mm_per_year."
        )
    return StationValues(
        station_id=station_ids,
        longitude=longitude,
        latitude=latitude,
        insar=insar,
        gnss=gnss,
        residual=residual,
        sigma=sigma,
        include=include,
    )


def parse_filters(items: Sequence[str]) -> dict[str, str]:
    filters: dict[str, str] = {}
    for item in items:
        if "=" not in item:
            raise ValidationError(
                f"Metric filter must have COLUMN=VALUE form, found {item!r}"
            )
        column, value = item.split("=", 1)
        column = column.strip()
        if not column:
            raise ValidationError(f"Metric filter has an empty column: {item!r}")
        filters[column] = value.strip()
    return filters


def load_metric_row(path: Path, filters: Sequence[str]) -> dict[str, float | int | str]:
    rows, fields = read_csv_rows(path)
    require_columns(path, fields, METRIC_COLUMNS)
    parsed_filters = parse_filters(filters)
    for column in parsed_filters:
        if column not in fields:
            raise ValidationError(f"Metric filter column {column!r} is absent from {path}")
    selected = [
        row
        for row in rows
        if all(row[column] == value for column, value in parsed_filters.items())
    ]
    if len(selected) != 1:
        detail = "no filter was supplied" if not parsed_filters else f"filters={parsed_filters}"
        raise ValidationError(
            f"Expected exactly one metrics row in {path} after {detail}; found {len(selected)}"
        )
    row = selected[0]
    result: dict[str, float | int | str] = dict(row)
    n_float = parse_float(row["n"], context=f"{path.name} n")
    if not n_float.is_integer() or n_float < 1:
        raise ValidationError(f"{path}: n must be a positive integer")
    result["n"] = int(n_float)
    for column in METRIC_COLUMNS - {"n"}:
        result[column] = parse_float(
            row[column], context=f"{path.name} {column}", allow_nan=(column == "correlation")
        )
    return result


def compute_metrics(stations: StationValues, np: Any) -> dict[str, float | int]:
    mask = stations.include
    residual = stations.residual[mask]
    insar = stations.insar[mask]
    gnss = stations.gnss[mask]
    bias = float(np.mean(residual))
    if len(residual) >= 2 and np.std(insar) > 0 and np.std(gnss) > 0:
        correlation = float(np.corrcoef(insar, gnss)[0, 1])
    else:
        correlation = math.nan
    return {
        "n": int(len(residual)),
        "bias_mm_per_year": bias,
        "rms_mm_per_year": float(np.sqrt(np.mean(residual**2))),
        "centered_rms_mm_per_year": float(np.sqrt(np.mean((residual - bias) ** 2))),
        "mae_mm_per_year": float(np.mean(np.abs(residual))),
        "correlation": correlation,
    }


def verify_metrics(
    stations: StationValues,
    supplied: dict[str, float | int | str],
    np: Any,
    *,
    tolerance: float = 1e-6,
) -> dict[str, float | int]:
    computed = compute_metrics(stations, np)
    if int(supplied["n"]) != computed["n"]:
        raise ValidationError(
            f"Metrics n={supplied['n']} but station_values contains "
            f"{computed['n']} included finite observations"
        )
    for column in (
        "bias_mm_per_year",
        "rms_mm_per_year",
        "centered_rms_mm_per_year",
        "mae_mm_per_year",
        "correlation",
    ):
        expected = float(computed[column])
        actual = float(supplied[column])
        if math.isnan(expected) and math.isnan(actual):
            continue
        if not math.isclose(actual, expected, rel_tol=tolerance, abs_tol=tolerance):
            raise ValidationError(
                f"Metrics mismatch for {column}: CSV={actual:.9g}, recomputed={expected:.9g}. "
                "The plots require raw, non-offset-removed metrics."
            )
    rms = float(supplied["rms_mm_per_year"])
    bias = float(supplied["bias_mm_per_year"])
    centered = float(supplied["centered_rms_mm_per_year"])
    if not math.isclose(rms * rms, bias * bias + centered * centered, rel_tol=2e-6, abs_tol=2e-6):
        raise ValidationError(
            "Metrics violate RMS^2 = bias^2 + centered_RMS^2; check whether an offset "
            "was removed or the statistics use inconsistent station sets."
        )
    return computed


def inspect_raster(path: Path, deps: Dependencies) -> RasterInfo:
    if deps.rasterio is None or deps.resampling is None:
        raise DependencyError("Raster inspection was requested without rasterio")
    if not path.is_file():
        raise ValidationError(f"Raster does not exist: {path}")
    np = deps.np
    with deps.rasterio.open(path) as dataset:
        if dataset.count != 1:
            raise ValidationError(f"Expected a single-band velocity raster: {path}")
        if dataset.crs is None:
            raise ValidationError(f"Raster has no CRS: {path}")
        if not dataset.crs.is_geographic:
            raise ValidationError(
                f"Map plotting currently requires a geographic lon/lat raster; found {dataset.crs}"
            )
        target_rows = min(dataset.height, 1024)
        target_cols = min(dataset.width, 1024)
        sampled = dataset.read(
            1,
            out_shape=(target_rows, target_cols),
            masked=True,
            resampling=deps.resampling.bilinear,
        )
        finite = np.asarray(sampled.compressed(), dtype=float)
        finite = finite[np.isfinite(finite)]
        if finite.size == 0:
            raise ValidationError(f"Raster contains no finite, unmasked values: {path}")
        bounds = dataset.bounds
        return RasterInfo(
            path=path.resolve(),
            region=(float(bounds.left), float(bounds.right), float(bounds.bottom), float(bounds.top)),
            finite_values=finite,
            crs=str(dataset.crs),
            width=int(dataset.width),
            height=int(dataset.height),
        )


def verify_raster_samples(
    path: Path,
    stations: StationValues,
    deps: Dependencies,
    *,
    tolerance: float,
) -> None:
    if tolerance < 0:
        raise ValidationError("--sample-tolerance must be non-negative")
    coordinate_finite = deps.np.isfinite(stations.longitude) & deps.np.isfinite(
        stations.latitude
    )
    valid_indices = deps.np.flatnonzero(coordinate_finite)
    coordinates = [
        (float(stations.longitude[index]), float(stations.latitude[index]))
        for index in valid_indices
    ]
    samples = [math.nan] * stations.count
    with deps.rasterio.open(path) as dataset:
        for index, value in zip(
            valid_indices,
            dataset.sample(coordinates, indexes=1, masked=True),
            strict=True,
        ):
            element = value[0]
            if deps.np.ma.is_masked(element):
                samples[int(index)] = math.nan
            else:
                samples[int(index)] = float(element)
    samples_array = deps.np.asarray(samples, dtype=float)
    required = stations.include
    if deps.np.any(required & ~deps.np.isfinite(samples_array)):
        indices = deps.np.flatnonzero(required & ~deps.np.isfinite(samples_array))
        names = ", ".join(stations.station_id[int(index)] for index in indices[:8])
        raise ValidationError(
            f"Included station(s) do not fall on a finite raster pixel: {names}. "
            "No nearest-pixel fallback is performed."
        )
    mismatch = required & deps.np.isfinite(samples_array) & (
        deps.np.abs(samples_array - stations.insar) > tolerance
    )
    if deps.np.any(mismatch):
        index = int(deps.np.flatnonzero(mismatch)[0])
        raise ValidationError(
            f"Station {stations.station_id[index]} records InSAR={stations.insar[index]:.9g} "
            f"but the exact raster pixel is {samples_array[index]:.9g} mm/yr; difference "
            f"exceeds --sample-tolerance={tolerance}."
        )


def read_raster_grid(path: Path, region: Sequence[float], deps: Dependencies) -> Any:
    """Read a geographic crop as an ascending-coordinate xarray grid.

    Passing a hand-written GeoTIFF directly to GMT can be fragile on some GDAL
    builds.  A DataArray also preserves the geographic coordinate contract
    explicitly and makes the plotted crop deterministic.
    """
    if deps.rasterio is None or deps.xarray is None:
        raise DependencyError("Raster plotting requires rasterio and xarray")
    west, east, south, north = map(float, region)
    rasterio = deps.rasterio
    np = deps.np
    with rasterio.open(path) as dataset:
        transform = dataset.transform
        if not math.isclose(transform.b, 0.0, abs_tol=1e-12) or not math.isclose(
            transform.d, 0.0, abs_tol=1e-12
        ):
            raise ValidationError("Rotated/sheared rasters are not supported by the map renderer")
        window = rasterio.windows.from_bounds(
            west,
            south,
            east,
            north,
            transform=transform,
        ).round_offsets().round_lengths()
        try:
            window = window.intersection(
                rasterio.windows.Window(0, 0, dataset.width, dataset.height)
            )
        except Exception as exc:
            raise ValidationError("Requested map region does not overlap the raster") from exc
        values = dataset.read(1, window=window, masked=True).filled(np.nan)
        crop_transform = dataset.window_transform(window)
    if values.size == 0 or not np.isfinite(values).any():
        raise ValidationError("Requested map region contains no finite raster pixels")
    x = crop_transform.c + crop_transform.a * (np.arange(values.shape[1]) + 0.5)
    y = crop_transform.f + crop_transform.e * (np.arange(values.shape[0]) + 0.5)
    if x[0] > x[-1]:
        x = x[::-1]
        values = values[:, ::-1]
    if y[0] > y[-1]:
        y = y[::-1]
        values = values[::-1, :]
    grid = deps.xarray.DataArray(
        values,
        coords={"lat": y, "lon": x},
        dims=("lat", "lon"),
        name="insar_velocity",
        attrs={"units": "mm yr-1"},
    )
    # These coordinates are pixel centers derived from the rasterio affine
    # transform.  Make the registration and geographic type explicit so GMT
    # does not guess gridline registration or shift the image by half a cell.
    grid.gmt.gtype = 1
    grid.gmt.registration = 1
    return grid


def parse_region(values: Sequence[float] | None, fallback: Sequence[float]) -> list[float]:
    region = list(fallback if values is None else values)
    if len(region) != 4 or not all(math.isfinite(float(value)) for value in region):
        raise ValidationError("Region must be four finite values: WEST EAST SOUTH NORTH")
    west, east, south, north = map(float, region)
    if not west < east or not south < north:
        raise ValidationError("Region must satisfy WEST < EAST and SOUTH < NORTH")
    return [west, east, south, north]


def robust_symmetric_limit(values: Any, np: Any, *, quantile: float, explicit: float | None) -> float:
    if explicit is not None:
        if not math.isfinite(explicit) or explicit <= 0:
            raise ValidationError("Color limit must be a positive finite number")
        return float(explicit)
    if not 50 <= quantile <= 100:
        raise ValidationError("--color-quantile must be between 50 and 100")
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        raise ValidationError("Cannot derive a color scale from an empty finite array")
    limit = float(np.percentile(np.abs(finite), quantile))
    if not math.isfinite(limit) or limit <= 0:
        limit = float(np.max(np.abs(finite)))
    if not math.isfinite(limit) or limit <= 0:
        limit = 1.0
    exponent = math.floor(math.log10(limit))
    scale = 10.0**exponent
    return math.ceil(limit / scale * 5.0) / 5.0 * scale


def make_cpt(deps: Dependencies, *, path: Path, limit: float) -> None:
    step = limit / 100.0
    deps.pygmt.makecpt(
        cmap="vik",
        series=[-limit, limit, step],
        continuous=True,
        output=str(path),
    )


def map_header(figure: Any, text: str) -> None:
    figure.text(
        position="TC",
        text=text,
        font="11p,Helvetica-Bold,black",
        justify="BC",
        offset="j0c/0.25c",
        no_clip=True,
    )


def plot_coast_after_raster(figure: Any, args: argparse.Namespace) -> None:
    try:
        coast_options = {
            "shorelines": "0.65p,black",
            "resolution": args.coast_resolution,
            "area_thresh": args.coast_area_threshold,
        }
        if args.water_fill.lower() != "none":
            coast_options["water"] = args.water_fill
        figure.coast(**coast_options)
    except Exception as exc:
        raise DependencyError(
            "GMT could not draw the requested coastline. Confirm that the GSHHG "
            "coastline database is installed, or choose a coarser "
            f"--coast-resolution. Original error: {exc}"
        ) from exc


def uncertainty_sizes(stations: StationValues, scale: float, np: Any) -> tuple[Any, Any]:
    if not math.isfinite(scale) or scale <= 0:
        raise ValidationError("--uncertainty-radius-scale must be positive")
    finite = np.isfinite(stations.sigma)
    # GMT circle size is diameter.  The requested scientific encoding is radius.
    diameters = 2.0 * stations.sigma[finite] * scale
    return finite, diameters


def auto_legend_values(sigma: Any, np: Any) -> list[float]:
    finite = np.asarray(sigma, dtype=float)
    finite = finite[np.isfinite(finite) & (finite > 0)]
    if finite.size == 0:
        return []
    candidates = np.percentile(finite, [25, 50, 75])
    rounded: list[float] = []
    for value in candidates:
        exponent = math.floor(math.log10(float(value))) if value > 0 else 0
        unit = 10.0 ** (exponent - 1)
        item = round(float(value) / unit) * unit
        if item > 0 and not any(math.isclose(item, existing) for existing in rounded):
            rounded.append(item)
    return rounded


def uncertainty_legend(
    figure: Any,
    values: Sequence[float],
    *,
    scale: float,
    position: str = "jTR+o0.12c",
) -> None:
    if not values:
        return
    lines = ["H 8.5p Helvetica-Bold GNSS rate uncertainty"]
    for value in values:
        diameter = 2.0 * value * scale
        lines.append(
            f"S 0.12c c {diameter:.4f}c - 0.65p,black 0.25c {value:g} mm/yr"
        )
    figure.legend(
        spec=io.StringIO("\n".join(lines) + "\n"),
        position=position,
        box="+gwhite@10+p0.45p,gray35",
    )


def plot_station_squares(
    figure: Any,
    stations: StationValues,
    values: Any,
    *,
    cpt: Path,
    np: Any,
) -> None:
    finite = (
        np.isfinite(stations.longitude)
        & np.isfinite(stations.latitude)
        & np.isfinite(values)
    )
    excluded = ~stations.include & finite
    if np.any(excluded):
        figure.plot(
            x=stations.longitude[excluded],
            y=stations.latitude[excluded],
            style="s0.30c",
            fill="white",
            pen=EXCLUDED_HALO_PEN,
        )
    figure.plot(
        x=stations.longitude[finite],
        y=stations.latitude[finite],
        style="s0.22c",
        fill=values[finite],
        cmap=str(cpt),
        pen=INCLUDED_PEN,
    )


def plot_uncertainty_circles(
    figure: Any,
    stations: StationValues,
    *,
    scale: float,
    np: Any,
) -> None:
    finite, diameters = uncertainty_sizes(stations, scale, np)
    finite = finite & np.isfinite(stations.longitude) & np.isfinite(stations.latitude)
    # Recompute after applying the coordinate mask so size and coordinate arrays align.
    diameters = 2.0 * stations.sigma[finite] * scale
    if np.any(finite):
        figure.plot(
            x=stations.longitude[finite],
            y=stations.latitude[finite],
            size=diameters,
            style="c",
            pen="0.65p,black",
        )


def annotate_raw_metrics(
    figure: Any,
    metrics: dict[str, float | int],
    region: Sequence[float],
    *,
    include_centered: bool = True,
) -> None:
    west, east, south, north = map(float, region)
    x = west + 0.025 * (east - west)
    y0 = north - 0.045 * (north - south)
    spacing = 0.047 * (north - south)
    lines = [
        f"N = {metrics['n']}",
        f"Bias = {float(metrics['bias_mm_per_year']):+.2f} mm/yr",
        f"Raw RMS = {float(metrics['rms_mm_per_year']):.2f} mm/yr",
    ]
    if include_centered:
        lines.append(
            f"cRMS (diagnostic) = {float(metrics['centered_rms_mm_per_year']):.2f} mm/yr"
        )
    correlation = float(metrics["correlation"])
    if math.isfinite(correlation):
        lines.append(f"r = {correlation:.2f}")
    for index, line in enumerate(lines):
        figure.text(
            x=x,
            y=y0 - spacing * index,
            text=line,
            font="7.4p,Helvetica,black",
            justify="TL",
            fill="white@12",
            clearance="1.2p/0.8p",
        )


def ensure_output(path: Path, *, overwrite: bool) -> Path:
    suffix = path.suffix.lower()
    if suffix not in {".png", ".pdf", ".svg", ".eps", ".tif", ".tiff"}:
        raise ValidationError(
            "Output suffix must be .png, .pdf, .svg, .eps, .tif, or .tiff"
        )
    output = path.expanduser().resolve()
    if output.exists() and not overwrite:
        raise ValidationError(f"Output already exists; use --overwrite to replace it: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


def save_figure(figure: Any, output: Path, *, dpi: int) -> None:
    if dpi < 72:
        raise ValidationError("--dpi must be at least 72")
    try:
        figure.savefig(str(output), dpi=dpi, crop=True)
    except Exception as exc:
        raise DependencyError(
            "PyGMT/GMT failed while rendering the figure. Confirm that GMT and "
            f"Ghostscript are available and the output directory is writable. Original error: {exc}"
        ) from exc


def nice_tick_increment(span: float, *, target_intervals: int) -> float:
    raw = span / max(1, target_intervals)
    exponent = math.floor(math.log10(raw))
    scale = 10.0**exponent
    for multiplier in (1.0, 2.0, 2.5, 5.0, 10.0):
        candidate = multiplier * scale
        if candidate >= raw * (1.0 - 1e-12):
            return candidate
    return 10.0 * scale


def common_map_frame(region: Sequence[float], *, target_intervals: int) -> list[str]:
    west, east, south, north = map(float, region)
    x_annot = nice_tick_increment(east - west, target_intervals=target_intervals)
    y_annot = nice_tick_increment(north - south, target_intervals=target_intervals)
    return [
        "WSne",
        f"xa{x_annot:g}f{x_annot / 2.0:g}",
        f"ya{y_annot:g}f{y_annot / 2.0:g}",
    ]


def prepare_map_inputs(
    args: argparse.Namespace,
    *,
    mode: str,
    deps: Dependencies,
) -> tuple[StationValues, dict[str, float | int], RasterInfo, list[float]]:
    stations = load_station_values(args.station_values, mode=mode, np=deps.np)
    metric_row = load_metric_row(args.metrics, args.metric_filter)
    computed = verify_metrics(stations, metric_row, deps.np)
    raster = inspect_raster(args.raster, deps)
    verify_raster_samples(
        args.raster,
        stations,
        deps,
        tolerance=args.sample_tolerance,
    )
    region = parse_region(args.region, raster.region)
    return stations, computed, raster, region


def plot_vlm(args: argparse.Namespace, deps: Dependencies) -> dict[str, Any]:
    np = deps.np
    pygmt = deps.pygmt
    stations, metrics, raster, region = prepare_map_inputs(args, mode="vlm", deps=deps)
    output = ensure_output(args.output, overwrite=args.overwrite)

    vlm_limit = robust_symmetric_limit(
        np.concatenate((raster.finite_values, stations.gnss)),
        np,
        quantile=args.color_quantile,
        explicit=args.vlm_limit,
    )
    residual_limit = robust_symmetric_limit(
        stations.residual,
        np,
        quantile=args.color_quantile,
        explicit=args.residual_limit,
    )
    legend_values = auto_legend_values(stations.sigma, np)

    cpt_dir = Path(os.environ["GMT_TMPDIR"])
    vlm_cpt = cpt_dir / f"vlm_{os.getpid()}.cpt"
    residual_cpt = cpt_dir / f"vlm_residual_{os.getpid()}.cpt"
    make_cpt(deps, path=vlm_cpt, limit=vlm_limit)
    make_cpt(deps, path=residual_cpt, limit=residual_limit)

    included = stations.include
    residual_included = stations.residual[included]
    bin_width = args.hist_bin_width or histogram_bin_width(residual_included, np)
    if not math.isfinite(bin_width) or bin_width <= 0:
        raise ValidationError("--hist-bin-width must be a positive finite number")
    histogram_limit = max(
        residual_limit,
        float(np.max(np.abs(residual_included))) * 1.05,
    )
    hist_edges = np.arange(
        -histogram_limit,
        histogram_limit + bin_width * 1.001,
        bin_width,
    )
    if hist_edges.size < 3:
        hist_edges = np.linspace(-histogram_limit, histogram_limit, 8)
        bin_width = float(hist_edges[1] - hist_edges[0])
    counts, _ = np.histogram(residual_included, bins=hist_edges)
    hist_ymax = max(1.0, float(np.max(counts)) * 1.20)
    scatter_region = scatter_axis_region(stations, np, explicit=args.scatter_range)
    grid = read_raster_grid(raster.path, region, deps)

    figure = pygmt.Figure()
    with pygmt.config(
        MAP_FRAME_TYPE="plain",
        MAP_FRAME_PEN="1p,black",
        FONT_ANNOT_PRIMARY="8p,Helvetica,black",
        FONT_LABEL="9p,Helvetica,black",
        MAP_TICK_LENGTH_PRIMARY="3p",
        MAP_GRID_PEN_PRIMARY="0.3p,gray85",
        FORMAT_GEO_MAP="ddd.xxF",
    ):
        with figure.subplot(
            nrows=2,
            ncols=2,
            figsize=("19.0c", "16.4c"),
            margins=["0.8c", "1.8c"],
        ):
            with figure.set_panel(panel=[0, 0]):
                figure.basemap(
                    region=region,
                    projection="M?",
                    frame=common_map_frame(region, target_intervals=3),
                )
                # Painter order is an invariant: velocity raster, then ocean/coast, then stations.
                figure.grdimage(
                    grid=grid,
                    cmap=str(vlm_cpt),
                    nan_transparent=True,
                    interpolation="n",
                )
                plot_coast_after_raster(figure, args)
                plot_station_squares(
                    figure,
                    stations,
                    stations.gnss,
                    cpt=vlm_cpt,
                    np=np,
                )
                plot_uncertainty_circles(
                    figure,
                    stations,
                    scale=args.uncertainty_radius_scale,
                    np=np,
                )
                uncertainty_legend(
                    figure,
                    legend_values,
                    scale=args.uncertainty_radius_scale,
                )
                map_header(figure, "A  InSAR VLM and GNSS Up")
                figure.colorbar(
                    cmap=str(vlm_cpt),
                    position="JBC+w6.2c/0.30c+h+o0c/0.75c",
                    frame="x+lVLM / GNSS Up (mm yr@+-1@+)",
                )

            with figure.set_panel(panel=[0, 1]):
                figure.basemap(
                    region=region,
                    projection="M?",
                    frame=common_map_frame(region, target_intervals=3),
                )
                figure.grdimage(
                    grid=grid,
                    cmap=str(vlm_cpt),
                    nan_transparent=True,
                    interpolation="n",
                )
                plot_coast_after_raster(figure, args)
                plot_station_squares(
                    figure,
                    stations,
                    stations.residual,
                    cpt=residual_cpt,
                    np=np,
                )
                plot_uncertainty_circles(
                    figure,
                    stations,
                    scale=args.uncertainty_radius_scale,
                    np=np,
                )
                annotate_raw_metrics(figure, metrics, region)
                map_header(figure, "B  Raw InSAR minus GNSS residual")
                figure.colorbar(
                    cmap=str(residual_cpt),
                    position="JBC+w6.2c/0.30c+h+o0c/0.75c",
                    frame="x+lResidual (mm yr@+-1@+)",
                )

            with figure.set_panel(panel=[1, 0]):
                figure.basemap(
                    region=[-histogram_limit, histogram_limit, 0, hist_ymax],
                    projection="X?/?",
                    frame=[
                        "WSne",
                        "xafg+lResidual (mm yr@+-1@+)",
                        "yafg+lCount",
                    ],
                )
                figure.histogram(
                    data=residual_included,
                    series=bin_width,
                    histtype=0,
                    fill="#5F8FBF",
                    pen="0.45p,black",
                )
                figure.plot(x=[0, 0], y=[0, hist_ymax], pen="0.8p,black,--")
                bias = float(metrics["bias_mm_per_year"])
                figure.plot(x=[bias, bias], y=[0, hist_ymax], pen="1.0p,#D55E00")
                map_header(figure, f"C  Raw residual distribution (N={metrics['n']})")

            with figure.set_panel(panel=[1, 1]):
                figure.basemap(
                    region=scatter_region,
                    projection="X?/?",
                    frame=[
                        "WSne",
                        "xafg+lInSAR VLM (mm yr@+-1@+)",
                        "yafg+lGNSS Up (mm yr@+-1@+)",
                    ],
                )
                bounds = np.asarray(scatter_region[:2], dtype=float)
                figure.plot(x=bounds, y=bounds, pen="1.0p,black")
                error_mask = included & np.isfinite(stations.sigma)
                if np.any(error_mask):
                    figure.plot(
                        data=np.column_stack(
                            (
                                stations.insar[error_mask],
                                stations.gnss[error_mask],
                                stations.sigma[error_mask],
                            )
                        ),
                        style="c0.001c",
                        pen="0p,white@100",
                        error_bar=f"y+w3p+p{ERROR_PEN}",
                    )
                figure.plot(
                    x=stations.insar[included],
                    y=stations.gnss[included],
                    style="s0.18c",
                    fill=stations.residual[included],
                    cmap=str(residual_cpt),
                    pen="0.45p,black",
                )
                map_header(figure, "D  Pointwise comparison with GNSS uncertainty")

    save_figure(figure, output, dpi=args.dpi)
    return {
        "command": "vlm",
        "output": str(output),
        "station_rows": stations.count,
        "included_in_metrics": int(np.count_nonzero(stations.include)),
        "vlm_color_limit_mm_per_year": vlm_limit,
        "residual_color_limit_mm_per_year": residual_limit,
        "raw_metrics": metrics,
        "offset_removed": False,
    }


def histogram_bin_width(values: Any, np: Any) -> float:
    values = np.asarray(values, dtype=float)
    if values.size < 2:
        return 1.0
    q25, q75 = np.percentile(values, [25, 75])
    iqr = float(q75 - q25)
    if iqr > 0:
        width = 2.0 * iqr / (values.size ** (1.0 / 3.0))
    else:
        span = float(np.max(values) - np.min(values))
        width = span / max(5.0, math.sqrt(values.size)) if span > 0 else 1.0
    return max(width, 1e-9)


def scatter_axis_region(
    stations: StationValues,
    np: Any,
    *,
    explicit: Sequence[float] | None,
) -> list[float]:
    if explicit is not None:
        if len(explicit) != 2 or not float(explicit[0]) < float(explicit[1]):
            raise ValidationError("--scatter-range must be MIN MAX with MIN < MAX")
        low, high = map(float, explicit)
    else:
        mask = stations.include
        values = np.concatenate((stations.insar[mask], stations.gnss[mask]))
        sigma = stations.sigma[mask]
        finite_sigma = np.where(np.isfinite(sigma), sigma, 0.0)
        lower = np.concatenate((stations.insar[mask], stations.gnss[mask] - finite_sigma))
        upper = np.concatenate((stations.insar[mask], stations.gnss[mask] + finite_sigma))
        low = float(np.min(lower))
        high = float(np.max(upper))
        span = high - low
        padding = 0.08 * span if span > 0 else 1.0
        low -= padding
        high += padding
    return [low, high, low, high]


def plot_los(args: argparse.Namespace, deps: Dependencies) -> dict[str, Any]:
    np = deps.np
    pygmt = deps.pygmt
    stations, metrics, raster, region = prepare_map_inputs(args, mode="los", deps=deps)
    output = ensure_output(args.output, overwrite=args.overwrite)
    los_limit = robust_symmetric_limit(
        raster.finite_values,
        np,
        quantile=args.color_quantile,
        explicit=args.los_limit,
    )
    residual_limit = robust_symmetric_limit(
        stations.residual,
        np,
        quantile=args.color_quantile,
        explicit=args.residual_limit,
    )
    legend_values = auto_legend_values(stations.sigma, np)
    grid = read_raster_grid(raster.path, region, deps)

    cpt_dir = Path(os.environ["GMT_TMPDIR"])
    los_cpt = cpt_dir / f"los_{os.getpid()}.cpt"
    residual_cpt = cpt_dir / f"los_residual_{os.getpid()}.cpt"
    make_cpt(deps, path=los_cpt, limit=los_limit)
    make_cpt(deps, path=residual_cpt, limit=residual_limit)

    figure = pygmt.Figure()
    with pygmt.config(
        MAP_FRAME_TYPE="plain",
        MAP_FRAME_PEN="1p,black",
        FONT_ANNOT_PRIMARY="8p,Helvetica,black",
        FONT_LABEL="9p,Helvetica,black",
        MAP_TICK_LENGTH_PRIMARY="3p",
        FORMAT_GEO_MAP="ddd.xxF",
    ):
        figure.basemap(
            region=region,
            projection="M16c",
            frame=common_map_frame(region, target_intervals=5),
        )
        figure.grdimage(
            grid=grid,
            cmap=str(los_cpt),
            nan_transparent=True,
            interpolation="n",
        )
        plot_coast_after_raster(figure, args)
        plot_station_squares(
            figure,
            stations,
            stations.residual,
            cpt=residual_cpt,
            np=np,
        )
        plot_uncertainty_circles(
            figure,
            stations,
            scale=args.uncertainty_radius_scale,
            np=np,
        )
        uncertainty_legend(
            figure,
            legend_values,
            scale=args.uncertainty_radius_scale,
        )
        annotate_raw_metrics(figure, metrics, region)
        map_header(figure, f"A  {args.title}")
        figure.colorbar(
            cmap=str(los_cpt),
            position="JBC+w6.2c/0.30c+h+o-3.65c/0.85c",
            frame="x+lInSAR LOS (mm yr@+-1@+)",
        )
        figure.colorbar(
            cmap=str(residual_cpt),
            position="JBC+w6.2c/0.30c+h+o3.65c/0.85c",
            frame="x+lRaw residual (mm yr@+-1@+)",
        )

    save_figure(figure, output, dpi=args.dpi)
    return {
        "command": "los",
        "output": str(output),
        "station_rows": stations.count,
        "included_in_metrics": int(np.count_nonzero(stations.include)),
        "los_color_limit_mm_per_year": los_limit,
        "residual_color_limit_mm_per_year": residual_limit,
        "raw_metrics": metrics,
        "offset_removed": False,
    }


def load_trials(path: Path, *, group_column: str, np: Any) -> dict[str, dict[str, Any]]:
    rows, fields = read_csv_rows(path)
    require_columns(path, fields, TRIAL_COLUMNS | {group_column})
    groups: dict[str, dict[str, list[float]]] = {}
    seen_trials: set[tuple[str, str]] = set()
    for index, row in enumerate(rows):
        group = row[group_column]
        if not group:
            raise ValidationError(f"{path.name} row {index + 2}: {group_column} is empty")
        trial = row["trial"]
        key = (group, trial)
        if key in seen_trials:
            raise ValidationError(f"Duplicate trial {trial!r} for group {group!r}")
        seen_trials.add(key)
        all_rms = parse_float(
            row["all_rms_mm_per_year"],
            context=f"{path.name} row {index + 2} all_rms_mm_per_year",
        )
        holdout_rms = parse_float(
            row["holdout_rms_mm_per_year"],
            context=f"{path.name} row {index + 2} holdout_rms_mm_per_year",
        )
        if all_rms < 0 or holdout_rms < 0:
            raise ValidationError("RMS values cannot be negative")
        group_values = groups.setdefault(group, {"all": [], "holdout": [], "trial": []})
        group_values["trial"].append(trial)
        group_values["all"].append(all_rms)
        group_values["holdout"].append(holdout_rms)
    return {
        group: {
            "trial": values["trial"],
            "all": np.asarray(values["all"], dtype=float),
            "holdout": np.asarray(values["holdout"], dtype=float),
        }
        for group, values in groups.items()
    }


def load_loocv_metrics(
    path: Path | None,
    *,
    group_column: str,
) -> dict[str, float]:
    if path is None:
        return {}
    rows, fields = read_csv_rows(path)
    required = {group_column, "rms_mm_per_year"}
    require_columns(path, fields, required)
    result: dict[str, float] = {}
    for index, row in enumerate(rows):
        group = row[group_column]
        if group in result:
            raise ValidationError(f"Duplicate LOOCV metrics row for group {group!r}")
        result[group] = parse_float(
            row["rms_mm_per_year"],
            context=f"{path.name} row {index + 2} rms_mm_per_year",
        )
    return result


def plot_cv_rms(args: argparse.Namespace, deps: Dependencies) -> dict[str, Any]:
    np = deps.np
    pygmt = deps.pygmt
    groups = load_trials(args.trials, group_column=args.group_column, np=np)
    loocv = load_loocv_metrics(args.metrics, group_column=args.group_column)
    unknown = sorted(set(loocv) - set(groups))
    if unknown:
        raise ValidationError(
            f"LOOCV metrics contain group(s) absent from trials: {', '.join(unknown)}"
        )
    output = ensure_output(args.output, overwrite=args.overwrite)

    group_names = list(groups)
    all_values = np.concatenate(
        [np.concatenate((values["all"], values["holdout"])) for values in groups.values()]
    )
    if loocv:
        all_values = np.concatenate((all_values, np.asarray(list(loocv.values()), dtype=float)))
    ymax = args.rms_limit or float(np.max(all_values)) * 1.18
    if not math.isfinite(ymax) or ymax <= 0:
        raise ValidationError("--rms-limit must be a positive finite number")
    rng = np.random.default_rng(args.seed)

    figure = pygmt.Figure()
    with pygmt.config(
        MAP_FRAME_TYPE="plain",
        MAP_FRAME_PEN="1p,black",
        FONT_ANNOT_PRIMARY="8p,Helvetica,black",
        FONT_LABEL="9p,Helvetica,black",
        MAP_GRID_PEN_PRIMARY="0.3p,gray85",
        MAP_TICK_LENGTH_PRIMARY="3p",
    ):
        figure.basemap(
            region=[-0.6, len(group_names) - 0.4, 0, ymax],
            projection="X16c/8.5c",
            frame=["WSne", "x0", "yafg+lRMS (mm yr@+-1@+)"],
        )
        for group_index, group in enumerate(group_names):
            values = groups[group]
            for offset, key, color in (
                (-0.16, "all", "#0072B2"),
                (+0.16, "holdout", "#D55E00"),
            ):
                rms = values[key]
                jitter = rng.uniform(-0.07, 0.07, len(rms))
                x = group_index + offset + jitter
                figure.plot(
                    x=x,
                    y=rms,
                    style="c0.085c",
                    fill=f"{color}@38",
                    pen=f"0.15p,{color}",
                )
                median = float(np.median(rms))
                figure.plot(
                    x=[group_index + offset - 0.12, group_index + offset + 0.12],
                    y=[median, median],
                    pen="1.4p,black",
                )
            if group in loocv:
                value = loocv[group]
                figure.plot(
                    x=[group_index - 0.42, group_index + 0.42],
                    y=[value, value],
                    pen="1.0p,magenta4,--",
                )
            figure.plot(
                x=[group_index, group_index],
                y=[0, -0.025 * ymax],
                pen="0.8p,black",
                no_clip=True,
            )
            figure.text(
                x=group_index,
                y=-0.065 * ymax,
                text=group,
                font="8p,Helvetica,black",
                justify="TC",
                no_clip=True,
            )
        legend_lines = (
            "S 0.10c c 0.18c #0072B2 0.2p,#0072B2 0.28c All stations\n"
            "S 0.10c c 0.18c #D55E00 0.2p,#D55E00 0.28c Held-out stations\n"
            "S 0.10c - 0.28c - 1.0p,black 0.28c Median\n"
        )
        if loocv:
            legend_lines += "S 0.10c - 0.28c - 1.0p,magenta4,-- 0.28c Group-wise LOOCV\n"
        figure.legend(
            spec=io.StringIO(legend_lines),
            position="jTR+o0.12c",
            box="+gwhite@10+p0.45p,gray35",
        )
        map_header(figure, "A  Repeated spatial holdout RMS distributions")

    save_figure(figure, output, dpi=args.dpi)
    return {
        "command": "cv-rms",
        "output": str(output),
        "groups": {
            group: {
                "trials": len(values["trial"]),
                "all_median_rms_mm_per_year": float(np.median(values["all"])),
                "holdout_median_rms_mm_per_year": float(np.median(values["holdout"])),
                "loocv_rms_mm_per_year": loocv.get(group),
            }
            for group, values in groups.items()
        },
        "note": "Only held-out RMS is fully out of sample.",
    }


def preflight(args: argparse.Namespace, deps: Dependencies) -> dict[str, Any]:
    report: dict[str, Any] = {
        "status": "ok",
        "python": sys.version.split()[0],
        "numpy": getattr(deps.np, "__version__", "unknown"),
        "pygmt": getattr(deps.pygmt, "__version__", "unknown"),
        "GMT_USERDIR": os.environ.get("GMT_USERDIR"),
        "GMT_TMPDIR": os.environ.get("GMT_TMPDIR"),
        "mode": args.mode,
    }
    if deps.rasterio is not None:
        report["rasterio"] = getattr(deps.rasterio, "__version__", "unknown")
    if deps.xarray is not None:
        report["xarray"] = getattr(deps.xarray, "__version__", "unknown")
    if args.mode in {"vlm", "los"}:
        if args.station_values is None or args.metrics is None or args.raster is None:
            raise ValidationError(
                f"preflight --mode {args.mode} requires --station-values, --metrics, and --raster"
            )
        stations = load_station_values(args.station_values, mode=args.mode, np=deps.np)
        metric_row = load_metric_row(args.metrics, args.metric_filter)
        metrics = verify_metrics(stations, metric_row, deps.np)
        raster = inspect_raster(args.raster, deps)
        verify_raster_samples(
            args.raster,
            stations,
            deps,
            tolerance=args.sample_tolerance,
        )
        report.update(
            {
                "station_rows": stations.count,
                "included_in_metrics": int(deps.np.count_nonzero(stations.include)),
                "raw_metrics": metrics,
                "raster": {
                    "path": str(raster.path),
                    "crs": raster.crs,
                    "width": raster.width,
                    "height": raster.height,
                    "region": list(raster.region),
                },
                "exact_pixel_sampling_verified": True,
                "offset_removed": False,
            }
        )
    elif args.mode == "cv-rms":
        if args.trials is None:
            raise ValidationError("preflight --mode cv-rms requires --trials")
        trials = load_trials(args.trials, group_column=args.group_column, np=deps.np)
        loocv = load_loocv_metrics(args.metrics, group_column=args.group_column)
        report["groups"] = {
            group: {
                "trials": len(values["trial"]),
                "has_loocv_metric": group in loocv,
            }
            for group, values in trials.items()
        }
    return report


def add_common_map_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--raster", type=Path, required=True, help="Geographic single-band InSAR GeoTIFF")
    parser.add_argument("--station-values", type=Path, required=True, help="Standardized station_values CSV")
    parser.add_argument("--metrics", type=Path, required=True, help="Standardized raw metrics CSV")
    parser.add_argument(
        "--metric-filter",
        action="append",
        default=[],
        metavar="COLUMN=VALUE",
        help="Select one metrics row; repeat for multiple equality filters",
    )
    parser.add_argument("--output", type=Path, required=True, help="Output .png/.pdf/.svg/.eps/.tif path")
    parser.add_argument("--region", nargs=4, type=float, metavar=("W", "E", "S", "N"), help="Optional map region; defaults to raster bounds")
    parser.add_argument("--residual-limit", type=float, help="Positive symmetric residual color limit")
    parser.add_argument("--color-quantile", type=float, default=98.0, help="Percentile used for automatic symmetric limits (default: 98)")
    parser.add_argument("--sample-tolerance", type=float, default=1e-4, help="Maximum station-table versus exact-raster-pixel mismatch in mm/yr")
    parser.add_argument("--uncertainty-radius-scale", type=float, default=0.055, help="Circle radius in cm per mm/yr of GNSS rate uncertainty")
    parser.add_argument("--coast-resolution", choices=("c", "l", "i", "h", "f"), default="f", help="GSHHG coastline resolution (default: f)")
    parser.add_argument("--coast-area-threshold", type=float, default=0.0, help="Minimum coastline polygon area in km^2")
    parser.add_argument(
        "--water-fill",
        default=OCEAN_COLOR,
        help=(
            "GMT color used to fill water after the raster (default: light blue); "
            "use 'none' for diagnostic rasters intentionally extending offshore"
        ),
    )
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--overwrite", action="store_true", help="Allow replacement of an existing output figure")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Create optional PyGMT validation figures from standardized, already "
            "computed InSAR-GNSS tables. No offset is fitted or removed."
        )
    )
    parser.add_argument(
        "--gmt-workdir",
        type=Path,
        default=default_gmt_workdir(),
        help="Writable directory for GMT state/cache (default: system temporary directory)",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    preflight_parser = subparsers.add_parser(
        "preflight", help="Check dependencies and optionally validate standardized inputs without plotting"
    )
    preflight_parser.add_argument("--mode", choices=("dependencies", "vlm", "los", "cv-rms"), default="dependencies")
    preflight_parser.add_argument("--station-values", type=Path)
    preflight_parser.add_argument("--metrics", type=Path)
    preflight_parser.add_argument("--raster", type=Path)
    preflight_parser.add_argument("--trials", type=Path)
    preflight_parser.add_argument("--group-column", default="orbit")
    preflight_parser.add_argument("--metric-filter", action="append", default=[], metavar="COLUMN=VALUE")
    preflight_parser.add_argument("--sample-tolerance", type=float, default=1e-4)

    vlm_parser = subparsers.add_parser("vlm", help="Plot 2x2 VLM maps, histogram, and scatter comparison")
    add_common_map_arguments(vlm_parser)
    vlm_parser.add_argument("--vlm-limit", type=float, help="Positive symmetric VLM color limit")
    vlm_parser.add_argument("--hist-bin-width", type=float, help="Residual histogram bin width in mm/yr")
    vlm_parser.add_argument("--scatter-range", nargs=2, type=float, metavar=("MIN", "MAX"), help="Common InSAR/GNSS scatter-axis limits")

    los_parser = subparsers.add_parser("los", help="Plot an InSAR LOS map with raw GNSS residuals")
    add_common_map_arguments(los_parser)
    los_parser.add_argument("--los-limit", type=float, help="Positive symmetric InSAR LOS color limit")
    los_parser.add_argument("--title", default="LOS validation", help="Panel subtitle placed outside the map frame")

    cv_parser = subparsers.add_parser("cv-rms", help="Plot all-station and held-out RMS distributions")
    cv_parser.add_argument("--trials", type=Path, required=True, help="Standardized trials CSV")
    cv_parser.add_argument("--metrics", type=Path, help="Optional group-wise LOOCV metrics CSV")
    cv_parser.add_argument("--group-column", default="orbit", help="Grouping column present in both CSVs (default: orbit)")
    cv_parser.add_argument("--output", type=Path, required=True)
    cv_parser.add_argument("--rms-limit", type=float, help="Positive y-axis maximum")
    cv_parser.add_argument("--seed", type=int, default=0, help="Seed used only for deterministic point jitter")
    cv_parser.add_argument("--dpi", type=int, default=300)
    cv_parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        require_raster = args.command in {"vlm", "los"} or (
            args.command == "preflight" and args.mode in {"vlm", "los"}
        )
        deps = load_dependencies(require_raster=require_raster, gmt_workdir=args.gmt_workdir)
        if args.command == "preflight":
            report = preflight(args, deps)
        elif args.command == "vlm":
            report = plot_vlm(args, deps)
        elif args.command == "los":
            report = plot_los(args, deps)
        elif args.command == "cv-rms":
            report = plot_cv_rms(args, deps)
        else:  # pragma: no cover
            parser.error(f"unknown command {args.command!r}")
            return 2
    except (ValidationError, DependencyError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

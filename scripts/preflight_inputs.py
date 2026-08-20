#!/usr/bin/env python3
"""Read-only schema and raster preflight for InSAR–GNSS validation."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
from typing import Any


COMMON_STATION_COLUMNS = {
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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def csv_header(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.reader(stream)
        try:
            return next(reader)
        except StopIteration as exc:
            raise ValueError(f"CSV is empty: {path}") from exc


def row_count(path: Path) -> int:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return max(0, sum(1 for _ in csv.reader(stream)) - 1)


def require_columns(path: Path, required: set[str]) -> list[str]:
    header = csv_header(path)
    missing = sorted(required - set(header))
    if missing:
        raise ValueError(f"{path.name} is missing columns: {', '.join(missing)}")
    return header


def load_metadata(path: Path, mode: str) -> dict[str, Any]:
    with path.open(encoding="utf-8") as stream:
        metadata = json.load(stream)
    common = {
        "velocity_units",
        "coordinate_crs",
        "vertical_positive",
        "reference_frame",
        "observation_window",
        "comparison_role",
    }
    missing = sorted(common - set(metadata))
    if missing:
        raise ValueError(f"Metadata is missing keys: {', '.join(missing)}")
    if metadata["velocity_units"] not in {"mm/year", "mm/yr"}:
        raise ValueError("velocity_units must be 'mm/year' or 'mm/yr'")
    if metadata["vertical_positive"] != "up":
        raise ValueError("This implementation requires vertical_positive='up'")
    if metadata["comparison_role"] not in {"internal", "independent"}:
        raise ValueError("comparison_role must be 'internal' or 'independent'")
    if mode == "los":
        required_los = {
            "los_sign",
            "incidence_definition",
            "heading_convention",
            "angle_units",
            "insar_product_stage",
        }
        missing = sorted(required_los - set(metadata))
        if missing:
            raise ValueError(f"LOS metadata is missing keys: {', '.join(missing)}")
        if metadata["los_sign"] not in {"toward_satellite", "away_from_satellite"}:
            raise ValueError("los_sign must be toward_satellite or away_from_satellite")
        if metadata["incidence_definition"] != "from_vertical":
            raise ValueError("Only incidence_definition='from_vertical' is implemented")
        expected_heading = "satellite_flight_heading_clockwise_from_north"
        if metadata["heading_convention"] != expected_heading:
            raise ValueError(
                f"Only heading_convention={expected_heading!r} is implemented; "
                "a radar look azimuth is not interchangeable with flight heading"
            )
        if metadata["angle_units"] != "degrees":
            raise ValueError("Only angle_units='degrees' is implemented")
    return metadata


def import_rasterio() -> Any:
    try:
        import rasterio
    except ImportError as exc:
        raise RuntimeError(
            "rasterio is required for raster preflight. Use a Python environment "
            "that provides rasterio/GDAL; this script does not install packages."
        ) from exc
    return rasterio


def raster_summary(path: Path) -> dict[str, Any]:
    rasterio = import_rasterio()
    with rasterio.open(path) as dataset:
        return {
            "file": path.name,
            "sha256": sha256_file(path),
            "shape_rows_cols": [dataset.height, dataset.width],
            "count": dataset.count,
            "dtype": dataset.dtypes[0] if dataset.count else None,
            "crs": str(dataset.crs),
            "transform": list(dataset.transform)[:6],
            "bounds": list(dataset.bounds),
            "nodata": dataset.nodata,
            "has_mask": bool(dataset.count and dataset.mask_flag_enums[0]),
        }


def audit_coregistration(paths: list[Path]) -> None:
    rasterio = import_rasterio()
    with rasterio.open(paths[0]) as reference:
        expected = (reference.shape, reference.crs, reference.transform)
        expected_count = reference.count
    if expected_count != 1:
        raise ValueError(f"Expected a single-band raster: {paths[0]}")
    for path in paths[1:]:
        with rasterio.open(path) as dataset:
            actual = (dataset.shape, dataset.crs, dataset.transform)
            if dataset.count != 1:
                raise ValueError(f"Expected a single-band raster: {path}")
            if actual != expected:
                raise ValueError(
                    f"Raster is not co-registered with {paths[0].name}: {path.name}"
                )


def ensure_files(paths: list[Path]) -> None:
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Read-only audit of station schemas and raster alignment."
    )
    subparsers = parser.add_subparsers(dest="mode", required=True)

    los = subparsers.add_parser("los", help="Audit an ENU-to-LOS comparison bundle")
    los.add_argument("--stations", type=Path, required=True)
    los.add_argument("--los", type=Path, required=True)
    los.add_argument("--incidence", type=Path, required=True)
    los.add_argument("--heading", type=Path, required=True)
    los.add_argument("--metadata", type=Path, required=True)
    los.add_argument("--require-rate-uncertainties", action="store_true")

    vlm = subparsers.add_parser("vlm", help="Audit a VLM-versus-Up comparison bundle")
    vlm.add_argument("--stations", type=Path, required=True)
    vlm.add_argument("--vlm", type=Path, required=True)
    vlm.add_argument("--metadata", type=Path, required=True)
    vlm.add_argument("--require-rate-uncertainties", action="store_true")

    cv = subparsers.add_parser("cross-validate", help="Audit grouped-CV observations")
    cv.add_argument("--observations", type=Path, required=True)
    cv.add_argument("--metadata", type=Path, required=True)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.mode == "los":
        paths = [args.stations, args.los, args.incidence, args.heading, args.metadata]
        ensure_files(paths)
        required = set(COMMON_STATION_COLUMNS)
        if args.require_rate_uncertainties:
            required |= SIGMA_COLUMNS
        header = require_columns(args.stations, required)
        metadata = load_metadata(args.metadata, "los")
        if args.require_rate_uncertainties and metadata.get(
            "los_uncertainty_model"
        ) != "diagonal_enu_rate_covariance_zero_cross_terms":
            raise ValueError(
                "Rate-uncertainty projection requires metadata "
                "los_uncertainty_model='diagonal_enu_rate_covariance_zero_cross_terms'"
            )
        audit_coregistration([args.los, args.incidence, args.heading])
        rasters = [raster_summary(path) for path in (args.los, args.incidence, args.heading)]
        station_file = args.stations
    elif args.mode == "vlm":
        paths = [args.stations, args.vlm, args.metadata]
        ensure_files(paths)
        required = set(COMMON_STATION_COLUMNS)
        if args.require_rate_uncertainties:
            required |= SIGMA_COLUMNS
        header = require_columns(args.stations, required)
        metadata = load_metadata(args.metadata, "vlm")
        rasters = [raster_summary(args.vlm)]
        station_file = args.stations
    else:
        paths = [args.observations, args.metadata]
        ensure_files(paths)
        header = require_columns(args.observations, CV_COLUMNS)
        metadata = load_metadata(args.metadata, "cross-validate")
        rasters = []
        station_file = args.observations

    result = {
        "status": "ok",
        "mode": args.mode,
        "table": {
            "file": station_file.name,
            "sha256": sha256_file(station_file),
            "rows": row_count(station_file),
            "columns": header,
        },
        "metadata": metadata,
        "rasters": rasters,
        "warnings": [
            "Preflight confirms structure, not the scientific correctness of declared conventions."
        ],
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc

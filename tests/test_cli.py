"""End-to-end tests for the portable validation and preflight CLIs.

All inputs are generated inside :class:`tempfile.TemporaryDirectory`; the
tests neither depend on nor disclose a project data path.  Raster cases use
small, single-band EPSG:4326 GeoTIFFs and are skipped with an actionable
message when rasterio/GDAL is unavailable.
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Iterable

try:
    import numpy as np
    import rasterio
    from rasterio.transform import from_origin
except ImportError as exc:  # pragma: no cover - environment-dependent skip
    np = None  # type: ignore[assignment]
    rasterio = None  # type: ignore[assignment]
    from_origin = None  # type: ignore[assignment]
    RASTER_IMPORT_ERROR = str(exc)
else:
    RASTER_IMPORT_ERROR = ""


SKILL_ROOT = Path(__file__).resolve().parents[1]
VALIDATE = SKILL_ROOT / "scripts" / "validate.py"
PREFLIGHT = SKILL_ROOT / "scripts" / "preflight_inputs.py"


def write_csv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, Any]]) -> None:
    """Write a deterministic UTF-8 CSV fixture."""

    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def common_metadata(*, comparison_role: str = "independent") -> dict[str, Any]:
    return {
        "velocity_units": "mm/year",
        "coordinate_crs": "EPSG:4326",
        "vertical_positive": "up",
        "reference_frame": "synthetic-test-frame",
        "observation_window": "2020-01-01/2026-01-01",
        "comparison_role": comparison_role,
    }


def los_metadata() -> dict[str, Any]:
    result = common_metadata(comparison_role="internal")
    result.update(
        {
            "los_sign": "toward_satellite",
            "incidence_definition": "from_vertical",
            "heading_convention": "satellite_flight_heading_clockwise_from_north",
            "angle_units": "degrees",
            "insar_product_stage": "synthetic-test-stage",
            "uncertainty_semantics": "one_standard_deviation_rate_uncertainty",
            "los_uncertainty_model": (
                "diagonal_enu_rate_covariance_zero_cross_terms"
            ),
        }
    )
    return result


def vlm_metadata() -> dict[str, Any]:
    result = common_metadata()
    result["uncertainty_semantics"] = (
        "one_standard_deviation_rate_uncertainty"
    )
    return result


def run_cli(script: Path, *arguments: object) -> subprocess.CompletedProcess[str]:
    """Run a skill CLI under the test interpreter and surface useful failures."""

    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    completed = subprocess.run(
        [sys.executable, str(script), *(str(value) for value in arguments)],
        cwd=SKILL_ROOT,
        env=environment,
        capture_output=True,
        text=True,
        timeout=90,
        check=False,
    )
    if completed.returncode != 0:
        raise AssertionError(
            f"Command failed ({completed.returncode}): {script.name} "
            f"{' '.join(str(value) for value in arguments)}\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
    return completed


def as_bool(value: str) -> bool:
    return value.strip().lower() in {"true", "1"}


RASTER_SKIP_REASON = (
    "rasterio/GDAL is required to generate temporary GeoTIFF fixtures; "
    "install rasterio in the test interpreter. Import error: " + RASTER_IMPORT_ERROR
)


@unittest.skipUnless(rasterio is not None, RASTER_SKIP_REASON)
class TestRasterValidationCLI(unittest.TestCase):
    station_fields = [
        "station_id",
        "longitude_deg",
        "latitude_deg",
        "east_mm_per_year",
        "north_mm_per_year",
        "up_mm_per_year",
        "sigma_east_mm_per_year",
        "sigma_north_mm_per_year",
        "sigma_up_mm_per_year",
    ]

    def write_raster(
        self,
        path: Path,
        values: list[list[float]],
        *,
        nodata: float = -9999.0,
    ) -> None:
        assert np is not None and rasterio is not None and from_origin is not None
        array = np.asarray(values, dtype=np.float32)
        transform = from_origin(-124.0, 48.0, 0.1, 0.1)
        with rasterio.open(
            path,
            "w",
            driver="GTiff",
            height=array.shape[0],
            width=array.shape[1],
            count=1,
            dtype="float32",
            crs="EPSG:4326",
            transform=transform,
            nodata=nodata,
        ) as dataset:
            dataset.write(array, 1)

    def test_los_zero_is_valid_exact_pixel_and_provenance_is_path_private(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            stations = root / "stations.csv"
            metadata = root / "metadata.json"
            los = root / "los.tif"
            incidence = root / "incidence.tif"
            heading = root / "heading.tif"
            output = root / "los-output"

            write_csv(
                stations,
                self.station_fields,
                [
                    {
                        "station_id": "ZERO",
                        "longitude_deg": -123.95,
                        "latitude_deg": 47.95,
                        "east_mm_per_year": 0.0,
                        "north_mm_per_year": 0.0,
                        "up_mm_per_year": 0.0,
                        "sigma_east_mm_per_year": 0.2,
                        "sigma_north_mm_per_year": 0.3,
                        "sigma_up_mm_per_year": 0.4,
                    },
                    {
                        "station_id": "ONE",
                        "longitude_deg": -123.85,
                        "latitude_deg": 47.95,
                        "east_mm_per_year": 0.0,
                        "north_mm_per_year": 0.0,
                        "up_mm_per_year": 1.0,
                        "sigma_east_mm_per_year": 0.2,
                        "sigma_north_mm_per_year": 0.3,
                        "sigma_up_mm_per_year": 0.4,
                    },
                ],
            )
            write_json(metadata, los_metadata())
            self.write_raster(los, [[0.0, 1.0], [5.0, 6.0]])
            self.write_raster(incidence, [[0.0, 0.0], [0.0, 0.0]])
            self.write_raster(heading, [[0.0, 0.0], [0.0, 0.0]])

            preflight = run_cli(
                PREFLIGHT,
                "los",
                "--stations",
                stations,
                "--los",
                los,
                "--incidence",
                incidence,
                "--heading",
                heading,
                "--metadata",
                metadata,
                "--require-rate-uncertainties",
            )
            audit = json.loads(preflight.stdout)
            self.assertEqual(audit["status"], "ok")
            self.assertEqual(audit["mode"], "los")
            self.assertEqual(audit["table"]["rows"], 2)
            self.assertEqual(len(audit["rasters"]), 3)

            run_cli(
                VALIDATE,
                "los",
                "--stations",
                stations,
                "--metadata",
                metadata,
                "--los",
                los,
                "--incidence",
                incidence,
                "--heading",
                heading,
                "--output-dir",
                output,
            )
            station_values = read_csv(output / "station_values.csv")
            by_id = {row["station_id"]: row for row in station_values}
            zero = by_id["ZERO"]
            self.assertEqual(zero["pixel_status"], "finite")
            self.assertTrue(as_bool(zero["include_in_metrics"]))
            self.assertEqual(float(zero["insar_los_mm_per_year"]), 0.0)
            self.assertEqual(float(zero["gnss_los_mm_per_year"]), 0.0)
            self.assertEqual(float(zero["residual_mm_per_year"]), 0.0)
            self.assertEqual(zero["selection_reason"], "included")
            self.assertEqual(zero["group_id"], "r0_c0")

            required_output_columns = {
                "pixel_row_zero_based",
                "pixel_column_zero_based",
                "pixel_center_x",
                "pixel_center_y",
                "pixel_status",
                "incidence_deg",
                "heading_deg",
                "insar_los_mm_per_year",
                "gnss_los_mm_per_year",
                "sigma_los_mm_per_year",
                "residual_mm_per_year",
                "selection_reason",
                "include_in_metrics",
                "group_id",
            }
            self.assertTrue(required_output_columns.issubset(station_values[0]))
            metrics = read_csv(output / "metrics.csv")[0]
            self.assertEqual(int(metrics["n"]), 2)
            self.assertAlmostEqual(float(metrics["rms_mm_per_year"]), 0.0)

            provenance_text = (output / "provenance.json").read_text(
                encoding="utf-8"
            )
            self.assertNotIn(str(root), provenance_text)
            provenance = json.loads(provenance_text)
            self.assertEqual(
                {entry["file"] for entry in provenance["inputs"].values()},
                {"stations.csv", "metadata.json", "los.tif", "incidence.tif", "heading.tif"},
            )
            for entry in provenance["inputs"].values():
                self.assertRegex(entry["sha256"], r"^[0-9a-f]{64}$")
            self.assertEqual(
                provenance["los_projection"]["heading_convention"],
                "satellite_flight_heading_clockwise_from_north",
            )
            self.assertEqual(
                provenance["los_projection"]["uncertainty_model"],
                "diagonal_enu_rate_covariance_zero_cross_terms",
            )
            self.assertFalse(
                provenance["los_projection"]["enu_cross_covariances_used"]
            )

    def test_vlm_strict_nodata_qc_and_high_sigma_exclusions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            stations = root / "stations.csv"
            metadata = root / "metadata.json"
            vlm = root / "vlm.tif"
            output = root / "vlm-output"
            fields = [*self.station_fields, "quality"]
            base_sigma = {
                "sigma_east_mm_per_year": 0.2,
                "sigma_north_mm_per_year": 0.3,
                "sigma_up_mm_per_year": 0.4,
            }
            write_csv(
                stations,
                fields,
                [
                    {
                        "station_id": "VALID",
                        "longitude_deg": -123.95,
                        "latitude_deg": 47.95,
                        "east_mm_per_year": 0.0,
                        "north_mm_per_year": 0.0,
                        "up_mm_per_year": 1.0,
                        "quality": "accepted",
                        **base_sigma,
                    },
                    {
                        "station_id": "NODATA",
                        "longitude_deg": -123.85,
                        "latitude_deg": 47.95,
                        "east_mm_per_year": 0.0,
                        "north_mm_per_year": 0.0,
                        "up_mm_per_year": 2.0,
                        "quality": "accepted",
                        **base_sigma,
                    },
                    {
                        "station_id": "HIGH_SIGMA",
                        "longitude_deg": -123.95,
                        "latitude_deg": 47.85,
                        "east_mm_per_year": 0.0,
                        "north_mm_per_year": 0.0,
                        "up_mm_per_year": 3.0,
                        "sigma_east_mm_per_year": 0.2,
                        "sigma_north_mm_per_year": 0.3,
                        "sigma_up_mm_per_year": 3.0,
                        "quality": "accepted",
                    },
                    {
                        "station_id": "BAD_QC",
                        "longitude_deg": -123.85,
                        "latitude_deg": 47.85,
                        "east_mm_per_year": 0.0,
                        "north_mm_per_year": 0.0,
                        "up_mm_per_year": 4.0,
                        "quality": "rejected",
                        **base_sigma,
                    },
                ],
            )
            write_json(metadata, vlm_metadata())
            self.write_raster(vlm, [[1.0, -9999.0], [3.0, 4.0]])

            run_cli(
                VALIDATE,
                "vlm",
                "--stations",
                stations,
                "--metadata",
                metadata,
                "--vlm",
                vlm,
                "--require-status",
                "quality=accepted",
                "--max-component-sigma",
                2.0,
                "--output-dir",
                output,
            )
            station_values = read_csv(output / "station_values.csv")
            by_id = {row["station_id"]: row for row in station_values}

            self.assertTrue(as_bool(by_id["VALID"]["include_in_metrics"]))
            self.assertEqual(by_id["VALID"]["pixel_status"], "finite")
            nodata = by_id["NODATA"]
            self.assertFalse(as_bool(nodata["include_in_metrics"]))
            self.assertEqual(nodata["pixel_status"], "no_data")
            self.assertEqual(nodata["insar_vlm_mm_per_year"], "")
            self.assertIn("raster_no_data", nodata["selection_reason"])
            self.assertNotEqual(nodata["insar_vlm_mm_per_year"], "3.0")

            high_sigma = by_id["HIGH_SIGMA"]
            self.assertFalse(as_bool(high_sigma["include_in_metrics"]))
            self.assertIn(
                "component_rate_uncertainty_above_threshold",
                high_sigma["selection_reason"],
            )
            bad_qc = by_id["BAD_QC"]
            self.assertFalse(as_bool(bad_qc["include_in_metrics"]))
            self.assertIn(
                "status_quality_not_accepted", bad_qc["selection_reason"]
            )

            metrics = read_csv(output / "metrics.csv")[0]
            self.assertEqual(int(metrics["n"]), 1)
            self.assertAlmostEqual(float(metrics["bias_mm_per_year"]), 0.0)
            exclusion_ids = {
                row["station_id"] for row in read_csv(output / "exclusions.csv")
            }
            self.assertEqual(exclusion_ids, {"NODATA", "HIGH_SIGMA", "BAD_QC"})

    def test_sigma_qc_precedes_physical_site_deduplication(self) -> None:
        """A rejected long record must not displace a qualified co-site record."""

        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            stations = root / "stations.csv"
            metadata = root / "metadata.json"
            vlm = root / "vlm.tif"
            output = root / "vlm-output"
            fields = [*self.station_fields, "duration_years"]
            write_csv(
                stations,
                fields,
                [
                    {
                        "station_id": "QUALIFIED",
                        "longitude_deg": -123.9500,
                        "latitude_deg": 47.9500,
                        "east_mm_per_year": 0.0,
                        "north_mm_per_year": 0.0,
                        "up_mm_per_year": 1.0,
                        "sigma_east_mm_per_year": 0.2,
                        "sigma_north_mm_per_year": 0.3,
                        "sigma_up_mm_per_year": 0.4,
                        "duration_years": 15.0,
                    },
                    {
                        "station_id": "LONG_HIGH_SIGMA",
                        "longitude_deg": -123.9502,
                        "latitude_deg": 47.9502,
                        "east_mm_per_year": 0.0,
                        "north_mm_per_year": 0.0,
                        "up_mm_per_year": 1.0,
                        "sigma_east_mm_per_year": 0.2,
                        "sigma_north_mm_per_year": 0.3,
                        "sigma_up_mm_per_year": 3.0,
                        "duration_years": 25.0,
                    },
                    {
                        "station_id": "CONTROL",
                        "longitude_deg": -123.8500,
                        "latitude_deg": 47.9500,
                        "east_mm_per_year": 0.0,
                        "north_mm_per_year": 0.0,
                        "up_mm_per_year": 2.0,
                        "sigma_east_mm_per_year": 0.2,
                        "sigma_north_mm_per_year": 0.3,
                        "sigma_up_mm_per_year": 0.4,
                        "duration_years": 16.0,
                    },
                ],
            )
            write_json(metadata, vlm_metadata())
            self.write_raster(vlm, [[1.0, 2.0], [3.0, 4.0]])

            run_cli(
                VALIDATE,
                "vlm",
                "--stations",
                stations,
                "--metadata",
                metadata,
                "--vlm",
                vlm,
                "--max-component-sigma",
                2.0,
                "--deduplicate-within-m",
                100.0,
                "--output-dir",
                output,
            )
            by_id = {
                row["station_id"]: row
                for row in read_csv(output / "station_values.csv")
            }
            self.assertTrue(as_bool(by_id["QUALIFIED"]["include_in_metrics"]))
            self.assertNotIn(
                "nonrepresentative_physical_site_record",
                by_id["QUALIFIED"]["selection_reason"],
            )
            self.assertFalse(
                as_bool(by_id["LONG_HIGH_SIGMA"]["include_in_metrics"])
            )
            self.assertIn(
                "component_rate_uncertainty_above_threshold",
                by_id["LONG_HIGH_SIGMA"]["selection_reason"],
            )
            self.assertTrue(as_bool(by_id["CONTROL"]["include_in_metrics"]))
            self.assertEqual(int(read_csv(output / "metrics.csv")[0]["n"]), 2)


class TestCrossValidationCLI(unittest.TestCase):
    def make_observations(self, path: Path) -> int:
        fields = [
            "orbit",
            "station_id",
            "group_id",
            "normalized_pixel_x",
            "normalized_pixel_y",
            "residual_before_fit_mm_per_year",
        ]
        rows: list[dict[str, Any]] = []
        coordinates = [-0.9, -0.54, -0.18, 0.18, 0.54, 0.9]
        for y in coordinates:
            for x in coordinates:
                index = len(rows)
                response = (
                    1.2 * x
                    - 0.7 * y
                    + 0.3
                    + 0.8 * x * y
                    - 0.4 * x * x
                    + 0.6 * y * y
                )
                rows.append(
                    {
                        "orbit": "ascending",
                        "station_id": f"S{index:03d}",
                        "group_id": f"G{index:03d}",
                        "normalized_pixel_x": x,
                        "normalized_pixel_y": y,
                        "residual_before_fit_mm_per_year": response,
                    }
                )
        duplicate = dict(rows[0])
        duplicate["station_id"] = "S036_DUPLICATE"
        rows.append(duplicate)
        write_csv(path, fields, rows)
        return len(rows)

    def test_cross_validation_outputs_are_group_safe_and_seed_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            observations = root / "observations.csv"
            metadata = root / "metadata.json"
            output_one = root / "cv-one"
            output_two = root / "cv-two"
            row_count = self.make_observations(observations)
            write_json(metadata, common_metadata(comparison_role="internal"))

            preflight = run_cli(
                PREFLIGHT,
                "cross-validate",
                "--observations",
                observations,
                "--metadata",
                metadata,
            )
            audit = json.loads(preflight.stdout)
            self.assertEqual(audit["status"], "ok")
            self.assertEqual(audit["mode"], "cross-validate")
            self.assertEqual(audit["table"]["rows"], row_count)

            common_arguments: tuple[object, ...] = (
                "cross-validate",
                "--observations",
                observations,
                "--metadata",
                metadata,
                "--trials",
                6,
                "--holdout-fraction",
                0.20,
                "--seed",
                314159,
            )
            run_cli(VALIDATE, *common_arguments, "--output-dir", output_one)
            run_cli(VALIDATE, *common_arguments, "--output-dir", output_two)

            expected_files = {
                "full_fit_metrics.csv",
                "loocv_station_values.csv",
                "loocv_metrics.csv",
                "loocv_folds.csv",
                "cv_trials.csv",
                "cv_assignments.csv",
                "cv_distribution.csv",
                "provenance.json",
            }
            self.assertTrue(expected_files.issubset({path.name for path in output_one.iterdir()}))
            for filename in expected_files:
                self.assertEqual(
                    (output_one / filename).read_bytes(),
                    (output_two / filename).read_bytes(),
                    msg=f"Seeded output differs for {filename}",
                )

            folds = read_csv(output_one / "loocv_folds.csv")
            self.assertEqual(len(folds), 36)
            duplicate_fold = next(row for row in folds if row["heldout_groups"] == "G000")
            self.assertEqual(int(duplicate_fold["n_heldout_observations"]), 2)
            self.assertEqual(int(duplicate_fold["n_fitting_observations"]), 35)

            trials = read_csv(output_one / "cv_trials.csv")
            self.assertEqual(len(trials), 6)
            self.assertTrue(
                {
                    "all_rms_mm_per_year",
                    "fitting_rms_mm_per_year",
                    "holdout_rms_mm_per_year",
                    "heldout_groups",
                }.issubset(trials[0])
            )
            assignments = read_csv(output_one / "cv_assignments.csv")
            self.assertEqual(len(assignments), 6 * 7)
            for trial_number in range(1, 7):
                heldout = [
                    row["group_id"]
                    for row in assignments
                    if int(row["trial"]) == trial_number
                ]
                self.assertEqual(len(heldout), 7)
                self.assertEqual(len(set(heldout)), 7)

            loocv_values = read_csv(output_one / "loocv_station_values.csv")
            self.assertEqual(len(loocv_values), row_count)
            self.assertTrue(
                {
                    "predicted_correction_mm_per_year",
                    "residual_mm_per_year",
                }.issubset(loocv_values[0])
            )
            loocv_metrics = read_csv(output_one / "loocv_metrics.csv")[0]
            self.assertEqual(int(loocv_metrics["n"]), row_count)
            self.assertLess(float(loocv_metrics["rms_mm_per_year"]), 1.0e-10)

            provenance_text = (output_one / "provenance.json").read_text(
                encoding="utf-8"
            )
            self.assertNotIn(str(root), provenance_text)
            provenance = json.loads(provenance_text)
            self.assertEqual(provenance["model"]["seed"], 314159)
            self.assertEqual(provenance["model"]["trials"], 6)
            self.assertEqual(provenance["model"]["split_unit"], "group_id")
            self.assertEqual(provenance["counts"]["unique_groups"], 36)


if __name__ == "__main__":
    unittest.main()

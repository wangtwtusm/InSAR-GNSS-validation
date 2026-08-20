"""Deterministic unit tests for the reusable validation numerical core."""

from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

import numpy as np


SKILL_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKILL_ROOT / "scripts"))

import validation_core as core  # noqa: E402


CONVENTION = core.HEADING_CONVENTION


def exact_quadratic_sample(
    side: int = 5, *, duplicate_first_group: bool = False
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    coordinates = np.linspace(-0.9, 0.9, side)
    x, y = np.meshgrid(coordinates, coordinates)
    x = x.ravel()
    y = y.ravel()
    coefficients = np.array([1.2, -0.7, 0.3, 0.8, -0.4, 0.6])
    response = core.quadratic_basis(x, y) @ coefficients
    groups = [f"site-{index:03d}" for index in range(len(x))]
    if duplicate_first_group:
        x = np.append(x, x[0])
        y = np.append(y, y[0])
        response = np.append(response, response[0])
        groups.append(groups[0])
    return x, y, response, groups


class TestLosProjection(unittest.TestCase):
    def test_cardinal_heading_and_positive_direction_signs(self) -> None:
        # Incidence 90, northward flight heading: toward-satellite LOS is West.
        west = core.enu_to_los(
            2.0,
            0.0,
            0.0,
            90.0,
            0.0,
            heading_convention=CONVENTION,
            positive_direction="toward_satellite",
        )
        self.assertAlmostEqual(float(west), -2.0, places=12)

        # Eastward flight heading: the same convention selects positive North.
        north = core.enu_to_los(
            0.0,
            3.0,
            0.0,
            90.0,
            90.0,
            heading_convention=CONVENTION,
            positive_direction="toward_satellite",
        )
        self.assertAlmostEqual(float(north), 3.0, places=12)

        away = core.enu_to_los(
            2.0,
            0.0,
            0.0,
            90.0,
            0.0,
            heading_convention=CONVENTION,
            positive_direction="away_from_satellite",
        )
        self.assertAlmostEqual(float(away), 2.0, places=12)

    def test_vertical_and_broadcast_projection(self) -> None:
        result = core.enu_to_los(
            np.array([5.0, 8.0]),
            np.array([3.0, 4.0]),
            np.array([2.0, -1.0]),
            0.0,
            np.array([0.0, 180.0]),
            heading_convention=CONVENTION,
            positive_direction="toward_satellite",
        )
        np.testing.assert_allclose(result, [2.0, -1.0], atol=1e-12)

    def test_explicit_convention_is_enforced(self) -> None:
        with self.assertRaises(core.ValidationCoreError):
            core.enu_to_los(
                1.0,
                2.0,
                3.0,
                30.0,
                10.0,
                heading_convention="look_azimuth",
                positive_direction="toward_satellite",
            )

    def test_diagonal_sigma_projection_and_sign_invariance(self) -> None:
        toward = core.project_diagonal_sigma(
            2.0,
            3.0,
            4.0,
            0.0,
            123.0,
            heading_convention=CONVENTION,
            positive_direction="toward_satellite",
        )
        away = core.project_diagonal_sigma(
            2.0,
            3.0,
            4.0,
            0.0,
            123.0,
            heading_convention=CONVENTION,
            positive_direction="away_from_satellite",
        )
        self.assertAlmostEqual(float(toward), 4.0, places=12)
        self.assertAlmostEqual(float(away), 4.0, places=12)
        with self.assertRaises(core.ValidationCoreError):
            core.project_diagonal_sigma(
                -1.0,
                1.0,
                1.0,
                30.0,
                0.0,
                heading_convention=CONVENTION,
                positive_direction="toward_satellite",
            )


class TestMetricsAndDistance(unittest.TestCase):
    def test_raw_and_centered_metric_identity(self) -> None:
        metrics = core.accuracy_metrics([1.0, 4.0], [0.0, 1.0])
        self.assertEqual(metrics["n"], 2)
        self.assertAlmostEqual(metrics["bias_mm_per_year"], 2.0)
        self.assertAlmostEqual(metrics["rms_mm_per_year"], math.sqrt(5.0))
        self.assertAlmostEqual(metrics["centered_rms_mm_per_year"], 1.0)
        self.assertAlmostEqual(
            metrics["rms_mm_per_year"] ** 2,
            metrics["bias_mm_per_year"] ** 2
            + metrics["centered_rms_mm_per_year"] ** 2,
            places=12,
        )

    def test_pairwise_nonfinite_filtering(self) -> None:
        metrics = core.accuracy_metrics(
            [1.0, np.nan, 5.0], [0.0, 2.0, 2.0]
        )
        self.assertEqual(metrics["n"], 2)
        self.assertAlmostEqual(metrics["bias_mm_per_year"], 2.0)
        with self.assertRaises(core.ValidationCoreError):
            core.accuracy_metrics(
                [1.0, np.nan], [0.0, 2.0], drop_nonfinite=False
            )

    def test_haversine_one_degree_at_equator(self) -> None:
        distance = core.haversine_m(0.0, 0.0, 1.0, 0.0)
        self.assertAlmostEqual(float(distance), 111_195.08, delta=0.5)
        vector = core.haversine_m([0.0, 0.0], [0.0, 0.0], [0.0, 1.0], [0.0, 0.0])
        np.testing.assert_allclose(vector, [0.0, distance], atol=1e-9)


class TestConnectedSiteDeduplication(unittest.TestCase):
    def test_single_link_chain_and_representative_rule(self) -> None:
        # Adjacent 0.0007-degree steps are about 78 m at the equator.  A--C is
        # about 156 m, but A--B--C must still form one connected component.
        result = core.deduplicate_connected_sites(
            [0.0, 0.0007, 0.0014, 0.01],
            [0.0, 0.0, 0.0, 0.0],
            ["A", "B", "C", "D"],
            [10.0, 20.0, 20.0, 5.0],
            [0.2, 0.5, 0.3, 0.1],
            threshold_m=100.0,
        )
        self.assertEqual(result.components, ((0, 1, 2), (3,)))
        # B and C tie in duration; C wins on lower uncertainty.
        self.assertEqual(result.representative_indices, (2, 3))
        self.assertEqual(result.component_by_row, (0, 0, 0, 1))


class TestQuadraticModel(unittest.TestCase):
    def test_pixel_normalization(self) -> None:
        x, y = core.normalize_pixel_coordinates(
            [0, 3], [0, 3], raster_height=4, raster_width=4
        )
        np.testing.assert_allclose(x, [-0.75, 0.75])
        np.testing.assert_allclose(y, [-0.75, 0.75])

    def test_exact_six_term_fit_and_prediction(self) -> None:
        x, y, response, _ = exact_quadratic_sample(side=5)
        model = core.fit_quadratic(x, y, response)
        self.assertEqual(model.rank, 6)
        self.assertTrue(math.isfinite(model.condition_number))
        np.testing.assert_allclose(
            core.predict_quadratic(model, x, y), response, atol=1e-12
        )

    def test_rank_and_condition_guards(self) -> None:
        with self.assertRaises(core.QuadraticFitError):
            core.fit_quadratic(
                np.zeros(8), np.linspace(-1.0, 1.0, 8), np.arange(8.0)
            )
        x, y, response, _ = exact_quadratic_sample(side=5)
        with self.assertRaises(core.QuadraticFitError):
            core.fit_quadratic(
                x, y, response, maximum_condition_number=1.01
            )


class TestGroupedCrossValidation(unittest.TestCase):
    def test_loocv_holds_duplicate_group_together_and_aggregates(self) -> None:
        x, y, response, groups = exact_quadratic_sample(
            side=5, duplicate_first_group=True
        )
        result = core.grouped_loocv(x, y, response, groups)
        self.assertEqual(len(result.folds), 25)
        self.assertEqual(result.metrics["n"], 26)
        self.assertLess(result.metrics["rms_mm_per_year"], 1e-11)

        duplicate_fold = next(
            fold for fold in result.folds if fold.heldout_groups == ("site-000",)
        )
        self.assertEqual(set(duplicate_fold.heldout_indices), {0, 25})
        self.assertTrue(
            set(duplicate_fold.fitting_indices).isdisjoint(
                duplicate_fold.heldout_indices
            )
        )
        self.assertNotIn(0, duplicate_fold.fitting_indices)
        self.assertNotIn(25, duplicate_fold.fitting_indices)

    def test_seeded_repeated_holdout_is_deterministic_and_group_safe(self) -> None:
        x, y, response, groups = exact_quadratic_sample(
            side=6, duplicate_first_group=True
        )
        first = core.repeated_group_holdout(
            x, y, response, groups, trials=100, holdout_fraction=0.20, seed=17
        )
        second = core.repeated_group_holdout(
            x, y, response, groups, trials=100, holdout_fraction=0.20, seed=17
        )
        third = core.repeated_group_holdout(
            x, y, response, groups, trials=5, holdout_fraction=0.20, seed=18
        )
        self.assertEqual(len(first.trials), 100)
        self.assertEqual(
            [trial.heldout_groups for trial in first.trials],
            [trial.heldout_groups for trial in second.trials],
        )
        self.assertNotEqual(
            first.trials[0].heldout_groups, third.trials[0].heldout_groups
        )

        for trial in first.trials:
            heldout_indices = set(trial.heldout_indices)
            self.assertTrue(heldout_indices.isdisjoint(trial.fitting_indices))
            expected = {
                index
                for index, group in enumerate(groups)
                if group in set(trial.heldout_groups)
            }
            self.assertEqual(heldout_indices, expected)
            self.assertLess(
                trial.heldout_metrics["rms_mm_per_year"], 1e-10
            )


if __name__ == "__main__":
    unittest.main()

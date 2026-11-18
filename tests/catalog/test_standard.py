"""
Unit tests for the validation standard, and how it determines validation tiers.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import date, timedelta

from pydantic import ValidationError

from tatc.catalog import STANDARD_VERSION, Grade, InstrumentClass, InstrumentEntry, Tier
from tatc.catalog.schemas import get_period_days
from tatc.catalog.standard import (
    COMPATIBLE_VERSIONS,
    MIN_HOLDOUT_GAP_DAYS,
    Gate,
    get_gates,
    get_required_grades,
    is_compatible,
)

# (metric, limit, at most) by class and grade. These thresholds are frozen: they
# must test new instruments rather than follow their results. A change needs a
# new STANDARD_VERSION, with its reason, in docs/validation/standard.rst.
FROZEN_VERSION = 6
FROZEN_GATES = {
    # version 6 added the profiler class, with the gates of the altimeter class,
    # and changed no other gate
    ("profiler", "coverage"): [
        ("segment_recall", 99, False),
        ("track_offset_median", 0.5, True),
        ("epoch_p95", 3, True),
    ],
    ("profiler", "timing"): [
        ("segment_recall", 99, False),
        ("track_offset_median", 0.5, True),
        ("epoch_p95", 1, True),
    ],
    # version 5 added the regional geostationary class, and changed no other gate
    ("geostationary_regional", "coverage"): [("field_of_regard_iou", 0.95, False)],
    # version 4 added the scatterometer class, with the gates of the sar class, and
    # changed no other gate
    ("scatterometer", "coverage"): [
        ("pass_recall", 99, False),
        ("swath_edge_p95", 10, True),
        ("footprint_iou_median", 0.9, False),
        ("footprint_iou_p5", 0.9, False),
        ("epoch_p95", 3, True),
    ],
    ("scatterometer", "timing"): [
        ("pass_recall", 99, False),
        ("swath_edge_p95", 10, True),
        ("footprint_iou_median", 0.9, False),
        ("footprint_iou_p5", 0.9, False),
        ("epoch_p95", 1, True),
    ],
    # version 3 added the altimeter class, and changed no other gate
    ("altimeter", "coverage"): [
        ("segment_recall", 99, False),
        ("track_offset_median", 0.5, True),
        ("epoch_p95", 3, True),
    ],
    ("altimeter", "timing"): [
        ("segment_recall", 99, False),
        ("track_offset_median", 0.5, True),
        ("epoch_p95", 1, True),
    ],
    # version 2 added the lidar class, and changed no other gate
    ("lidar", "coverage"): [
        ("segment_recall", 99, False),
        ("beam_spacing_error_median", 0.2, True),
        ("epoch_p95", 3, True),
    ],
    ("lidar", "timing"): [
        ("segment_recall", 99, False),
        ("beam_spacing_error_median", 0.2, True),
        ("epoch_p95", 1, True),
    ],
    ("cross_track", "coverage"): [
        ("pass_recall", 99, False),
        ("pass_precision", 95, False),
        ("swath_width_median", 5, True),
        ("footprint_iou_median", 0.95, False),
        ("epoch_p95", 3, True),
    ],
    ("cross_track", "timing"): [
        ("pass_recall", 99, False),
        ("pass_precision", 95, False),
        ("swath_width_median", 5, True),
        ("footprint_iou_median", 0.95, False),
        ("epoch_p95", 1, True),
    ],
    ("cross_track", "sample_size"): [
        ("sample_size_within_3pct", 90, False),
        ("pixel_position_median", 6, True),
    ],
    ("pushbroom", "coverage"): [
        ("pass_recall", 99, False),
        ("swath_width_median", 5, True),
        ("footprint_iou_median", 0.95, False),
        ("epoch_p95", 3, True),
    ],
    ("pushbroom", "timing"): [
        ("pass_recall", 99, False),
        ("swath_width_median", 5, True),
        ("footprint_iou_median", 0.95, False),
        ("epoch_p95", 1, True),
    ],
    ("conical", "coverage"): [
        ("pass_recall", 99, False),
        ("pass_precision", 95, False),
        ("swath_width_median", 5, True),
        ("footprint_iou_median", 0.95, False),
        ("epoch_p95", 3, True),
    ],
    ("conical", "timing"): [
        ("pass_recall", 99, False),
        ("pass_precision", 95, False),
        ("swath_width_median", 5, True),
        ("footprint_iou_median", 0.95, False),
        ("epoch_p95", 1, True),
    ],
    ("conical", "sample_size"): [("sample_size_within_3pct", 90, False)],
    ("sar", "coverage"): [
        ("pass_recall", 99, False),
        ("swath_edge_p95", 10, True),
        ("footprint_iou_median", 0.90, False),
        ("footprint_iou_p5", 0.90, False),
        ("epoch_p95", 3, True),
    ],
    ("sar", "timing"): [
        ("pass_recall", 99, False),
        ("swath_edge_p95", 10, True),
        ("footprint_iou_median", 0.90, False),
        ("footprint_iou_p5", 0.90, False),
        ("epoch_p95", 1, True),
    ],
    ("geostationary", "coverage"): [("product_limit", 0.5, True)],
}


class TestFrozenGates(unittest.TestCase):
    """
    Unit tests that hold the thresholds fixed for a version of the standard.
    """

    def test_version(self):
        """
        Test that the expected thresholds below are of the current version.
        """
        self.assertEqual(STANDARD_VERSION, FROZEN_VERSION)

    def test_gates_unchanged(self):
        """
        Test each combination of class and grade: gates are unchanged, and
        combinations without gates are errors.
        """
        for instrument_class in InstrumentClass:
            for grade in Grade:
                key = (instrument_class.value, grade.value)
                if key in FROZEN_GATES:
                    self.assertEqual(
                        [
                            (g.metric, g.limit, g.at_most)
                            for g in get_gates(instrument_class, grade)
                        ],
                        FROZEN_GATES[key],
                        key,
                    )
                else:
                    with self.assertRaises(ValueError, msg=key):
                        get_gates(instrument_class, grade)

    def test_required_grades(self):
        """
        Test that classes with sample size gates need configurations for them.
        """
        self.assertEqual(
            get_required_grades(InstrumentClass.CONICAL),
            (Grade.COVERAGE, Grade.SAMPLE_SIZE),
        )
        self.assertEqual(
            get_required_grades(InstrumentClass.GEOSTATIONARY), (Grade.COVERAGE,)
        )

    def test_compatible_versions(self):
        """
        Test that results apply under the compatible versions that have their
        class's gates, as each version after the first only added a class.
        """
        self.assertEqual(COMPATIBLE_VERSIONS, tuple(range(1, STANDARD_VERSION + 1)))
        self.assertTrue(is_compatible(InstrumentClass.CROSS_TRACK, 1))
        self.assertFalse(is_compatible(InstrumentClass.LIDAR, 1))
        self.assertTrue(is_compatible(InstrumentClass.LIDAR, 2))
        self.assertFalse(is_compatible(InstrumentClass.PROFILER, 5))
        self.assertFalse(
            is_compatible(InstrumentClass.CROSS_TRACK, STANDARD_VERSION + 1)
        )

    def test_every_class_is_defined(self):
        """
        Test that every class has the version of the standard that added it and
        gates at its required grades, so that adding a class misses neither.
        """
        for instrument_class in InstrumentClass:
            with self.subTest(instrument_class.value):
                self.assertTrue(
                    any(
                        is_compatible(instrument_class, version)
                        for version in COMPATIBLE_VERSIONS
                    )
                )
                grades = get_required_grades(instrument_class)
                self.assertIn(Grade.COVERAGE, grades)
                for grade in grades:
                    gates = get_gates(instrument_class, grade)
                    self.assertTrue(gates)
                    for gate in gates:
                        # labels are shown on the catalog page
                        self.assertNotEqual(gate.label, gate.metric)

    def test_rounding_keeps_outcome(self):
        """
        Test that a rounded value never appears to meet a gate it fails.
        """
        gate = Gate("footprint_iou_median", 0.95, False, "")
        self.assertEqual(gate.format_value(0.9486), "0.9486")
        self.assertEqual(gate.format_value(0.949996), "0.949996")
        gate = Gate("pass_recall", 99, False, "%")
        self.assertEqual(gate.format_value(98.9996), "98.9996%")
        self.assertEqual(gate.format_limit(), "≥ 99%")
        gate = Gate("epoch_p95", 1, True, "s")
        self.assertEqual(gate.format_value(1.001), "1.001 s")
        self.assertEqual(gate.format_value(0.43), "0.43 s")
        self.assertEqual(gate.format_limit(), "≤ 1 s")
        self.assertEqual(gate.label, "epoch_p95")

    def test_tier_rank(self):
        """
        Test that tiers rank from weakest to strongest.
        """
        self.assertEqual(
            sorted(Tier, key=lambda tier: tier.rank),
            [Tier.NOMINAL, Tier.PARTIAL, Tier.VALIDATED],
        )


class TestPeriods(unittest.TestCase):
    """
    Unit tests for the periods of fitting and test data.
    """

    def test_forms(self):
        """
        Test a date, a time with a zone, a range, and a list.
        """
        day = date(2026, 10, 2)
        self.assertEqual(get_period_days("2026-10-02"), ((day, day),))
        self.assertEqual(get_period_days("2026-10-02T23:30Z"), ((day, day),))
        self.assertEqual(
            get_period_days("2026-10-02T23:30-02:00"),
            ((day + timedelta(days=1), day + timedelta(days=1)),),
        )
        self.assertEqual(
            get_period_days("2026-08-15 to 2026-10-02"), ((date(2026, 8, 15), day),)
        )
        self.assertEqual(
            get_period_days("2026-10-01, 2026-10-02"),
            ((day - timedelta(days=1),) * 2, (day, day)),
        )

    def test_invalid(self):
        """
        Test that periods of other forms are rejected.
        """
        for period in (
            "October 2026",
            "2026-10-02 and 2026-10-03",
            "2026-10-02T12:00",
            "2026-10-03 to 2026-10-02",
        ):
            with self.assertRaises(ValueError, msg=period):
                get_period_days(period)


def _measured():
    return {
        "kind": "measured",
        "source": "test",
        "platforms": ["test"],
        "period": "2025-12-01 to 2026-01-01",
    }


def _entry(
    *, kind="nominal", metrics=None, record=None, grade="coverage", cls="geostationary"
):
    """
    Builds an entry with one configuration, and a record of it if it has metrics.
    """
    provenance = _measured() if kind == "measured" else {"kind": kind, "source": "test"}
    data = {
        "id": "test",
        "name": "Test",
        "instrument_class": cls,
        "selection": "test",
        "platforms": [
            {
                "id": platform,
                "name": platform,
                "norad_id": 1,
                "status": "operational",
                "as_of": "2026-01-01",
            }
            for platform in ("test", "other")
        ],
        "data_products": [
            {"name": "T", "provider": "T", "url": "http://t", "access": "open"}
        ],
        "default_configuration": "c",
        "configurations": {
            "c": {
                "grade": grade,
                "description": "test",
                "platforms": ["test"],
                "velocity_frame": "earth_fixed",
                "instruments": [
                    {
                        "type": "instrument",
                        "parameters": {"name": "T", "field_of_regard": 17},
                        "provenance": {"field_of_regard": provenance},
                    }
                ],
            }
        },
    }
    if metrics is None:
        return InstrumentEntry(**data)
    digest = InstrumentEntry(**data).configurations["c"].digest
    data["record"] = {
        "entry_id": "test",
        "configurations": {
            "c": {
                "digest": digest,
                "standard_version": 1,
                "tatc_version": "0",
                "source": "test",
                "reference": "test",
                "period": "2026-03-01",
                "platforms": ["test"],
                "metrics": metrics,
                **(record or {}),
            }
        },
    }
    return InstrumentEntry(**data)


class TestTiers(unittest.TestCase):
    """
    Unit tests for the validation tiers computed from records and the standard.
    """

    good = [{"name": "product_limit", "value": 0.3}]

    def test_no_record_is_nominal(self):
        """
        Test that an entry without results is nominal.
        """
        self.assertEqual(_entry().validation.tier, Tier.NOMINAL)

    def test_passing_rerun_is_validated(self):
        """
        Test that passing all gates with reproduced evidence is validated.
        """
        validation = _entry(metrics=self.good).validation
        self.assertEqual(validation.tier, Tier.VALIDATED)
        self.assertEqual(validation.configurations["c"].gaps, ())
        self.assertTrue(validation.configurations["c"].gates[0].passed)

    def test_failed_gate_is_partial(self):
        """
        Test that a failed gate gives partial, reporting the value and limit.
        """
        validation = _entry(
            metrics=[{"name": "product_limit", "value": 0.6}]
        ).validation
        self.assertEqual(validation.tier, Tier.PARTIAL)
        self.assertIs(validation.configurations["c"].gates[0].passed, False)
        self.assertIn("0.6", validation.configurations["c"].gaps[0])
        self.assertEqual(validation.configurations["c"].evidence_gaps, ())
        self.assertEqual(
            validation.gaps, ("c: " + validation.configurations["c"].gaps[0],)
        )

    def test_gate_limit_inclusive(self):
        """
        Test that a value at the limit meets the gate.
        """
        validation = _entry(
            metrics=[{"name": "product_limit", "value": 0.5}]
        ).validation
        self.assertEqual(validation.tier, Tier.VALIDATED)

    def test_untested_gate_is_partial(self):
        """
        Test that metrics other than the gates' do not satisfy them.
        """
        validation = _entry(metrics=[{"name": "other", "value": 0}]).validation
        self.assertEqual(validation.tier, Tier.PARTIAL)
        self.assertIsNone(validation.configurations["c"].gates[0].passed)

    def test_worst_scope_decides(self):
        """
        Test that one failing part of a configuration fails its gate.
        """
        metrics = [
            {"name": "product_limit", "value": 0.1, "scope": "A"},
            {"name": "product_limit", "value": 0.9, "scope": "B"},
        ]
        result = _entry(metrics=metrics).validation.configurations["c"].gates[0]
        self.assertEqual(result.value, 0.9)
        self.assertIs(result.passed, False)

    def test_incompatible_version_is_partial(self):
        """
        Test that results under a version of the standard that does not apply are
        not validated.
        """
        validation = _entry(
            metrics=self.good, record={"standard_version": STANDARD_VERSION + 1}
        ).validation
        self.assertEqual(validation.tier, Tier.PARTIAL)
        self.assertIn("does not apply", validation.configurations["c"].gaps[0])

    def test_record_of_untested_platform_invalid(self):
        """
        Test that a record is of the platforms on which its configuration is
        tested, and its transfers of other platforms of the entry.
        """
        with self.assertRaisesRegex(ValidationError, "not those tested"):
            _entry(metrics=self.good, record={"platforms": ["other"]})
        transfer = {"period": "2026-03-02", "metrics": self.good}
        validation = _entry(
            metrics=self.good, record={"transfers": {"other": transfer}}
        ).validation
        self.assertEqual(validation.tier, Tier.VALIDATED)
        for platform in ("test", "unknown"):
            with self.assertRaisesRegex(ValidationError, "transfers", msg=platform):
                _entry(metrics=self.good, record={"transfers": {platform: transfer}})

    def test_changed_configuration_is_partial(self):
        """
        Test that results of another configuration are not validated.
        """
        validation = _entry(metrics=self.good, record={"digest": "other"}).validation
        self.assertEqual(validation.tier, Tier.PARTIAL)

    def test_fitted_values_need_holdout(self):
        """
        Test the holdout rule for fitted values, at the limit of its gap, which is
        computed from the periods of the fitting data and the test, before or
        after the fit.
        """
        fit_first, fit_last = date(2025, 12, 1), date(2026, 1, 1)
        for period, gap, expected in (
            ("2025-12-15", 0, Tier.PARTIAL),
            (
                str(fit_last + timedelta(days=MIN_HOLDOUT_GAP_DAYS - 1)),
                29,
                Tier.PARTIAL,
            ),
            (str(fit_last + timedelta(days=MIN_HOLDOUT_GAP_DAYS)), 30, Tier.VALIDATED),
            (str(fit_first - timedelta(days=MIN_HOLDOUT_GAP_DAYS)), 30, Tier.VALIDATED),
            # the nearest day of a test of several days decides
            ("2025-12-31, 2026-03-01", 0, Tier.PARTIAL),
        ):
            validation = _entry(
                kind="measured", metrics=self.good, record={"period": period}
            ).validation
            self.assertEqual(validation.configurations["c"].holdout_gap_days, gap)
            self.assertEqual(validation.tier, expected, period)

    def test_derived_values_need_no_holdout(self):
        """
        Test that values computed from published bounds need no holdout.
        """
        validation = _entry(kind="derived", metrics=self.good).validation
        self.assertEqual(validation.tier, Tier.VALIDATED)
        self.assertIsNone(validation.configurations["c"].holdout_gap_days)

    def test_missing_required_grade_limits_entry(self):
        """
        Test that an entry without a sample size configuration is not validated,
        even if all gates of its other configurations pass.
        """
        metrics = [
            {"name": "pass_recall", "value": 100},
            {"name": "pass_precision", "value": 100},
            {"name": "swath_width_median", "value": 0},
            {"name": "footprint_iou_median", "value": 1},
            {"name": "epoch_p95", "value": 0},
        ]
        validation = _entry(metrics=metrics, cls="conical", grade="timing").validation
        self.assertEqual(validation.configurations["c"].tier, Tier.VALIDATED)
        self.assertEqual(validation.tier, Tier.PARTIAL)
        self.assertEqual(
            validation.gaps, ("no configuration at the sample_size grade",)
        )

    def test_class_added_after_version_is_partial(self):
        """
        Test that results under a version before a class's gates are not
        validated.
        """
        metrics = [
            {"name": "segment_recall", "value": 100},
            {"name": "track_offset_median", "value": 0},
            {"name": "epoch_p95", "value": 0},
        ]
        for version, expected in ((5, Tier.PARTIAL), (6, Tier.VALIDATED)):
            validation = _entry(
                metrics=metrics,
                cls="profiler",
                record={"standard_version": version},
            ).validation
            self.assertEqual(validation.tier, expected, version)

    def test_grade_without_gates_invalid(self):
        """
        Test that a configuration cannot use a grade its class does not have.
        """
        with self.assertRaises(ValueError):
            _entry(grade="timing")


if __name__ == "__main__":
    unittest.main()

"""
Unit tests for the validation standard, and how it determines validation tiers.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

from tatc.catalog import (
    MIN_HOLDOUT_GAP_DAYS,
    STANDARD_VERSION,
    Grade,
    InstrumentClass,
    InstrumentEntry,
    Tier,
    get_gates,
    get_required_grades,
)

# (metric, limit, at most) by class and grade. These thresholds are frozen: they
# must test new instruments rather than follow their results. A change needs a
# new STANDARD_VERSION, with its reason, in docs/validation/standard.rst.
FROZEN_VERSION = 1
FROZEN_GATES = {
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


def _measured():
    return {
        "kind": "measured",
        "source": "test",
        "platform": "test",
        "period": "2026-01-01",
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
                "id": "test",
                "name": "Test",
                "norad_id": 1,
                "status": "operational",
                "as_of": "2026-01-01",
            }
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
        "standard_version": 1,
        "configurations": {
            "c": {
                "digest": digest,
                "evidence": "rerun",
                "tatc_version": "0",
                "source": "test",
                "reference": "test",
                "period": "2026-01-01",
                "platform": "test",
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

    def test_reported_evidence_is_partial(self):
        """
        Test that results not reproduced with exact values are not validated.
        """
        validation = _entry(
            metrics=self.good, record={"evidence": "reported"}
        ).validation
        self.assertEqual(validation.tier, Tier.PARTIAL)

    def test_changed_configuration_is_partial(self):
        """
        Test that results of another configuration are not validated.
        """
        validation = _entry(metrics=self.good, record={"digest": "other"}).validation
        self.assertEqual(validation.tier, Tier.PARTIAL)

    def test_fitted_values_need_holdout(self):
        """
        Test the holdout rule for fitted values, at the limit of its gap.
        """
        for gap, expected in (
            (None, Tier.PARTIAL),
            (MIN_HOLDOUT_GAP_DAYS - 1, Tier.PARTIAL),
            (MIN_HOLDOUT_GAP_DAYS, Tier.VALIDATED),
        ):
            validation = _entry(
                kind="measured", metrics=self.good, record={"holdout_gap_days": gap}
            ).validation
            self.assertEqual(validation.tier, expected, gap)

    def test_derived_values_need_no_holdout(self):
        """
        Test that values computed from published bounds need no holdout.
        """
        validation = _entry(kind="derived", metrics=self.good).validation
        self.assertEqual(validation.tier, Tier.VALIDATED)

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

    def test_grade_without_gates_invalid(self):
        """
        Test that a configuration cannot use a grade its class does not have.
        """
        with self.assertRaises(ValueError):
            _entry(grade="timing")


if __name__ == "__main__":
    unittest.main()

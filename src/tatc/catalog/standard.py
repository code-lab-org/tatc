"""
The standard against which catalog entries are validated.

The acceptance thresholds are frozen so that they test new instruments rather
than following their results: changing one requires a new `STANDARD_VERSION`.
See `docs/validation/standard.rst` for their back-test against earlier results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

STANDARD_VERSION = 6
"""The version of the thresholds in `get_gates`."""

_LAST_THRESHOLD_CHANGE = 1
# the version of the standard that added each class's gates
_CLASS_VERSIONS = {
    "cross_track": 1,
    "pushbroom": 1,
    "conical": 1,
    "sar": 1,
    "geostationary": 1,
    "lidar": 2,
    "altimeter": 3,
    "scatterometer": 4,
    "geostationary_regional": 5,
    "profiler": 6,
}

COMPATIBLE_VERSIONS = tuple(range(_LAST_THRESHOLD_CHANGE, STANDARD_VERSION + 1))
"""
Versions whose results still apply, as no threshold changed after the first:
version 2 only added the lidar class, version 3 the altimeter class, version 4
the scatterometer class, version 5 the regional geostationary class, and
version 6 the profiler class.
"""

MIN_HOLDOUT_GAP_DAYS = 30
"""Minimum days between the fitting data of measured values and their test data."""


class InstrumentClass(str, Enum):
    """
    Enumeration of instrument classes, which select the applicable gates.
    """

    CROSS_TRACK = "cross_track"
    PUSHBROOM = "pushbroom"
    CONICAL = "conical"
    SAR = "sar"
    GEOSTATIONARY = "geostationary"
    GEOSTATIONARY_REGIONAL = "geostationary_regional"
    LIDAR = "lidar"
    ALTIMETER = "altimeter"
    SCATTEROMETER = "scatterometer"
    PROFILER = "profiler"


class Grade(str, Enum):
    """
    Enumeration of the uses for which a configuration is validated.
    """

    COVERAGE = "coverage"
    """Which points and regions are observed, with epochs to a few seconds."""
    TIMING = "timing"
    """As for coverage, with epochs to about a second."""
    SAMPLE_SIZE = "sample_size"
    """The size, spacing, and position of samples on the ground."""


class Tier(str, Enum):
    """
    Enumeration of validation tiers, from weakest to strongest.
    """

    NOMINAL = "nominal"
    """Published values, without results compared with reference data."""
    PARTIAL = "partial"
    """Some gates are unmet or untested, or the evidence is incomplete."""
    VALIDATED = "validated"
    """All applicable gates pass, with results of the configuration as cataloged."""

    @property
    def rank(self) -> int:
        """
        The strength of the tier, from 0 (`nominal`) to 2 (`validated`).
        """
        return list(Tier).index(self)


@dataclass(frozen=True)
class Gate:
    """
    An acceptance threshold on a metric.
    """

    metric: str
    limit: float
    at_most: bool
    unit: str

    def is_met(self, value: float) -> bool:
        """
        Checks whether a value of the metric meets this gate.

        Args:
            value (float): The value.

        Returns:
            bool: `True`, if the value meets the threshold.
        """
        return value <= self.limit if self.at_most else value >= self.limit


_PASS_RECALL = Gate("pass_recall", 99, False, "%")
_PASS_PRECISION = Gate("pass_precision", 95, False, "%")
_EPOCH_COVERAGE = Gate("epoch_p95", 3, True, "s")
_EPOCH_TIMING = Gate("epoch_p95", 1, True, "s")
_SWATH_WIDTH = Gate("swath_width_median", 5, True, "km")
_SWATH_EDGE = Gate("swath_edge_p95", 10, True, "km")
_IOU = Gate("footprint_iou_median", 0.95, False, "")
_IOU_SAR = Gate("footprint_iou_median", 0.90, False, "")
_IOU_SAR_LOW = Gate("footprint_iou_p5", 0.90, False, "")
_PIXEL_POSITION = Gate("pixel_position_median", 6, True, "km")
_SAMPLE_SIZE = Gate("sample_size_within_3pct", 90, False, "%")
_PRODUCT_LIMIT = Gate("product_limit", 0.5, True, "deg")
# version 2: the share of beam segments within the views at their times, where
# the beams follow their reference tracks, and the error of the spacing of beams
_SEGMENT_RECALL = Gate("segment_recall", 99, False, "%")
_BEAM_SPACING = Gate("beam_spacing_error_median", 0.2, True, "km")
# version 3: the distance across track of a nadir altimeter's measurements from the
# predicted ground track
_TRACK_OFFSET = Gate("track_offset_median", 0.5, True, "km")
# version 5: the intersection over union of a geostationary instrument's predicted
# field of regard with the area imaged in a scan of a region
_FIELD_OF_REGARD_IOU = Gate("field_of_regard_iou", 0.95, False, "")

# gates of the coverage and timing grades that do not depend on epoch accuracy
_SWEEP_GATES: dict[InstrumentClass, tuple[Gate, ...]] = {
    InstrumentClass.CROSS_TRACK: (_PASS_RECALL, _PASS_PRECISION, _SWATH_WIDTH, _IOU),
    # the acquisitions of pushbroom imagers and SAR follow plans, so predicted
    # passes that were not acquired are not errors (no precision gate)
    InstrumentClass.PUSHBROOM: (_PASS_RECALL, _SWATH_WIDTH, _IOU),
    InstrumentClass.CONICAL: (_PASS_RECALL, _PASS_PRECISION, _SWATH_WIDTH, _IOU),
    InstrumentClass.SAR: (_PASS_RECALL, _SWATH_EDGE, _IOU_SAR, _IOU_SAR_LOW),
    # lidar beams are pointed off their reference tracks according to plans, so
    # predicted observations that were not made are not errors (no precision gate)
    InstrumentClass.LIDAR: (_SEGMENT_RECALL, _BEAM_SPACING),
    # products of nadir altimeters omit measurements by surface and quality, so
    # predicted measurements that were not made are not errors (no precision gate)
    InstrumentClass.ALTIMETER: (_SEGMENT_RECALL, _TRACK_OFFSET),
    # version 6: single-beam atmospheric profilers (lidar or radar), near nadir,
    # with the altimeter's gates; products omit profiles by quality and mode, so
    # predicted profiles that were not measured are not errors (no precision gate)
    InstrumentClass.PROFILER: (_SEGMENT_RECALL, _TRACK_OFFSET),
    # version 4: side-looking swaths, as for SAR; wind products omit land and sea
    # ice, so predicted passes that were not measured are not errors (no precision)
    InstrumentClass.SCATTEROMETER: (
        _PASS_RECALL,
        _SWATH_EDGE,
        _IOU_SAR,
        _IOU_SAR_LOW,
    ),
}

# the gates of each class at each of its grades
_GATES: dict[tuple[InstrumentClass, Grade], tuple[Gate, ...]] = {
    **{
        (instrument_class, grade): (*gates, epoch)
        for instrument_class, gates in _SWEEP_GATES.items()
        for grade, epoch in (
            (Grade.COVERAGE, _EPOCH_COVERAGE),
            (Grade.TIMING, _EPOCH_TIMING),
        )
    },
    # geostationary instruments have only a coverage grade
    (InstrumentClass.GEOSTATIONARY, Grade.COVERAGE): (_PRODUCT_LIMIT,),
    (InstrumentClass.GEOSTATIONARY_REGIONAL, Grade.COVERAGE): (_FIELD_OF_REGARD_IOU,),
    (InstrumentClass.CROSS_TRACK, Grade.SAMPLE_SIZE): (_SAMPLE_SIZE, _PIXEL_POSITION),
    (InstrumentClass.CONICAL, Grade.SAMPLE_SIZE): (_SAMPLE_SIZE,),
}


def get_gates(instrument_class: InstrumentClass, grade: Grade) -> tuple[Gate, ...]:
    """
    Gets the gates of a configuration.

    Args:
        instrument_class (InstrumentClass): The class of the instrument.
        grade (Grade): The grade of the configuration.

    Returns:
        tuple[Gate, ...]: The applicable gates.

    Raises:
        ValueError: If the class is not validated at the grade.
    """
    instrument_class, grade = InstrumentClass(instrument_class), Grade(grade)
    if (instrument_class, grade) not in _GATES:
        grades = [g.value for c, g in _GATES if c == instrument_class]
        raise ValueError(
            f"{instrument_class.value} has no {grade.value} grade, only {grades}"
        )
    return _GATES[instrument_class, grade]


def get_required_grades(instrument_class: InstrumentClass) -> tuple[Grade, ...]:
    """
    Gets the grades that an entry needs to be validated: the coverage grade, and
    the sample size grade if its class has one. A configuration at the timing
    grade also meets the coverage grade, as timing also tests coverage.

    Args:
        instrument_class (InstrumentClass): The class of the instrument.

    Returns:
        tuple[Grade, ...]: The required grades.
    """
    return tuple(
        grade
        for grade in (Grade.COVERAGE, Grade.SAMPLE_SIZE)
        if (InstrumentClass(instrument_class), grade) in _GATES
    )


def is_compatible(instrument_class: InstrumentClass, version: int) -> bool:
    """
    Checks whether results tested under a version of the standard apply to a
    class: the version is compatible and has the class's gates.

    Args:
        instrument_class (InstrumentClass): The class of the instrument.
        version (int): The version of the standard.

    Returns:
        bool: `True`, if the results apply.
    """
    return (
        version in COMPATIBLE_VERSIONS
        and version >= _CLASS_VERSIONS[InstrumentClass(instrument_class).value]
    )

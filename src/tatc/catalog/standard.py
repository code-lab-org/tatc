"""
The standard against which catalog entries are validated, and its evaluation of
their records (`evaluate_entry`).

The acceptance thresholds are frozen so that they test new instruments rather
than following their results: changing one requires a new `STANDARD_VERSION`.
See `docs/validation/standard.rst` for their back-test against earlier results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from collections.abc import Collection, Mapping
from typing import TYPE_CHECKING

import pandas as pd

if TYPE_CHECKING:
    from .records import ConfigurationRecord
    from .schemas import InstrumentEntry

STANDARD_VERSION = 8
"""The version of the thresholds in `get_gates`."""

_LAST_THRESHOLD_CHANGE = 1

COMPATIBLE_VERSIONS = tuple(range(_LAST_THRESHOLD_CHANGE, STANDARD_VERSION + 1))
"""
Versions whose results still apply, as no threshold changed after the first:
later versions only added classes (see `is_compatible`).
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
    LIMB = "limb"
    OCCULTATION = "occultation"


# the version of the standard that added each class's gates
_CLASS_VERSIONS: dict[InstrumentClass, int] = {
    InstrumentClass.CROSS_TRACK: 1,
    InstrumentClass.PUSHBROOM: 1,
    InstrumentClass.CONICAL: 1,
    InstrumentClass.SAR: 1,
    InstrumentClass.GEOSTATIONARY: 1,
    InstrumentClass.LIDAR: 2,
    InstrumentClass.ALTIMETER: 3,
    InstrumentClass.SCATTEROMETER: 4,
    InstrumentClass.GEOSTATIONARY_REGIONAL: 5,
    InstrumentClass.PROFILER: 6,
    InstrumentClass.LIMB: 7,
    InstrumentClass.OCCULTATION: 7,
}


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


# the decimals with which values of a unit are shown, at least
_DECIMALS = {"": 4}


@dataclass(frozen=True)
class Gate:
    """
    An acceptance threshold on a metric, with the metric's `label` for display
    (by default, its name).
    """

    metric: str
    limit: float
    at_most: bool
    unit: str
    label: str = ""

    def __post_init__(self):
        if not self.label:
            object.__setattr__(self, "label", self.metric)

    def is_met(self, value: float) -> bool:
        """
        Checks whether a value of the metric meets this gate.

        Args:
            value (float): The value.

        Returns:
            bool: `True`, if the value meets the threshold.
        """
        return value <= self.limit if self.at_most else value >= self.limit

    def _with_unit(self, text: str) -> str:
        return text + (
            "%" if self.unit == "%" else f" {self.unit}" if self.unit else ""
        )

    def format_limit(self) -> str:
        """
        Formats the threshold, such as `"≤ 5 km"`.

        Returns:
            str: The threshold, with its unit.
        """
        return self._with_unit(f"{'≤' if self.at_most else '≥'} {self.limit:g}")

    def format_value(self, value: float) -> str:
        """
        Formats a value of the metric, with enough decimals that its rounded
        value meets this gate exactly when the value does.

        Args:
            value (float): The value.

        Returns:
            str: The value, with its unit.
        """
        decimals = _DECIMALS.get(self.unit, 2)
        text = f"{value:.{decimals}f}"
        while decimals < 8 and self.is_met(float(text)) != self.is_met(value):
            decimals += 1
            text = f"{value:.{decimals}f}"
        return self._with_unit(text)


_PASS_RECALL = Gate("pass_recall", 99, False, "%", "Pass recall")
_PASS_PRECISION = Gate("pass_precision", 95, False, "%", "Pass precision")
_EPOCH_COVERAGE = Gate("epoch_p95", 3, True, "s", "Epoch, 95th percentile")
_EPOCH_TIMING = Gate("epoch_p95", 1, True, "s", "Epoch, 95th percentile")
_SWATH_WIDTH = Gate("swath_width_median", 5, True, "km", "Swath width, median")
_SWATH_EDGE = Gate("swath_edge_p95", 10, True, "km", "Swath edges, 95th percentile")
_IOU_LABEL = "Footprint or swath IoU, median"
_IOU = Gate("footprint_iou_median", 0.95, False, "", _IOU_LABEL)
_IOU_SAR = Gate("footprint_iou_median", 0.90, False, "", _IOU_LABEL)
_IOU_SAR_LOW = Gate(
    "footprint_iou_p5", 0.90, False, "", "Footprint or swath IoU, 5th percentile"
)
_PIXEL_POSITION = Gate("pixel_position_median", 6, True, "km", "Pixel position, median")
_SAMPLE_SIZE = Gate("sample_size_within_3pct", 90, False, "%", "Dimensions within 3%")
_PRODUCT_LIMIT = Gate(
    "product_limit", 0.5, True, "deg", "Product limit of local zenith angle"
)
# version 2: the share of beam segments within the views at their times, where
# the beams follow their reference tracks, and the error of the spacing of beams
_SEGMENT_RECALL = Gate("segment_recall", 99, False, "%", "Segment recall")
_BEAM_SPACING = Gate(
    "beam_spacing_error_median", 0.2, True, "km", "Beam spacing error, median"
)
# version 3: the distance across track of a nadir altimeter's measurements from the
# predicted ground track
_TRACK_OFFSET = Gate("track_offset_median", 0.5, True, "km", "Track offset, median")
# version 5: the intersection over union of a geostationary instrument's predicted
# field of regard with the area imaged in a scan of a region
_FIELD_OF_REGARD_IOU = Gate(
    "field_of_regard_iou", 0.95, False, "", "Field of regard IoU"
)
# version 7: limb sounders, whose scans are started at the times of the reference's
# scans and whose representative points (at the sample elevation) are compared
# with the reference's tangent points at that elevation; the period of their cycle
# of scans; and the extent in latitude of a day of representative points
_LIMB_POSITION = Gate(
    "tangent_position_median", 10, True, "km", "Tangent point position, median"
)
_SCAN_PERIOD = Gate("scan_period_error", 1, True, "%", "Scan cycle period, error")
_LATITUDE_EXTENT = Gate(
    "latitude_extent_error", 0.5, True, "deg", "Latitude extent, error"
)
# version 7: radio occultation receivers, whose observed occultations are
# predicted, with tangent points (at the sample elevation) compared with those of
# the straight line between the recorded positions of receiver and transmitter
_OCCULTATION_RECALL = Gate("occultation_recall", 99, False, "%", "Occultation recall")
_OCCULTATION_POSITION = Gate(
    "tangent_position_median", 2, True, "km", "Tangent point position, median"
)

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
    # version 7: limb sounders profile the atmosphere beside or ahead of their
    # track, in scans whose phase depends on the instrument's clock
    InstrumentClass.LIMB: (_LIMB_POSITION, _SCAN_PERIOD, _LATITUDE_EXTENT),
    # version 7: receivers track transmitters selectively, so predicted
    # occultations that were not recorded are not errors (no precision gate)
    InstrumentClass.OCCULTATION: (_OCCULTATION_RECALL, _OCCULTATION_POSITION),
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
    # version 8: the sample sizes of pushbroom imagers, with the gates of
    # cross-track scanners
    (InstrumentClass.PUSHBROOM, Grade.SAMPLE_SIZE): (_SAMPLE_SIZE, _PIXEL_POSITION),
}
# the versions of the standard that added a grade to a class that had others
_GRADE_VERSIONS: dict[tuple[InstrumentClass, Grade], int] = {
    (InstrumentClass.PUSHBROOM, Grade.SAMPLE_SIZE): 8,
}
# grades whose configurations are gated but that an entry does not need to be
# validated: the sample sizes of pushbroom imagers, whose published pixel sizes are
# rarely given away from nadir
_OPTIONAL_GRADES: frozenset[tuple[InstrumentClass, Grade]] = frozenset(
    {(InstrumentClass.PUSHBROOM, Grade.SAMPLE_SIZE)}
)


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
    the sample size grade if its class has one that is not optional (the
    pushbroom class's is optional). A configuration at the timing grade also meets
    the coverage grade, as timing also tests coverage.

    Args:
        instrument_class (InstrumentClass): The class of the instrument.

    Returns:
        tuple[Grade, ...]: The required grades.
    """
    instrument_class = InstrumentClass(instrument_class)
    return tuple(
        grade
        for grade in (Grade.COVERAGE, Grade.SAMPLE_SIZE)
        if (instrument_class, grade) in _GATES
        and (instrument_class, grade) not in _OPTIONAL_GRADES
    )


def is_compatible(
    instrument_class: InstrumentClass, version: int, grade: Grade | None = None
) -> bool:
    """
    Checks whether results tested under a version of the standard apply to a
    class (and, if given, a grade): the version is compatible and has the gates
    of the class at the grade. Version 2 added the lidar class, 3 the altimeter
    class, 4 the scatterometer class, 5 the regional geostationary class, 6 the
    profiler class, 7 the limb and occultation classes, and 8 the sample size
    grade of the pushbroom class.

    Args:
        instrument_class (InstrumentClass): The class of the instrument.
        version (int): The version of the standard.
        grade (Grade | None): The grade of the results, or `None` for any grade
            of the class.

    Returns:
        bool: `True`, if the results apply.
    """
    instrument_class = InstrumentClass(instrument_class)
    added = _CLASS_VERSIONS[instrument_class]
    if grade is not None:
        added = max(added, _GRADE_VERSIONS.get((instrument_class, Grade(grade)), 0))
    return version in COMPATIBLE_VERSIONS and version >= added


def serves(configuration_grade: Grade | str, grade: Grade | str) -> bool:
    """
    Checks whether a configuration at a grade serves a use: its own grade, or
    the coverage grade for a configuration at the timing grade, which also tests
    coverage.

    Args:
        configuration_grade (Grade | str): The grade of the configuration.
        grade (Grade | str): The grade of the use.

    Returns:
        bool: `True`, if it serves the use.
    """
    configuration_grade, grade = Grade(configuration_grade), Grade(grade)
    return configuration_grade == grade or (
        grade == Grade.COVERAGE and configuration_grade == Grade.TIMING
    )


@dataclass(frozen=True)
class GateResult:
    """
    The outcome of a gate: `None` if the metric was not tested.
    """

    gate: Gate
    value: float | None

    @property
    def passed(self) -> bool | None:
        """
        Whether the gate is met, or `None` if it was not tested.
        """
        return None if self.value is None else self.gate.is_met(self.value)

    @property
    def gap(self) -> str | None:
        """
        The gap, if the gate failed or was not tested.
        """
        if self.value is None:
            return f"{self.gate.metric}: not tested"
        if not self.passed:
            return (
                f"{self.gate.metric}: {self.value} fails the limit of "
                f"{self.gate.limit}"
            )
        return None


@dataclass(frozen=True)
class ConfigurationValidation:
    """
    The validation status of a configuration: the results of its gates, and its
    `evidence_gaps`, of its test rather than of a gate (such as a configuration
    that changed after it was tested). Its `holdout_gap_days` are the days
    between the data to which its values were fitted and the data of its test,
    or `None` if it has no fitted values or no test.
    """

    tier: Tier
    grade: Enum
    """Its grade (a `Grade`, or a `SystemGrade` of a system)."""
    gates: tuple[GateResult, ...]
    evidence_gaps: tuple[str, ...]
    record: ConfigurationRecord | None
    holdout_gap_days: int | None = None

    @property
    def gaps(self) -> tuple[str, ...]:
        """
        Every gap: of failed gates, then of untested gates, then of evidence.
        """
        return (
            *(r.gap for r in self.gates if r.passed is False),
            *(r.gap for r in self.gates if r.passed is None),
            *self.evidence_gaps,
        )


@dataclass(frozen=True)
class EntryValidation:
    """
    The validation status of an entry, computed from its records and the
    standard: of each configuration, and the grades of its class that no
    configuration serves (`missing_grades`).
    """

    tier: Tier
    configurations: dict[str, ConfigurationValidation]
    missing_grades: tuple[Enum, ...]

    @property
    def gaps(self) -> tuple[str, ...]:
        """
        Every gap: of missing grades, then of each configuration, after its
        identifier (such as `"swath: epoch_p95: not tested"`).
        """
        return (
            *(f"no configuration at the {g.value} grade" for g in self.missing_grades),
            *(
                f"{key}: {gap}"
                for key, validation in self.configurations.items()
                for gap in validation.gaps
            ),
        )

    def to_frame(self) -> pd.DataFrame:
        """
        Tabulates the gates of every configuration.

        Returns:
            pandas.DataFrame: A row for each gate of each configuration, with
            its `configuration`, `grade`, `tier`, metric `label`, `value` (`NaN`
            if not tested), formatted `limit`, and whether it `passed` (`None`
            if not tested).
        """
        return pd.DataFrame(
            [
                {
                    "configuration": key,
                    "grade": validation.grade.value,
                    "tier": validation.tier.value,
                    "metric": result.gate.label,
                    "value": result.value,
                    "limit": result.gate.format_limit(),
                    "passed": result.passed,
                }
                for key, validation in self.configurations.items()
                for result in validation.gates
            ],
            columns=[
                "configuration",
                "grade",
                "tier",
                "metric",
                "value",
                "limit",
                "passed",
            ],
        )


def _evaluate_record(  # pylint: disable=too-many-arguments
    grade: Enum,
    record: ConfigurationRecord | None,
    gates: tuple[Gate, ...],
    *,
    digest: str,
    holdout: int | None,
    evidence_gaps: tuple[str, ...] = (),
    scopes: Mapping[str, Collection[str]] | None = None,
) -> ConfigurationValidation:
    """
    Applies gates to the record of a configuration (of an instrument or a
    system): each gate takes the worst value of its metric over scopes, and the
    evidence has gaps if the configuration changed after its test (its `digest`
    differs from the record's), fitted values were tested on data too close to
    their fitting data (`holdout`), or the caller finds other gaps
    (`evidence_gaps`, such as a version of the standard that does not apply).

    Args:
        grade (Enum): The grade of the configuration.
        record (ConfigurationRecord | None): Its record, if any.
        gates (tuple[Gate, ...]): The gates of its grade.
        digest (str): The digest of the configuration as cataloged.
        holdout (int | None): The days between the data to which its values were
            fitted and the data of the record's test, or `None` if no value was
            fitted (or it has no record).
        evidence_gaps (tuple[str, ...]): Other gaps of the evidence.
        scopes (Mapping[str, Collection[str]] | None): The scopes (such as
            spacecraft) of which a gate's metric needs a value, by metric, if
            any: a gate with values of some scopes but not of these has a gap.

    Returns:
        ConfigurationValidation: Its validation status.
    """
    if record is None or not record.metrics:
        return ConfigurationValidation(
            Tier.NOMINAL,
            grade,
            tuple(GateResult(gate, None) for gate in gates),
            ("no comparison with reference data",),
            record,
        )
    results = []
    for gate in gates:
        values = [m.value for m in record.metrics if m.name == gate.metric]
        worst = (max if gate.at_most else min)(values) if values else None
        results.append(GateResult(gate, worst))
    gaps = []
    for gate in gates:
        present = {m.scope for m in record.metrics if m.name == gate.metric}
        missing = [s for s in (scopes or {}).get(gate.metric, ()) if s not in present]
        if present and missing:
            gaps.append(f"{gate.metric}: not tested on {', '.join(missing)}")
    if record.digest != digest:
        gaps.append("the configuration changed after it was tested")
    gaps.extend(evidence_gaps)
    if holdout is not None and holdout < MIN_HOLDOUT_GAP_DAYS:
        gaps.append(
            "fitted values are not tested on data at least "
            f"{MIN_HOLDOUT_GAP_DAYS} days from their fitting data"
        )
    passed = not gaps and all(result.passed for result in results)
    return ConfigurationValidation(
        Tier.VALIDATED if passed else Tier.PARTIAL,
        grade,
        tuple(results),
        tuple(gaps),
        record,
        holdout,
    )


def _combine(
    configurations: dict[str, ConfigurationValidation], missing: tuple[Enum, ...]
) -> EntryValidation:
    """
    Combines the validations of an entry's configurations: an entry is only as
    strong as its weakest configuration, is partial rather than validated if no
    configuration serves a grade it needs (`missing`), and is partial rather
    than nominal if any configuration was compared with reference data.

    Args:
        configurations (dict[str, ConfigurationValidation]): The validations.
        missing (tuple[Enum, ...]): The grades it needs that none serves.

    Returns:
        EntryValidation: Its validation status.
    """
    tier = min((c.tier for c in configurations.values()), key=lambda t: t.rank)
    if tier == Tier.VALIDATED and missing:
        tier = Tier.PARTIAL
    if tier == Tier.NOMINAL and any(
        c.tier != Tier.NOMINAL for c in configurations.values()
    ):
        tier = Tier.PARTIAL
    return EntryValidation(tier, configurations, missing)


def _evaluate_configuration(
    entry: InstrumentEntry, key: str
) -> ConfigurationValidation:
    """
    Applies the standard to a configuration of an entry (see `_evaluate_record`).

    Args:
        entry (InstrumentEntry): The entry.
        key (str): The identifier of the configuration.

    Returns:
        ConfigurationValidation: Its validation status.
    """
    configuration = entry.configurations[key]
    record = None if entry.record is None else entry.record.configurations.get(key)
    gaps = ()
    if record is not None and not is_compatible(
        entry.instrument_class, record.standard_version, configuration.grade
    ):
        gaps = (
            f"version {record.standard_version} of the standard, of the test, "
            f"does not apply to the {entry.instrument_class.value} class",
        )
    return _evaluate_record(
        configuration.grade,
        record,
        get_gates(entry.instrument_class, configuration.grade),
        digest=configuration.digest,
        holdout=(
            None
            if record is None
            else configuration.get_holdout_gap(record.period, record.platform_periods)
        ),
        evidence_gaps=gaps,
    )


def evaluate_entry(entry: InstrumentEntry) -> EntryValidation:
    """
    Applies the standard to an entry's records. An entry is only as strong as its
    weakest configuration, needs configurations that serve the grades of its
    class (see `get_required_grades`), and is partial rather than nominal if any
    configuration was compared with reference data. Configurations at a grade
    that is optional for the class (the pushbroom class's sample size grade) have
    their own tiers, which do not set the entry's.

    Args:
        entry (InstrumentEntry): The entry.

    Returns:
        EntryValidation: Its validation status.
    """
    configurations = {
        key: _evaluate_configuration(entry, key) for key in entry.configurations
    }
    missing = tuple(
        grade
        for grade in get_required_grades(entry.instrument_class)
        if not any(c.serves(grade) for c in entry.configurations.values())
    )
    instrument_class = InstrumentClass(entry.instrument_class)
    counted = {
        key: validation
        for key, validation in configurations.items()
        if (instrument_class, Grade(validation.grade)) not in _OPTIONAL_GRADES
    }
    return EntryValidation(_combine(counted, missing).tier, configurations, missing)

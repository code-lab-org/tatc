"""
Object schemas for catalog entries of validated instruments.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import copy
import hashlib
import json
import warnings
from dataclasses import dataclass
from datetime import date, datetime, time, timedelta, timezone
from typing import Annotated, Any, Literal

import numpy as np
from pydantic import AfterValidator, BaseModel, ConfigDict, Field, model_validator

from ..schemas import (
    Attitude,
    AttitudeSchedule,
    ConicalInstrument,
    GlintAttitude,
    Instrument,
    PointedInstrument,
    Satellite,
    TargetAttitude,
)
from ..schemas.orbit import AllOrbits
from ..schemas.space.base import AllAttitudes
from ..utils import (
    VelocityFrame,
    compute_along_track_field_of_view,
    compute_ascending_argument_of_latitude,
    compute_subsolar_latitude,
)
from .standard import (
    MIN_HOLDOUT_GAP_DAYS,
    Gate,
    Grade,
    InstrumentClass,
    Tier,
    get_gates,
    get_required_grades,
    is_compatible,
)

# the instrument schemas, by their `type`
_INSTRUMENT_CLASSES: dict[str, type[Instrument]] = {
    cls.model_fields["type"].default: cls  # pylint: disable=unsubscriptable-object
    for cls in (Instrument, PointedInstrument, ConicalInstrument)
}
_TAILORED = {"tailor": "time_step"}
_TAILORED_FIELD = "along_track_field_of_view"
# the provenance of parameters whose origin is the same in every configuration
_DEFAULT_PROVENANCE = {
    _TAILORED_FIELD: {
        "kind": "derived",
        "source": "compute_along_track_field_of_view, from the orbit and time step",
    },
    "is_swept": {
        "kind": "derived",
        "source": "Models the view as sweeping along track, so that its swaths are "
        "the strips its cross-track line sweeps, as observations of points do",
    },
}
# the attitudes that a configuration may give, by their `type`
_ATTITUDE_TYPES = {
    "attitude": Attitude,
    "glint": GlintAttitude,
    "target": TargetAttitude,
}


def _parse_day(text: str) -> date:
    """
    Parses an ISO 8601 date, or the UTC date of an ISO 8601 time with a zone.
    """
    try:
        return date.fromisoformat(text)
    except ValueError:
        moment = datetime.fromisoformat(text.replace("Z", "+00:00"))
        if moment.tzinfo is None:
            raise ValueError(f"the time {text!r} needs a time zone") from None
        return moment.astimezone(timezone.utc).date()


def get_period_days(period: str) -> tuple[tuple[date, date], ...]:
    """
    Gets the days of a period: one or more items separated by commas, each a
    date or a time with a zone in ISO 8601 (such as `2026-10-02` or
    `2026-10-03T12:00Z`), or a range of dates (`2026-08-15 to 2026-10-08`).

    Args:
        period (str): The period.

    Returns:
        tuple[tuple[datetime.date, datetime.date], ...]: The first and last day
        (UTC) of each item.

    Raises:
        ValueError: If the period is not of this form.
    """
    intervals = []
    for item in period.split(","):
        first, separator, last = item.strip().partition(" to ")
        start = _parse_day(first)
        end = _parse_day(last) if separator else start
        if end < start:
            raise ValueError(f"the range {item.strip()!r} ends before it starts")
        intervals.append((start, end))
    return tuple(intervals)


def _check_period(period: str) -> str:
    get_period_days(period)
    return period


Period = Annotated[str, AfterValidator(_check_period)]
"""
A period of dates (see `get_period_days`), such as `2026-10-02`,
`2026-10-03T12:00Z`, `2026-08-15 to 2026-10-08`, or `2026-11-01, 2026-11-03`.
"""


def _get_days_apart(first: str, second: str) -> int:
    """
    Gets the days between the nearest days of two periods (0 if they share one).
    """
    return min(
        max(0, (b_start - a_end).days, (a_start - b_end).days)
        for a_start, a_end in get_period_days(first)
        for b_start, b_end in get_period_days(second)
    )


class _CatalogModel(BaseModel):
    """
    Base class that rejects unknown fields, as catalog files are written by hand.
    """

    model_config = ConfigDict(extra="forbid", frozen=True)


class Provenance(_CatalogModel):
    """
    The origin of a parameter value.
    """

    kind: Literal["nominal", "measured", "derived", "placeholder"] = Field(
        ...,
        description="Origin of the value: `nominal`, a published specification; "
        + "`measured`, fitted to reference data; `derived`, computed from other "
        + "values; or `placeholder`, without physical meaning.",
    )
    source: str = Field(..., description="Where the value comes from.")
    platforms: tuple[str, ...] | None = Field(
        default=None,
        min_length=1,
        description="The platforms to whose data a measured value was fitted.",
    )
    period: Period | None = Field(
        default=None,
        description="The period of the data to which a measured value was fitted, "
        + "from which the gap to the data of its tests is computed.",
    )

    @model_validator(mode="after")
    def _check_measured(self) -> "Provenance":
        if self.kind == "measured" and (self.platforms is None or self.period is None):
            raise ValueError("a measured value needs its platforms and period")
        return self


def _check_provenance(described: set[str], provenance: dict[str, Any]) -> None:
    """
    Checks that each value has a provenance.
    """
    if described != set(provenance):
        raise ValueError(
            "each value needs provenance: "
            f"{sorted(described.symmetric_difference(provenance))}"
        )


class InstrumentSpecification(_CatalogModel):
    """
    The class and parameters of an instrument, from which objects are constructed.
    The provenance of an along-track field of view tailored to the time step and
    of `is_swept` may be omitted, as it is the same in every configuration.
    """

    type: str = Field(
        ...,
        description="The instrument class: `instrument`, `pointed`, or `conical` "
        + "(the `type` of its schema).",
    )
    parameters: dict[str, Any] = Field(..., description="The instrument's fields.")
    provenance: dict[str, Provenance] = Field(
        default_factory=dict,
        description="The origin of each parameter, except for its name.",
    )

    @model_validator(mode="before")
    @classmethod
    def _add_default_provenance(cls, data: Any) -> Any:
        if not isinstance(data, dict) or not isinstance(data.get("parameters"), dict):
            return data
        parameters, provenance = data["parameters"], dict(data.get("provenance", {}))
        for field, default in _DEFAULT_PROVENANCE.items():
            applies = (
                parameters.get(field) == _TAILORED
                if field == _TAILORED_FIELD
                else field in parameters
            )
            if applies and field not in provenance:
                provenance[field] = default
        return {**data, "provenance": provenance}

    @model_validator(mode="after")
    def _check_parameters(self) -> "InstrumentSpecification":
        if self.type not in _INSTRUMENT_CLASSES:
            raise ValueError(
                f"unknown instrument type {self.type!r}: {sorted(_INSTRUMENT_CLASSES)}"
            )
        fields = set(_INSTRUMENT_CLASSES[self.type].model_fields) - {"type"}
        if "velocity_frame" in self.parameters:
            raise ValueError(
                "velocity_frame belongs to the configuration, as it sets the "
                "attitude of the space system"
            )
        unknown = set(self.parameters) - fields
        if unknown:
            raise ValueError(f"unknown {self.type} parameters: {sorted(unknown)}")
        _check_provenance(set(self.parameters) - {"name"}, self.provenance)
        if "name" not in self.parameters:
            raise ValueError("an instrument specification needs a name")
        if (
            self.type == "pointed"
            and self.is_tailored
            and not self.parameters.get("is_swept")
        ):
            raise ValueError(
                "a pointed instrument whose along-track field of view is tailored "
                "to a time step must be swept (is_swept), so that its swaths are "
                "the strips its cross-track line sweeps rather than its view "
                "tailored to the time step"
            )
        # constructing an example checks the values of the parameters
        self.construct_instrument(1.0)
        return self

    @property
    def is_tailored(self) -> bool:
        """
        Whether the along-track field of view is tailored to a time step.
        """
        return self.parameters.get(_TAILORED_FIELD) == _TAILORED

    def construct_instrument(self, tailored_field_of_view: float | None) -> Instrument:
        """
        Constructs a new instrument that shares no state with this specification.

        Args:
            tailored_field_of_view (float | None): The along-track field of view
                (degrees), if this specification tailors it to a time step.

        Returns:
            Instrument: The instrument.
        """
        parameters = copy.deepcopy(self.parameters)
        if self.is_tailored:
            parameters[_TAILORED_FIELD] = tailored_field_of_view
        return _INSTRUMENT_CLASSES[self.type](**parameters)


class TiltSchedule(_CatalogModel):
    """
    A pitch of the body that changes between aft and fore twice an orbit, on a
    schedule that depends on the day, as OCI's tilt does to avoid sun glint. After
    the ascending node, the pitch changes from aft to fore where the satellite's
    geodetic latitude reaches the subsolar latitude at noon (UTC) of the day plus an
    offset that alternates with the parity of the day of the year (UTC); it changes
    back at a fixed argument of latitude. Each change is linear over an arc of the
    orbit, within an interval in which the instruments do not observe.
    """

    fore_pitch_angle: float = Field(
        ..., ge=-90, le=90, description="The pitch angle (degrees) of the fore tilt."
    )
    aft_pitch_angle: float = Field(
        ..., ge=-90, le=90, description="The pitch angle (degrees) of the aft tilt."
    )
    odd_day_offset: float = Field(
        ...,
        ge=-45,
        le=45,
        description="On days of odd day of the year, the latitude (degrees) of the "
        + "change to fore, from the subsolar latitude.",
    )
    even_day_offset: float = Field(
        ...,
        ge=-45,
        le=45,
        description="On days of even day of the year, the latitude (degrees) of the "
        + "change to fore, from the subsolar latitude.",
    )
    transition: float = Field(
        ...,
        gt=0,
        lt=90,
        description="The arc of argument of latitude (degrees) over which each "
        + "change is linear.",
    )
    return_argument_of_latitude: float = Field(
        ...,
        ge=90,
        le=270,
        description="The argument of latitude (degrees) of the middle of the "
        + "change back to aft.",
    )
    blind_before: float = Field(
        ...,
        ge=0,
        le=600,
        description="The time (seconds) before the middle of each change from "
        + "which the instruments do not observe.",
    )
    blind_after: float = Field(
        ...,
        ge=0,
        le=600,
        description="The time (seconds) after the middle of each change at which "
        + "the instruments observe again.",
    )
    provenance: dict[str, Provenance] = Field(
        ..., description="The origin of each value."
    )

    @model_validator(mode="after")
    def _check_provenance(self) -> "TiltSchedule":
        _check_provenance(
            set(type(self).model_fields) - {"provenance"}, self.provenance
        )
        return self

    def get_switch_latitude(self, day: date) -> float:
        """
        Gets the geodetic latitude (degrees) at which the pitch changes from aft to
        fore on a day.

        Args:
            day (datetime.date): The day (UTC).

        Returns:
            float: The latitude (degrees).
        """
        noon = datetime.combine(day, time(12), tzinfo=timezone.utc)
        offset = (
            self.odd_day_offset
            if day.timetuple().tm_yday % 2 == 1
            else self.even_day_offset
        )
        return compute_subsolar_latitude(noon) + offset

    def get_pitch_angle_profile(
        self, day: date, orbit: AllOrbits
    ) -> list[tuple[float, float]]:
        """
        Gets the pitch angle profile of a day, as (argument of latitude, pitch angle)
        pairs.

        Args:
            day (datetime.date): The day (UTC).
            orbit (AllOrbits): The orbit.

        Returns:
            list[tuple[float, float]]: The profile.
        """
        switch = self._get_switch_argument_of_latitude(day, orbit)
        half = self.transition / 2
        back = self.return_argument_of_latitude
        return sorted(
            (float(np.mod(argument, 360)), angle)
            for argument, angle in (
                (switch - half, self.aft_pitch_angle),
                (switch + half, self.fore_pitch_angle),
                (back - half, self.fore_pitch_angle),
                (back + half, self.aft_pitch_angle),
            )
        )

    def get_blind_intervals(
        self, day: date, orbit: AllOrbits
    ) -> list[tuple[float, float]] | None:
        """
        Gets the intervals of argument of latitude (degrees) of a day in which the
        instruments do not observe, around each change, at the orbit's mean rate.

        Args:
            day (datetime.date): The day (UTC).
            orbit (AllOrbits): The orbit.

        Returns:
            list[tuple[float, float]] | None: The (start, end) pairs, or `None` if
                the instruments always observe.
        """
        if self.blind_before + self.blind_after == 0:
            return None
        rate = 360 / orbit.get_orbit_period().total_seconds()
        return sorted(
            (
                float(np.mod(middle - self.blind_before * rate, 360)),
                float(np.mod(middle + self.blind_after * rate, 360)),
            )
            for middle in (
                self._get_switch_argument_of_latitude(day, orbit),
                self.return_argument_of_latitude,
            )
        )

    def _get_switch_argument_of_latitude(self, day: date, orbit: AllOrbits) -> float:
        """
        Gets the argument of latitude (degrees) of the middle of the change from aft
        to fore on a day.
        """
        return compute_ascending_argument_of_latitude(
            self.get_switch_latitude(day),
            orbit.get_mean_altitude(),
            orbit.get_inclination(),
        )

    def build_schedule(
        self,
        orbit: AllOrbits,
        start: datetime,
        end: datetime,
        velocity_frame: VelocityFrame,
        attitude: dict[str, Any],
    ) -> AttitudeSchedule:
        """
        Builds a schedule of attitudes, one for each day (UTC) from the start to the
        end of an analysis.

        Args:
            orbit (AllOrbits): The orbit.
            start (datetime.datetime): The start of the analysis (timezone-aware).
            end (datetime.datetime): The end of the analysis (timezone-aware).
            velocity_frame (VelocityFrame): The velocity frame of the attitude.
            attitude (dict[str, Any]): Other fields of the attitude.

        Returns:
            AttitudeSchedule: The schedule.
        """
        first, last = (t.astimezone(timezone.utc).date() for t in (start, end))
        if last < first:
            raise ValueError("the end of the analysis precedes its start")
        days = [first + timedelta(days=i) for i in range((last - first).days + 1)]
        return AttitudeSchedule.from_attitudes(
            [datetime.combine(day, time(0), tzinfo=timezone.utc) for day in days],
            [
                Attitude(
                    velocity_frame=velocity_frame,
                    pitch_angle_profile=self.get_pitch_angle_profile(day, orbit),
                    blind_intervals=self.get_blind_intervals(day, orbit),
                    **copy.deepcopy(attitude),
                )
                for day in days
            ],
        )


class Configuration(_CatalogModel):
    """
    A way of modeling an instrument for a use, as it was validated.
    """

    grade: Grade = Field(..., description="The use for which it is validated.")
    description: str = Field(..., description="What the configuration represents.")
    platforms: tuple[str, ...] = Field(
        ..., min_length=1, description="The platforms on which it was tested."
    )
    velocity_frame: VelocityFrame = Field(
        ..., description="The velocity frame of the attitude of the space system."
    )
    attitude: dict[str, Any] = Field(
        default_factory=dict,
        description="Other fields of the attitude of the space system.",
    )
    tilt_schedule: TiltSchedule | None = Field(
        default=None,
        description="A pitch of the body that changes by day, from which a "
        + "schedule of attitudes is built for the period of an analysis.",
    )
    instruments: tuple[InstrumentSpecification, ...] = Field(
        ..., min_length=1, description="The instruments."
    )
    instrument_index: int = Field(
        default=0,
        ge=0,
        description="The index of the instrument to analyze by default, of "
        + "`alternatives`.",
    )
    grouping: Literal["union", "alternatives"] | None = Field(
        default=None,
        description="How multiple instruments relate: `union`, parts of one "
        + "view that are analyzed together; or `alternatives`, separately valid "
        + "views that are analyzed one at a time. Only for multiple instruments.",
    )
    limitations: tuple[str, ...] = Field(
        default=(), description="What this configuration does not represent."
    )

    @model_validator(mode="after")
    def _check_instruments(self) -> "Configuration":
        names = [i.parameters["name"] for i in self.instruments]
        if len(set(names)) != len(names):
            # compute_sample_sizes finds instruments by name
            raise ValueError(f"instrument names must be unique: {names}")
        if self.instrument_index >= len(self.instruments):
            raise ValueError("instrument_index is out of range")
        if self.instrument_index and self.grouping != "alternatives":
            raise ValueError("only alternatives have an instrument_index")
        if (self.grouping is None) != (len(self.instruments) == 1):
            raise ValueError("multiple instruments, and only they, need a grouping")
        self.make_attitude()
        if self.tilt_schedule is not None and {
            "pitch_angle",
            "pitch_angle_profile",
            "blind_intervals",
        } & set(self.attitude):
            raise ValueError(
                "a tilt schedule sets the pitch and blind intervals of the attitude"
            )
        return self

    def make_attitude(self) -> Attitude | GlintAttitude | TargetAttitude:
        """
        Makes the attitude of the configuration: of the type its `type` names
        (`attitude` by default, `glint`, or `target`), in its velocity frame.

        Returns:
            Attitude | GlintAttitude | TargetAttitude: The attitude.
        """
        fields = copy.deepcopy(self.attitude)
        kind = _ATTITUDE_TYPES[fields.pop("type", "attitude")]
        return kind(velocity_frame=self.velocity_frame, **fields)

    def _get_measured(self) -> list[Provenance]:
        """
        Gets the provenance of the values fitted to reference data.
        """
        provenances = [
            provenance
            for instrument in self.instruments
            for provenance in instrument.provenance.values()
        ]
        if self.tilt_schedule is not None:
            provenances += list(self.tilt_schedule.provenance.values())
        return [
            provenance for provenance in provenances if provenance.kind == "measured"
        ]

    def get_holdout_gap(self, period: str) -> int | None:
        """
        Gets the days between the data to which values were fitted and the data
        of a test: between the nearest days of the periods.

        Args:
            period (str): The period of the data of the test (see `Period`).

        Returns:
            int | None: The days, or `None` if no value was fitted.
        """
        measured = self._get_measured()
        if not measured:
            return None
        return min(_get_days_apart(p.period, period) for p in measured)

    @property
    def digest(self) -> str:
        """
        A digest of the values that affect analysis results, with which a
        validation record identifies the configuration it tested.
        """
        content = {
            "velocity_frame": self.velocity_frame.value,
            "attitude": self.attitude,
            "instruments": [
                {"type": i.type, "parameters": i.parameters} for i in self.instruments
            ],
        }
        if self.tilt_schedule is not None:
            # keyed as when the field was named `attitude_schedule`, so that the
            # digests of the records stay those of the configurations they tested
            content["attitude_schedule"] = self.tilt_schedule.model_dump(
                mode="json", exclude={"provenance"}
            )
        encoded = json.dumps(content, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(encoded.encode("utf-8")).hexdigest()

    def serves(self, grade: Grade | str) -> bool:
        """
        Checks whether this configuration serves a use: its grade, or the
        coverage grade for a configuration at the timing grade, which also tests
        coverage.

        Args:
            grade (Grade | str): The grade of the use.

        Returns:
            bool: `True`, if it serves the use.
        """
        grade = Grade(grade)
        return self.grade == grade or (
            grade == Grade.COVERAGE and self.grade == Grade.TIMING
        )


class Platform(_CatalogModel):
    """
    A platform that carries an instrument, with its status on a date.
    """

    id: str = Field(..., description="Platform identifier.")
    name: str = Field(..., description="Platform name.")
    norad_id: int = Field(
        ...,
        gt=0,
        description="NORAD catalog number, from which to obtain current orbital "
        + "elements, which are not cataloged because they age.",
    )
    status: Literal["operational", "extended", "retired"] = Field(
        ..., description="The status on the date `as_of`."
    )
    as_of: date = Field(..., description="The date of the status.")


class DataProduct(_CatalogModel):
    """
    A data product of an instrument.
    """

    name: str = Field(..., description="Product name.")
    provider: str = Field(..., description="The organization that distributes it.")
    url: str = Field(..., description="Where to find it.")
    access: Literal["open", "registration"] = Field(
        ..., description="`open`, or free with `registration`."
    )


class Metric(_CatalogModel):
    """
    A value of a metric from a comparison with reference data.
    """

    name: str = Field(..., description="The metric (see `get_gates`).")
    value: float = Field(..., description="The value.")
    scope: str | None = Field(
        default=None,
        description="The part or condition of the test it applies to (such as an "
        + "instrument or a time step), if not all of it. A gate takes the worst "
        + "value of its metric over scopes.",
    )


class Transfer(_CatalogModel):
    """
    The results of a configuration on a platform on which it was not tested,
    which do not count toward its validation.
    """

    period: Period = Field(..., description="The period of the reference data.")
    metrics: tuple[Metric, ...] = Field(..., min_length=1, description="The results.")


class ConfigurationRecord(_CatalogModel):
    """
    The results of validating a configuration, written by a validation script.
    """

    digest: str = Field(..., description="The `Configuration.digest` tested.")
    standard_version: int = Field(
        ..., description="The `STANDARD_VERSION` of the test."
    )
    tatc_version: str = Field(..., description="The version of TAT-C tested.")
    source: str = Field(..., description="The script of the results.")
    reference: str = Field(..., description="The reference data.")
    period: Period = Field(..., description="The period of the reference data.")
    platforms: tuple[str, ...] = Field(
        ..., min_length=1, description="The platforms of the reference data."
    )
    metrics: tuple[Metric, ...] = Field(..., description="The results.")
    transfers: dict[str, Transfer] = Field(
        default_factory=dict,
        description="Results on platforms on which the configuration was not "
        + "tested, by platform.",
    )
    fixture_metrics: dict[str, float] = Field(
        default_factory=dict,
        description="Results from the reduced test fixture, which detect "
        + "regressions but do not themselves establish agreement with the "
        + "reference data.",
    )


class ValidationRecord(_CatalogModel):
    """
    The results of validating the configurations of an entry.
    """

    entry_id: str = Field(..., description="The identifier of the entry.")
    configurations: dict[str, ConfigurationRecord] = Field(
        ..., description="The results of each configuration."
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
    missing_grades: tuple[Grade, ...]

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


@dataclass(frozen=True)
class InstrumentPayload:
    """
    Instruments with the attitude that they were validated with, which must be
    used together. Analyze them with the `instrument_index`: of the instrument
    to analyze, or `None` to analyze the instruments together (as the parts of
    one view).
    """

    instruments: list[Instrument]
    attitude: AllAttitudes | AttitudeSchedule | None
    instrument_index: int | None
    references: tuple[str, ...]
    """The configurations, as `entry:configuration`."""
    transferred: bool
    """Whether a configuration was not tested on the platform modeled."""
    limitations: tuple[str, ...]

    def to_satellite(self, name: str, orbit: AllOrbits) -> Satellite:
        """
        Builds a satellite that carries the instruments with their attitude.
        Analyze it with the `instrument_index` of this payload.

        Args:
            name (str): The satellite's name.
            orbit (AllOrbits): The orbit.

        Returns:
            Satellite: The satellite.
        """
        return Satellite(
            name=name, orbit=orbit, instruments=self.instruments, attitude=self.attitude
        )


class InstrumentEntry(_CatalogModel):
    """
    An instrument of a past or present mission, in configurations validated
    against its data products.
    """

    id: str = Field(..., description="Entry identifier.")
    name: str = Field(..., description="Instrument name.")
    aliases: tuple[str, ...] = Field(default=(), description="Other names.")
    instrument_class: InstrumentClass = Field(..., description="Instrument class.")
    selection: str = Field(..., description="Why the instrument is cataloged.")
    platforms: tuple[Platform, ...] = Field(
        ..., min_length=1, description="The platforms that carry it."
    )
    data_products: tuple[DataProduct, ...] = Field(
        ..., min_length=1, description="Its principal data products."
    )
    default_configuration: str = Field(
        ..., description="The configuration used if none is specified."
    )
    configurations: dict[str, Configuration] = Field(
        ..., min_length=1, description="The ways of modeling it, by identifier."
    )
    limitations: tuple[str, ...] = Field(
        default=(), description="What no configuration represents."
    )
    record: ValidationRecord | None = Field(
        default=None, description="The results of its validation."
    )

    @model_validator(mode="after")
    def _check_references(self) -> "InstrumentEntry":
        if self.default_configuration not in self.configurations:
            raise ValueError("default_configuration is not a configuration")
        platforms = {platform.id for platform in self.platforms}
        for key, configuration in self.configurations.items():
            unknown = set(configuration.platforms) - platforms
            if unknown:
                raise ValueError(f"{key}: unknown platforms {sorted(unknown)}")
            # checks that the class has gates at this grade
            get_gates(self.instrument_class, configuration.grade)
        if self.record is not None:
            if self.record.entry_id != self.id:
                raise ValueError("the record belongs to a different entry")
            unknown = set(self.record.configurations) - set(self.configurations)
            if unknown:
                raise ValueError(f"records of unknown {sorted(unknown)}")
            for key, record in self.record.configurations.items():
                tested = set(self.configurations[key].platforms)
                if not set(record.platforms) <= tested:
                    raise ValueError(
                        f"{key}: the record's platforms {list(record.platforms)} "
                        f"are not those tested, {sorted(tested)}"
                    )
                transferred = set(record.transfers)
                if transferred & tested or not transferred <= platforms:
                    raise ValueError(
                        f"{key}: transfers must be to other platforms of the entry"
                    )
        return self

    def _evaluate_configuration(self, key: str) -> ConfigurationValidation:
        """
        Applies the standard to a configuration.

        Args:
            key (str): The identifier of the configuration.

        Returns:
            ConfigurationValidation: Its validation status.
        """
        configuration = self.configurations[key]
        gates = get_gates(self.instrument_class, configuration.grade)
        record = None if self.record is None else self.record.configurations.get(key)
        if record is None or not record.metrics:
            return ConfigurationValidation(
                Tier.NOMINAL,
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
        if record.digest != configuration.digest:
            gaps.append("the configuration changed after it was tested")
        if not is_compatible(self.instrument_class, record.standard_version):
            gaps.append(
                f"version {record.standard_version} of the standard, of the test, "
                f"does not apply to the {self.instrument_class.value} class"
            )
        holdout = configuration.get_holdout_gap(record.period)
        if holdout is not None and holdout < MIN_HOLDOUT_GAP_DAYS:
            gaps.append(
                "fitted values are not tested on data at least "
                f"{MIN_HOLDOUT_GAP_DAYS} days from their fitting data"
            )
        passed = not gaps and all(result.passed for result in results)
        return ConfigurationValidation(
            Tier.VALIDATED if passed else Tier.PARTIAL,
            tuple(results),
            tuple(gaps),
            record,
            holdout,
        )

    @property
    def validation(self) -> EntryValidation:
        """
        The validation status of this entry, computed from its records and the
        standard rather than declared. An entry is only as strong as its weakest
        configuration, and needs configurations at the grades of its class.
        """
        configurations = {
            key: self._evaluate_configuration(key) for key in self.configurations
        }
        missing = tuple(
            grade
            for grade in get_required_grades(self.instrument_class)
            if not any(c.serves(grade) for c in self.configurations.values())
        )
        tier = min((c.tier for c in configurations.values()), key=lambda t: t.rank)
        if tier == Tier.VALIDATED and missing:
            tier = Tier.PARTIAL
        if tier == Tier.NOMINAL and any(
            c.tier != Tier.NOMINAL for c in configurations.values()
        ):
            tier = Tier.PARTIAL
        return EntryValidation(tier, configurations, missing)

    def select_configuration(
        self,
        configuration: str | None = None,
        *,
        grade: Grade | str | None = None,
        platform: str | None = None,
    ) -> str:
        """
        Selects a configuration: the one named, or else one that serves a use
        (see `Configuration.serves`) and was tested on a platform, preferring
        the default configuration and then those at the grade of the default.
        Without a configuration tested on the platform, it selects among the
        others (a transfer, on which `build` warns).

        Args:
            configuration (str | None): The identifier of the configuration, or
                `None` to select one.
            grade (Grade | str | None): The grade of the use, if any.
            platform (str | None): The identifier of the platform modeled, if
                any.

        Returns:
            str: The identifier of the configuration.

        Raises:
            ValueError: If the configuration is unknown, or none serves the use.
        """
        if configuration is not None:
            if configuration not in self.configurations:
                raise ValueError(
                    f"unknown configuration {configuration!r} of {self.id!r}: "
                    f"{sorted(self.configurations)}"
                )
            if grade is not None and not self.configurations[configuration].serves(
                grade
            ):
                raise ValueError(
                    f"configuration {configuration!r} of {self.id!r} is at the "
                    f"{self.configurations[configuration].grade.value} grade"
                )
            return configuration
        candidates = [
            key
            for key, c in self.configurations.items()
            if grade is None or c.serves(grade)
        ]
        if not candidates:
            raise ValueError(
                f"{self.id!r} has no configuration at the {Grade(grade).value} grade"
            )
        tested = [k for k in candidates if platform in self.configurations[k].platforms]
        candidates = tested or candidates
        default = self.configurations[self.default_configuration]
        return min(
            candidates,
            key=lambda k: (
                k != self.default_configuration,
                self.configurations[k].grade != default.grade,
            ),
        )

    def build(
        self,
        *,
        configuration: str | None = None,
        grade: Grade | str | None = None,
        orbit: AllOrbits | None = None,
        time_step: float | None = None,
        attitude: AllAttitudes | AttitudeSchedule | None = None,
        platform: str | None = None,
        start: datetime | None = None,
        end: datetime | None = None,
    ) -> InstrumentPayload:
        """
        Builds the instruments of a configuration with the attitude they were
        validated with. Each call builds new objects.

        Args:
            configuration (str | None): The identifier of the configuration, or
                `None` to select one by grade and platform (see
                `select_configuration`), by default the default configuration.
            grade (Grade | str | None): The grade of the use, if any.
            orbit (AllOrbits | None): The orbit, if the configuration tailors a
                field of view to a time step or schedules its tilt.
            time_step (float | None): The analysis time step (seconds), if the
                configuration tailors a field of view to it.
            attitude (AllAttitudes | AttitudeSchedule | None): An attitude to use
                instead of the validated one (for example, a schedule of yaw
                flips), in the same velocity frame.
            platform (str | None): The identifier of the platform modeled, if
                any, which selects the configuration tested on it. Using a
                configuration on a platform on which it was not tested gives a
                warning.
            start (datetime.datetime | None): The start of the analysis
                (timezone-aware), if the configuration schedules its tilt by
                day.
            end (datetime.datetime | None): The end of the analysis
                (timezone-aware), if the configuration schedules its tilt by
                day.

        Returns:
            InstrumentPayload: The instruments, attitude, and analysis settings.

        Raises:
            ValueError: If the configuration is unknown or does not serve the
                grade, tailoring lacks an orbit or time step, a schedule lacks
                an orbit or the period of the analysis, or an attitude has
                another velocity frame.
        """
        key = self.select_configuration(configuration, grade=grade, platform=platform)
        selected = self.configurations[key]
        tailored = None
        if any(instrument.is_tailored for instrument in selected.instruments):
            if orbit is None or time_step is None:
                raise ValueError(
                    f"configuration {key!r} of {self.id!r} tailors its along-track "
                    "field of view, which needs an orbit and a time step"
                )
            tailored = compute_along_track_field_of_view(
                orbit.get_mean_altitude(),
                time_step,
                orbit.get_inclination(),
                eccentricity=orbit.get_eccentricity(),
            )
        transferred = platform is not None and platform not in selected.platforms
        if transferred:
            warnings.warn(
                f"configuration {key!r} of {self.id!r} was tested on "
                f"{list(selected.platforms)}, not on {platform!r}",
                stacklevel=2,
            )
        if attitude is None and selected.tilt_schedule is not None:
            if orbit is None or start is None or end is None:
                raise ValueError(
                    f"configuration {key!r} of {self.id!r} schedules its tilt by "
                    "day, which needs an orbit and the start and end of the analysis"
                )
            attitude = selected.tilt_schedule.build_schedule(
                orbit, start, end, selected.velocity_frame, selected.attitude
            )
        return InstrumentPayload(
            instruments=[
                instrument.construct_instrument(tailored)
                for instrument in selected.instruments
            ],
            attitude=self._build_attitude(selected, attitude),
            instrument_index=(
                None if selected.grouping == "union" else selected.instrument_index
            ),
            references=(f"{self.id}:{key}",),
            transferred=transferred,
            limitations=(*self.limitations, *selected.limitations),
        )

    @staticmethod
    def _build_attitude(
        selected: Configuration, attitude: AllAttitudes | AttitudeSchedule | None
    ) -> AllAttitudes | AttitudeSchedule | None:
        """
        Gets the attitude of a configuration, or checks one that replaces it.

        Args:
            selected (Configuration): The configuration.
            attitude (AllAttitudes | AttitudeSchedule | None): The replacement,
                or `None`.

        Returns:
            AllAttitudes | AttitudeSchedule | None: The attitude, or `None` for a
            body aligned with the nadir frame (as for an Earth-fixed frame
            without other fields).

        Raises:
            ValueError: If the replacement has another velocity frame.
        """
        if attitude is not None:
            attitudes = (
                [segment.attitude for segment in attitude.segments]
                if isinstance(attitude, AttitudeSchedule)
                else [attitude]
            )
            frames = {a.velocity_frame for a in attitudes}
            if frames != {selected.velocity_frame}:
                raise ValueError(
                    "the attitude must have the validated velocity frame: "
                    f"{selected.velocity_frame.value}"
                )
            return attitude
        if (
            selected.velocity_frame == VelocityFrame.EARTH_FIXED
            and not selected.attitude
        ):
            return None
        return selected.make_attitude()

    def to_satellite(
        self,
        orbit: AllOrbits,
        *,
        configuration: str | None = None,
        grade: Grade | str | None = None,
        time_step: float | None = None,
        attitude: AllAttitudes | AttitudeSchedule | None = None,
        platform: str | None = None,
        name: str | None = None,
        start: datetime | None = None,
        end: datetime | None = None,
    ) -> Satellite:
        """
        Builds a satellite that carries the instruments of a configuration (see
        `build`). For several instruments, `build(...).to_satellite(name, orbit)`
        also gives the `instrument_index` to analyze them with, which is 0 for
        one instrument.

        Args:
            orbit (AllOrbits): The orbit.
            configuration (str | None): The identifier of the configuration, or
                `None` to select one by grade and platform.
            grade (Grade | str | None): The grade of the use, if any.
            time_step (float | None): The analysis time step (seconds), if the
                configuration tailors a field of view to it.
            attitude (AllAttitudes | AttitudeSchedule | None): An attitude to use
                instead of the validated one.
            platform (str | None): The identifier of the platform modeled.
            name (str | None): The satellite's name, or `None` for the
                instrument's.
            start (datetime.datetime | None): The start of the analysis, if the
                configuration schedules its tilt by day.
            end (datetime.datetime | None): The end of the analysis, if the
                configuration schedules its tilt by day.

        Returns:
            Satellite: The satellite.
        """
        return self.build(
            configuration=configuration,
            grade=grade,
            orbit=orbit,
            time_step=time_step,
            attitude=attitude,
            platform=platform,
            start=start,
            end=end,
        ).to_satellite(self.name if name is None else name, orbit)

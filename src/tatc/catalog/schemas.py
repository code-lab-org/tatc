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
from typing import Any, Literal

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, model_validator
from skyfield.framelib import itrs

from .. import constants
from ..schemas import (
    Attitude,
    AttitudeSchedule,
    ConicalInstrument,
    Instrument,
    PointedInstrument,
    Satellite,
)
from ..schemas.space.base import AllAttitudes
from ..utils import VelocityFrame, compute_along_track_field_of_view
from .standard import (
    MIN_HOLDOUT_GAP_DAYS,
    Gate,
    Grade,
    InstrumentClass,
    Tier,
    get_gates,
    get_required_grades,
)

_INSTRUMENT_CLASSES: dict[str, type[Instrument]] = {
    "instrument": Instrument,
    "pointed": PointedInstrument,
    "conical": ConicalInstrument,
}
_TAILORED = {"tailor": "time_step"}
_TAILORED_FIELD = "along_track_field_of_view"


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
    platform: str | None = Field(
        default=None, description="The platform to which a measured value was fitted."
    )
    period: str | None = Field(
        default=None, description="The period to which a measured value was fitted."
    )

    @model_validator(mode="after")
    def _check_measured(self) -> "Provenance":
        if self.kind == "measured" and (self.platform is None or self.period is None):
            raise ValueError("a measured value needs its platform and period")
        return self


class InstrumentSpecification(_CatalogModel):
    """
    The class and parameters of an instrument, from which objects are constructed.
    """

    type: Literal["instrument", "pointed", "conical"] = Field(
        ..., description="The instrument class."
    )
    parameters: dict[str, Any] = Field(..., description="The instrument's fields.")
    provenance: dict[str, Provenance] = Field(
        default_factory=dict,
        description="The origin of each parameter, except for its name.",
    )

    @model_validator(mode="after")
    def _check_parameters(self) -> "InstrumentSpecification":
        fields = set(_INSTRUMENT_CLASSES[self.type].model_fields) - {"type"}
        if "velocity_frame" in self.parameters:
            raise ValueError(
                "velocity_frame belongs to the configuration, as it sets the "
                "attitude of the space system"
            )
        unknown = set(self.parameters) - fields
        if unknown:
            raise ValueError(f"unknown {self.type} parameters: {sorted(unknown)}")
        described = set(self.parameters) - {"name"}
        if described != set(self.provenance):
            raise ValueError(
                "each parameter needs provenance: "
                f"{sorted(described.symmetric_difference(self.provenance))}"
            )
        if "name" not in self.parameters:
            raise ValueError("an instrument specification needs a name")
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


def get_subsolar_latitude(at: datetime) -> float:
    """
    Gets the subsolar latitude (degrees): the geocentric latitude of the Sun's
    apparent direction, or its declination of date.

    Args:
        at (datetime.datetime): The time (timezone-aware).

    Returns:
        float: The subsolar latitude (degrees).
    """
    t = constants.timescale.from_datetime(at)
    x, y, z = (
        constants.de421["earth"]
        .at(t)
        .observe(constants.de421["sun"])
        .apparent()
        .frame_xyz(itrs)
        .m
    )
    return float(np.degrees(np.arctan2(z, np.hypot(x, y))))


def _get_ascending_argument_of_latitude(
    latitude: float, altitude: float, inclination: float
) -> float:
    """
    Gets the argument of latitude (degrees) at which a circular orbit, ascending,
    passes over a geodetic latitude.

    Args:
        latitude (float): The geodetic latitude (degrees) of the satellite.
        altitude (float): The satellite's altitude (meters).
        inclination (float): The orbit's inclination (degrees).

    Returns:
        float: The argument of latitude (degrees, from 0 to 360).
    """
    phi = np.radians(latitude)
    e2 = constants.EARTH_ECCENTRICITY**2
    normal = constants.EARTH_EQUATORIAL_RADIUS / np.sqrt(1 - e2 * np.sin(phi) ** 2)
    geocentric = np.arctan2(
        (normal * (1 - e2) + altitude) * np.sin(phi), (normal + altitude) * np.cos(phi)
    )
    ratio = np.sin(geocentric) / np.sin(np.radians(inclination))
    if abs(ratio) > 1:
        raise ValueError(
            f"an orbit inclined by {inclination} deg does not reach latitude {latitude}"
        )
    return float(np.mod(np.degrees(np.arcsin(ratio)), 360))


class TiltSchedule(_CatalogModel):
    """
    A pitch of the body that changes between aft and fore twice an orbit, on a
    schedule that depends on the day, as OCI's tilt does to avoid sun glint. After
    the ascending node, the pitch changes from aft to fore where the satellite's
    geodetic latitude reaches the subsolar latitude at noon (UTC) of the day plus an
    offset that alternates with the parity of the day of the year (UTC); it changes
    back at a fixed argument of latitude. Each change is linear over an arc of the
    orbit.
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
    provenance: dict[str, Provenance] = Field(
        ..., description="The origin of each value."
    )

    @model_validator(mode="after")
    def _check_provenance(self) -> "TiltSchedule":
        described = set(type(self).model_fields) - {"provenance"}
        if described != set(self.provenance):
            raise ValueError(
                "each value needs provenance: "
                f"{sorted(described.symmetric_difference(self.provenance))}"
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
        return get_subsolar_latitude(noon) + offset

    def get_pitch_angle_profile(
        self, day: date, orbit: Any
    ) -> list[tuple[float, float]]:
        """
        Gets the pitch angle profile of a day, as (argument of latitude, pitch angle)
        pairs.

        Args:
            day (datetime.date): The day (UTC).
            orbit (Any): The orbit.

        Returns:
            list[tuple[float, float]]: The profile.
        """
        switch = _get_ascending_argument_of_latitude(
            self.get_switch_latitude(day),
            orbit.get_mean_altitude(),
            orbit.get_inclination(),
        )
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

    def build_schedule(
        self,
        orbit: Any,
        start: datetime,
        end: datetime,
        velocity_frame: VelocityFrame,
        attitude: dict[str, Any],
    ) -> AttitudeSchedule:
        """
        Builds a schedule of attitudes, one for each day (UTC) from the start to the
        end of an analysis.

        Args:
            orbit (Any): The orbit.
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
    attitude_schedule: TiltSchedule | None = Field(
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
        description="The index of the instrument to analyze by default.",
    )
    grouping: Literal["single", "union", "alternatives"] = Field(
        default="single",
        description="How multiple instruments relate: `union`, parts of one "
        + "view that are analyzed together; or `alternatives`, separately valid "
        + "views that are analyzed one at a time.",
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
        if (self.grouping == "single") != (len(self.instruments) == 1):
            raise ValueError("grouping is `single` only for one instrument")
        Attitude(velocity_frame=self.velocity_frame, **self.attitude)
        if self.attitude_schedule is not None and {
            "pitch_angle",
            "pitch_angle_profile",
        } & set(self.attitude):
            raise ValueError("an attitude schedule sets the pitch of the attitude")
        return self

    @property
    def is_fitted(self) -> bool:
        """
        Whether any parameter was fitted to reference data.
        """
        provenances = [
            provenance
            for instrument in self.instruments
            for provenance in instrument.provenance.values()
        ]
        if self.attitude_schedule is not None:
            provenances += list(self.attitude_schedule.provenance.values())
        return any(provenance.kind == "measured" for provenance in provenances)

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
        if self.attitude_schedule is not None:
            content["attitude_schedule"] = self.attitude_schedule.model_dump(
                mode="json", exclude={"provenance"}
            )
        encoded = json.dumps(content, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


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
        default=None, description="What it applies to, if not the whole configuration."
    )


class ConfigurationRecord(_CatalogModel):
    """
    The results of validating a configuration, written by a validation script.
    """

    digest: str = Field(..., description="The `Configuration.digest` tested.")
    evidence: Literal["rerun", "reported"] = Field(
        ...,
        description="`rerun`, computed with the configuration's exact values "
        + "through public methods; or `reported`, taken from the outputs of a "
        + "notebook that did not use them exactly.",
    )
    tatc_version: str = Field(..., description="The version of TAT-C tested.")
    source: str = Field(..., description="The notebook or script of the results.")
    reference: str = Field(..., description="The reference data.")
    period: str = Field(..., description="The period of the reference data.")
    platform: str = Field(..., description="The platform of the reference data.")
    holdout_gap_days: float | None = Field(
        default=None,
        description="For fitted values, the days between their fitting data and "
        + "the reference data of these metrics, or `None` if tested on the "
        + "fitting data.",
    )
    metrics: tuple[Metric, ...] = Field(..., description="The results.")
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
    standard_version: int = Field(..., description="The `STANDARD_VERSION` applied.")
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


@dataclass(frozen=True)
class ConfigurationValidation:
    """
    The validation status of a configuration.
    """

    tier: Tier
    gates: tuple[GateResult, ...]
    gaps: tuple[str, ...]
    record: ConfigurationRecord | None


@dataclass(frozen=True)
class EntryValidation:
    """
    The validation status of an entry, computed from its records and the standard.
    """

    tier: Tier
    configurations: dict[str, ConfigurationValidation]
    gaps: tuple[str, ...]


@dataclass(frozen=True)
class InstrumentPayload:
    """
    Instruments with the attitude and analysis settings that they were validated
    with, which must be used together.
    """

    instruments: list[Instrument]
    attitude: Any
    instrument_index: int
    grouping: str
    transferred: bool
    limitations: tuple[str, ...]


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
        gaps = [
            f"{result.gate.metric}: {result.value} fails the limit of "
            f"{result.gate.limit}"
            for result in results
            if result.passed is False
        ]
        gaps += [
            f"{result.gate.metric}: not tested"
            for result in results
            if result.passed is None
        ]
        if record.digest != configuration.digest:
            gaps.append("the configuration changed after it was tested")
        if record.evidence != "rerun":
            gaps.append("results are not reproduced with its exact values")
        if configuration.is_fitted and (
            record.holdout_gap_days is None
            or record.holdout_gap_days < MIN_HOLDOUT_GAP_DAYS
        ):
            gaps.append(
                "fitted values are not tested on data at least "
                f"{MIN_HOLDOUT_GAP_DAYS} days from their fitting data"
            )
        return ConfigurationValidation(
            Tier.PARTIAL if gaps else Tier.VALIDATED,
            tuple(results),
            tuple(gaps),
            record,
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
        grades = {configuration.grade for configuration in self.configurations.values()}
        if Grade.TIMING in grades:
            grades.add(Grade.COVERAGE)
        gaps = [
            f"no configuration at the {grade.value} grade"
            for grade in get_required_grades(self.instrument_class)
            if grade not in grades
        ]
        order = [Tier.NOMINAL, Tier.PARTIAL, Tier.VALIDATED]
        tier = min((c.tier for c in configurations.values()), key=order.index)
        if tier == Tier.VALIDATED and gaps:
            tier = Tier.PARTIAL
        if tier == Tier.NOMINAL and any(
            c.tier != Tier.NOMINAL for c in configurations.values()
        ):
            tier = Tier.PARTIAL
        return EntryValidation(tier, configurations, tuple(gaps))

    def build(
        self,
        configuration: str | None = None,
        *,
        orbit: Any = None,
        time_step: float | None = None,
        attitude: Any = None,
        platform: str | None = None,
        start: datetime | None = None,
        end: datetime | None = None,
    ) -> InstrumentPayload:
        """
        Builds the instruments of a configuration with the attitude they were
        validated with. Each call builds new objects.

        Args:
            configuration (str | None): The identifier of the configuration, or
                `None` for the default.
            orbit (Any): The orbit, if the configuration tailors a field of view
                to a time step or schedules its attitude.
            time_step (float | None): The analysis time step (seconds), if the
                configuration tailors a field of view to it.
            attitude (Any): An attitude to use instead of the validated one (for
                example, a schedule of yaw flips), in the same velocity frame.
            platform (str | None): The identifier of the platform being modeled,
                if any. Using a configuration on a platform on which it was not
                tested gives a warning.
            start (datetime.datetime | None): The start of the analysis
                (timezone-aware), if the configuration schedules its attitude
                by day.
            end (datetime.datetime | None): The end of the analysis
                (timezone-aware), if the configuration schedules its attitude
                by day.

        Returns:
            InstrumentPayload: The instruments, attitude, and analysis settings.

        Raises:
            ValueError: If the configuration is unknown, tailoring lacks an
                orbit or time step, a schedule lacks an orbit or the period of
                the analysis, or an attitude has another velocity frame.
        """
        key = self.default_configuration if configuration is None else configuration
        if key not in self.configurations:
            raise ValueError(
                f"unknown configuration {key!r} of {self.id!r}: "
                f"{sorted(self.configurations)}"
            )
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
        if (
            selected.attitude_schedule is not None
            and attitude is None
            and (orbit is None or start is None or end is None)
        ):
            raise ValueError(
                f"configuration {key!r} of {self.id!r} schedules its attitude by "
                "day, which needs an orbit and the start and end of the analysis"
            )
        if attitude is None and selected.attitude_schedule is not None:
            attitude = selected.attitude_schedule.build_schedule(
                orbit, start, end, selected.velocity_frame, selected.attitude
            )
        return InstrumentPayload(
            instruments=[
                instrument.construct_instrument(tailored)
                for instrument in selected.instruments
            ],
            attitude=self._build_attitude(selected, attitude),
            instrument_index=selected.instrument_index,
            grouping=selected.grouping,
            transferred=transferred,
            limitations=(*self.limitations, *selected.limitations),
        )

    @staticmethod
    def _build_attitude(selected: Configuration, attitude: Any) -> AllAttitudes | None:
        """
        Gets the attitude of a configuration, or checks one that replaces it.

        Args:
            selected (Configuration): The configuration.
            attitude (Any): The replacement, or `None`.

        Returns:
            AllAttitudes | None: The attitude, or `None` for a body aligned with
            the nadir frame (as for an Earth-fixed frame without other fields).

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
        return Attitude(
            velocity_frame=selected.velocity_frame,
            **copy.deepcopy(selected.attitude),
        )

    def to_satellite(
        self,
        orbit: Any,
        configuration: str | None = None,
        *,
        time_step: float | None = None,
        attitude: Any = None,
        platform: str | None = None,
        name: str | None = None,
        start: datetime | None = None,
        end: datetime | None = None,
    ) -> Satellite:
        """
        Builds a satellite that carries the instruments of a configuration. To
        analyze an instrument other than the first, use the `instrument_index`
        of `build` in analysis methods.

        Args:
            orbit (Any): The orbit.
            configuration (str | None): The identifier of the configuration, or
                `None` for the default.
            time_step (float | None): The analysis time step (seconds), if the
                configuration tailors a field of view to it.
            attitude (Any): An attitude to use instead of the validated one.
            platform (str | None): The identifier of the platform being modeled.
            name (str | None): The satellite's name, or `None` for the
                instrument's.
            start (datetime.datetime | None): The start of the analysis, if the
                configuration schedules its attitude by day.
            end (datetime.datetime | None): The end of the analysis, if the
                configuration schedules its attitude by day.

        Returns:
            Satellite: The satellite.
        """
        payload = self.build(
            configuration,
            orbit=orbit,
            time_step=time_step,
            attitude=attitude,
            platform=platform,
            start=start,
            end=end,
        )
        return Satellite(
            name=self.name if name is None else name,
            orbit=orbit,
            instruments=payload.instruments,
            attitude=payload.attitude,
        )

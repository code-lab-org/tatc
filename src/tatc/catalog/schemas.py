"""
Object schemas for catalog entries of validated instruments, and the payloads
built from them. The schemas of their validation records are in
`tatc.catalog.records`.

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
from pydantic import Field, model_validator

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
from .records import Period, ValidationRecord, _CatalogModel, _get_days_apart
from .standard import (
    EntryValidation,
    Grade,
    InstrumentClass,
    evaluate_entry,
    get_gates,
    serves,
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
        description="Other fields of the attitude of the space system: that of "
        + "its body, which every instrument it carries shares, so that "
        + "configurations carried together (see `build_payload`) must have the "
        + "same. Angles that orient one instrument on the body, such as the yaw "
        + "of its scan line, are its instrument's parameters (`yaw_angle`, "
        + "`roll_angle`, and `pitch_angle` of a pointed instrument).",
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
        Checks whether this configuration serves a use (see
        `tatc.catalog.standard.serves`).

        Args:
            grade (Grade | str): The grade of the use.

        Returns:
            bool: `True`, if it serves the use.
        """
        return serves(self.grade, grade)


def _get_platform_suffix(configuration: Configuration) -> str | None:
    """
    Gets the suffix of the platform of a configuration tested on one platform
    (such as `_noaa_20`), or `None` for several platforms.
    """
    if len(configuration.platforms) != 1:
        return None
    return "_" + configuration.platforms[0].replace("-", "_")


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

    def to_satellite(self, orbit: AllOrbits, *, name: str) -> Satellite:
        """
        Builds a satellite that carries the instruments with their attitude.
        Analyze it with the `instrument_index` of this payload.

        Args:
            orbit (AllOrbits): The orbit.
            name (str): The satellite's name.

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
        if len({c.platforms for c in self.configurations.values()}) > 1:
            # references may name a configuration without its platform
            for key, configuration in self.configurations.items():
                suffix = _get_platform_suffix(configuration)
                if suffix is None or not key.endswith(suffix):
                    raise ValueError(
                        f"{key}: in an entry with configurations tested on "
                        "different platforms, each is tested on one platform and "
                        "its name ends with it (such as swath_noaa_20)"
                    )
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

    @property
    def validation(self) -> EntryValidation:
        """
        The validation status of this entry, computed from its records and the
        standard rather than declared (see `tatc.catalog.standard.evaluate_entry`).
        """
        return evaluate_entry(self)

    def get_variants(self, name: str) -> list[str]:
        """
        Gets the configurations that a name refers to: the configuration of that
        identifier, or else those whose identifier is the name followed by their
        platform (such as `swath_noaa_20` and `swath_noaa_21` for `swath`), in
        an entry with configurations tested on different platforms.

        Args:
            name (str): The name.

        Returns:
            list[str]: The identifiers of the configurations, empty if none.
        """
        if name in self.configurations:
            return [name]
        if len({c.platforms for c in self.configurations.values()}) == 1:
            return []
        return [
            key
            for key, c in self.configurations.items()
            if key == f"{name}{_get_platform_suffix(c)}"
        ]

    def select_configuration(
        self,
        configuration: str | None = None,
        *,
        grade: Grade | str | None = None,
        platform: str | None = None,
    ) -> str:
        """
        Selects a configuration: one that serves a use (see `Configuration.serves`)
        and was tested on a platform, preferring the default configuration and then
        those at the grade of the default. A name selects among the configurations
        it refers to (see `get_variants`), so that `swath` selects `swath_noaa_21`
        for the platform `noaa-21`, and otherwise prefers the platforms of the
        default configuration. Without a configuration tested on the platform, it
        selects among the others (a transfer, on which `build` warns).

        Args:
            configuration (str | None): The identifier of the configuration, or a
                name of its variants (see `get_variants`), or `None` to select
                any configuration.
            grade (Grade | str | None): The grade of the use, if any.
            platform (str | None): The identifier of the platform modeled, if
                any.

        Returns:
            str: The identifier of the configuration.

        Raises:
            ValueError: If the configuration is unknown, or none serves the use.
        """
        default = self.configurations[self.default_configuration]
        if configuration is not None:
            variants = self.get_variants(configuration)
            if not variants:
                raise ValueError(
                    f"unknown configuration {configuration!r} of {self.id!r}: "
                    f"{sorted(self.configurations)}"
                )
            tested = [
                k for k in variants if platform in self.configurations[k].platforms
            ]
            key = min(
                tested or variants,
                key=lambda k: self.configurations[k].platforms != default.platforms,
            )
            if grade is not None and not self.configurations[key].serves(grade):
                raise ValueError(
                    f"configuration {key!r} of {self.id!r} is at the "
                    f"{self.configurations[key].grade.value} grade"
                )
            return key
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
            configuration (str | None): The identifier of the configuration, or a
                name of its variants (such as `swath` for `swath_noaa_20`), or
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
        `build`), which is analyzed with an `instrument_index` of 0, as the
        analysis functions do by default. A configuration analyzed otherwise (the
        parts of one view, or alternatives that analyze another instrument by
        default) needs `build(...).to_satellite(orbit, name=name)` and the
        payload's `instrument_index`.

        Args:
            orbit (AllOrbits): The orbit.
            configuration (str | None): The identifier of the configuration, or a
                name of its variants, or `None` to select one by grade and
                platform (see `select_configuration`).
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

        Raises:
            ValueError: If the configuration is not analyzed with an
                `instrument_index` of 0.
        """
        payload = self.build(
            configuration=configuration,
            grade=grade,
            orbit=orbit,
            time_step=time_step,
            attitude=attitude,
            platform=platform,
            start=start,
            end=end,
        )
        if payload.instrument_index != 0:
            raise ValueError(
                f"{payload.references[0]} is analyzed with instrument_index="
                f"{payload.instrument_index}, which a satellite does not carry: use "
                "build(...).to_satellite(orbit, name=name) and the payload's "
                "instrument_index"
            )
        return payload.to_satellite(orbit, name=self.name if name is None else name)

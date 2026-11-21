"""
Object schemas for catalog systems: satellite programs, radar networks, and
surface networks, in configurations validated against reference data (such as
the element sets of satellites, or the metadata of stations), and the platforms
built from them. Their records are those of instrument entries
(`tatc.catalog.records`), and their standard is `tatc.catalog.system_standard`.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import functools
import hashlib
import json
import pathlib
from collections.abc import Mapping, Sequence
from datetime import date, datetime, time
from typing import Annotated, Any, ClassVar, Literal

import numpy as np
import pandas as pd
import shapely
from pydantic import Field, TypeAdapter, model_validator
from shapely.geometry.base import BaseGeometry

from ..schemas import (
    CircularOrbit,
    GeosynchronousOrbit,
    InSituInstrument,
    RadarStation,
    RepeatGroundTrackOrbit,
    Satellite,
    SunSynchronousOrbit,
    SurfacePlatform,
)
from ..schemas.orbit import AllOrbits
from ..schemas.surface.radar import RadarBand
from ..utils import compute_geoid_height
from .records import ValidationRecord, _CatalogModel, _get_days_apart
from .schemas import InstrumentSpecification, Provenance, _check_provenance
from .standard import EntryValidation, InstrumentClass
from .system_standard import (
    SystemClass,
    SystemGrade,
    evaluate_system,
    get_system_gates,
    serves_system,
)

_TABLES = pathlib.Path(__file__).parents[1] / "resources" / "systems" / "tables"


class Spacecraft(_CatalogModel):
    """
    A satellite of a system, with its status on a date. Instrument entries refer
    to the satellites that carry them by identifier.
    """

    id: str = Field(..., description="Spacecraft identifier.")
    name: str = Field(..., description="Spacecraft name.")
    norad_id: int = Field(
        ...,
        gt=0,
        description="NORAD catalog number, from which to obtain current orbital "
        + "elements, which are not cataloged because they age (the catalog's "
        + "orbits are idealized).",
    )
    status: Literal["operational", "extended", "retired"] = Field(
        ..., description="The status on the date `as_of`."
    )
    as_of: date = Field(..., description="The date of the status.")


class SunSynchronousPlane(_CatalogModel):
    """
    The plane of a sun-synchronous orbit, at the plane grade: its mean altitude
    and the local time of its ascending node, and, for an orbit whose local time
    is no longer maintained (and drifts), its inclination. Built at an epoch as
    a `SunSynchronousOrbit` (whose inclination is sun-synchronous, so that its
    local time stays constant), or, with an inclination, as a `CircularOrbit`
    with the local time at the epoch, which drifts at the rate of its inclination
    during an analysis. Its phase is that of a constellation whose operator
    maintains the satellites' relative phases, or else is spaced evenly among
    the satellites built (see `OrbitConfiguration.build_orbits`).
    """

    type: Literal["sso"] = Field(default="sso", description="Orbit discriminator.")
    mean_altitude: float = Field(
        ..., gt=0, description="Mean altitude (meters), as TAT-C defines it."
    )
    equator_crossing_time: time = Field(
        ..., description="Local mean solar time of the ascending node."
    )
    inclination: float | None = Field(
        default=None,
        gt=0,
        lt=180,
        description="Inclination (degrees) of an orbit whose local time drifts, or "
        + "None (the default) for the sun-synchronous inclination.",
    )
    phase: float | None = Field(
        default=None,
        ge=0,
        lt=360,
        description="Argument of latitude (degrees) at the epoch, relative to the "
        + "other satellites of the configuration, for a constellation whose "
        + "operator maintains their relative phases, or None (the default) for an "
        + "arbitrary phase.",
    )

    def build(
        self, epoch: datetime, phase: float | None = None
    ) -> SunSynchronousOrbit | CircularOrbit:
        """
        Builds the orbit at an epoch.

        Args:
            epoch (datetime): The epoch.
            phase (float | None): The argument of latitude (degrees) at the
                epoch, if not the cataloged one (or zero).

        Returns:
            SunSynchronousOrbit | CircularOrbit: The orbit.
        """
        true_anomaly = _get_value(phase, self.phase)
        orbit = SunSynchronousOrbit(
            mean_altitude=self.mean_altitude,
            equator_crossing_time=self.equator_crossing_time,
            epoch=epoch,
            true_anomaly=true_anomaly,
        )
        if self.inclination is None:
            return orbit
        return CircularOrbit(
            mean_altitude=self.mean_altitude,
            inclination=self.inclination,
            right_ascension_ascending_node=orbit.get_right_ascension_ascending_node(),
            epoch=epoch,
            true_anomaly=true_anomaly,
        )

    def describe(self) -> str:
        """
        Describes the orbit.

        Returns:
            str: The description.
        """
        text = (
            f"Sun-synchronous, {self.mean_altitude / 1e3:.2f} km, ascending node "
            f"at {self.equator_crossing_time:%H:%M:%S}"
        )
        if self.inclination is not None:
            text += f", drifting at {self.inclination} deg"
        if self.phase is not None:
            text += f", phase {self.phase:g} deg"
        return text


def _get_value(given: float | None, cataloged: float | None) -> float:
    """
    Gets an angle given to build an orbit, or else its cataloged value, or else
    zero.
    """
    if given is not None:
        return given
    return 0 if cataloged is None else cataloged


class InclinedPlane(_CatalogModel):
    """
    The plane of an inclined circular orbit, at the plane grade: its mean
    altitude and inclination. Built as a `CircularOrbit` at an epoch. Its right
    ascension of ascending node and phase are those of a constellation whose
    operator maintains the satellites' relative planes and phases, or else are
    spaced evenly among the satellites built (see
    `OrbitConfiguration.build_orbits`).
    """

    type: Literal["circular"] = Field(
        default="circular", description="Orbit discriminator."
    )
    mean_altitude: float = Field(
        ..., gt=0, description="Mean altitude (meters), as TAT-C defines it."
    )
    inclination: float = Field(..., gt=0, lt=180, description="Inclination (degrees).")
    right_ascension_ascending_node: float | None = Field(
        default=None,
        ge=0,
        lt=360,
        description="Right ascension of ascending node (degrees) at the epoch, "
        + "relative to the other satellites of the configuration, or None (the "
        + "default) for an arbitrary one.",
    )
    phase: float | None = Field(
        default=None,
        ge=0,
        lt=360,
        description="Argument of latitude (degrees) at the epoch, relative to the "
        + "other satellites of the configuration, or None (the default) for an "
        + "arbitrary phase.",
    )

    def build(
        self,
        epoch: datetime,
        phase: float | None = None,
        right_ascension_ascending_node: float | None = None,
    ) -> CircularOrbit:
        """
        Builds the orbit at an epoch.

        Args:
            epoch (datetime): The epoch.
            phase (float | None): The argument of latitude (degrees) at the
                epoch, if not the cataloged one (or zero).
            right_ascension_ascending_node (float | None): The right ascension
                of ascending node (degrees), if not the cataloged one (or zero).

        Returns:
            CircularOrbit: The orbit.
        """
        return CircularOrbit(
            mean_altitude=self.mean_altitude,
            inclination=self.inclination,
            right_ascension_ascending_node=_get_value(
                right_ascension_ascending_node, self.right_ascension_ascending_node
            ),
            epoch=epoch,
            true_anomaly=_get_value(phase, self.phase),
        )

    def describe(self) -> str:
        """
        Describes the orbit.

        Returns:
            str: The description.
        """
        text = f"{self.mean_altitude / 1e3:.2f} km, {self.inclination} deg"
        if self.right_ascension_ascending_node is not None:
            text += f", node {self.right_ascension_ascending_node:g} deg"
        if self.phase is not None:
            text += f", phase {self.phase:g} deg"
        return text


class GeostationarySlot(_CatalogModel):
    """
    The longitude of a geostationary orbit, at the ground track grade. Built as
    a `GeosynchronousOrbit` at an epoch, maintained at its longitude.
    """

    type: Literal["geosynchronous"] = Field(
        default="geosynchronous", description="Orbit discriminator."
    )
    longitude: float = Field(
        ..., ge=-180, le=180, description="Longitude (degrees) of the subpoint."
    )
    inclination: float = Field(
        default=0, ge=0, lt=180, description="Inclination (degrees)."
    )

    def build(self, epoch: datetime) -> GeosynchronousOrbit:
        """
        Builds the orbit at an epoch.

        Args:
            epoch (datetime): The epoch.

        Returns:
            GeosynchronousOrbit: The orbit.
        """
        return GeosynchronousOrbit(
            longitude=self.longitude, inclination=self.inclination, epoch=epoch
        )

    def describe(self) -> str:
        """
        Describes the orbit.

        Returns:
            str: The description.
        """
        return f"Geostationary at {self.longitude} deg, {self.inclination} deg"


CatalogOrbit = Annotated[
    SunSynchronousPlane | InclinedPlane | GeostationarySlot | RepeatGroundTrackOrbit,
    Field(discriminator="type"),
]
"""An orbit of a catalog system, by its `type`."""

# the orbit types of each class at each grade
_ORBIT_TYPES = {
    (SystemClass.SUN_SYNCHRONOUS, SystemGrade.PLANE): "sso",
    (SystemClass.SUN_SYNCHRONOUS, SystemGrade.GROUND_TRACK): "repeat",
    (SystemClass.INCLINED, SystemGrade.PLANE): "circular",
    (SystemClass.INCLINED, SystemGrade.GROUND_TRACK): "repeat",
    (SystemClass.GEOSTATIONARY, SystemGrade.GROUND_TRACK): "geosynchronous",
}


# the optional values of orbits at the plane grade, which a configuration gives
# for every orbit or for none
_RELATIVE_VALUES = ("phase", "right_ascension_ascending_node")


def describe_orbit(orbit: CatalogOrbit) -> str:
    """
    Describes an orbit of a configuration, such as "Sun-synchronous, 834.10 km,
    ascending node at 13:30:03".

    Args:
        orbit (CatalogOrbit): The orbit.

    Returns:
        str: The description.
    """
    if isinstance(orbit, RepeatGroundTrackOrbit):
        node = orbit.ascending_node
        return (
            f"Repeat of {orbit.revolutions} orbits in {orbit.days} days "
            f"({orbit.get_mean_altitude() / 1e3:.2f} km), ascending node at "
            f"{node.longitude} deg on {node.time:%Y-%m-%d %H:%M:%S} UTC"
        )
    return orbit.describe()


def _get_values(orbit: Any) -> set[str]:
    """
    Gets the values of an orbit that need provenance: its fields, without its
    `type` or unset optional values (such as a derived inclination).
    """
    return {field for field, value in orbit if field != "type" and value is not None}


class _SystemConfiguration(_CatalogModel):
    """
    Base class of system configurations.
    """

    grade: SystemGrade = Field(..., description="The use for which it is validated.")
    description: str = Field(..., description="What the configuration represents.")
    limitations: tuple[str, ...] = Field(
        default=(), description="What it does not represent."
    )

    def get_measured(self) -> list[Provenance]:
        """
        Gets the provenance of the values fitted to reference data.
        """
        raise NotImplementedError

    def get_content(self) -> dict[str, Any]:
        """
        Gets the values that affect analysis results.
        """
        raise NotImplementedError

    def get_holdout_gap(
        self, period: str, platform_periods: Mapping[str, str] | None = None
    ) -> int | None:
        """
        Gets the days between the data to which values were fitted and the data
        of a test: between the nearest days of the periods.

        Args:
            period (str): The period of the data of the test (see `Period`).
            platform_periods (Mapping[str, str] | None): The period of each
                platform's data, if they differ, with which the values fitted to
                a spacecraft's own data are compared.

        Returns:
            int | None: The days, or `None` if no value was fitted.
        """
        del platform_periods
        return _get_gap(self.get_measured(), period)

    def get_members(self) -> list[str]:
        """
        Gets the platforms on which the configuration is tested: its spacecraft,
        or, for a network, nothing (its record names the network).

        Returns:
            list[str]: The identifiers.
        """
        return []

    def serves(self, grade: SystemGrade | str) -> bool:
        """
        Checks whether this configuration serves a use (see
        `tatc.catalog.system_standard.serves_system`).

        Args:
            grade (SystemGrade | str): The grade of the use.

        Returns:
            bool: `True`, if it serves the use.
        """
        return serves_system(self.grade, grade)


def _get_gap(measured: list[Provenance], period: str) -> int | None:
    """
    Gets the least days between the periods of fitted values and a period.
    """
    if not measured:
        return None
    return min(_get_days_apart(p.period, period) for p in measured)


def _get_digest(content: dict[str, Any]) -> str:
    """
    Gets the digest of the values that affect analysis results.
    """
    encoded = json.dumps(content, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


class OrbitConfiguration(_SystemConfiguration):
    """
    The orbits of a satellite system's spacecraft, as validated for a use.
    """

    orbits: dict[str, CatalogOrbit] = Field(
        ..., min_length=1, description="The orbit of each spacecraft, by identifier."
    )
    provenance: dict[str, dict[str, Provenance]] = Field(
        ...,
        description="The origin of each value of each orbit, by spacecraft and "
        + "field (`ascending_node` for both its time and longitude).",
    )

    @model_validator(mode="after")
    def _check_provenance(self) -> OrbitConfiguration:
        if set(self.provenance) != set(self.orbits):
            raise ValueError("each orbit needs provenance, by spacecraft")
        for spacecraft, orbit in self.orbits.items():
            _check_provenance(_get_values(orbit), self.provenance[spacecraft])
        for field in _RELATIVE_VALUES:
            given = {
                getattr(orbit, field) is not None
                for orbit in self.orbits.values()
                if field in type(orbit).model_fields
            }
            if len(given) > 1:
                raise ValueError(f"{field} is given for every orbit or for none")
        return self

    def build_orbits(
        self, members: Sequence[str], epoch: datetime | None = None
    ) -> list[AllOrbits]:
        """
        Builds the orbits of spacecraft. An orbit at the plane grade is built at
        the epoch with its cataloged phase (and, for an inclined orbit, right
        ascension of ascending node), if any, or else with phases (and nodes)
        spaced evenly among the orbits built at the plane grade, so that
        satellites that share a plane are not built at the same position.

        Args:
            members (Sequence[str]): The identifiers of the spacecraft.
            epoch (datetime | None): The epoch of orbits at the plane grade or of
                geostationary orbits.

        Returns:
            list[AllOrbits]: The orbits, in the order of the members.

        Raises:
            ValueError: If a spacecraft has no orbit, or an orbit needs an epoch.
        """
        unknown = [m for m in members if m not in self.orbits]
        if unknown:
            raise ValueError(f"no orbits of {unknown} in this configuration")
        planes = [
            m
            for m in members
            if isinstance(self.orbits[m], (SunSynchronousPlane, InclinedPlane))
        ]
        spread = {m: 360 * i / len(planes) for i, m in enumerate(planes)}
        orbits = []
        for member in members:
            orbit = self.orbits[member]
            if isinstance(orbit, RepeatGroundTrackOrbit):
                orbits.append(orbit)
                continue
            if epoch is None:
                raise ValueError(f"{member}'s orbit is built at an epoch (or start)")
            if isinstance(orbit, GeostationarySlot):
                orbits.append(orbit.build(epoch))
            elif isinstance(orbit, InclinedPlane):
                orbits.append(
                    orbit.build(
                        epoch,
                        phase=_get_value(orbit.phase, spread[member]),
                        right_ascension_ascending_node=_get_value(
                            orbit.right_ascension_ascending_node, spread[member]
                        ),
                    )
                )
            else:
                orbits.append(
                    orbit.build(epoch, phase=_get_value(orbit.phase, spread[member]))
                )
        return orbits

    def get_members(self) -> list[str]:
        return list(self.orbits)

    def get_holdout_gap(
        self, period: str, platform_periods: Mapping[str, str] | None = None
    ) -> int | None:
        """
        Gets the days between the data to which values were fitted and the data
        of a test: the least, over spacecraft, of the days between the nearest
        days of the periods of the data to which its values were fitted and of
        its own test data (from `platform_periods`, if given, or else `period`).

        Args:
            period (str): The period of the data of the test (see `Period`).
            platform_periods (Mapping[str, str] | None): The period of each
                tested spacecraft's data, if they differ (see
                `ConfigurationRecord.platform_periods`).

        Returns:
            int | None: The days, or `None` if no tested spacecraft has fitted
            values.
        """
        periods = platform_periods or dict.fromkeys(self.orbits, period)
        gaps = [
            _get_days_apart(provenance.period, test_period)
            for spacecraft, test_period in periods.items()
            for provenance in self.provenance.get(spacecraft, {}).values()
            if provenance.kind == "measured"
        ]
        return min(gaps) if gaps else None

    def get_measured(self) -> list[Provenance]:
        return [
            provenance
            for values in self.provenance.values()
            for provenance in values.values()
            if provenance.kind == "measured"
        ]

    def get_content(self) -> dict[str, Any]:
        return {
            "orbits": {
                key: orbit.model_dump(mode="json", exclude_none=True)
                for key, orbit in self.orbits.items()
            }
        }


class RadarConfiguration(_SystemConfiguration):
    """
    The scanning of a radar network's stations, as validated for a use, such as
    a volume coverage pattern of weather radars.
    """

    _VALUES: ClassVar[tuple[str, ...]] = (
        "band",
        "max_range",
        "min_elevation_angle",
        "max_elevation_angle",
        "beam_width",
    )

    band: RadarBand = Field(..., description="The frequency band.")
    max_range: float = Field(..., gt=0, description="Maximum slant range (meters).")
    min_elevation_angle: float = Field(
        ...,
        ge=-90,
        le=90,
        description="Lowest scanned elevation angle (degrees), unless a station "
        + "of the table gives its own.",
    )
    max_elevation_angle: float = Field(
        ..., ge=0, le=90, description="Highest scanned elevation angle (degrees)."
    )
    beam_width: float = Field(
        ..., ge=0, le=180, description="Half-power beam width (degrees)."
    )
    volume_coverage_patterns: tuple[int, ...] = Field(
        default=(),
        description="The volume coverage patterns that the configuration "
        + "represents, by number (such as 12 and 212), with which its coverage is "
        + "compared with the volumes the radars scanned.",
    )
    provenance: dict[str, Provenance] = Field(
        ..., description="The origin of each value, by field."
    )

    @model_validator(mode="after")
    def _check_values(self) -> RadarConfiguration:
        _check_provenance(set(self._VALUES), self.provenance)
        if self.max_elevation_angle < self.min_elevation_angle:
            raise ValueError("max_elevation_angle is below min_elevation_angle")
        return self

    def get_measured(self) -> list[Provenance]:
        return [p for p in self.provenance.values() if p.kind == "measured"]

    def get_content(self) -> dict[str, Any]:
        return {"radar": self.model_dump(mode="json", include=set(self._VALUES))}


class SurfaceConfiguration(_SystemConfiguration):
    """
    The instruments of a surface network's stations, as validated for a use.
    """

    height: float = Field(
        ..., gt=0, description="Height (meters) of the instruments above the ground."
    )
    instruments: tuple[InstrumentSpecification, ...] = Field(
        ..., min_length=1, description="The instruments of each station."
    )
    provenance: dict[str, Provenance] = Field(
        ..., description="The origin of the height (`height`)."
    )

    @model_validator(mode="after")
    def _check_values(self) -> SurfaceConfiguration:
        _check_provenance({"height"}, self.provenance)
        return self

    def get_measured(self) -> list[Provenance]:
        provenances = list(self.provenance.values()) + [
            provenance
            for instrument in self.instruments
            for provenance in instrument.provenance.values()
        ]
        return [p for p in provenances if p.kind == "measured"]

    def get_content(self) -> dict[str, Any]:
        return {
            "height": self.height,
            "instruments": [
                {"type": i.type, "parameters": i.parameters} for i in self.instruments
            ],
        }


@functools.lru_cache(maxsize=None)
def _read_table(file: str) -> tuple[pd.DataFrame, str]:
    """
    Reads a table of stations and its digest (of its bytes, with line endings
    normalized, as git may check the file out with either).
    """
    content = (_TABLES / file).read_bytes().replace(b"\r\n", b"\n")
    table = pd.read_csv(_TABLES / file, dtype={"id": str, "name": str, "wban": str})
    return table, hashlib.sha256(content).hexdigest()


class StationTable(_CatalogModel):
    """
    The stations of a network, in a table of comma-separated values: an
    identifier (`id`), `name`, `latitude` and `longitude` (decimal degrees),
    the elevation of the `ground` (meters above mean sea level, as published),
    and, for radars, the `height` (meters) of the antenna above the ground and
    an optional lowest elevation angle of the station's own
    (`min_elevation_angle`, degrees), or, for surface stations, an optional
    WBAN identifier (`wban`). Heights above mean sea level are converted
    to the WGS 84 ellipsoid with the EGM2008 geoid when stations are built (see
    `tatc.utils.compute_geoid_height`).
    """

    file: str = Field(..., description="The file, among the packaged tables.")
    as_of: date = Field(..., description="The date of the snapshot.")
    provenance: Provenance = Field(
        ...,
        description="Where the table comes from: `measured` (with the date of the "
        + "snapshot as its period), so that tests of its stations' positions and "
        + "membership are held out from it.",
    )

    @model_validator(mode="after")
    def _check_file(self) -> StationTable:
        if not (_TABLES / self.file).is_file():
            raise ValueError(f"no packaged table {self.file!r}")
        return self

    def load(self) -> pd.DataFrame:
        """
        Loads the table.

        Returns:
            pandas.DataFrame: A copy of the table.
        """
        return _read_table(self.file)[0].copy()

    @property
    def digest(self) -> str:
        """
        A digest of the table's contents.
        """
        return _read_table(self.file)[1]

    def get_period(self) -> str:
        """
        Gets the period of the snapshot, for its holdout from a test.

        Returns:
            str: The period (see `Period`).
        """
        return self.as_of.isoformat()


class SystemEntry(_CatalogModel):
    """
    A system of a past or present program, in configurations validated against
    reference data, from which to build the platforms of baseline systems in
    pre-phase A studies: a `SatelliteSystem`, a `RadarNetwork`, or a
    `SurfaceNetwork` (parsed by its `system_class` with `parse_system`).
    """

    id: str = Field(..., description="Entry identifier.")
    name: str = Field(..., description="System name.")
    aliases: tuple[str, ...] = Field(default=(), description="Other names.")
    system_class: SystemClass = Field(..., description="System class.")
    selection: str = Field(..., description="Why the system is cataloged.")
    operator: str = Field(..., description="The organizations that operate it.")
    default_configuration: str = Field(
        ..., description="The configuration used if none is specified."
    )
    configurations: dict[str, _SystemConfiguration] = Field(
        ..., min_length=1, description="The ways of modeling it, by identifier."
    )
    limitations: tuple[str, ...] = Field(
        default=(), description="What no configuration represents."
    )
    record: ValidationRecord | None = Field(
        default=None, description="The results of its validation."
    )

    @model_validator(mode="after")
    def _check_entry(self) -> SystemEntry:
        if self.default_configuration not in self.configurations:
            raise ValueError("default_configuration is not a configuration")
        for configuration in self.configurations.values():
            # checks that the class has gates at this grade
            get_system_gates(self.system_class, configuration.grade)
        if self.record is not None:
            self._check_record()
        return self

    def _get_platforms(self) -> set[str]:
        """
        Gets the platforms that a record may name.
        """
        return {self.id}

    def _check_record(self) -> None:
        """
        Checks that the record is of this system's configurations and platforms.
        """
        if self.record.entry_id != self.id:
            raise ValueError("the record belongs to a different system")
        unknown = set(self.record.configurations) - set(self.configurations)
        if unknown:
            raise ValueError(f"records of unknown {sorted(unknown)}")
        platforms = self._get_platforms()
        for key, record in self.record.configurations.items():
            if not set(record.platforms) <= platforms:
                raise ValueError(
                    f"{key}: the record's platforms {list(record.platforms)} are "
                    f"not of the system, {sorted(platforms)}"
                )
            if record.transfers:
                raise ValueError(f"{key}: a system's record has no transfers")

    @property
    def validation(self) -> EntryValidation:
        """
        The validation status of this system, computed from its records and the
        standard rather than declared (see
        `tatc.catalog.system_standard.evaluate_system`).
        """
        return evaluate_system(self)

    def get_digest(self, configuration: str) -> str:
        """
        Gets the digest of the values of a configuration that affect analysis
        results, with which a validation record identifies the configuration it
        tested.

        Args:
            configuration (str): The identifier of the configuration.

        Returns:
            str: The digest.
        """
        return _get_digest(self.configurations[configuration].get_content())

    def get_holdout_gap(
        self,
        configuration: str,
        period: str,
        platform_periods: Mapping[str, str] | None = None,
    ) -> int | None:
        """
        Gets the days between the data to which a configuration's values were
        fitted and the data of a test (see `get_holdout_gap` of its class).

        Args:
            configuration (str): The identifier of the configuration.
            period (str): The period of the data of the test (see `Period`).
            platform_periods (Mapping[str, str] | None): The period of each
                platform's data, if they differ.

        Returns:
            int | None: The days, or `None` if no value was fitted.
        """
        return self.configurations[configuration].get_holdout_gap(
            period, platform_periods
        )

    def get_members(self) -> list[str]:
        """
        Gets the identifiers of the system's members: its spacecraft, or its
        stations.

        Returns:
            list[str]: The identifiers.
        """
        raise NotImplementedError

    def summarize(self) -> dict[str, Any]:
        """
        Summarizes the system, as a row of a table of systems (such as
        `pandas.DataFrame([s.summarize() for s in list_systems()])`).

        Returns:
            dict[str, Any]: Its `system` (identifier), `name`, `class`,
            validation `tier`, `configurations` (comma-separated), and number
            of `members` (spacecraft or stations).
        """
        return {
            "system": self.id,
            "name": self.name,
            "class": self.system_class.value,
            "tier": self.validation.tier.value,
            "configurations": ", ".join(self.configurations),
            "members": len(self.get_members()),
        }

    def select_configuration(
        self,
        configuration: str | None = None,
        *,
        grade: SystemGrade | str | None = None,
    ) -> str:
        """
        Selects a configuration: the one named, or else the default if it
        serves the grade (if any), or else the first that serves the grade.

        Args:
            configuration (str | None): The identifier of a configuration.
            grade (SystemGrade | str | None): The grade of the use.

        Returns:
            str: The identifier of the configuration.

        Raises:
            ValueError: If the configuration is unknown or does not serve the
                grade, or no configuration serves the grade.
        """
        if configuration is not None:
            if configuration not in self.configurations:
                raise ValueError(
                    f"{self.id} has no configuration {configuration!r}: "
                    f"{sorted(self.configurations)}"
                )
            if grade is not None and not self.configurations[configuration].serves(
                grade
            ):
                raise ValueError(
                    f"{self.id}:{configuration} does not serve the "
                    f"{SystemGrade(grade).value} grade"
                )
            return configuration
        if grade is None or self.configurations[self.default_configuration].serves(
            grade
        ):
            return self.default_configuration
        for key, candidate in self.configurations.items():
            if candidate.serves(grade):
                return key
        raise ValueError(
            f"{self.id} has no configuration that serves the "
            f"{SystemGrade(grade).value} grade"
        )


class SatelliteSystem(SystemEntry):
    """
    A satellite program: its spacecraft, and the orbits of its configurations.
    """

    system_class: Literal[
        SystemClass.SUN_SYNCHRONOUS, SystemClass.INCLINED, SystemClass.GEOSTATIONARY
    ] = Field(..., description="System class.")
    spacecraft: tuple[Spacecraft, ...] = Field(
        ..., min_length=1, description="The satellites of the system."
    )
    configurations: dict[str, OrbitConfiguration] = Field(
        ..., min_length=1, description="The ways of modeling it, by identifier."
    )

    @model_validator(mode="after")
    def _check_orbits(self) -> SatelliteSystem:
        spacecraft = {s.id for s in self.spacecraft}
        if len(spacecraft) != len(self.spacecraft):
            raise ValueError("spacecraft identifiers must be distinct")
        for key, configuration in self.configurations.items():
            unknown = set(configuration.orbits) - spacecraft
            if unknown:
                raise ValueError(f"{key}: unknown spacecraft {sorted(unknown)}")
            wanted = _ORBIT_TYPES[self.system_class, configuration.grade]
            for name, orbit in configuration.orbits.items():
                if orbit.type != wanted:
                    raise ValueError(
                        f"{key}: {name}'s orbit at the "
                        f"{configuration.grade.value} grade of a "
                        f"{self.system_class.value} system is of type {wanted}"
                    )
                if isinstance(orbit, RepeatGroundTrackOrbit) and (
                    orbit.is_sun_synchronous()
                    != (self.system_class == SystemClass.SUN_SYNCHRONOUS)
                ):
                    raise ValueError(
                        f"{key}: {name}'s repeat orbit has an inclination only "
                        "if it is not sun-synchronous"
                    )
        return self

    def _get_platforms(self) -> set[str]:
        return {s.id for s in self.spacecraft}

    def get_members(self) -> list[str]:
        return [s.id for s in self.spacecraft]

    def get_spacecraft(self, spacecraft_id: str) -> Spacecraft:
        """
        Gets a spacecraft of this system.

        Args:
            spacecraft_id (str): The identifier of the spacecraft.

        Returns:
            Spacecraft: The spacecraft.

        Raises:
            KeyError: If the system has no such spacecraft.
        """
        for spacecraft in self.spacecraft:
            if spacecraft.id == spacecraft_id:
                return spacecraft
        raise KeyError(
            f"{self.id} has no spacecraft {spacecraft_id!r}: "
            f"{[s.id for s in self.spacecraft]}"
        )

    def get_payload_references(self, spacecraft_id: str) -> list[str]:
        """
        Gets the instrument catalog's entries that list a spacecraft of this
        system, except those of limb sounders and radio occultation receivers,
        whose soundings are analyzed with their own methods (see
        `tatc.catalog.InstrumentPayload.sounding`).

        Args:
            spacecraft_id (str): The identifier of the spacecraft.

        Returns:
            list[str]: The identifiers of the entries.
        """
        self.get_spacecraft(spacecraft_id)
        sounders = (InstrumentClass.LIMB, InstrumentClass.OCCULTATION)
        return [
            entry.id
            for entry in _loader().list_entries(platform=spacecraft_id)
            if entry.instrument_class not in sounders
        ]

    def build(  # pylint: disable=too-many-arguments,too-many-locals
        self,
        configuration: str | None = None,
        *,
        grade: SystemGrade | str | None = None,
        members: Sequence[str] | None = None,
        epoch: datetime | None = None,
        start: datetime | None = None,
        end: datetime | None = None,
        payload: bool | Sequence[str] = True,
        time_step: float | None = None,
    ) -> list[Satellite]:
        """
        Builds the satellites of a configuration (see `select_configuration`),
        named by their spacecraft's `name`, with the orbits of
        `OrbitConfiguration.build_orbits`: an orbit at the plane grade is built
        at the `epoch` (by default, the `start`), with the relative phases of a
        constellation whose operator maintains them, or else with phases (and,
        for inclined orbits, nodes) spaced evenly among the satellites built; an
        orbit maintained on a repeat ground track is pinned to its ascending
        node, at any epoch. Satellites carry the instruments of the instrument
        catalog's entries that list them (with `payload`, see
        `get_payload_references`), built with `tatc.catalog.build_payload` for
        the spacecraft.

        Args:
            configuration (str | None): The configuration, if not selected by
                grade.
            grade (SystemGrade | str | None): The grade of the use, if any.
            members (Sequence[str] | None): The identifiers of the spacecraft to
                build, in the order of the satellites. By default, the
                operational and extended spacecraft of the configuration.
            epoch (datetime | None): The epoch of orbits at the plane grade or of
                geostationary orbits.
            start (datetime | None): The start of the analysis, for instruments
                that schedule their tilt by day, and the default epoch.
            end (datetime | None): The end of the analysis, for instruments that
                schedule their tilt by day.
            payload (bool | Sequence[str]): Whether satellites carry their
                cataloged instruments, or the references of the instrument
                catalog to carry (such as `["viirs", "atms"]`).
            time_step (float | None): The analysis time step (seconds), for
                instruments that tailor a field of view to it.

        Returns:
            list[Satellite]: The satellites.

        Raises:
            ValueError: If a member has no orbit in the configuration, an orbit
                needs an epoch, or every spacecraft of the configuration is
                retired (and no members are given).
            KeyError: If a member is not a spacecraft of the system.
        """
        key = self.select_configuration(configuration, grade=grade)
        selected = self.configurations[key]
        if members is None:
            members = [
                s.id
                for s in self.spacecraft
                if s.id in selected.orbits and s.status != "retired"
            ]
            if not members:
                raise ValueError(
                    f"every spacecraft of {self.id}:{key} is retired; give members"
                )
        spacecraft = [self.get_spacecraft(member) for member in members]
        orbits = selected.build_orbits(members, epoch or start)
        satellites = []
        for member, orbit in zip(spacecraft, orbits):
            references = (
                self.get_payload_references(member.id)
                if payload is True
                else list(payload or [])
            )
            if not references:
                satellites.append(Satellite(name=member.name, orbit=orbit))
                continue
            satellites.append(
                _loader()
                .build_payload(
                    references,
                    orbit=orbit,
                    time_step=time_step,
                    platform=member.id,
                    start=start,
                    end=end,
                )
                .to_satellite(orbit, name=member.name)
            )
        return satellites


class _NetworkSystem(SystemEntry):
    """
    Base class of networks of stations, whose positions are in a table.
    """

    # the columns of the table that every station gives, and those it may give
    _REQUIRED_COLUMNS: ClassVar[tuple[str, ...]] = (
        "id",
        "name",
        "latitude",
        "longitude",
        "ground",
    )
    _OPTIONAL_COLUMNS: ClassVar[tuple[str, ...]] = ()

    stations: StationTable = Field(..., description="The stations of the network.")

    @model_validator(mode="after")
    def _check_table(self) -> _NetworkSystem:
        table = self.stations.load()
        required = self._REQUIRED_COLUMNS
        missing = [c for c in required if c not in table.columns]
        unknown = [
            c for c in table.columns if c not in required + self._OPTIONAL_COLUMNS
        ]
        if missing or unknown:
            raise ValueError(
                f"the table of stations has columns {list(table.columns)}: "
                f"missing {missing}, unknown {unknown}"
            )
        if table["id"].duplicated().any() or table["id"].isna().any():
            raise ValueError("stations need distinct identifiers")
        if table[list(required)].isna().any().any():
            raise ValueError(f"stations need values of {list(required)}")
        if (
            not (table.latitude.abs() <= 90).all()
            or not (table.longitude.abs() <= 180).all()
        ):
            raise ValueError("stations need latitudes and longitudes in degrees")
        return self

    def get_members(self) -> list[str]:
        return list(self.stations.load()["id"])

    def get_station_names(self) -> dict[str, str]:
        """
        Gets the names of the stations (such as "Stillwater 2 W") by their
        identifiers, by which they are built.

        Returns:
            dict[str, str]: The names, by identifier, in the order of the table.
        """
        table = self.stations.load()
        return dict(zip(table["id"], table["name"]))

    def get_digest(self, configuration: str) -> str:
        return _get_digest(
            {
                **self.configurations[configuration].get_content(),
                "table": self.stations.digest,
            }
        )

    def get_holdout_gap(
        self,
        configuration: str,
        period: str,
        platform_periods: Mapping[str, str] | None = None,
    ) -> int | None:
        del platform_periods
        measured = self.configurations[configuration].get_measured()
        if self.stations.provenance.kind == "measured":
            measured.append(self.stations.provenance)
        return _get_gap(measured, period)

    def _select_stations(
        self, members: Sequence[str] | None, region: BaseGeometry | None
    ) -> tuple[pd.DataFrame, np.ndarray]:
        """
        Selects the rows of stations to build, with the elevations (meters) of
        their ground above the WGS 84 ellipsoid.
        """
        table = self.stations.load()
        if members is not None:
            unknown = set(members) - set(table["id"])
            if unknown:
                raise ValueError(f"{self.id} has no stations {sorted(unknown)}")
            # in the order of the members
            table = table.set_index("id", drop=False).loc[list(members)]
        if region is not None:
            table = table[
                shapely.intersects_xy(
                    region, table.longitude.to_numpy(), table.latitude.to_numpy()
                )
            ]
        ground = table.ground.to_numpy() + compute_geoid_height(
            table.longitude.to_numpy(), table.latitude.to_numpy()
        )
        return table, np.atleast_1d(ground)


class RadarNetwork(_NetworkSystem):
    """
    A network of weather radars: its stations, and their scanning in each
    configuration.
    """

    _REQUIRED_COLUMNS: ClassVar[tuple[str, ...]] = (
        *_NetworkSystem._REQUIRED_COLUMNS,
        "height",
    )
    _OPTIONAL_COLUMNS: ClassVar[tuple[str, ...]] = ("min_elevation_angle",)

    system_class: Literal[SystemClass.RADAR_NETWORK] = Field(
        ..., description="System class."
    )
    configurations: dict[str, RadarConfiguration] = Field(
        ..., min_length=1, description="The ways of modeling it, by identifier."
    )

    def build(
        self,
        configuration: str | None = None,
        *,
        grade: SystemGrade | str | None = None,
        members: Sequence[str] | None = None,
        region: BaseGeometry | None = None,
    ) -> list[RadarStation]:
        """
        Builds the radar stations of a configuration (see
        `select_configuration`), named by their identifiers, at the elevations
        of their ground above the WGS 84 ellipsoid (the table's height above
        mean sea level plus the EGM2008 geoid height). Terrain masks are not
        built (see `tatc.preprocess.compute_terrain_mask_for_station`).

        Args:
            configuration (str | None): The configuration, if not selected by
                grade.
            grade (SystemGrade | str | None): The grade of the use, if any.
            members (Sequence[str] | None): The identifiers of the stations to
                build, in the order of the stations. By default, every station
                (in the order of the table).
            region (shapely.geometry.base.BaseGeometry | None): A region
                (longitude and latitude) within which to build stations.

        Returns:
            list[RadarStation]: The stations.

        Raises:
            ValueError: If a member is unknown.
        """
        selected = self.configurations[
            self.select_configuration(configuration, grade=grade)
        ]
        table, ground = self._select_stations(members, region)
        return [
            _build_radar(selected, row, elevation)
            for row, elevation in zip(table.itertuples(), ground)
        ]


class SurfaceNetwork(_NetworkSystem):
    """
    A network of surface stations: its stations, and their instruments in each
    configuration.
    """

    _OPTIONAL_COLUMNS: ClassVar[tuple[str, ...]] = ("wban",)

    system_class: Literal[SystemClass.SURFACE_NETWORK] = Field(
        ..., description="System class."
    )
    configurations: dict[str, SurfaceConfiguration] = Field(
        ..., min_length=1, description="The ways of modeling it, by identifier."
    )

    def build(
        self,
        configuration: str | None = None,
        *,
        grade: SystemGrade | str | None = None,
        members: Sequence[str] | None = None,
        region: BaseGeometry | None = None,
        radius: float | None = None,
    ) -> list[SurfacePlatform]:
        """
        Builds the surface stations of a configuration (see
        `select_configuration`), named by their identifiers (such as ICAO
        codes; see `get_station_names`), at the elevations of their ground above
        the WGS 84 ellipsoid (the table's height above mean sea level plus the
        EGM2008 geoid height).

        Args:
            configuration (str | None): The configuration, if not selected by
                grade.
            grade (SystemGrade | str | None): The grade of the use, if any.
            members (Sequence[str] | None): The identifiers of the stations to
                build, in the order of the stations. By default, every station
                (in the order of the table).
            region (shapely.geometry.base.BaseGeometry | None): A region
                (longitude and latitude) within which to build stations.
            radius (float | None): The radius (meters) of the area that the
                in-situ instruments represent, for the variable of a study,
                instead of the cataloged assumption.

        Returns:
            list[SurfacePlatform]: The stations.

        Raises:
            ValueError: If a member is unknown.
        """
        selected = self.configurations[
            self.select_configuration(configuration, grade=grade)
        ]
        table, ground = self._select_stations(members, region)
        instruments = [spec.construct_instrument(None) for spec in selected.instruments]
        if radius is not None:
            instruments = [
                (
                    InSituInstrument(**{**instrument.model_dump(), "radius": radius})
                    if isinstance(instrument, InSituInstrument)
                    else instrument
                )
                for instrument in instruments
            ]
        return [
            SurfacePlatform(
                name=row.id,
                latitude=row.latitude,
                longitude=row.longitude,
                elevation=float(elevation),
                height=selected.height,
                instruments=[instrument.model_copy() for instrument in instruments],
            )
            for row, elevation in zip(table.itertuples(), ground)
        ]


AnySystem = Annotated[
    SatelliteSystem | RadarNetwork | SurfaceNetwork,
    Field(discriminator="system_class"),
]
"""A catalog system, by its `system_class`."""

_SYSTEM_ADAPTER: TypeAdapter[SatelliteSystem | RadarNetwork | SurfaceNetwork] = (
    TypeAdapter(AnySystem)
)


def parse_system(data: Mapping[str, Any]) -> SystemEntry:
    """
    Parses a system as the schema of its class: a `SatelliteSystem`, a
    `RadarNetwork`, or a `SurfaceNetwork`.

    Args:
        data (Mapping[str, Any]): The system's fields.

    Returns:
        SystemEntry: The system.

    Raises:
        pydantic.ValidationError: If the data are not of a valid system.
    """
    return _SYSTEM_ADAPTER.validate_python(data)


def _loader() -> Any:
    """
    Gets the catalog's loader, whose entries satellites carry; deferred, as the
    loader imports this module.
    """
    # pylint: disable-next=import-outside-toplevel,cyclic-import
    from . import loader

    return loader


def _build_radar(
    configuration: RadarConfiguration, row: Any, elevation: float
) -> RadarStation:
    """
    Builds a radar station of a table's row.
    """
    own = getattr(row, "min_elevation_angle", np.nan)
    return RadarStation(
        name=row.id,
        latitude=row.latitude,
        longitude=row.longitude,
        elevation=float(elevation),
        height=float(row.height),
        band=configuration.band,
        max_range=configuration.max_range,
        min_elevation_angle=(
            configuration.min_elevation_angle if pd.isna(own) else float(own)
        ),
        max_elevation_angle=configuration.max_elevation_angle,
        beam_width=configuration.beam_width,
    )

"""
Object schemas for catalog systems: satellite programs, radar networks, and
surface networks, in configurations validated against reference data (such as
the element sets of satellites, or the metadata of stations), and the platforms
built from them. The orbits of satellite systems are in
`tatc.catalog.system_orbits`. Their records are those of instrument entries
(`tatc.catalog.records`), and their standard is `tatc.catalog.system_standard`.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import functools
import hashlib
import json
import pathlib
import warnings
from collections.abc import Mapping, Sequence
from datetime import date, datetime, timedelta
from typing import Annotated, Any, ClassVar, Literal

import numpy as np
import pandas as pd
import shapely
from pydantic import Field, TypeAdapter, model_validator
from shapely.geometry.base import BaseGeometry

from ..schemas import (
    AttitudeSchedule,
    GroundStation,
    InSituInstrument,
    RadarStation,
    RepeatGroundTrackOrbit,
    Satellite,
    SurfacePlatform,
)
from ..schemas.architecture import Relay
from ..schemas.orbit import AllOrbits
from ..schemas.space.base import AllAttitudes
from ..schemas.surface.radar import RadarBand
from ..utils import compute_geoid_height
from .records import ValidationRecord, _CatalogModel, _get_days_apart
from .schemas import InstrumentSpecification, Provenance, _check_provenance
from .standard import EntryValidation, Grade, InstrumentClass, Tier
from .system_orbits import (
    RELATIVE_VALUES,
    CatalogOrbit,
    InclinedPlane,
    SunSynchronousPlane,
    get_orbit_values,
)
from .system_standard import (
    SystemClass,
    SystemGrade,
    evaluate_system,
    get_orbit_type,
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
    status: Literal["operational", "extended", "storage", "retired"] = Field(
        ...,
        description="The status on the date `as_of` (`storage` for a spare kept "
        + "in orbit, out of service), which selects the satellites built by "
        + "default (those neither retired nor in storage) but is not a value of "
        + "a configuration (its tests include every spacecraft).",
    )
    as_of: date = Field(..., description="The date of the status.")


def _as_list(values: str | Sequence[str]) -> list[str]:
    """Gets identifiers given as one string or a sequence of strings."""
    return [values] if isinstance(values, str) else list(values)


def _get_digest(content: dict[str, Any]) -> str:
    """
    Gets the digest of the values that affect analysis results.
    """
    encoded = json.dumps(content, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


class _SystemConfiguration(_CatalogModel):
    """
    Base class of system configurations.
    """

    grade: SystemGrade = Field(..., description="The use for which it is validated.")
    description: str = Field(..., description="What the configuration represents.")
    limitations: tuple[str, ...] = Field(
        default=(), description="What it does not represent."
    )

    def get_fitted(self, platform: str | None = None) -> list[Provenance]:
        """
        Gets the provenance of the values fitted to reference data: of a
        platform's values, for a configuration whose values are each
        platform's, or else of all of them.

        Args:
            platform (str | None): The identifier of a tested platform, if any.

        Returns:
            list[Provenance]: The provenance of the fitted values.
        """
        raise NotImplementedError

    def get_content(self) -> dict[str, Any]:
        """
        Gets the values that affect analysis results.

        Returns:
            dict[str, Any]: The values.
        """
        raise NotImplementedError

    def get_tested_platforms(self) -> list[str]:
        """
        Gets the platforms on which the configuration is tested: its spacecraft,
        or, for a network, nothing (its record names the network).

        Returns:
            list[str]: The identifiers.
        """
        return []

    def get_relative_values(self) -> tuple[str, ...]:
        """
        Gets the values of the orbits relative to each other that the
        configuration gives (`phase`, `right_ascension_ascending_node`), whose
        gates the standard applies too: none, but for orbits.

        Returns:
            tuple[str, ...]: The names of the values.
        """
        return ()

    def get_holdout_gap(
        self,
        period: str,
        platform_periods: Mapping[str, str] | None = None,
        *,
        platforms: Sequence[str] | None = None,
        others: Sequence[Provenance] = (),
    ) -> int | None:
        """
        Gets the days between the data to which values were fitted and the data
        of a test: the least days between the nearest days of the periods of the
        data to which values were fitted and of the data with which they were
        tested. A spacecraft's values are compared with its own test data (from
        `platform_periods`, if given, or else `period`).

        Args:
            period (str): The period of the data of the test (see `Period`).
            platform_periods (Mapping[str, str] | None): The period of each
                tested platform's data, if they differ (see
                `ConfigurationRecord.platform_periods`).
            platforms (Sequence[str] | None): The platforms tested, by default
                those of `platform_periods`, or else every platform of the
                configuration.
            others (Sequence[Provenance]): The provenance of other values
                tested with the configuration (such as a network's table).

        Returns:
            int | None: The days, or `None` if no value was fitted.
        """
        if platform_periods:
            periods = dict(platform_periods)
        else:
            tested = (
                self.get_tested_platforms() if platforms is None else list(platforms)
            )
            periods = dict.fromkeys(tested or [None], period)
        fitted = [
            (provenance, test_period)
            for platform, test_period in periods.items()
            for provenance in self.get_fitted(platform)
        ]
        fitted += [(p, period) for p in others if p.kind == "measured"]
        if not fitted:
            return None
        return min(_get_days_apart(p.period, test) for p, test in fitted)

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
            _check_provenance(get_orbit_values(orbit), self.provenance[spacecraft])
        for field in RELATIVE_VALUES:
            given = {
                getattr(orbit, field) is not None
                for orbit in self.orbits.values()
                if field in type(orbit).model_fields
            }
            if len(given) > 1:
                raise ValueError(f"{field} is given for every orbit or for none")
        return self

    def get_relative_values(self) -> tuple[str, ...]:
        """
        Gets the values of the orbits relative to each other that the
        configuration gives (`phase`, `right_ascension_ascending_node`), which
        the standard tests.

        Returns:
            tuple[str, ...]: The names of the values.
        """
        return tuple(
            field
            for field in RELATIVE_VALUES
            if any(getattr(o, field, None) is not None for o in self.orbits.values())
        )

    def build_orbits(
        self, members: Sequence[str], epoch: datetime | None = None
    ) -> list[AllOrbits]:
        """
        Builds the orbits of spacecraft. An orbit at the plane grade is built at
        the epoch with its cataloged phase (and, for an inclined orbit, right
        ascension of ascending node), if any, or else with phases (and nodes)
        spaced evenly among the orbits built at the plane grade, so that
        satellites that share a plane are not built at the same position. A
        geostationary orbit, maintained at its longitude, is built at the epoch
        (or at `SLOT_EPOCH`), and an orbit maintained on a repeat ground track is
        pinned to its ascending node.

        Args:
            members (Sequence[str]): The identifiers of the spacecraft.
            epoch (datetime | None): The epoch of orbits at the plane grade.

        Returns:
            list[AllOrbits]: The orbits, in the order of the members.

        Raises:
            ValueError: If a spacecraft has no orbit, or an orbit at the plane
                grade has no epoch.
        """
        unknown = [m for m in members if m not in self.orbits]
        if unknown:
            raise ValueError(f"no orbits of {unknown} in this configuration")
        planes = [
            m
            for m in members
            if isinstance(self.orbits[m], (SunSynchronousPlane, InclinedPlane))
        ]
        if planes and epoch is None:
            raise ValueError(f"{planes[0]}'s orbit is built at an epoch (or start)")
        spread = {m: 360 * i / len(planes) for i, m in enumerate(planes)}
        orbits = []
        for member in members:
            orbit = self.orbits[member]
            phase = getattr(orbit, "phase", None)
            node = getattr(orbit, "right_ascension_ascending_node", None)
            orbits.append(
                orbit.build(
                    epoch,
                    phase=spread.get(member, 0) if phase is None else phase,
                    right_ascension_ascending_node=(
                        spread.get(member, 0) if node is None else node
                    ),
                )
            )
        return orbits

    def get_tested_platforms(self) -> list[str]:
        return list(self.orbits)

    def get_fitted(self, platform: str | None = None) -> list[Provenance]:
        values = (
            self.provenance.values()
            if platform is None
            else [self.provenance.get(platform, {})]
        )
        return [p for value in values for p in value.values() if p.kind == "measured"]

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
        + "compared with the volumes the radars scanned (part of its digest, as "
        + "they select the volumes of its tests).",
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

    def get_fitted(self, platform: str | None = None) -> list[Provenance]:
        del platform
        return [p for p in self.provenance.values() if p.kind == "measured"]

    def get_content(self) -> dict[str, Any]:
        return {
            "radar": self.model_dump(mode="json", include=set(self._VALUES)),
            "volume_coverage_patterns": list(self.volume_coverage_patterns),
        }


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

    def get_fitted(self, platform: str | None = None) -> list[Provenance]:
        del platform
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
        description="Where the table comes from: `measured`, with the date of the "
        + "snapshot (`as_of`) as its period, so that tests of its stations' "
        + "positions and membership are held out from it.",
    )

    @model_validator(mode="after")
    def _check_table(self) -> StationTable:
        if not (_TABLES / self.file).is_file():
            raise ValueError(f"no packaged table {self.file!r}")
        if (
            self.provenance.kind == "measured"
            and self.provenance.period != self.as_of.isoformat()
        ):
            raise ValueError(
                f"the table's provenance has the period {self.provenance.period!r}, "
                f"not the date of its snapshot ({self.as_of})"
            )
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


class SystemEntry(_CatalogModel):
    """
    A system of a past or present program, in configurations validated against
    reference data, from which to build the platforms of baseline systems in
    pre-phase A studies: a `SatelliteSystem`, a `RadarNetwork`, or a
    `SurfaceNetwork`, parsed by its `system_class` with `parse_system`.
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

    @model_validator(mode="before")
    @classmethod
    def _check_class(cls, data: Any) -> Any:
        if cls is SystemEntry:
            raise TypeError(
                "a SystemEntry is a SatelliteSystem, RadarNetwork, or "
                "SurfaceNetwork: parse one with parse_system"
            )
        return data

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
        return _get_digest(
            {
                **self.configurations[configuration].get_content(),
                **self._get_shared_content(),
            }
        )

    def _get_shared_content(self) -> dict[str, Any]:
        """
        Gets the values of the entry that every configuration uses (such as a
        network's table of stations), which are part of their digests.
        """
        return {}

    def _get_shared_provenance(self) -> list[Provenance]:
        """
        Gets the provenance of the entry's values that every configuration uses,
        from which their tests are held out.
        """
        return []

    def get_holdout_gap(
        self,
        configuration: str,
        period: str,
        platform_periods: Mapping[str, str] | None = None,
        *,
        platforms: Sequence[str] | None = None,
    ) -> int | None:
        """
        Gets the days between the data to which a configuration's values were
        fitted and the data of a test (see
        `_SystemConfiguration.get_holdout_gap`).

        Args:
            configuration (str): The identifier of the configuration.
            period (str): The period of the data of the test (see `Period`).
            platform_periods (Mapping[str, str] | None): The period of each
                tested platform's data, if they differ.
            platforms (Sequence[str] | None): The platforms tested, if not
                every platform of the configuration.

        Returns:
            int | None: The days, or `None` if no value was fitted.
        """
        return self.configurations[configuration].get_holdout_gap(
            period,
            platform_periods,
            platforms=platforms,
            others=self._get_shared_provenance(),
        )

    def get_member_names(self) -> dict[str, str]:
        """
        Gets the names of the system's members (its spacecraft, or its
        stations) by their identifiers. Satellites are built under their names
        (such as "NOAA-20"), and stations under their identifiers (such as
        "KTLX"), as they are known.

        Returns:
            dict[str, str]: The names, by identifier.
        """
        raise NotImplementedError

    def get_member_id(self, member: str) -> str:
        """
        Gets the identifier of a member given by its identifier or name
        (ignoring case), such as the name of a built satellite.

        Args:
            member (str): The identifier or name.

        Returns:
            str: The identifier.

        Raises:
            KeyError: If the system has no such member.
        """
        names = self.get_member_names()
        if member in names:
            return member
        wanted = member.casefold()
        for identifier, name in names.items():
            if wanted in (identifier.casefold(), str(name).casefold()):
                return identifier
        raise KeyError(f"{self.id} has no member {member!r}")

    def _select(
        self, configuration: str | None, grade: SystemGrade | str | None
    ) -> str:
        """
        Selects the configuration to build (see `select_configuration`), with a
        warning if it is not validated.
        """
        key = self.select_configuration(configuration, grade=grade)
        validation = self.validation.configurations[key]
        if validation.tier != Tier.VALIDATED:
            warnings.warn(
                f"{self.id}:{key} is {validation.tier.value}, not validated: "
                f"{'; '.join(validation.gaps)}",
                stacklevel=3,
            )
        return key

    def get_members(self) -> list[str]:
        """
        Gets the identifiers of the system's members: its spacecraft, or its
        stations.

        Returns:
            list[str]: The identifiers.
        """
        return list(self.get_member_names())

    def get_default_members(self) -> list[str]:
        """
        Gets the identifiers of the members built by default: the spacecraft
        that are not retired, or every station.

        Returns:
            list[str]: The identifiers.
        """
        return self.get_members()

    def summarize(self) -> dict[str, Any]:
        """
        Summarizes the system, as a row of a table of systems (such as
        `pandas.DataFrame([s.summarize() for s in list_systems()])`).

        Returns:
            dict[str, Any]: Its `system` (identifier), `name`, `class`,
            validation `tier`, `configurations` (comma-separated), number of
            `members` (spacecraft or stations), number of those built by
            default (`active`: spacecraft neither retired nor in storage, or
            every station), and `role` (`observing`, `transmitting` for
            navigation satellites, or `relaying` for data relay satellites).
        """
        return {
            "system": self.id,
            "name": self.name,
            "class": self.system_class.value,
            "tier": self.validation.tier.value,
            "configurations": ", ".join(self.configurations),
            "members": len(self.get_members()),
            "active": len(self.get_default_members()),
            "role": getattr(self, "role", "observing"),
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
    role: Literal["observing", "transmitting", "relaying"] = Field(
        default="observing",
        description="What the satellites do: `observing`, carrying the instruments "
        + "of the instrument catalog that list them; `transmitting`, as the "
        + "navigation satellites whose signals radio occultation receivers sound "
        + "(see `tatc.analysis.collect_ro_observations`), which carry none; or "
        + "`relaying`, as data relay satellites, which carry none and relay the "
        + "data of other satellites to the ground (see `build_relay`).",
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
            wanted = get_orbit_type(self.system_class, configuration.grade)
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

    def get_member_names(self) -> dict[str, str]:
        return {s.id: s.name for s in self.spacecraft}

    def get_default_members(self) -> list[str]:
        return [s.id for s in self.spacecraft if s.status not in ("retired", "storage")]

    def get_spacecraft(self, spacecraft: str) -> Spacecraft:
        """
        Gets a spacecraft of this system by its identifier or name (ignoring
        case).

        Args:
            spacecraft (str): The identifier or name of the spacecraft.

        Returns:
            Spacecraft: The spacecraft.

        Raises:
            KeyError: If the system has no such spacecraft.
        """
        wanted = spacecraft.casefold()
        for candidate in self.spacecraft:
            if wanted in (candidate.id.casefold(), candidate.name.casefold()):
                return candidate
        raise KeyError(
            f"{self.id} has no spacecraft {spacecraft!r}: "
            f"{[s.id for s in self.spacecraft]}"
        )

    def get_payload_references(self, spacecraft: str) -> list[str]:
        """
        Gets the instrument catalog's entries that list a spacecraft of this
        system, except those of limb sounders and radio occultation receivers,
        whose soundings are analyzed with their own methods (see
        `tatc.catalog.InstrumentPayload.sounding`).

        Args:
            spacecraft (str): The identifier or name of the spacecraft.

        Returns:
            list[str]: The identifiers of the entries.
        """
        spacecraft_id = self.get_spacecraft(spacecraft).id
        sounders = (InstrumentClass.LIMB, InstrumentClass.OCCULTATION)
        return [
            entry.id
            for entry in _loader().list_entries(platform=spacecraft_id)
            if entry.instrument_class not in sounders
        ]

    def _get_references(
        self, payload: bool | str | Sequence[str], members: list[str]
    ) -> dict[str, list[str]]:
        """
        Gets the references of the instrument catalog that each member carries:
        its cataloged instruments, or those given that list it.
        """
        if payload is False:
            return {member: [] for member in members}
        if payload is True:
            return {member: self.get_payload_references(member) for member in members}
        given = {
            reference: _loader().get_entry(reference) for reference in _as_list(payload)
        }
        carried = {
            member: [ref for ref, entry in given.items() if member in entry.platforms]
            for member in members
        }
        uncarried = [
            ref for ref in given if not any(ref in refs for refs in carried.values())
        ]
        if uncarried:
            raise ValueError(
                f"no satellite built of {self.id} carries {uncarried} (see "
                "get_payload_references)"
            )
        return carried

    def build(  # pylint: disable=too-many-arguments,too-many-locals
        self,
        configuration: str | None = None,
        *,
        grade: SystemGrade | str | None = None,
        members: str | Sequence[str] | None = None,
        epoch: datetime | None = None,
        start: datetime | None = None,
        end: datetime | None = None,
        payload: bool | str | Sequence[str] | None = None,
        time_step: float | None = None,
        instrument_grade: Grade | str | None = None,
        attitude: AllAttitudes | AttitudeSchedule | None = None,
    ) -> list[Satellite]:
        """
        Builds the satellites of a configuration (see `select_configuration`),
        named by their spacecraft's `name`, with the orbits of
        `OrbitConfiguration.build_orbits`: an orbit at the plane grade is built
        at the `epoch` (by default, the `start`), with the relative phases of a
        constellation whose operator maintains them, or else with phases (and,
        for inclined orbits, nodes) spaced evenly among the satellites built; a
        geostationary orbit or one maintained on a repeat ground track is built
        for any epoch. Satellites carry the instruments of the instrument
        catalog's entries that list them (see `get_payload_references`), built
        with `tatc.catalog.build_payload` for the spacecraft, when a `time_step`
        is given (most imagers and sounders tailor a field of view to it). A
        configuration that is not validated warns. A satellite that
        carries none (such as a radio occultation satellite, whose soundings are
        analyzed with their own methods) has TAT-C's default instrument, which
        observes within a field of regard of 180 degrees, with a warning, as
        do the satellites of a system whose `role` is `transmitting` (such as
        those of a navigation system, for `tatc.analysis.collect_ro_observations`)
        or `relaying` (see `build_relay`), without one.

        Args:
            configuration (str | None): The configuration, if not selected by
                grade.
            grade (SystemGrade | str | None): The grade of the use, if any.
            members (str | Sequence[str] | None): The identifiers (or names) of
                the spacecraft to build, in the order of the satellites. By
                default, the spacecraft of the configuration that are neither
                retired nor in storage.
            epoch (datetime | None): The epoch of orbits at the plane grade (or
                of geostationary orbits).
            start (datetime | None): The start of the analysis, for instruments
                that schedule their tilt by day, and the default epoch.
            end (datetime | None): The end of the analysis, for instruments that
                schedule their tilt by day.
            payload (bool | str | Sequence[str] | None): Whether satellites
                carry their cataloged instruments, or the references of the
                instrument catalog to carry (such as `["viirs", "atms"]`), each
                on the satellites it lists. By default, their cataloged
                instruments if a `time_step` is given, or else none (with a
                warning).
            time_step (float | None): The analysis time step (seconds), for
                instruments that tailor a field of view to it (most imagers and
                sounders that scan across track).
            instrument_grade (Grade | str | None): The grade of the instruments'
                configurations, if any (see `tatc.catalog.build_payload`).
            attitude (AllAttitudes | AttitudeSchedule | None): The attitude of
                the satellites, if not their instruments' (see
                `tatc.catalog.build_payload`).

        Returns:
            list[Satellite]: The satellites.

        Raises:
            ValueError: If a member has no orbit in the configuration, an orbit
                needs an epoch, every spacecraft of the configuration is retired
                or in storage (and no members are given), a given reference is
                carried by no
                satellite built, or an instrument needs a time step.
            KeyError: If a member is not a spacecraft of the system.
        """
        key = self._select(configuration, grade)
        selected = self.configurations[key]
        if self.role in ("transmitting", "relaying"):
            if payload not in (None, False):
                kind = "transmitters" if self.role == "transmitting" else "relays"
                raise ValueError(
                    f"{self.id}'s satellites are {kind}, which carry no "
                    "cataloged instruments"
                )
            payload = False
        if payload is None:
            payload = time_step is not None
            if not payload:
                warnings.warn(
                    "the satellites carry no cataloged instruments without a "
                    "time_step (give one, or payload=False)",
                    stacklevel=2,
                )
        if members is None:
            members = [m for m in self.get_default_members() if m in selected.orbits]
            if not members:
                raise ValueError(
                    f"every spacecraft of {self.id}:{key} is retired or in "
                    "storage; give members"
                )
        spacecraft = [self.get_spacecraft(member) for member in _as_list(members)]
        ids = [s.id for s in spacecraft]
        references = self._get_references(payload, ids)
        orbits = selected.build_orbits(ids, epoch or start)
        satellites = []
        for member, orbit in zip(spacecraft, orbits):
            if not references[member.id]:
                if payload is not False:
                    carried = (
                        "no cataloged instrument other than sounders"
                        if payload is True
                        else "none of the instruments given"
                    )
                    warnings.warn(
                        f"{member.name} carries {carried}, so it has TAT-C's "
                        "default instrument (a field of regard of 180 deg)",
                        stacklevel=2,
                    )
                satellites.append(Satellite(name=member.name, orbit=orbit))
                continue
            try:
                built = _loader().build_payload(
                    references[member.id],
                    orbit=orbit,
                    grade=instrument_grade,
                    time_step=time_step,
                    attitude=attitude,
                    platform=member.id,
                    start=start,
                    end=end,
                )
            except ValueError as error:
                if time_step is None and "time step" in str(error):
                    raise ValueError(
                        f"{member.name}'s instruments need build(time_step=...), the "
                        f"time step of the analysis: {error}"
                    ) from error
                raise
            satellites.append(built.to_satellite(orbit, name=member.name))
        return satellites

    def build_relay(
        self,
        configuration: str | None = None,
        *,
        grade: SystemGrade | str | None = None,
        members: str | Sequence[str] | None = None,
        epoch: datetime | None = None,
        name: str | None = None,
        ground_stations: Sequence[GroundStation] = (),
        min_grazing_altitude: float = 0,
        max_range: float | None = None,
        min_access_time: timedelta = timedelta(0),
    ) -> Relay:
        """
        Builds the relay of a system whose `role` is `relaying` (see
        `tatc.schemas.Relay`): its satellites (see `build`), which satellites of
        an architecture contact by line of sight (see
        `tatc.analysis.collect_relay_contacts`). Without `ground_stations`, each
        relay satellite is assumed to be always in contact with the ground (as
        those of a geostationary relay are with the terminals that serve them).

        Args:
            configuration (str | None): The configuration, if not selected by
                grade.
            grade (SystemGrade | str | None): The grade of the use, if any.
            members (str | Sequence[str] | None): The identifiers (or names) of
                the spacecraft to build, by default those neither retired nor
                in storage.
            epoch (datetime | None): The epoch of the orbits, if needed.
            name (str | None): The name of the relay, by default the system's
                identifier.
            ground_stations (Sequence[GroundStation]): The ground stations that
                receive the relay's data, if any.
            min_grazing_altitude (float): The minimum altitude (meters) of the
                line of sight between a satellite and a relay satellite.
            max_range (float | None): The maximum range (meters), if any.
            min_access_time (timedelta): The minimum duration of a contact.

        Returns:
            Relay: The relay.

        Raises:
            ValueError: If the system's role is not `relaying`.
        """
        # pylint: disable=too-many-arguments
        if self.role != "relaying":
            raise ValueError(
                f"{self.id}'s satellites are not relays (its role is {self.role})"
            )
        satellites = self.build(
            configuration, grade=grade, members=members, epoch=epoch, payload=False
        )
        return Relay(
            name=self.id if name is None else name,
            satellites=satellites,
            ground_stations=list(ground_stations),
            min_grazing_altitude=min_grazing_altitude,
            max_range=max_range,
            min_access_time=min_access_time,
        )


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
    def _check_stations(self) -> _NetworkSystem:
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

    def get_member_names(self) -> dict[str, str]:
        table = self.stations.load()
        return dict(zip(table["id"], table["name"]))

    def _get_shared_content(self) -> dict[str, Any]:
        return {"table": self.stations.digest}

    def _get_shared_provenance(self) -> list[Provenance]:
        return [self.stations.provenance]

    def _select_stations(
        self, members: str | Sequence[str] | None, region: BaseGeometry | None
    ) -> tuple[pd.DataFrame, np.ndarray]:
        """
        Selects the rows of stations to build, with the elevations (meters) of
        their ground above the WGS 84 ellipsoid.
        """
        table = self.stations.load()
        if members is not None:
            identifiers, unknown = [], []
            for member in _as_list(members):
                try:
                    identifiers.append(self.get_member_id(member))
                except KeyError:
                    unknown.append(member)
            if unknown:
                raise ValueError(f"{self.id} has no stations {sorted(unknown)}")
            # in the order of the members
            table = table.set_index("id", drop=False).loc[identifiers]
        if region is not None:
            within = shapely.intersects_xy(
                region, table.longitude.to_numpy(), table.latitude.to_numpy()
            )
            if members is not None and not within.all():
                warnings.warn(
                    f"stations {list(table.id[~within])} are not within the region",
                    stacklevel=3,
                )
            table = table[within]
            if table.empty:
                warnings.warn(
                    f"no stations of {self.id} are within the region", stacklevel=3
                )
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
        members: str | Sequence[str] | None = None,
        region: BaseGeometry | None = None,
    ) -> list[RadarStation]:
        """
        Builds the radar stations of a configuration (see
        `select_configuration`), named by their identifiers (see
        `get_member_names`), at the elevations of their ground above the WGS 84
        ellipsoid (the table's height above mean sea level plus the EGM2008
        geoid height). Terrain masks are not built (see
        `tatc.preprocess.compute_terrain_mask_for_station`).

        Args:
            configuration (str | None): The configuration, if not selected by
                grade.
            grade (SystemGrade | str | None): The grade of the use, if any.
            members (str | Sequence[str] | None): The identifiers (or names)
                of the stations to build, in the order of the stations. By
                default, every station (in the order of the table).
            region (shapely.geometry.base.BaseGeometry | None): A region
                (longitude and latitude) within which to build stations.

        Returns:
            list[RadarStation]: The stations.

        Raises:
            ValueError: If a member is unknown.
        """
        selected = self.configurations[self._select(configuration, grade)]
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
        members: str | Sequence[str] | None = None,
        region: BaseGeometry | None = None,
        radius: float | None = None,
    ) -> list[SurfacePlatform]:
        """
        Builds the surface stations of a configuration (see
        `select_configuration`), named by their identifiers (such as ICAO
        codes; see `get_member_names`), at the elevations of their ground above
        the WGS 84 ellipsoid (the table's height above mean sea level plus the
        EGM2008 geoid height).

        Args:
            configuration (str | None): The configuration, if not selected by
                grade.
            grade (SystemGrade | str | None): The grade of the use, if any.
            members (str | Sequence[str] | None): The identifiers (or names)
                of the stations to build, in the order of the stations. By
                default, every station (in the order of the table).
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
        selected = self.configurations[self._select(configuration, grade)]
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

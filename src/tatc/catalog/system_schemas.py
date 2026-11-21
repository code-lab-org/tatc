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
from collections.abc import Sequence
from datetime import date, datetime, time
from typing import Annotated, Any, Literal

import numpy as np
import pandas as pd
from pydantic import Field, model_validator
from shapely.geometry import Point as ShapelyPoint
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
from ..schemas.surface.radar import RadarBand
from ..utils import compute_geoid_height
from .records import ValidationRecord, _CatalogModel, _get_days_apart
from .schemas import InstrumentSpecification, Provenance, _check_provenance
from .standard import EntryValidation
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
    is no longer maintained (and drifts), its inclination. Built at an epoch,
    with an arbitrary phase, as a `SunSynchronousOrbit` (whose inclination is
    sun-synchronous, so that its local time stays constant), or, with an
    inclination, as a `CircularOrbit` with the local time at the epoch, which
    drifts at the rate of its inclination during an analysis.
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

    def build(
        self, epoch: datetime, true_anomaly: float = 0
    ) -> SunSynchronousOrbit | CircularOrbit:
        """
        Builds the orbit at an epoch.

        Args:
            epoch (datetime): The epoch.
            true_anomaly (float): The true anomaly (degrees) at the epoch.

        Returns:
            SunSynchronousOrbit | CircularOrbit: The orbit.
        """
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


class InclinedPlane(_CatalogModel):
    """
    The plane of an inclined circular orbit, at the plane grade: its mean
    altitude and inclination. Built as a `CircularOrbit` at an epoch, with an
    arbitrary right ascension of ascending node and phase.
    """

    type: Literal["circular"] = Field(
        default="circular", description="Orbit discriminator."
    )
    mean_altitude: float = Field(
        ..., gt=0, description="Mean altitude (meters), as TAT-C defines it."
    )
    inclination: float = Field(..., gt=0, lt=180, description="Inclination (degrees).")

    def build(self, epoch: datetime, true_anomaly: float = 0) -> CircularOrbit:
        """
        Builds the orbit at an epoch, with a right ascension of ascending node of
        zero.

        Args:
            epoch (datetime): The epoch.
            true_anomaly (float): The true anomaly (degrees) at the epoch.

        Returns:
            CircularOrbit: The orbit.
        """
        return CircularOrbit(
            mean_altitude=self.mean_altitude,
            inclination=self.inclination,
            right_ascension_ascending_node=0,
            epoch=epoch,
            true_anomaly=true_anomaly,
        )


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

    def build(self, epoch: datetime, true_anomaly: float = 0) -> GeosynchronousOrbit:
        """
        Builds the orbit at an epoch.

        Args:
            epoch (datetime): The epoch.
            true_anomaly (float): Not used (the longitude sets the position).

        Returns:
            GeosynchronousOrbit: The orbit.
        """
        del true_anomaly
        return GeosynchronousOrbit(
            longitude=self.longitude, inclination=self.inclination, epoch=epoch
        )


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


def _get_values(orbit: Any) -> set[str]:
    """
    Gets the values of an orbit that need provenance: its fields, without its
    `type` or an unset (derived) inclination.
    """
    return {
        field
        for field, value in orbit
        if field != "type" and not (field == "inclination" and value is None)
    }


class _SystemConfiguration(_CatalogModel):
    """
    Base class of system configurations.
    """

    grade: SystemGrade = Field(..., description="The use for which it is validated.")
    description: str = Field(..., description="What the configuration represents.")
    limitations: tuple[str, ...] = Field(
        default=(), description="What it does not represent."
    )

    def _get_measured(self) -> list[Provenance]:
        """
        Gets the provenance of the values fitted to reference data.
        """
        raise NotImplementedError

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

    def _get_content(self) -> dict[str, Any]:
        """
        Gets the values that affect analysis results.
        """
        raise NotImplementedError

    @property
    def digest(self) -> str:
        """
        A digest of the values that affect analysis results, with which a
        validation record identifies the configuration it tested.
        """
        encoded = json.dumps(self._get_content(), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(encoded.encode("utf-8")).hexdigest()

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
        return self

    def get_members(self) -> list[str]:
        return list(self.orbits)

    def get_platform_holdout_gap(self, periods: dict[str, str]) -> int | None:
        """
        Gets the days between the data to which each spacecraft's values were
        fitted and the data of its test: the least, over spacecraft, of the days
        between the nearest days of its periods.

        Args:
            periods (dict[str, str]): The period of the data of each spacecraft's
                test (see `Period`), by identifier.

        Returns:
            int | None: The days, or `None` if no tested spacecraft has fitted
            values.
        """
        gaps = [
            _get_days_apart(provenance.period, period)
            for spacecraft, period in periods.items()
            for provenance in self.provenance.get(spacecraft, {}).values()
            if provenance.kind == "measured"
        ]
        return min(gaps) if gaps else None

    def _get_measured(self) -> list[Provenance]:
        return [
            provenance
            for values in self.provenance.values()
            for provenance in values.values()
            if provenance.kind == "measured"
        ]

    def _get_content(self) -> dict[str, Any]:
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
    provenance: dict[str, Provenance] = Field(
        ..., description="The origin of each value, by field."
    )
    table_digest: str = Field(
        default="",
        description="The digest of the stations' table, set by the system.",
        exclude=True,
    )
    table_provenance: Provenance | None = Field(
        default=None,
        description="The provenance of the stations' table, set by the system.",
        exclude=True,
    )

    @model_validator(mode="after")
    def _check_values(self) -> RadarConfiguration:
        _check_provenance(
            {
                "band",
                "max_range",
                "min_elevation_angle",
                "max_elevation_angle",
                "beam_width",
            },
            self.provenance,
        )
        if self.max_elevation_angle < self.min_elevation_angle:
            raise ValueError("max_elevation_angle is below min_elevation_angle")
        return self

    def _get_measured(self) -> list[Provenance]:
        provenances = list(self.provenance.values())
        if self.table_provenance is not None:
            provenances.append(self.table_provenance)
        return [p for p in provenances if p.kind == "measured"]

    def _get_content(self) -> dict[str, Any]:
        return {
            "radar": self.model_dump(
                mode="json",
                include={
                    "band",
                    "max_range",
                    "min_elevation_angle",
                    "max_elevation_angle",
                    "beam_width",
                },
            ),
            "table": self.table_digest,
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
    table_digest: str = Field(
        default="",
        description="The digest of the stations' table, set by the system.",
        exclude=True,
    )
    table_provenance: Provenance | None = Field(
        default=None,
        description="The provenance of the stations' table, set by the system.",
        exclude=True,
    )

    @model_validator(mode="after")
    def _check_values(self) -> SurfaceConfiguration:
        _check_provenance({"height"}, self.provenance)
        return self

    def _get_measured(self) -> list[Provenance]:
        provenances = list(self.provenance.values()) + [
            provenance
            for instrument in self.instruments
            for provenance in instrument.provenance.values()
        ]
        if self.table_provenance is not None:
            provenances.append(self.table_provenance)
        return [p for p in provenances if p.kind == "measured"]

    def _get_content(self) -> dict[str, Any]:
        return {
            "height": self.height,
            "instruments": [
                {"type": i.type, "parameters": i.parameters} for i in self.instruments
            ],
            "table": self.table_digest,
        }


# the columns of the tables of stations, by class (those after the first five
# are optional)
_TABLE_COLUMNS = {
    SystemClass.RADAR_NETWORK: (
        "id",
        "name",
        "latitude",
        "longitude",
        "ground",
        "height",
        "min_elevation_angle",
    ),
    SystemClass.SURFACE_NETWORK: (
        "id",
        "name",
        "latitude",
        "longitude",
        "ground",
        "wban",
    ),
}
_REQUIRED_COLUMNS = {
    SystemClass.RADAR_NETWORK: 6,
    SystemClass.SURFACE_NETWORK: 5,
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
    A system of a past or present program (a satellite program, a radar network,
    or a surface network), in configurations validated against reference data,
    from which to build the platforms of baseline systems in pre-phase A studies.
    """

    id: str = Field(..., description="Entry identifier.")
    name: str = Field(..., description="System name.")
    aliases: tuple[str, ...] = Field(default=(), description="Other names.")
    system_class: SystemClass = Field(..., description="System class.")
    selection: str = Field(..., description="Why the system is cataloged.")
    operator: str = Field(..., description="The organizations that operate it.")
    spacecraft: tuple[Spacecraft, ...] = Field(
        default=(), description="The satellites of a satellite system."
    )
    stations: StationTable | None = Field(
        default=None, description="The stations of a network."
    )
    default_configuration: str = Field(
        ..., description="The configuration used if none is specified."
    )
    configurations: dict[
        str, OrbitConfiguration | RadarConfiguration | SurfaceConfiguration
    ] = Field(..., min_length=1, description="The ways of modeling it, by identifier.")
    limitations: tuple[str, ...] = Field(
        default=(), description="What no configuration represents."
    )
    record: ValidationRecord | None = Field(
        default=None, description="The results of its validation."
    )

    @model_validator(mode="before")
    @classmethod
    def _parse_configurations(cls, data: Any) -> Any:
        """
        Parses the configurations as those of the class, with the digest of the
        stations' table of a network.
        """
        if not isinstance(data, dict) or not isinstance(
            data.get("configurations"), dict
        ):
            return data
        system_class = SystemClass(data.get("system_class"))
        if system_class.is_satellite():
            schema: type[_SystemConfiguration] = OrbitConfiguration
        elif system_class == SystemClass.RADAR_NETWORK:
            schema = RadarConfiguration
        else:
            schema = SurfaceConfiguration
        stations = data.get("stations")
        table = {}
        if stations is not None:
            stations = (
                stations
                if isinstance(stations, StationTable)
                else StationTable(**stations)
            )
            try:
                digest = stations.digest
            except FileNotFoundError:
                raise ValueError(f"no packaged table {stations.file!r}") from None
            table = {"table_digest": digest, "table_provenance": stations.provenance}
        configurations = {}
        for key, value in data["configurations"].items():
            if isinstance(value, _SystemConfiguration):
                value = value.model_dump()
            if not isinstance(value, dict):
                raise ValueError(f"{key}: a configuration is a mapping")
            value = {**value, **table}
            configurations[key] = schema(**value)
        return {**data, "configurations": configurations}

    @model_validator(mode="after")
    def _check_references(self) -> SystemEntry:
        if self.default_configuration not in self.configurations:
            raise ValueError("default_configuration is not a configuration")
        satellite = self.system_class.is_satellite()
        if satellite and (not self.spacecraft or self.stations is not None):
            raise ValueError("a satellite system has spacecraft and no stations")
        if not satellite and (self.spacecraft or self.stations is None):
            raise ValueError("a network has stations and no spacecraft")
        spacecraft = {s.id for s in self.spacecraft}
        if len(spacecraft) != len(self.spacecraft):
            raise ValueError("spacecraft identifiers must be distinct")
        for key, configuration in self.configurations.items():
            # checks that the class has gates at this grade
            get_system_gates(self.system_class, configuration.grade)
            if isinstance(configuration, OrbitConfiguration):
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
        if self.stations is not None:
            self._check_table()
        if self.record is not None:
            self._check_record(spacecraft)
        return self

    def _check_table(self) -> None:
        """
        Checks the columns and values of the table of stations.
        """
        table = self.stations.load()
        columns = _TABLE_COLUMNS[self.system_class]
        required = columns[: _REQUIRED_COLUMNS[self.system_class]]
        missing = [c for c in required if c not in table.columns]
        unknown = [c for c in table.columns if c not in columns]
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

    def _check_record(self, spacecraft: set[str]) -> None:
        """
        Checks that the record is of this system's configurations and platforms.
        """
        if self.record.entry_id != self.id:
            raise ValueError("the record belongs to a different system")
        unknown = set(self.record.configurations) - set(self.configurations)
        if unknown:
            raise ValueError(f"records of unknown {sorted(unknown)}")
        platforms = spacecraft if spacecraft else {self.id}
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

    def build(  # pylint: disable=too-many-arguments,too-many-locals
        self,
        configuration: str | None = None,
        *,
        grade: SystemGrade | str | None = None,
        members: Sequence[str] | None = None,
        region: BaseGeometry | None = None,
        epoch: datetime | None = None,
        payload: bool | Sequence[str] = True,
        time_step: float | None = None,
        start: datetime | None = None,
        end: datetime | None = None,
        radius: float | None = None,
    ) -> list[Satellite] | list[RadarStation] | list[SurfacePlatform]:
        """
        Builds the platforms of a configuration (see `select_configuration`):
        the satellites of a satellite system, or the stations of a network.

        Satellites are named by their spacecraft's `name` and carry the
        instruments of the instrument catalog's entries that list them (with
        `payload`, see `get_payload_references`), built with
        `tatc.catalog.build_payload` for the spacecraft. An orbit at the plane
        grade is built at the `epoch` (by default, the `start`), with an
        arbitrary phase (and, for an inclined orbit, node); an orbit maintained
        on a repeat ground track is pinned to its ascending node, at any epoch.

        Stations are named by their identifier (such as an ICAO code), at the
        elevation of their ground above the WGS 84 ellipsoid (the table's height
        above mean sea level plus the EGM2008 geoid height). Terrain masks are
        not built (see `tatc.preprocess.compute_terrain_mask_for_station`).

        Args:
            configuration (str | None): The configuration, if not selected by
                grade.
            grade (SystemGrade | str | None): The grade of the use, if any.
            members (Sequence[str] | None): The identifiers of the spacecraft or
                stations to build, in the order of the platforms. By default,
                the operational and extended spacecraft of the configuration, or
                every station (in the order of the table).
            region (shapely.geometry.base.BaseGeometry | None): A region
                (longitude and latitude) within which to build stations.
            epoch (datetime | None): The epoch of orbits at the plane grade or
                of geostationary orbits.
            payload (bool | Sequence[str]): Whether satellites carry their
                cataloged instruments, or the references of the instrument
                catalog to carry (such as `["viirs", "atms"]`).
            time_step (float | None): The analysis time step (seconds), for
                instruments that tailor a field of view to it.
            start (datetime | None): The start of the analysis, for instruments
                that schedule their tilt by day, and the default epoch.
            end (datetime | None): The end of the analysis, for instruments that
                schedule their tilt by day.
            radius (float | None): The radius (meters) of the area that the
                in-situ instruments of a surface network represent, for the
                variable of a study, instead of the cataloged assumption.

        Returns:
            list[Satellite] | list[RadarStation] | list[SurfacePlatform]: The
            platforms.

        Raises:
            ValueError: If an argument does not apply to the system, a member is
                unknown, or an orbit needs an epoch.
        """
        key = self.select_configuration(configuration, grade=grade)
        selected = self.configurations[key]
        if isinstance(selected, OrbitConfiguration):
            if region is not None or radius is not None:
                raise ValueError("region and radius apply to the stations of a network")
            return self._build_satellites(
                selected, members, epoch or start, payload, time_step, start, end
            )
        if epoch is not None or payload is not True or time_step is not None:
            raise ValueError(
                "epoch, payload, and time_step apply to the satellites of a "
                "satellite system"
            )
        if radius is not None and not isinstance(selected, SurfaceConfiguration):
            raise ValueError(
                "radius applies to the in-situ instruments of a surface network"
            )
        table = self.stations.load()
        if members is not None:
            unknown = set(members) - set(table["id"])
            if unknown:
                raise ValueError(f"{self.id} has no stations {sorted(unknown)}")
            # in the order of the members
            table = table.set_index("id", drop=False).loc[list(members)]
        if region is not None:
            table = table[
                [
                    region.intersects(ShapelyPoint(lon, lat))
                    for lon, lat in zip(table.longitude, table.latitude)
                ]
            ]
        ground = table.ground.to_numpy() + compute_geoid_height(
            table.longitude.to_numpy(), table.latitude.to_numpy()
        )
        if isinstance(selected, RadarConfiguration):
            return [
                _build_radar(selected, row, elevation)
                for row, elevation in zip(table.itertuples(), ground)
            ]
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

    def _build_satellites(  # pylint: disable=too-many-arguments,too-many-positional-arguments
        self,
        configuration: OrbitConfiguration,
        members: Sequence[str] | None,
        epoch: datetime | None,
        payload: bool | Sequence[str],
        time_step: float | None,
        start: datetime | None,
        end: datetime | None,
    ) -> list[Satellite]:
        """
        Builds the satellites of a configuration (see `build`).
        """
        if members is None:
            members = [
                s.id
                for s in self.spacecraft
                if s.id in configuration.orbits and s.status != "retired"
            ]
        satellites = []
        for spacecraft_id in members:
            spacecraft = self.get_spacecraft(spacecraft_id)
            if spacecraft_id not in configuration.orbits:
                raise ValueError(
                    f"{self.id} has no orbit of {spacecraft_id} in this configuration"
                )
            orbit = configuration.orbits[spacecraft_id]
            if not isinstance(orbit, RepeatGroundTrackOrbit):
                if epoch is None:
                    raise ValueError(
                        f"{spacecraft.name}'s orbit is built at an epoch (or start)"
                    )
                orbit = orbit.build(epoch)
            references = (
                self.get_payload_references(spacecraft_id)
                if payload is True
                else list(payload or [])
            )
            if not references:
                satellites.append(Satellite(name=spacecraft.name, orbit=orbit))
                continue
            # deferred import, as the loader imports this module
            # pylint: disable-next=import-outside-toplevel,cyclic-import
            from .loader import build_payload

            satellites.append(
                build_payload(
                    references,
                    orbit=orbit,
                    time_step=time_step,
                    platform=spacecraft_id,
                    start=start,
                    end=end,
                ).to_satellite(orbit, name=spacecraft.name)
            )
        return satellites

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
        # deferred import, as the loader imports this module
        # pylint: disable-next=import-outside-toplevel,cyclic-import
        from .loader import get_instrument_ids

        return get_instrument_ids(spacecraft_id, soundings=False)


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

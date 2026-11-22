"""
The standard against which catalog systems (satellite programs, radar networks,
and surface networks) are validated, and its evaluation of their records
(`evaluate_system`).

As for instruments (see `tatc.catalog.standard`), the acceptance thresholds are
frozen so that they test new systems rather than following their results:
changing one requires a new `SYSTEM_STANDARD_VERSION`. The grades are the uses a
pre-phase A study makes of a baseline system: the planes of its orbits, for the
statistics of revisits and coverage; the ground tracks of maintained orbits
(and the longitudes of geostationary ones), for where and roughly when passes
occur; and the siting and membership of stations, and the coverage of radars.
See `docs/validation/standard.rst`.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Collection
from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING

from .standard import (
    EntryValidation,
    Gate,
    _combine,
    _evaluate_record,
)

if TYPE_CHECKING:
    from .system_schemas import SystemEntry

SYSTEM_STANDARD_VERSION = 3
"""The version of the thresholds in `get_system_gates`."""

COMPATIBLE_SYSTEM_VERSIONS = (3,)
"""Versions whose results still apply."""


class SystemClass(str, Enum):
    """
    Enumeration of system classes, which select the applicable gates.
    """

    SUN_SYNCHRONOUS = "sun_synchronous"
    """Satellites in sun-synchronous orbits."""
    INCLINED = "inclined"
    """Satellites in other (non-sun-synchronous) low Earth orbits."""
    GEOSTATIONARY = "geostationary"
    """Satellites in geostationary orbits."""
    RADAR_NETWORK = "radar_network"
    """A network of weather radars."""
    SURFACE_NETWORK = "surface_network"
    """A network of in-situ surface stations."""

    def is_satellite(self) -> bool:
        """
        Checks whether systems of this class are satellites.

        Returns:
            bool: `True`, for a class of satellites.
        """
        return self in (
            SystemClass.SUN_SYNCHRONOUS,
            SystemClass.INCLINED,
            SystemClass.GEOSTATIONARY,
        )


class SystemGrade(str, Enum):
    """
    Enumeration of the uses for which a system configuration is validated.
    """

    PLANE = "plane"
    """
    The planes of the orbits (altitude, inclination, and the local time of the
    ascending node of sun-synchronous orbits), for statistics of revisits and
    coverage, with the relative phases (and nodes) of a constellation whose
    operator maintains them, which are tested if cataloged; other phases are
    arbitrary.
    """
    GROUND_TRACK = "ground_track"
    """
    As for the plane, with the ground track of an orbit maintained on a repeat
    ground track (or the longitude of a geostationary orbit), for where and, to
    about a minute, when passes occur, as tested over the weeks to months of
    its element sets (the operator's maintenance keeps it so for longer, but
    the tests do not show it).
    """
    SITING = "siting"
    """The positions of a network's stations, and its current membership."""
    COVERAGE = "coverage"
    """As for siting, with the volume a network's radars observe."""


# orbits: errors of mean elements (from element sets of the test period)
_ALTITUDE = Gate("altitude_error_p95", 10, True, "km", "Mean altitude error (p95)")
_INCLINATION = Gate(
    "inclination_error_p95", 0.1, True, "deg", "Inclination error (p95)"
)
_LOCAL_TIME = Gate(
    "local_time_error_p95", 10, True, "min", "Local time of ascending node error (p95)"
)
# maintained orbits: the ascending nodes of the element sets of the test period,
# compared with those of the cataloged orbit
_NODE_LONGITUDE = Gate(
    "node_longitude_error_p95", 10, True, "km", "Node longitude error (p95)"
)
_NODE_TIME = Gate("node_time_error_p95", 60, True, "s", "Node time error (p95)")
# constellations: the relative phases (and nodes) of the satellites, from the
# element sets of the days on which each satellite and the first have sets
_RELATIVE_PHASE = Gate(
    "relative_phase_error_p95", 5, True, "deg", "Relative phase error (p95)"
)
_RELATIVE_NODE = Gate(
    "relative_node_error_p95", 1, True, "deg", "Relative node error (p95)"
)
# orbits: the element sets tested of each spacecraft, enough for a 95th
# percentile and spanning at least a week
_TEST_SETS = Gate("test_sets", 3, False, "", "Element sets tested")
_TEST_SPAN = Gate("test_span_days", 7, False, "days", "Days spanned by the sets")
_SLOT_LONGITUDE = Gate(
    "slot_longitude_error_p95", 0.1, True, "deg", "Subsatellite longitude error (p95)"
)
# networks: positions compared with an independent source, and the stations
# reporting during the test period compared with the table
_RADAR_HORIZONTAL = Gate(
    "horizontal_error_p95", 50, True, "m", "Horizontal position error (p95)"
)
_RADAR_VERTICAL = Gate(
    "antenna_altitude_error_p95", 5, True, "m", "Antenna altitude error (p95)"
)
_SURFACE_HORIZONTAL = Gate(
    "horizontal_error_p95", 100, True, "m", "Horizontal position error (p95)"
)
_SURFACE_VERTICAL = Gate(
    "ground_elevation_error_p95", 5, True, "m", "Ground elevation error (p95)"
)
_MEMBERSHIP_RECALL = Gate(
    "membership_recall", 98, False, "%", "Reporting stations in the table"
)
_MEMBERSHIP_PRECISION = Gate(
    "membership_precision", 98, False, "%", "Stations of the table reporting"
)
# the stations of the table found in the source with which positions are
# compared, so that errors are not of a few stations
_MATCHED = Gate(
    "reference_matched", 90, False, "%", "Stations compared with the reference"
)

# radars: the volume coverage pattern of a configuration compared with the
# volumes the network's radars scanned with its patterns
_LOWEST_TILT = Gate(
    "min_elevation_error_p95", 0.2, True, "deg", "Lowest elevation angle error (p95)"
)
# one-sided: a radar that ends a volume below its highest cut (as AVSET does
# when few echoes are present) observes less than the configuration, as its
# limitations state, while one that scans above it observes more
_HIGHEST_TILT = Gate(
    "max_elevation_excess_p95",
    1.0,
    True,
    "deg",
    "Highest elevation angle above the configuration's (p95)",
)
# two-sided: the highest cut that each sampled volume's pattern plans
_PLANNED_HIGHEST_TILT = Gate(
    "planned_max_elevation_error_p95",
    1.0,
    True,
    "deg",
    "Highest planned elevation angle error (p95)",
)
_RANGE = Gate("range_error_p95", 5, True, "%", "Maximum range error (p95)")
_VOLUMES = Gate("coverage_volumes", 20, False, "", "Volumes sampled")

_SITING = {
    SystemClass.RADAR_NETWORK: (
        _RADAR_HORIZONTAL,
        _RADAR_VERTICAL,
        _MATCHED,
        _MEMBERSHIP_RECALL,
        _MEMBERSHIP_PRECISION,
    ),
    SystemClass.SURFACE_NETWORK: (
        _SURFACE_HORIZONTAL,
        _SURFACE_VERTICAL,
        _MATCHED,
        _MEMBERSHIP_RECALL,
        _MEMBERSHIP_PRECISION,
    ),
}
_SAMPLES = (_TEST_SETS, _TEST_SPAN)


@dataclass(frozen=True)
class _Grade:
    """
    A grade of a class: its gates, and the type of the orbits of a satellite
    system's configurations at it.
    """

    gates: tuple[Gate, ...]
    orbit_type: str | None = None


# the grades of each class
_GRADES: dict[tuple[SystemClass, SystemGrade], _Grade] = {
    (SystemClass.SUN_SYNCHRONOUS, SystemGrade.PLANE): _Grade(
        (_ALTITUDE, _INCLINATION, _LOCAL_TIME, *_SAMPLES), "sso"
    ),
    (SystemClass.SUN_SYNCHRONOUS, SystemGrade.GROUND_TRACK): _Grade(
        (
            _ALTITUDE,
            _INCLINATION,
            _LOCAL_TIME,
            _NODE_LONGITUDE,
            _NODE_TIME,
            *_SAMPLES,
        ),
        "repeat",
    ),
    (SystemClass.INCLINED, SystemGrade.PLANE): _Grade(
        (_ALTITUDE, _INCLINATION, *_SAMPLES), "circular"
    ),
    (SystemClass.INCLINED, SystemGrade.GROUND_TRACK): _Grade(
        (_ALTITUDE, _INCLINATION, _NODE_LONGITUDE, _NODE_TIME, *_SAMPLES), "repeat"
    ),
    (SystemClass.GEOSTATIONARY, SystemGrade.GROUND_TRACK): _Grade(
        (_SLOT_LONGITUDE, _INCLINATION, *_SAMPLES), "geosynchronous"
    ),
    (SystemClass.RADAR_NETWORK, SystemGrade.SITING): _Grade(
        _SITING[SystemClass.RADAR_NETWORK]
    ),
    (SystemClass.RADAR_NETWORK, SystemGrade.COVERAGE): _Grade(
        (
            *_SITING[SystemClass.RADAR_NETWORK],
            _LOWEST_TILT,
            _HIGHEST_TILT,
            _PLANNED_HIGHEST_TILT,
            _RANGE,
            _VOLUMES,
        )
    ),
    (SystemClass.SURFACE_NETWORK, SystemGrade.SITING): _Grade(
        _SITING[SystemClass.SURFACE_NETWORK]
    ),
}

# the grades each class needs to be validated
_REQUIRED_GRADES: dict[SystemClass, tuple[SystemGrade, ...]] = {
    SystemClass.SUN_SYNCHRONOUS: (SystemGrade.PLANE,),
    SystemClass.INCLINED: (SystemGrade.PLANE,),
    SystemClass.GEOSTATIONARY: (SystemGrade.GROUND_TRACK,),
    SystemClass.RADAR_NETWORK: (SystemGrade.COVERAGE,),
    SystemClass.SURFACE_NETWORK: (SystemGrade.SITING,),
}

# the grade that each grade also tests
_INCLUDES = {
    SystemGrade.GROUND_TRACK: SystemGrade.PLANE,
    SystemGrade.COVERAGE: SystemGrade.SITING,
}


RELATIVE_GATES: dict[str, Gate] = {
    "phase": _RELATIVE_PHASE,
    "right_ascension_ascending_node": _RELATIVE_NODE,
}
"""
The gates of the values of orbits relative to each other, by field, whose
metrics are of each spacecraft but the first tested of a configuration.
"""


def _get_grade(system_class: SystemClass | str, grade: SystemGrade | str) -> _Grade:
    """
    Gets a grade of a class.

    Raises:
        ValueError: If the class is not validated at the grade.
    """
    system_class, grade = SystemClass(system_class), SystemGrade(grade)
    if (system_class, grade) not in _GRADES:
        grades = [g.value for c, g in _GRADES if c == system_class]
        raise ValueError(
            f"{system_class.value} has no {grade.value} grade, only {grades}"
        )
    return _GRADES[system_class, grade]


def get_system_gates(
    system_class: SystemClass | str,
    grade: SystemGrade | str,
    relative: Collection[str] = (),
) -> tuple[Gate, ...]:
    """
    Gets the gates of a system configuration.

    Args:
        system_class (SystemClass | str): The class of the system.
        grade (SystemGrade | str): The grade of the configuration.
        relative (Collection[str]): The values of its orbits relative to each
            other that it gives (`phase`, `right_ascension_ascending_node`; see
            `OrbitConfiguration.get_relative_values`), whose gates apply too.

    Returns:
        tuple[Gate, ...]: The applicable gates.

    Raises:
        ValueError: If the class is not validated at the grade.
    """
    return _get_grade(system_class, grade).gates + tuple(
        gate for field, gate in RELATIVE_GATES.items() if field in relative
    )


def get_orbit_type(system_class: SystemClass | str, grade: SystemGrade | str) -> str:
    """
    Gets the type of the orbits of a satellite system's configurations at a
    grade (such as `sso` for the planes of sun-synchronous orbits).

    Args:
        system_class (SystemClass | str): The class of the system.
        grade (SystemGrade | str): The grade of the configuration.

    Returns:
        str: The `type` of the orbits.

    Raises:
        ValueError: If the class is not validated at the grade, or is not of
            satellites.
    """
    orbit_type = _get_grade(system_class, grade).orbit_type
    if orbit_type is None:
        raise ValueError(f"{SystemClass(system_class).value} has no orbits")
    return orbit_type


def get_required_system_grades(
    system_class: SystemClass | str,
) -> tuple[SystemGrade, ...]:
    """
    Gets the grades that a system needs to be validated: the plane grade for
    satellites in low Earth orbits (whose ground tracks are validated only for
    orbits maintained on one), the ground track grade for geostationary
    satellites, the coverage grade for radar networks, and the siting grade for
    surface networks.

    Args:
        system_class (SystemClass | str): The class of the system.

    Returns:
        tuple[SystemGrade, ...]: The required grades.
    """
    return _REQUIRED_GRADES[SystemClass(system_class)]


def serves_system(
    configuration_grade: SystemGrade | str, grade: SystemGrade | str
) -> bool:
    """
    Checks whether a system configuration at a grade serves a use: its own
    grade, or the grade it also tests (the plane, for a ground track; the
    siting, for coverage).

    Args:
        configuration_grade (SystemGrade | str): The grade of the configuration.
        grade (SystemGrade | str): The grade of the use.

    Returns:
        bool: `True`, if it serves the use.
    """
    configuration_grade, grade = SystemGrade(configuration_grade), SystemGrade(grade)
    return configuration_grade == grade or _INCLUDES.get(configuration_grade) == grade


def evaluate_system(entry: SystemEntry) -> EntryValidation:
    """
    Applies the standard to a system's records, as `evaluate_entry` does to an
    instrument's: each gate takes the worst value of its metric over the
    spacecraft (each of which needs a value of it: every spacecraft of the
    configuration, but the first tested for the gates of relative values), a
    configuration is validated only if every gate passes with evidence of the
    configuration as cataloged (its digest, a compatible version of the
    standard, and fitted values tested on data at least `MIN_HOLDOUT_GAP_DAYS`
    from their fitting data), and a system is only as strong as its weakest
    configuration and needs one that serves its required grade.

    Args:
        entry (SystemEntry): The system.

    Returns:
        EntryValidation: Its validation status.
    """
    configurations = {}
    for key, configuration in entry.configurations.items():
        record = None if entry.record is None else entry.record.configurations.get(key)
        gates = get_system_gates(
            entry.system_class,
            configuration.grade,
            configuration.get_relative_values(),
        )
        gaps, scopes = [], None
        if record is not None and configuration.get_tested_platforms():
            tested = [
                p for p in configuration.get_tested_platforms() if p in record.platforms
            ]
            relative = {gate.metric for gate in RELATIVE_GATES.values()}
            scopes = {
                gate.metric: tested[1:] if gate.metric in relative else tested
                for gate in gates
            }
        if record is not None:
            if record.standard_version not in COMPATIBLE_SYSTEM_VERSIONS:
                gaps.append(
                    f"version {record.standard_version} of the system standard, "
                    "of the test, does not apply"
                )
            untested = sorted(
                set(configuration.get_tested_platforms()) - set(record.platforms)
            )
            if untested:
                gaps.append(f"not tested on {', '.join(untested)}")
        configurations[key] = _evaluate_record(
            configuration.grade,
            record,
            gates,
            digest=entry.get_digest(key),
            holdout=(
                None
                if record is None
                else entry.get_holdout_gap(
                    key,
                    record.period,
                    record.platform_periods,
                    platforms=record.platforms,
                )
            ),
            evidence_gaps=tuple(gaps),
            scopes=scopes,
        )
    missing = tuple(
        grade
        for grade in get_required_system_grades(entry.system_class)
        if not any(c.serves(grade) for c in entry.configurations.values())
    )
    return _combine(configurations, missing)

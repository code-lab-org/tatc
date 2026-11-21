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

SYSTEM_STANDARD_VERSION = 1
"""The version of the thresholds in `get_system_gates`."""

COMPATIBLE_SYSTEM_VERSIONS = (1,)
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
    coverage; the phase of the satellites in their planes is arbitrary.
    """
    GROUND_TRACK = "ground_track"
    """
    As for the plane, with the ground track of an orbit maintained on a repeat
    ground track (or the longitude of a geostationary orbit), for where and, to
    about a minute, when passes occur, for years.
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

# radars: the volume coverage pattern of a configuration compared with the
# volumes the network's radars scanned with its patterns
_LOWEST_TILT = Gate(
    "min_elevation_error_p95", 0.2, True, "deg", "Lowest elevation angle error (p95)"
)
_HIGHEST_TILT = Gate(
    "max_elevation_error_p95", 1.0, True, "deg", "Highest elevation angle error (p95)"
)
_RANGE = Gate("range_error_p95", 5, True, "%", "Maximum range error (p95)")

_SITING = {
    SystemClass.RADAR_NETWORK: (
        _RADAR_HORIZONTAL,
        _RADAR_VERTICAL,
        _MEMBERSHIP_RECALL,
        _MEMBERSHIP_PRECISION,
    ),
    SystemClass.SURFACE_NETWORK: (
        _SURFACE_HORIZONTAL,
        _SURFACE_VERTICAL,
        _MEMBERSHIP_RECALL,
        _MEMBERSHIP_PRECISION,
    ),
}

# the gates of each class at each of its grades
_GATES: dict[tuple[SystemClass, SystemGrade], tuple[Gate, ...]] = {
    (SystemClass.SUN_SYNCHRONOUS, SystemGrade.PLANE): (
        _ALTITUDE,
        _INCLINATION,
        _LOCAL_TIME,
    ),
    (SystemClass.SUN_SYNCHRONOUS, SystemGrade.GROUND_TRACK): (
        _ALTITUDE,
        _INCLINATION,
        _LOCAL_TIME,
        _NODE_LONGITUDE,
        _NODE_TIME,
    ),
    (SystemClass.INCLINED, SystemGrade.PLANE): (_ALTITUDE, _INCLINATION),
    (SystemClass.INCLINED, SystemGrade.GROUND_TRACK): (
        _ALTITUDE,
        _INCLINATION,
        _NODE_LONGITUDE,
        _NODE_TIME,
    ),
    (SystemClass.GEOSTATIONARY, SystemGrade.GROUND_TRACK): (
        _SLOT_LONGITUDE,
        _INCLINATION,
    ),
    (SystemClass.RADAR_NETWORK, SystemGrade.SITING): _SITING[SystemClass.RADAR_NETWORK],
    (SystemClass.RADAR_NETWORK, SystemGrade.COVERAGE): (
        *_SITING[SystemClass.RADAR_NETWORK],
        _LOWEST_TILT,
        _HIGHEST_TILT,
        _RANGE,
    ),
    (SystemClass.SURFACE_NETWORK, SystemGrade.SITING): _SITING[
        SystemClass.SURFACE_NETWORK
    ],
}

# the grade each class needs to be validated
_REQUIRED_GRADES: dict[SystemClass, SystemGrade] = {
    SystemClass.SUN_SYNCHRONOUS: SystemGrade.PLANE,
    SystemClass.INCLINED: SystemGrade.PLANE,
    SystemClass.GEOSTATIONARY: SystemGrade.GROUND_TRACK,
    SystemClass.RADAR_NETWORK: SystemGrade.COVERAGE,
    SystemClass.SURFACE_NETWORK: SystemGrade.SITING,
}

# the grade that each grade also tests
_INCLUDES = {
    SystemGrade.GROUND_TRACK: SystemGrade.PLANE,
    SystemGrade.COVERAGE: SystemGrade.SITING,
}


def get_system_gates(
    system_class: SystemClass | str, grade: SystemGrade | str
) -> tuple[Gate, ...]:
    """
    Gets the gates of a system configuration.

    Args:
        system_class (SystemClass | str): The class of the system.
        grade (SystemGrade | str): The grade of the configuration.

    Returns:
        tuple[Gate, ...]: The applicable gates.

    Raises:
        ValueError: If the class is not validated at the grade.
    """
    system_class, grade = SystemClass(system_class), SystemGrade(grade)
    if (system_class, grade) not in _GATES:
        grades = [g.value for c, g in _GATES if c == system_class]
        raise ValueError(
            f"{system_class.value} has no {grade.value} grade, only {grades}"
        )
    return _GATES[system_class, grade]


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
    return (_REQUIRED_GRADES[SystemClass(system_class)],)


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
    instrument's: each gate takes the worst value of its metric over scopes
    (such as satellites), a configuration is validated only if every gate
    passes with evidence of the configuration as cataloged (its digest, a
    compatible version of the standard, and fitted values tested on data at
    least `MIN_HOLDOUT_GAP_DAYS` from their fitting data, tested on every
    spacecraft of the configuration), and a system is only as strong as its
    weakest configuration and needs one that serves its required grade.

    Args:
        entry (SystemEntry): The system.

    Returns:
        EntryValidation: Its validation status.
    """
    configurations = {}
    for key, configuration in entry.configurations.items():
        record = None if entry.record is None else entry.record.configurations.get(key)
        gaps = []
        if record is not None:
            if record.standard_version not in COMPATIBLE_SYSTEM_VERSIONS:
                gaps.append(
                    f"version {record.standard_version} of the system standard, "
                    "of the test, does not apply"
                )
            untested = sorted(set(configuration.get_members()) - set(record.platforms))
            if untested:
                gaps.append(f"not tested on {', '.join(untested)}")
        configurations[key] = _evaluate_record(
            configuration,
            record,
            get_system_gates(entry.system_class, configuration.grade),
            tuple(gaps),
        )
    missing = tuple(
        grade
        for grade in get_required_system_grades(entry.system_class)
        if not any(c.serves(grade) for c in entry.configurations.values())
    )
    return _combine(configurations, missing)

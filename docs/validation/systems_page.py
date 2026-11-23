"""
Writes the documentation page of the catalog's systems (satellite programs,
radar networks, and surface networks) from their entries and validation records,
so that the page cannot differ from them. `docs/conf.py` runs it on each build.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import pathlib

from catalog_page import REPOSITORY, escape, heading, list_table, render_sourced

from tatc.catalog import (
    SYSTEM_STANDARD_VERSION,
    SatelliteSystem,
    get_system,
    get_system_ids,
)
from tatc.catalog.standard import MIN_HOLDOUT_GAP_DAYS
from tatc.catalog.system_orbits import describe_orbit
from tatc.catalog.system_schemas import (
    GroundConfiguration,
    OrbitConfiguration,
    RadarConfiguration,
)


def render_configuration(entry, key: str) -> str:
    """
    Writes the section of a configuration of a system.

    Args:
        entry (SystemEntry): The system.
        key (str): The identifier of the configuration.

    Returns:
        str: The section.
    """
    configuration = entry.configurations[key]
    validation = entry.validation.configurations[key]
    record = validation.record
    default = " (default)" if key == entry.default_configuration else ""
    parts = [
        heading(f"``{key}``{default}", "-"),
        f"**Tier:** `{validation.tier.value}`. **Grade:** "
        f"`{configuration.grade.value}`.\n",
        escape(configuration.description) + "\n",
    ]
    if isinstance(configuration, OrbitConfiguration):
        names = {s.id: s.name for s in entry.spacecraft}
        parts.append(
            list_table(
                ["Spacecraft", "Orbit"],
                [
                    [escape(names[sid]), escape(describe_orbit(orbit))]
                    for sid, orbit in configuration.orbits.items()
                ],
                [25, 75],
            )
        )
    elif isinstance(configuration, RadarConfiguration):
        parts.append(
            f"{configuration.band.value} band, {configuration.max_range / 1e3:g} km, "
            f"{configuration.min_elevation_angle:g} to "
            f"{configuration.max_elevation_angle:g} deg, beam width "
            f"{configuration.beam_width:g} deg.\n"
        )
    elif isinstance(configuration, GroundConfiguration):
        members = ", ".join(configuration.members) or "every station"
        parts.append(
            f"Stations: {escape(members)}; minimum elevation angle "
            f"{configuration.min_elevation_angle:g} deg.\n"
        )
    else:
        parts.append(
            f"Instruments {configuration.height:g} m above the ground: "
            + ", ".join(
                f"{escape(spec.parameters['name'])} (``{spec.type}``)"
                for spec in configuration.instruments
            )
            + ".\n"
        )
    if record is not None:
        if validation.holdout_gap_days is None:
            gap = "no fitted values"
        elif validation.holdout_gap_days == 0:
            gap = "no holdout"
        else:
            gap = f"{validation.holdout_gap_days} days from the fitting data"
        source = f"`{record.source} <{REPOSITORY}{record.source}>`_"
        parts.append(
            f"**Reference:** {escape(record.reference)}, {record.period} "
            f"({gap}). **Results:** {source}, TAT-C {record.tatc_version}, version "
            f"{record.standard_version} of the system standard.\n"
        )
    rows = []
    for result in validation.gates:
        gate = result.gate
        if result.value is None:
            value, outcome = "not tested", "untested"
        else:
            value = gate.format_value(result.value)
            scopes = [
                m.scope
                for m in record.metrics
                if m.name == gate.metric and m.value == result.value and m.scope
            ]
            if scopes:
                value += f" ({escape(scopes[0])})"
            outcome = "pass" if result.passed else "**fail**"
        rows.append([gate.label, value, gate.format_limit(), outcome])
    parts.append("Gates, with the worst result over the spacecraft tested:\n")
    parts.append(
        list_table(["Metric", "Result", "Limit", "Outcome"], rows, [40, 30, 15, 15])
    )
    if validation.evidence_gaps:
        parts.append("Other gaps:\n")
        parts.append(
            "\n".join(f"* {escape(text)}" for text in validation.evidence_gaps) + "\n"
        )
    if configuration.limitations:
        parts.append("Limitations:\n")
        parts.append(
            "\n".join(f"* {escape(text)}" for text in configuration.limitations) + "\n"
        )
    return "\n".join(parts)


def render_entry(entry) -> str:
    """
    Writes the section of a system.

    Args:
        entry (SystemEntry): The system.

    Returns:
        str: The section.
    """
    validation = entry.validation
    parts = [
        f".. _system-{entry.id}:\n",
        heading(entry.name, "="),
        f'**System:** ``get_system("{entry.id}")``. **Class:** '
        f"``{entry.system_class.value}``. **Tier:** `{validation.tier.value}`. "
        f"**Operator:** {escape(entry.operator)}."
        + {
            "transmitting": " **Role:** transmitting (navigation satellites, the "
            "transmitters of radio occultations).\n",
            "relaying": " **Role:** relaying (data relay satellites, built as a "
            "relay with ``build_relay``).\n",
        }.get(getattr(entry, "role", "observing"), "\n"),
    ]
    if entry.aliases:
        parts.append("**Also known as:** " + escape(", ".join(entry.aliases)) + ".\n")
    parts.append(escape(entry.selection) + "\n")
    if isinstance(entry, SatelliteSystem):
        parts.append(
            list_table(
                ["Spacecraft", "NORAD ID", "Status"],
                [
                    [
                        f"{escape(s.name)} (``{s.id}``)",
                        str(s.norad_id),
                        f"{s.status} as of {s.as_of.isoformat()}",
                    ]
                    for s in entry.spacecraft
                ],
                [40, 20, 40],
            )
        )
        groups = []  # spacecraft of the same resources, in order
        for spacecraft in entry.spacecraft:
            if spacecraft.resources is None:
                continue
            for resources, names in groups:
                if resources == spacecraft.resources:
                    names.append(spacecraft.name)
                    break
            else:
                groups.append((spacecraft.resources, [spacecraft.name]))
        for resources, names in groups:
            parts.append(
                f"Resources of {escape(', '.join(names))} (nominal, and not "
                "validated):\n"
            )
            parts.append(render_sourced([resources]))
    else:
        table = entry.stations
        path = f"src/tatc/resources/systems/tables/{table.file}"
        parts.append(
            f"**Stations:** {len(entry.get_members())}, as of {table.as_of.isoformat()} "
            f"(`{table.file} <{REPOSITORY}{path}>`_), from: "
            f"{escape(table.provenance.source)}.\n"
        )
    if validation.missing_grades:
        parts.append(
            "Gaps of the system:\n\n"
            + "\n".join(
                f"* no configuration at the {grade.value} grade"
                for grade in validation.missing_grades
            )
            + "\n"
        )
    if entry.limitations:
        parts.append(
            "Limitations of every configuration:\n\n"
            + "\n".join(f"* {escape(text)}" for text in entry.limitations)
            + "\n"
        )
    for key in entry.configurations:
        parts.append(render_configuration(entry, key))
    return "\n".join(parts)


def render() -> str:
    """
    Writes the systems page.

    Returns:
        str: The page, as reStructuredText.
    """
    entries = sorted(
        (get_system(system_id) for system_id in get_system_ids()),
        key=lambda e: (e.system_class.value, e.name.lower()),
    )
    summary = list_table(
        ["System", "Class", "Tier", "Configurations", "Members"],
        [
            [
                f":ref:`{escape(e.name)} <system-{e.id}>`",
                f"``{e.system_class.value}``",
                f"`{e.validation.tier.value}`",
                ", ".join(
                    f"``{key}`` ({c.grade.value}, "
                    f"{e.validation.configurations[key].tier.value})"
                    for key, c in e.configurations.items()
                ),
                (
                    ", ".join(escape(s.name) for s in e.spacecraft)
                    if isinstance(e, SatelliteSystem)
                    else f"{len(e.get_members())} stations"
                ),
            ]
            for e in entries
        ],
        [20, 14, 10, 30, 26],
    )
    intro = f"""\
.. This page is written by systems_page.py on each build. Do not edit it.

==============
System Catalog
==============

The systems of :mod:`tatc.catalog` are satellite programs, radar networks, and
surface networks, from which to build the platforms of baseline systems in
pre-phase A studies, such as the programs of record that a proposed system adds
to. Each system states what was tested, against what reference data, and how
well. This page is written from the systems and their records when the
documentation is built.

A system's tier is computed from its results and the thresholds of the
:ref:`system standard <system-standard>` (version {SYSTEM_STANDARD_VERSION}),
and cannot be declared. Values fitted to reference data are tested on data at
least {MIN_HOLDOUT_GAP_DAYS} days from their fitting data.

.. code-block:: python

    from tatc.catalog import get_system, list_systems

    jpss = get_system("jpss")
    jpss.validation  # tier, gate results, and gaps of each configuration
    # the planes of the orbits, at an epoch, with their cataloged instruments
    satellites = jpss.build(epoch=start, time_step=10)
    # orbits maintained on JPSS's repeat ground track, at any epoch
    satellites = jpss.build("ground_track", time_step=10)
    radars = get_system("nexrad").build(region=region)
    stations = get_system("uscrn").build(radius=10000)

Satellites carry the instruments of the instrument catalog (:doc:`catalog`) that
list them, and instrument entries list their satellites by the systems'
spacecraft identifiers. Orbits are idealized, so that they do not age as
element sets do: an orbit at the plane grade has the relative phases of a
constellation whose operator maintains them, or else phases (and, if inclined,
nodes) spaced evenly among the satellites built, so it represents the
statistics of revisits and coverage; an orbit at the ground track grade follows the maintained repeat ground track, so
it also represents where, and to about a minute when, passes occur. Build orbits
from current element sets (by NORAD catalog number) for passes to seconds.

Stations are built from dated tables, at their ground elevations above the WGS
84 ellipsoid (heights above mean sea level plus the EGM2008 geoid height), with
the scanning of a radar network's configuration or the instruments of a surface
network's. The footprint of an in-situ instrument (``InSituInstrument``) is a
disk of an assumed radius, which a study should set for its variable. Aircraft
and their sensors are not cataloged here: HyPlan catalogs them, and
:mod:`tatc.preprocess` adapts them (``aircraft_from_hyplan``,
``trajectory_from_waypoints``, and ``instrument_from_sensor``).

"""
    body = "\n".join(render_entry(entry) for entry in entries)
    return intro + summary + "\n" + body


def write_systems_page(path: pathlib.Path) -> None:
    """
    Writes the systems page, if it changed, so that unchanged builds stay cached.

    Args:
        path (pathlib.Path): The file to write.
    """
    text = render()
    if not path.exists() or path.read_text(encoding="utf-8") != text:
        path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    write_systems_page(pathlib.Path(__file__).parent / "systems.rst")

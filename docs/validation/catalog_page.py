"""
Writes the documentation page of the instrument catalog from its entries and
validation records, so that the page cannot differ from them. `docs/conf.py` runs
it on each build.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import pathlib
import re

from tatc.catalog import STANDARD_VERSION, get_entry, get_entry_ids
from tatc.catalog.standard import MIN_HOLDOUT_GAP_DAYS

REPOSITORY = "https://github.com/code-lab-org/tatc/blob/main/"


def escape(text: str) -> str:
    """
    Escapes reStructuredText markup in plain text.

    Args:
        text (str): The text.

    Returns:
        str: The text, shown as written.
    """
    text = re.sub(r"([\\*`|])", r"\\\1", text)
    return re.sub(r"_(?=\s|$|[.,;:)])", r"\\_", text)


def list_table(header: list[str], rows: list[list[str]], widths: list[int]) -> str:
    """
    Writes a reStructuredText list table.

    Args:
        header (list[str]): The column headings.
        rows (list[list[str]]): The cells, as reStructuredText.
        widths (list[int]): The relative widths of the columns.

    Returns:
        str: The table.
    """
    lines = [
        ".. list-table::",
        "   :header-rows: 1",
        f"   :widths: {' '.join(str(w) for w in widths)}",
        "",
    ]
    for row in [header, *rows]:
        for index, cell in enumerate(row):
            prefix = "   * - " if index == 0 else "     - "
            lines.append(prefix + cell if cell else prefix.rstrip())
    return "\n".join(lines) + "\n"


# the label and unit of each sourced value, in the order they are shown
SOURCED_VALUES = {
    "reference_altitude": ("Reference altitude", 1e3, "km"),
    "ground_sample_distance": ("Ground sample distance", 1, "m"),
    "ground_sample_distance_exponent": ("Exponent of ground sample distance", 1, ""),
    "signal_to_noise_ratio": ("Signal-to-noise ratio", 1, ""),
    "signal_to_noise_exponent": ("Exponent of signal-to-noise ratio", 1, ""),
    "mass": ("Mass", 1, "kg"),
    "power": ("Power", 1, "W"),
    "data_rate": ("Data rate", 1e6, "Mbit/s"),
    "launch_mass": ("Launch mass", 1, "kg"),
    "dry_mass": ("Dry mass", 1, "kg"),
    "design_life": ("Design life", 1, "years"),
    "storage": ("Storage", 1e9, "Gbit"),
    "downlink_rate": ("Downlink rate", 1e6, "Mbit/s"),
}


def render_sourced(sourced: list) -> str:
    """
    Writes a table of sourced values (performance and resources), with the
    origin of each.

    Args:
        sourced (list): The models of sourced values, or `None`.

    Returns:
        str: The table, or an empty string without values.
    """
    rows = []
    for model in sourced:
        if model is None:
            continue
        for field, provenance in model.provenance.items():
            label, scale, unit = SOURCED_VALUES[field]
            value = f"{getattr(model, field) / scale:g}" + (f" {unit}" if unit else "")
            rows.append(
                [label, value, f"{provenance.kind}: {escape(provenance.source)}"]
            )
    if not rows:
        return ""
    return list_table(["Value", "Amount", "Source"], rows, [25, 15, 60])


def heading(text: str, character: str) -> str:
    """
    Writes a section heading.

    Args:
        text (str): The title.
        character (str): The underline character.

    Returns:
        str: The heading.
    """
    return f"{text}\n{character * len(text)}\n"


def render_configuration(entry, key: str) -> str:
    """
    Writes the section of a configuration of an entry.

    Args:
        entry (InstrumentEntry): The entry.
        key (str): The identifier of the configuration.

    Returns:
        str: The section.
    """
    configuration = entry.configurations[key]
    validation = entry.validation.configurations[key]
    record = validation.record
    names = {spacecraft.id: spacecraft.name for spacecraft in entry.get_spacecraft()}
    default = " (default)" if key == entry.default_configuration else ""
    parts = [
        heading(f"``{key}``{default}", "-"),
        f"**Tier:** `{validation.tier.value}`. **Grade:** "
        f"`{configuration.grade.value}`. **Tested on:** "
        + ", ".join(names[p] for p in configuration.platforms)
        + ".\n",
        escape(configuration.description) + "\n",
    ]
    if configuration.grouping == "union":
        parts.append(
            "Its instruments are a union: analyze them with "
            "``instrument_index=None``.\n"
        )
    elif len(configuration.instruments) > 1:
        parts.append(
            f"Its instruments are {configuration.grouping}; analyses use "
            f"``instrument_index={configuration.instrument_index}`` unless "
            "another is given.\n"
        )
    if configuration.tilt_schedule is not None:
        parts.append(
            "Its tilt is scheduled by day: build it with the start and end of the "
            "analysis (``start`` and ``end``).\n"
        )
    if configuration.sounding is not None and configuration.sounding.type == "limb":
        parts.append(
            "It is a limb sounding: collect its observations with "
            "``build(...).sounding.collect_observations(satellite, start, end)``.\n"
        )
    elif configuration.sounding is not None:
        parts.append(
            "It is a radio occultation sounding: collect its observations with "
            "``build(...).sounding.collect_observations(receivers, transmitters, "
            "start, end)``, with the satellites of the constellations it tracks as "
            "transmitters (from the catalog's navigation systems, "
            "``sounding.get_transmitters(start)``).\n"
        )
    if record is not None:
        if validation.holdout_gap_days is None:
            gap = "no fitted values"
        elif validation.holdout_gap_days == 0:
            gap = "no holdout"
        else:
            gap = f"{validation.holdout_gap_days} days from the fitting data"
        source = f"`{record.source} <{REPOSITORY}{record.source}>`_"
        tested = ", ".join(names.get(p, p) for p in record.platforms)
        parts.append(
            f"**Reference:** {escape(record.reference)}, {record.period} "
            f"({escape(tested)}; {gap}). **Results:** {source}, "
            f"TAT-C {record.tatc_version}, version {record.standard_version} of "
            "the standard.\n"
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
        rows.append(
            [
                gate.label,
                value,
                gate.format_limit(),
                outcome,
            ]
        )
    parts.append(
        "Gates, with the worst result over the scopes tested (such as "
        "instruments and time steps):\n"
    )
    parts.append(
        list_table(["Metric", "Result", "Limit", "Outcome"], rows, [34, 36, 15, 15])
    )
    if record is not None and record.transfers:
        parts.append(
            "**Transfers**, results on other platforms, which do not count toward "
            "its tier: "
            + ", ".join(
                f"{escape(names.get(p, p))} ({transfer.period})"
                for p, transfer in record.transfers.items()
            )
            + ".\n"
        )
    if record is not None:
        path = f"src/tatc/resources/validation/{entry.id}.json"
        parts.append(
            f"Every result, including those not gated, is in the record "
            f"`{entry.id}.json <{REPOSITORY}{path}>`_.\n"
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
    table = render_sourced([configuration.performance, configuration.resources])
    if table:
        parts.append("Performance and resources, which replace the entry's:\n")
        parts.append(table)
    return "\n".join(parts)


def render_entry(entry) -> str:
    """
    Writes the section of an entry.

    Args:
        entry (InstrumentEntry): The entry.

    Returns:
        str: The section.
    """
    validation = entry.validation
    parts = [
        f".. _catalog-{entry.id}:\n",
        heading(entry.name, "="),
        f'**Entry:** ``get_entry("{entry.id}")``. **Class:** '
        f"``{entry.instrument_class.value}``. **Tier:** `{validation.tier.value}`.\n",
    ]
    if entry.aliases:
        parts.append("**Also known as:** " + escape(", ".join(entry.aliases)) + ".\n")
    parts.append(escape(entry.selection) + "\n")
    parts.append(
        list_table(
            ["Platform", "NORAD ID", "Status"],
            [
                [
                    f"{p.name} (``{p.id}``)",
                    str(p.norad_id),
                    f"{p.status} as of {p.as_of.isoformat()}",
                ]
                for p in entry.get_spacecraft()
            ],
            [40, 20, 40],
        )
    )
    if entry.data_products:
        parts.append(
            "**Data products:** "
            + "; ".join(
                f"`{escape(d.name)} <{d.url}>`_ ({escape(d.provider)}, {d.access})"
                for d in entry.data_products
            )
            + ".\n"
        )
    if validation.missing_grades:
        parts.append(
            "Gaps of the entry:\n\n"
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
    table = render_sourced([entry.performance, entry.resources])
    if table:
        parts.append(
            "Performance and resources (nominal, and not validated), unless a "
            "configuration replaces them:\n"
        )
        parts.append(table)
    for key in entry.configurations:
        parts.append(render_configuration(entry, key))
    return "\n".join(parts)


def render() -> str:
    """
    Writes the catalog page.

    Returns:
        str: The page, as reStructuredText.
    """
    entries = sorted(
        (get_entry(entry_id) for entry_id in get_entry_ids()),
        key=lambda e: (-e.validation.tier.rank, e.name.lower()),
    )
    summary = list_table(
        ["Instrument", "Class", "Tier", "Configurations", "Platforms"],
        [
            [
                f":ref:`{escape(e.name)} <catalog-{e.id}>`",
                f"``{e.instrument_class.value}``",
                f"`{e.validation.tier.value}`",
                ", ".join(
                    f"``{key}`` ({c.grade.value}, "
                    f"{e.validation.configurations[key].tier.value})"
                    for key, c in e.configurations.items()
                ),
                ", ".join(p.name for p in e.get_spacecraft()),
            ]
            for e in entries
        ],
        [16, 14, 12, 34, 24],
    )
    intro = f"""\
.. This page is written by catalog_page.py on each build. Do not edit it.

==================
Instrument Catalog
==================

The entries of :mod:`tatc.catalog` are instruments of past and present missions,
in configurations compared with their data products. Each entry gives its
instruments with the attitude they were tested with, and states what was
tested, against what data, and how well. This page is written from the entries
and their records when the documentation is built.

An entry's tier is computed from its results and the thresholds of
:doc:`standard` (version {STANDARD_VERSION}), and cannot be declared. Values fitted
to reference data are tested on data at least {MIN_HOLDOUT_GAP_DAYS} days from
their fitting data. A `partial` entry can still be useful within its
limitations, which each configuration lists below; a `validated` entry also has
limitations.

.. code-block:: python

    from tatc.catalog import get_entry, list_entries

    list_entries(tier="validated")
    entry = get_entry("atms")
    entry.validation  # tier, gate results, and gaps of each configuration
    satellite = entry.to_satellite(orbit, time_step=10, platform="noaa-21")
    beams = entry.build(orbit=orbit, grade="sample_size", platform="noaa-21")
    collect_observations(point, beams.to_satellite(orbit, name="NOAA-21"),
                         start, end, instrument_index=beams.instrument_index)

Name the platform analyzed: it selects the configuration tested on it, and
using a configuration on a platform on which it was not tested gives a
warning. Where an entry has configurations fitted on more than one platform
(ATMS, CrIS, VIIRS, the OMPS Nadir Mapper, CERES, and TROPICS), each configuration's
name ends with its platform, such as ``swath_noaa_21``, and the name without
it (``swath``) selects the one tested on the platform. A grade (``coverage``,
``timing``, or ``sample_size``) selects a configuration for that use. Analyze
a configuration of several instruments with the ``instrument_index`` that
``build`` gives, as ``to_satellite`` builds only configurations analyzed with
the default index of 0. Pass the analysis time
step to configurations whose along-track field of view is tailored to it, and
keep the attitude of the configuration: a scanner modeled in another velocity
frame can misplace its swath edges by tens of kilometers. Pass the start and
end of the analysis to configurations whose tilt is scheduled by day, such as
OCI's. The entries hold no element sets, as they age; supply a recent one for
the analysis period.

Limb sounders (MLS and SABER) and radio occultation receivers (TGRS, of
COSMIC-2) have no views to analyze: their configurations hold the sounding with
which their analysis is called, as ``build(...).sounding.collect_observations``,
and ``to_satellite`` builds the satellite that carries them.

"""
    body = "\n".join(render_entry(entry) for entry in entries)
    return intro + summary + "\n" + body


def write_catalog_page(path: pathlib.Path) -> None:
    """
    Writes the catalog page, if it changed, so that unchanged builds stay cached.

    Args:
        path (pathlib.Path): The file to write.
    """
    text = render()
    if not path.exists() or path.read_text(encoding="utf-8") != text:
        path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    write_catalog_page(pathlib.Path(__file__).parent / "catalog.rst")

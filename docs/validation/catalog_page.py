"""
Writes the documentation page of the instrument catalog from its entries and
validation records, so that the page cannot differ from them. `docs/conf.py` runs
it on each build.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import pathlib
import re

from tatc.catalog import (
    MIN_HOLDOUT_GAP_DAYS,
    STANDARD_VERSION,
    Tier,
    get_entry,
    get_entry_ids,
)

REPOSITORY = "https://github.com/code-lab-org/tatc/blob/main/"
METRIC_LABELS = {
    "pass_recall": "Pass recall",
    "pass_precision": "Pass precision",
    "epoch_p95": "Epoch, 95th percentile",
    "swath_width_median": "Swath width, median",
    "swath_edge_p95": "Swath edges, 95th percentile",
    "footprint_iou_median": "Footprint or swath IoU, median",
    "footprint_iou_p5": "Footprint or swath IoU, 5th percentile",
    "pixel_position_median": "Pixel position, median",
    "sample_size_within_3pct": "Dimensions within 3%",
    "product_limit": "Product limit of local zenith angle",
}
DECIMALS = {"%": 2, "s": 2, "km": 2, "deg": 2, "": 4}
ORDER = {Tier.VALIDATED: 0, Tier.PARTIAL: 1, Tier.NOMINAL: 2}


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


def format_value(value: float, gate) -> str:
    """
    Formats a metric value, with enough decimals that its rounded value meets
    its gate exactly when the value does.

    Args:
        value (float): The value.
        gate (Gate): The gate of the metric.

    Returns:
        str: The value with its unit.
    """
    unit = gate.unit
    decimals = DECIMALS.get(unit, 4)
    text = f"{value:.{decimals}f}"
    while decimals < 8:
        if gate.is_met(float(text)) == gate.is_met(value):
            break
        decimals += 1
        text = f"{value:.{decimals}f}"
    return text + ("%" if unit == "%" else f" {unit}" if unit else "")


def format_limit(gate) -> str:
    """
    Formats the threshold of a gate.

    Args:
        gate (Gate): The gate.

    Returns:
        str: The threshold.
    """
    sign = "≤" if gate.at_most else "≥"
    limit = f"{gate.limit:g}"
    unit = "%" if gate.unit == "%" else f" {gate.unit}" if gate.unit else ""
    return f"{sign} {limit}{unit}"


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
    names = {platform.id: platform.name for platform in entry.platforms}
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
    if configuration.attitude_schedule is not None:
        parts.append(
            "Its attitude is scheduled by day: build it with the start and end of "
            "the analysis (``start`` and ``end``).\n"
        )
    if record is not None:
        gap = (
            "no holdout"
            if record.holdout_gap_days is None
            else f"{record.holdout_gap_days:g} days from the fitting data"
        )
        if not configuration.is_fitted:
            gap = "no fitted values"
        source = f"`{record.source} <{REPOSITORY}{record.source}>`_"
        parts.append(
            f"**Reference:** {escape(record.reference)}, {record.period} "
            f"({escape(record.platform)}; {gap}). **Results:** {source}, "
            f"TAT-C {record.tatc_version}.\n"
        )
    rows = []
    for result in validation.gates:
        gate = result.gate
        if result.value is None:
            value, outcome = "not tested", "untested"
        else:
            value = format_value(result.value, gate)
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
                METRIC_LABELS.get(gate.metric, gate.metric),
                value,
                format_limit(gate),
                outcome,
            ]
        )
    parts.append(
        "Gates, with the worst result over the scopes tested (such as platforms "
        "and time steps):\n"
    )
    parts.append(
        list_table(["Metric", "Result", "Limit", "Outcome"], rows, [34, 36, 15, 15])
    )
    if record is not None:
        path = f"src/tatc/resources/validation/{entry.id}.json"
        parts.append(
            f"Every result, including those not gated, is in the record "
            f"`{entry.id}.json <{REPOSITORY}{path}>`_.\n"
        )
    gated = {result.gate.metric for result in validation.gates}
    gaps = [
        text
        for text in validation.gaps
        if not any(text.startswith(f"{metric}:") for metric in gated)
    ]
    if gaps:
        parts.append("Other gaps:\n")
        parts.append("\n".join(f"* {escape(text)}" for text in gaps) + "\n")
    if configuration.limitations:
        parts.append("Limitations:\n")
        parts.append(
            "\n".join(f"* {escape(text)}" for text in configuration.limitations) + "\n"
        )
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
                for p in entry.platforms
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
    if validation.gaps:
        parts.append(
            "Gaps of the entry:\n\n"
            + "\n".join(f"* {escape(text)}" for text in validation.gaps)
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
    Writes the catalog page.

    Returns:
        str: The page, as reStructuredText.
    """
    entries = sorted(
        (get_entry(entry_id) for entry_id in get_entry_ids()),
        key=lambda e: (ORDER[e.validation.tier], e.name.lower()),
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
                ", ".join(p.name for p in e.platforms),
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
    entry = get_entry("viirs")
    entry.validation  # tier, gate results, and gaps of each configuration
    satellite = entry.to_satellite(orbit, time_step=10, platform="noaa-20")

Pass the analysis time step to configurations whose along-track field of view
is tailored to it, and keep the attitude of the configuration: a scanner modeled
in another velocity frame can misplace its swath edges by tens of kilometers.
Pass the start and end of the analysis to configurations whose attitude is
scheduled by day, such as OCI's tilt. Where a platform has configurations of its
own (NOAA-21's ATMS and VIIRS), use them on that platform. The entries hold no
element sets, as they age; supply a recent one for the analysis period.

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

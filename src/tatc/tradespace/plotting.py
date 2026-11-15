"""
Plots of the results of trade studies (with matplotlib, an optional
dependency, imported when a plot is drawn): designs by two metrics, with the
efficient designs ringed and labeled; the progress of searches; the means
of metrics over replicates; the probabilities that the metrics of designs
exceed those of others; and the shares of the variance of metrics that
factors explain. Each plot draws on given axes (or new ones, with a
constrained layout) and returns them, for further labels and annotations.
`use_style` gives the figures of a session (such as of a notebook) the
style of these plots: titles and legends of their sizes, frameless legends,
light grids behind the data, no top or right spines, and constrained
layouts.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd

from .tables import feasible_records

COLORS = (
    "#2a78d6",
    "#eb6834",
    "#1baf7a",
    "#eda100",
    "#8f5bd6",
    "#d6457a",
    "#5aa8b0",
    "#8a6b4f",
)
"""The colors of groups, in order."""

INK = "#2b2b29"
"""The color of text and of the rings of efficient designs."""

GRAY = "#8a8a86"
"""The color of references (such as lines of equal values)."""

GRID = "#e6e6e3"
"""The color of grids."""

NEUTRAL = "#f0efec"
"""The color of the midpoint of diverging scales (such as of even odds)."""

STYLE = {
    "axes.titlesize": 10,
    "axes.labelcolor": INK,
    "axes.axisbelow": True,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.grid.axis": "y",
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "legend.fontsize": 8,
    "legend.frameon": False,
    "figure.constrained_layout.use": True,
}
"""The matplotlib settings of `use_style`."""


def use_style() -> None:
    """
    Gives the figures of the session the style of the plots of trade studies
    (see `STYLE`): titles of 10 points and legends of 8 without frames, light
    grids (of the vertical axis) behind the data, no top or right spines, and
    constrained layouts (which need no `tight_layout`).
    """
    try:
        # pylint: disable-next=import-outside-toplevel
        import matplotlib as mpl
    except ImportError as error:
        raise ImportError(
            "plots require matplotlib (pip install matplotlib)"
        ) from error
    mpl.rcParams.update(STYLE)


def _axes(ax: Any, size: tuple[float, float]) -> Any:
    """Gets the axes to draw on: those given, or those of a new figure."""
    if ax is not None:
        return ax
    try:
        # pylint: disable-next=import-outside-toplevel
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise ImportError(
            "plots require matplotlib (pip install matplotlib)"
        ) from error
    _, ax = plt.subplots(figsize=size, layout="constrained")
    return ax


def style_axes(ax: Any, grid: str | None = "y") -> Any:
    """
    Styles axes plainly: a light grid behind the data (of the given axes
    only), and no top or right spines.

    Args:
        ax (matplotlib.axes.Axes): The axes.
        grid (str | None): The axes of the grid (`"x"`, `"y"`, or `"both"`),
            or `None` for none.

    Returns:
        matplotlib.axes.Axes: the axes
    """
    ax.grid(False)
    if grid is not None:
        ax.grid(True, axis=grid, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    return ax


_JITTER_SEED = 0
"""The seed of the random shifts of `plot_designs` (see its `jitter`)."""

_OTHERS = object()
"""The group of the rows without a value of the column that groups them."""


def _color(colors: Sequence[str] | Mapping[Any, str], value: Any, k: int) -> str:
    if isinstance(colors, Mapping):
        if value is None and None not in colors:
            return next(iter(colors.values()))  # of rows without groups
        return colors[value]
    return colors[k % len(colors)]


def _mask(frame: pd.DataFrame, rows: str | Sequence[bool], name: str) -> np.ndarray:
    """Gets the rows of a frame that a boolean column (or values) selects."""
    if isinstance(rows, str):
        return frame[rows].to_numpy(dtype=bool)
    mask = np.asarray(rows, dtype=bool)
    if mask.shape != (len(frame),):
        raise ValueError(f"{name} needs a value for each row ({len(frame)})")
    return mask


def _name(names: Mapping[Any, str] | Callable[[Any], str] | None, value: Any) -> str:
    if names is None:
        return str(value)
    if callable(names):
        return names(value)
    return names.get(value, str(value))


def _front(frame: pd.DataFrame, front: str | Sequence[bool] | None) -> Any:
    """
    Gets the efficient rows of a frame to ring (a boolean column, or values),
    or `None` for none: without the column `"pareto"` (of a summary of
    criteria), none.
    """
    if front is None:
        return None
    if isinstance(front, str) and front not in frame.columns:
        if front == "pareto":
            return None
        raise ValueError(f"the frame has no column {front!r} of efficient designs")
    return _mask(frame, front, "front")


def _distinct(values: pd.Series) -> list:
    """Gets the distinct values (of a column), sorted if they can be."""
    distinct = list(pd.unique(values.dropna()))
    try:
        return sorted(distinct)
    except TypeError:
        return distinct


def _values_of(frame: pd.DataFrame, by: str | Sequence, name: str) -> pd.Series:
    """Gets the values of a column of a frame, or values given for its rows."""
    if isinstance(by, str):
        return frame[by]
    values = pd.Series(np.asarray(by), index=frame.index)
    if len(values) != len(frame):
        raise ValueError(f"{name} needs a value for each row ({len(frame)})")
    return values


def plot_designs(
    frame: pd.DataFrame,
    x: str,
    y: str,
    *,
    by: str | Sequence | None = None,
    front: str | Sequence[bool] | None = "pareto",
    annotate: str | Callable[[pd.Series], str] | None = None,
    names: Mapping[Any, str] | Callable[[Any], str] | None = None,
    highlight: str | Sequence[bool] | None = None,
    label: str | None = None,
    others: str | None = None,
    front_label: str = "Pareto-efficient",
    jitter: float = 0.0,
    colors: Sequence[str] | Mapping[Any, str] = COLORS,
    legend: str | None = "best",
    ax: Any = None,
) -> Any:
    """
    Plots designs (such as a summary of a study, see `summarize`) by two
    metrics: a point per row, colored by the value of a column `by` (if
    any), with the efficient designs (`front`, by default those of the
    column `pareto` of a summary of criteria, if it has one) ringed. Rows
    without a value of `by` (`NaN`, such as of a variable of a branch that a
    design does not choose) are drawn in gray, named `others` (or "other")
    in the legend. With `highlight`, only some designs are colored, over the
    others in gray (such as the designs a search evaluated, over every
    design)::

        plot_designs(summary, "satellites", "sufficient", by="layout")
        plot_designs(
            summary, "cost", "accurate", highlight=summary.policy == "degrade",
            label="degrade", others="other policies",
        )

    Args:
        frame (pandas.DataFrame): A row per design.
        x (str): The column of the horizontal axis.
        y (str): The column of the vertical axis.
        by (str | Sequence | None): A column whose values color the designs,
            values for the rows (such as of a column grouped, as
            `summary.planes.clip(upper=3)`), or `None`.
        front (str | Sequence[bool] | None): The efficient designs: a boolean
            column (by default `"pareto"`, of a summary of criteria, which a
            frame without the column does not ring) or values (one per row,
            such as of `pareto_front`), or `None`.
        annotate (str | Callable | None): The labels of the efficient
            designs, written next to them: a column, a function of a row, or
            `None`.
        names (Mapping | Callable | None): The names of the values of `by` in
            the legend (by default, the values).
        highlight (str | Sequence[bool] | None): The designs to color (a
            boolean column, or values), the others drawn in gray, or `None`
            to color every design.
        label (str | None): The name of the colored designs in the legend,
            without `by`.
        others (str | None): The name of the designs in gray in the legend:
            those not highlighted, or without a value of `by` ("other" if
            not given).
        front_label (str): The name of the efficient designs in the legend.
        jitter (float): The most that `x` is shifted at random (in its units),
            to separate designs of the same `x` (such as numbers of
            satellites).
        colors (Sequence[str] | Mapping): The colors of the values of `by`,
            in order (or by value).
        legend (str | None): The location of the legend, or `None` for none.
        ax (matplotlib.axes.Axes | None): The axes, or `None` for a new figure.

    Returns:
        matplotlib.axes.Axes: the axes

    Raises:
        ValueError: if `front` is a column (other than `"pareto"`) that the
            frame lacks, or values of a length other than the frame's.
    """
    ax = _axes(ax, (9, 4.5))
    shift = np.random.default_rng(_JITTER_SEED).uniform(-jitter, jitter, len(frame))
    xs = frame[x].to_numpy(dtype=float) + (shift if jitter else 0.0)
    ys = frame[y].to_numpy(dtype=float)
    shown = np.ones(len(frame), dtype=bool)
    if highlight is not None:
        shown = _mask(frame, highlight, "highlight")
        ax.scatter(
            xs[~shown], ys[~shown], s=10, color=GRAY, alpha=0.3, lw=0, label=others
        )
    if by is None:
        ax.scatter(
            xs[shown],
            ys[shown],
            s=14,
            color=_color(colors, None, 0),
            alpha=0.55 if highlight is None else 0.8,
            lw=0,
            label=label,
        )
    else:
        groups = _values_of(frame, by, "by")
        for k, value in enumerate(_distinct(groups)):
            selected = (groups == value).to_numpy() & shown
            ax.scatter(
                xs[selected],
                ys[selected],
                s=14,
                color=_color(colors, value, k),
                alpha=0.55,
                lw=0,
                label=_name(names, value),
            )
        missing = groups.isna().to_numpy() & shown
        if missing.any():
            # the designs not highlighted, if any, are named so already
            named = highlight is not None and others is not None
            ax.scatter(
                xs[missing],
                ys[missing],
                s=14,
                color=GRAY,
                alpha=0.55,
                lw=0,
                label=None if named else (others or "other"),
            )
    ringed = _front(frame, front)
    if ringed is not None:
        ax.scatter(
            xs[ringed],
            ys[ringed],
            s=70,
            facecolor="none",
            edgecolor=INK,
            lw=1.5,
            label=front_label,
        )
        if annotate is not None:
            rows = frame[ringed]
            for (_, row), point in zip(rows.iterrows(), zip(xs[ringed], ys[ringed])):
                ax.annotate(
                    annotate(row) if callable(annotate) else str(row[annotate]),
                    point,
                    xytext=(8, -3),
                    textcoords="offset points",
                    fontsize=8,
                    color=INK,
                )
    ax.set(xlabel=x, ylabel=y)
    style_axes(ax)
    named = by is not None or ringed is not None or label is not None
    if legend is not None and (named or others is not None):
        ax.legend(loc=legend, fontsize=8, frameon=False, markerscale=1.6)
    return ax


def plot_progress(
    progress: Mapping[str, pd.Series],
    *,
    colors: Sequence[str] | Mapping[str, str] = COLORS,
    legend: str | None = "lower right",
    ax: Any = None,
) -> Any:
    """
    Plots the progress of searches: a line of each (such as the hypervolume
    of the designs evaluated, see `hypervolume_progress`) by the number of
    designs evaluated, labeled at its end.

    Args:
        progress (Mapping[str, pandas.Series]): The progress of each search,
            by its name: a value by the number of designs evaluated.
        colors (Sequence[str] | Mapping): The colors of the searches, in
            order (or by name).
        legend (str | None): The location of the legend, or `None` for none.
        ax (matplotlib.axes.Axes | None): The axes, or `None` for a new figure.

    Returns:
        matplotlib.axes.Axes: the axes
    """
    ax = _axes(ax, (9, 4.2))
    for k, (name, series) in enumerate(progress.items()):
        ax.plot(
            series.index, series.values, color=_color(colors, name, k), lw=2, label=name
        )
        if len(series) > 0:
            ax.annotate(
                name,
                (series.index[-1], series.values[-1]),
                xytext=(6, 0),
                textcoords="offset points",
                va="center",
                fontsize=8,
                color=INK,
            )
    ax.set(xlabel="Designs evaluated", xlim=(0, None))
    style_axes(ax)
    if legend is not None and progress:
        ax.legend(loc=legend, fontsize=8, frameon=False)
    return ax


def plot_means(
    results: pd.DataFrame,
    x: str,
    y: str,
    *,
    by: str | None = None,
    names: Mapping[Any, str] | Callable[[Any], str] | None = None,
    others: str | None = None,
    colors: Sequence[str] | Mapping[Any, str] = COLORS,
    legend: str | None = "best",
    ax: Any = None,
) -> Any:
    """
    Plots the mean of a metric over the feasible records of each value of a
    column (such as over the replicates of the designs of each number of
    satellites), with a band of one standard deviation about it, a line for
    each value of a column `by` (if any), and a line in gray of the records
    without a value of `by` (`NaN`), named `others` (or "other").

    Args:
        results (pandas.DataFrame): The results of a study (see `Study.run`).
        x (str): The column of the horizontal axis.
        y (str): The metric.
        by (str | None): A column whose values have their own lines, or
            `None`.
        names (Mapping | Callable | None): The names of the values of `by` in
            the legend (by default, the values).
        others (str | None): The name of the records without a value of `by`
            in the legend ("other" if not given).
        colors (Sequence[str] | Mapping): The colors of the values of `by`,
            in order (or by value).
        legend (str | None): The location of the legend, or `None` for none.
        ax (matplotlib.axes.Axes | None): The axes, or `None` for a new figure.

    Returns:
        matplotlib.axes.Axes: the axes
    """
    ax = _axes(ax, (6.5, 4.2))
    results = feasible_records(results)[0]
    groups = (
        [(None, results)]
        if by is None
        else [
            (value, results[results[by] == value]) for value in _distinct(results[by])
        ]
    )
    missing = None if by is None else results[results[by].isna()]
    if missing is not None and len(missing) > 0:
        groups.append((_OTHERS, missing))
    for k, (value, group) in enumerate(groups):
        stats = group.groupby(x)[y].agg(["mean", "std"]).sort_index()
        spread = stats["std"].fillna(0.0)
        color = GRAY if value is _OTHERS else _color(colors, value, k)
        ax.fill_between(
            stats.index,
            stats["mean"] - spread,
            stats["mean"] + spread,
            color=color,
            alpha=0.12,
            lw=0,
        )
        ax.plot(
            stats.index,
            stats["mean"],
            marker="o",
            ms=5,
            lw=2,
            color=color,
            mec="white",
            mew=1.2,
            label=(
                None
                if value is None
                else (others or "other") if value is _OTHERS else _name(names, value)
            ),
        )
    ax.set(xlabel=x, ylabel=y)
    style_axes(ax)
    if legend is not None and by is not None:
        ax.legend(loc=legend, fontsize=8, frameon=False)
    return ax


def plot_exceedance(
    frame: pd.DataFrame,
    *,
    names: Mapping[str, str] | Callable[[str], str] | None = None,
    colorbar: bool = True,
    ax: Any = None,
) -> Any:
    """
    Plots the probabilities that the metric of each design exceeds that of
    each other design (see `exceedance`) as a matrix: a row per design and a
    column per other, in the order of the frame, colored from blue (0: the
    design's metric is below the other's in every replicate) through gray
    (0.5: as likely either way) to orange (1: above it in every replicate).
    Pairs whose order the replicates settle are blue or orange, and those
    whose order they leave open are gray: with designs in the order of the
    metric's mean, a band of gray along the diagonal is the designs that
    the replicates do not tell apart.

    Args:
        frame (pandas.DataFrame): The probabilities of a group (see
            `exceedance`).
        names (Mapping | Callable | None): The names of the designs (by
            key), or `None` for their keys.
        colorbar (bool): `True` to draw a color bar.
        ax (matplotlib.axes.Axes | None): The axes, or `None` for a new figure.

    Returns:
        matplotlib.axes.Axes: the axes

    Raises:
        ValueError: if the frame has several rows of a pair of designs (as of
            several groups, of which plot one).
    """
    if frame.duplicated(["design", "other"]).any():
        raise ValueError(
            "the frame has several rows of a pair of designs: plot one group"
        )
    ax = _axes(ax, (6.5, 5.5))
    # pylint: disable-next=import-outside-toplevel
    from matplotlib.colors import LinearSegmentedColormap

    rows = list(pd.unique(frame["design"]))
    columns = list(pd.unique(frame["other"]))
    matrix = frame.pivot(index="design", columns="other", values="exceedance")
    matrix = matrix.reindex(index=rows, columns=columns).to_numpy(dtype=float)
    scale = LinearSegmentedColormap.from_list(
        "exceedance", [COLORS[0], NEUTRAL, COLORS[1]]
    )
    image = ax.imshow(
        matrix, cmap=scale, vmin=0, vmax=1, interpolation="nearest", aspect="equal"
    )
    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels([_name(names, key) for key in columns], rotation=90)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([_name(names, key) for key in rows])
    ax.tick_params(length=0, labelsize=8, colors=INK)
    ax.set(xlabel="Other design", ylabel="Design")
    for spine in ax.spines.values():
        spine.set_visible(False)
    if colorbar:
        legend = ax.figure.colorbar(image, ax=ax, fraction=0.046, pad=0.03)
        legend.set_label("Probability that the design exceeds the other", fontsize=9)
        legend.set_ticks([0, 0.25, 0.5, 0.75, 1])
        legend.ax.tick_params(labelsize=8, length=0)
        legend.outline.set_visible(False)
    return ax


def plot_shares(
    table: pd.DataFrame,
    metrics: Sequence[str] | None = None,
    *,
    groups: str = "factors",
    names: Mapping[Any, str] | Callable[[Any], str] | None = None,
    labels: Mapping[Any, str] | Callable[[Any], str] | None = None,
    horizontal: bool = False,
    colors: Sequence[str] | Mapping[Any, str] = COLORS,
    legend: str | None = "best",
    ax: Any = None,
) -> Any:
    """
    Plots the shares of the variance of metrics that factors explain (such
    as those of `sensitivity`, by design variable, or of
    `Attribution.summary`, by uncertain parameter) as grouped bars: a group
    of bars for each factor (a row of the table), one bar for each metric
    (a column), or a group for each metric (`groups="metrics"`), one bar for
    each factor, with error bars of the standard errors of the shares
    (`<metric> se`) where the table gives them.

    Args:
        table (pandas.DataFrame): A row per factor and a column per metric
            (and optionally their standard errors, `<metric> se`).
        metrics (Sequence[str] | None): The metrics, or `None` for every
            numeric column but standard errors and `active`.
        groups (str): `"factors"` for a group of bars for each factor, or
            `"metrics"` for one for each metric.
        names (Mapping | Callable | None): The names of the factors (by the
            index of the table), by default their values.
        labels (Mapping | Callable | None): The names of the metrics, by
            default the metrics.
        horizontal (bool): `True` for horizontal bars (the groups down the
            vertical axis, the first at the top).
        colors (Sequence[str] | Mapping): The colors of the bars of each
            group, in order (or by metric or factor).
        legend (str | None): The location of the legend, or `None` for none.
        ax (matplotlib.axes.Axes | None): The axes, or `None` for a new figure.

    Returns:
        matplotlib.axes.Axes: the axes
    """
    # pylint: disable=too-many-locals
    if groups not in ("factors", "metrics"):
        raise ValueError(f"groups must be 'factors' or 'metrics', not {groups!r}")
    if metrics is None:
        metrics = [
            column
            for column in table.columns
            if not str(column).endswith(" se")
            and column != "active"
            and pd.api.types.is_numeric_dtype(table[column])
        ]
    metrics = list(metrics)
    factors = list(table.index)
    values = table[metrics].to_numpy(dtype=float)
    errors = np.column_stack(
        [
            (
                table[f"{metric} se"].to_numpy(dtype=float)
                if f"{metric} se" in table.columns
                else np.full(len(table), np.nan)
            )
            for metric in metrics
        ]
    )
    if groups == "metrics":
        values, errors = values.T, errors.T
        group_names = [_name(labels, metric) for metric in metrics]
        bar_keys, bar_names = factors, [_name(names, factor) for factor in factors]
    else:
        group_names = [_name(names, factor) for factor in factors]
        bar_keys, bar_names = metrics, [_name(labels, metric) for metric in metrics]
    ax = _axes(ax, (7.5, 4) if horizontal else (8, 4))
    positions = np.arange(len(group_names))
    width = 0.8 / max(len(bar_keys), 1)
    for k, (key, name) in enumerate(zip(bar_keys, bar_names)):
        offset = positions + (k - (len(bar_keys) - 1) / 2) * width
        error = None if np.all(np.isnan(errors[:, k])) else errors[:, k]
        options = {
            "color": _color(colors, key, k),
            "edgecolor": "white",
            "lw": 1,
            "label": name,
            "error_kw": {"ecolor": INK, "lw": 1},
        }
        if horizontal:
            ax.barh(offset, values[:, k], height=width, xerr=error, **options)
        else:
            ax.bar(offset, values[:, k], width=width, yerr=error, **options)
    if horizontal:
        ax.set_yticks(positions, group_names)
        ax.invert_yaxis()
        style_axes(ax, grid="x")
    else:
        ax.set_xticks(positions, group_names)
        style_axes(ax, grid="y")
    if legend is not None and len(bar_keys) > 1:
        ax.legend(loc=legend, fontsize=8, frameon=False)
    return ax

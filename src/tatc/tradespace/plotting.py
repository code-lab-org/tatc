"""
Plots of the results of trade studies (with matplotlib, an optional
dependency, imported when a plot is drawn): designs by two metrics, with the
efficient designs ringed and labeled; the progress of searches; the means
of metrics over replicates; and the probabilities that the metrics of
designs exceed those of others. Each plot draws on given axes (or new ones)
and returns them, for further labels and annotations.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd

from .criteria import Criteria
from .pareto import pareto_front

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
    _, ax = plt.subplots(figsize=size)
    return ax


def style_axes(ax: Any, grid: str | None = "y") -> Any:
    """
    Styles axes plainly: a light grid behind the data, and no top or right
    spines.

    Args:
        ax (matplotlib.axes.Axes): The axes.
        grid (str | None): The axes of the grid (`"x"`, `"y"`, or `"both"`),
            or `None` for none.

    Returns:
        matplotlib.axes.Axes: the axes
    """
    if grid is not None:
        ax.grid(axis=grid, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    return ax


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


def _groups(frame: pd.DataFrame, by: str) -> list:
    """Gets the values of a column, sorted if they can be."""
    values = list(pd.unique(frame[by].dropna()))
    try:
        return sorted(values)
    except TypeError:
        return values


def plot_designs(
    frame: pd.DataFrame,
    x: str,
    y: str,
    *,
    by: str | None = None,
    front: str | Sequence[bool] | Criteria | Mapping[str, str] | None = None,
    labels: str | Callable[[pd.Series], str] | None = None,
    names: Mapping[Any, str] | Callable[[Any], str] | None = None,
    highlight: str | Sequence[bool] | None = None,
    label: str | None = None,
    others: str | None = None,
    front_label: str = "Pareto-efficient",
    jitter: float = 0.0,
    colors: Sequence[str] | Mapping[Any, str] = COLORS,
    legend: str | None = "best",
    ax: Any = None,
    seed: int = 0,
) -> Any:
    """
    Plots designs (such as a summary of a study, see `aggregate`) by two
    metrics: a point per row, colored by the value of a column `by` (if
    any), with the efficient designs (`front`) ringed and labeled. With
    `highlight`, only some designs are colored, over the others in gray
    (such as the designs a search evaluated, over every design)::

        plot_designs(summary, "satellites", "sufficient", front="pareto")
        plot_designs(
            summary, "cost", "accurate", highlight=summary.policy == "degrade",
            label="degrade", others="other policies", front="pareto",
        )

    Args:
        frame (pandas.DataFrame): A row per design.
        x (str): The column of the horizontal axis.
        y (str): The column of the vertical axis.
        by (str | None): A column whose values color the designs, or `None`.
        front (str | Sequence[bool] | Criteria | Mapping | None): The
            efficient designs: a boolean column (such as `"pareto"` of a
            summary of criteria) or values (one per row), criteria (see
            `pareto_front`), or `None`.
        labels (str | Callable | None): The labels of the efficient designs:
            a column, a function of a row, or `None`.
        names (Mapping | Callable | None): The names of the values of `by` in
            the legend (by default, the values).
        highlight (str | Sequence[bool] | None): The designs to color (a
            boolean column, or values), the others drawn in gray, or `None`
            to color every design.
        label (str | None): The name of the colored designs in the legend,
            without `by`.
        others (str | None): The name of the designs in gray in the legend.
        front_label (str): The name of the efficient designs in the legend.
        jitter (float): The most that `x` is shifted at random (in its units),
            to separate designs of the same `x` (such as numbers of
            satellites).
        colors (Sequence[str] | Mapping): The colors of the values of `by`,
            in order (or by value).
        legend (str | None): The location of the legend, or `None` for none.
        ax (matplotlib.axes.Axes | None): The axes, or `None` for a new figure.
        seed (int): The seed of the jitter.

    Returns:
        matplotlib.axes.Axes: the axes
    """
    ax = _axes(ax, (9, 4.5))
    shift = np.random.default_rng(seed).uniform(-jitter, jitter, len(frame))
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
        for k, value in enumerate(_groups(frame, by)):
            selected = (frame[by] == value).to_numpy() & shown
            ax.scatter(
                xs[selected],
                ys[selected],
                s=14,
                color=_color(colors, value, k),
                alpha=0.55,
                lw=0,
                label=_name(names, value),
            )
    if front is not None:
        if isinstance(front, (Criteria, Mapping)):
            efficient = np.asarray(pareto_front(frame, front), dtype=bool)
        else:
            efficient = _mask(frame, front, "front")
        ax.scatter(
            xs[efficient],
            ys[efficient],
            s=70,
            facecolor="none",
            edgecolor=INK,
            lw=1.5,
            label=front_label,
        )
        if labels is not None:
            rows = frame[efficient]
            for (_, row), point in zip(
                rows.iterrows(), zip(xs[efficient], ys[efficient])
            ):
                ax.annotate(
                    labels(row) if callable(labels) else str(row[labels]),
                    point,
                    xytext=(8, -3),
                    textcoords="offset points",
                    fontsize=8,
                    color=INK,
                )
    ax.set(xlabel=x, ylabel=y)
    style_axes(ax)
    named = by is not None or front is not None or label is not None
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
    colors: Sequence[str] | Mapping[Any, str] = COLORS,
    legend: str | None = "best",
    ax: Any = None,
) -> Any:
    """
    Plots the mean of a metric over the feasible records of each value of a
    column (such as over the replicates of the designs of each number of
    satellites), with a band of one standard deviation about it, a line for
    each value of a column `by` (if any).

    Args:
        results (pandas.DataFrame): The results of a study (see `Study.run`).
        x (str): The column of the horizontal axis.
        y (str): The metric.
        by (str | None): A column whose values have their own lines, or
            `None`.
        names (Mapping | Callable | None): The names of the values of `by` in
            the legend (by default, the values).
        colors (Sequence[str] | Mapping): The colors of the values of `by`,
            in order (or by value).
        legend (str | None): The location of the legend, or `None` for none.
        ax (matplotlib.axes.Axes | None): The axes, or `None` for a new figure.

    Returns:
        matplotlib.axes.Axes: the axes
    """
    ax = _axes(ax, (6.5, 4.2))
    if "feasible" in results.columns:
        results = results[results["feasible"].astype(bool)]
    groups = (
        [(None, results)]
        if by is None
        else [(value, results[results[by] == value]) for value in _groups(results, by)]
    )
    for k, (value, group) in enumerate(groups):
        stats = group.groupby(x)[y].agg(["mean", "std"]).sort_index()
        spread = stats["std"].fillna(0.0)
        color = _color(colors, value, k)
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
            label=None if value is None else _name(names, value),
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

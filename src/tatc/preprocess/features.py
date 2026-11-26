"""
Methods to find features in the fields of a nature run: contiguous regions of
cells that meet a condition (such as storms, cold cloud tops, or cold pools)
at each time, as the features whose observations an observing system
simulation experiment (OSSE) collects (see
`tatc.analysis.collect_feature_observations`). Requires scipy and xarray
(`pip install tatc[osse]`).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from datetime import timedelta
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely

from ..utils._grids import FieldGrid, unit_vectors
from ..utils.geometry import split_polygon

_REDUCTIONS = ("min", "max", "mean", "integral")
"""The reductions of the values of a variable over a feature's cells."""


def _get_time_dimension(mask: Any, given: str | None) -> str:
    """Gets the time dimension of a mask."""
    name = given or "time"
    if name not in mask.dims:
        raise ValueError(f"the mask needs a time dimension (named {name!r})")
    return name


def _get_frame_duration(
    times: pd.DatetimeIndex, given: timedelta | None
) -> pd.Timedelta:
    """Gets the duration of each frame: given, or the median step of the times."""
    if given is not None:
        duration = pd.Timedelta(given)
    elif len(times) > 1:
        duration = pd.Timedelta(np.median(np.diff(times.as_unit("ns").asi8)), unit="ns")
    else:
        raise ValueError("a mask of one time needs a frame_duration")
    if duration <= pd.Timedelta(0):
        raise ValueError(f"frame_duration must be positive, not {duration}")
    return duration


def _to_degrees(vectors: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Gets the longitudes and latitudes (degrees) of Earth-centered vectors
    (with the components along the last axis)."""
    longitude = np.degrees(np.arctan2(vectors[..., 1], vectors[..., 0]))
    latitude = np.degrees(
        np.arctan2(vectors[..., 2], np.hypot(vectors[..., 0], vectors[..., 1]))
    )
    return longitude, latitude


def _edges(values: np.ndarray) -> np.ndarray:
    """Gets the edges of cells centered at sorted coordinate values: the
    midpoints between them, and half a spacing beyond the outer values."""
    if len(values) == 1:
        return np.array([values[0] - 0.5, values[0] + 0.5])
    middle = (values[1:] + values[:-1]) / 2
    return np.concatenate(
        [
            [values[0] - (values[1] - values[0]) / 2],
            middle,
            [values[-1] + (values[-1] - values[-2]) / 2],
        ]
    )


def _get_corners(grid: FieldGrid) -> tuple[np.ndarray, np.ndarray]:
    """
    Gets the longitudes and latitudes (degrees) of the corners of the cells
    of a grid (regular, with sorted coordinates, or curvilinear), with shape
    (rows + 1, columns + 1): of a regular grid, the edges between its
    latitudes and longitudes (latitudes within the poles); of a curvilinear
    grid, the means of the unit vectors of the four cells around each corner
    (extrapolated linearly beyond its outer cells).
    """
    if grid.kind == "regular":
        longitude, latitude = np.meshgrid(
            _edges(grid.longitude), np.clip(_edges(grid.latitude), -90, 90)
        )
        return longitude, latitude
    rows, columns = grid.shape
    vectors = unit_vectors(grid.longitude, grid.latitude).reshape(rows, columns, 3)
    padded = np.zeros((rows + 2, columns + 2, 3))
    padded[1:-1, 1:-1] = vectors
    if rows > 1:
        padded[0, 1:-1] = 2 * vectors[0] - vectors[1]
        padded[-1, 1:-1] = 2 * vectors[-1] - vectors[-2]
    else:
        padded[0, 1:-1] = padded[-1, 1:-1] = vectors[0]
    if columns > 1:
        padded[:, 0] = 2 * padded[:, 1] - padded[:, 2]
        padded[:, -1] = 2 * padded[:, -2] - padded[:, -3]
    else:
        padded[:, 0] = padded[:, -1] = padded[:, 1]
    corners = (
        padded[:-1, :-1] + padded[1:, :-1] + padded[:-1, 1:] + padded[1:, 1:]
    ) / 4
    corners /= np.linalg.norm(corners, axis=-1, keepdims=True)
    return _to_degrees(corners)


def _merge_around(labels: np.ndarray, count: int, connectivity: int) -> np.ndarray:
    """
    Merges the labels of cells that touch across the first and last columns
    of a grid that goes around the globe in longitude, relabeling them from 1.
    """
    # pylint: disable=import-outside-toplevel
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    pairs = [np.column_stack([labels[:, 0], labels[:, -1]])]
    if connectivity > 1:
        pairs.append(np.column_stack([labels[1:, 0], labels[:-1, -1]]))
        pairs.append(np.column_stack([labels[:-1, 0], labels[1:, -1]]))
    pairs = np.concatenate(pairs)
    pairs = pairs[(pairs[:, 0] > 0) & (pairs[:, 1] > 0)] - 1
    if len(pairs) == 0:
        return labels
    graph = coo_matrix(
        (np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(count, count)
    )
    _, component = connected_components(graph, directed=False)
    # renumber the merged labels by their first cell
    _, first = np.unique(component, return_index=True)
    order = np.empty(len(first), dtype=int)
    order[np.argsort(first)] = np.arange(len(first))
    mapping = np.concatenate([[0], order[component] + 1])
    return mapping[labels]


def _get_runs(
    labels: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Gets the runs of cells of the same label along each row of labels (with
    zero for no label): the row, first and last column, and label of each.
    """
    rows, columns = labels.shape
    # the boundaries before each column where the label changes, and at the
    # ends of each row
    change = np.ones((rows, columns + 1), dtype=bool)
    change[:, 1:-1] = labels[:, 1:] != labels[:, :-1]
    row, boundary = np.nonzero(change)
    same = row[1:] == row[:-1]
    row, start, end = row[:-1][same], boundary[:-1][same], boundary[1:][same] - 1
    label = labels[row, start]
    keep = label > 0
    return row[keep], start[keep], end[keep], label[keep]


def _get_polygons(
    labels: np.ndarray,
    keep: np.ndarray,
    corner_longitude: np.ndarray,
    corner_latitude: np.ndarray,
    reference: np.ndarray,
) -> list:
    """
    Gets the polygon of each kept label (by its index among the labels): the
    union of the strips of the runs of its cells along each row, each
    bounded by the corners of the cells, with longitudes unwrapped near the
    feature's reference longitude and split at the antimeridian (see
    `tatc.utils.split_polygon`).
    """
    row, start, end, label = _get_runs(labels)
    selected = keep[label - 1]
    row, start, end, label = (
        row[selected],
        start[selected],
        end[selected],
        label[selected],
    )
    polygons: list = [None] * len(keep)
    if len(row) == 0:
        return polygons
    # the corners of each strip: along the top edge, then back along the bottom
    length = end - start + 2
    total = 2 * length
    offset = np.repeat(np.cumsum(total) - total, total)
    position = np.arange(np.sum(total)) - offset
    run = np.repeat(np.arange(len(row)), total)
    top = position < np.repeat(length, total)
    column = np.where(
        top,
        np.repeat(start, total) + position,
        np.repeat(start, total) + 2 * np.repeat(length, total) - 1 - position,
    )
    corner_row = np.repeat(row, total) + np.where(top, 0, 1)
    longitude = corner_longitude[corner_row, column]
    latitude = corner_latitude[corner_row, column]
    center = reference[label[run] - 1]
    longitude = center + (longitude - center + 180) % 360 - 180
    strips = shapely.polygons(
        shapely.linearrings(np.column_stack([longitude, latitude]), indices=run)
    )
    order = np.argsort(label, kind="stable")
    bounds = np.flatnonzero(np.diff(label[order])) + 1
    for group in np.split(order, bounds):
        union = shapely.union_all(strips[group])
        polygons[label[group[0]] - 1] = split_polygon(shapely.make_valid(union))
    return polygons


def _reduce(
    labels: np.ndarray,
    count: int,
    areas: np.ndarray,
    variables: Mapping[str, np.ndarray],
    reductions: Sequence[str],
) -> dict[str, np.ndarray]:
    """
    Reduces the values of variables over the cells of each label (ignoring
    cells without values): the minimum, maximum, area-weighted mean, and
    integral (the sum of values times areas) of each.
    """
    flat = labels.ravel()
    labeled = flat > 0
    group = flat[labeled] - 1
    area = areas.ravel()[labeled]
    result = {}
    for name, values in variables.items():
        value = values.ravel()[labeled].astype(float)
        frame = pd.DataFrame({"group": group, "value": value})
        finite = np.isfinite(value)
        for reduction in reductions:
            column = f"{reduction}_{name}"
            if reduction in ("min", "max"):
                reduced = getattr(frame.groupby("group")["value"], reduction)()
            else:
                weighted = np.bincount(
                    group[finite], weights=value[finite] * area[finite], minlength=count
                )
                if reduction == "integral":
                    reduced = pd.Series(weighted)
                else:
                    covered = np.bincount(
                        group[finite], weights=area[finite], minlength=count
                    )
                    with np.errstate(invalid="ignore", divide="ignore"):
                        reduced = pd.Series(weighted / covered)
            result[column] = reduced.reindex(np.arange(count)).to_numpy(dtype=float)
    return result


def find_features(
    mask: Any,
    variables: Mapping[str, Any] | None = None,
    reductions: Sequence[str] = _REDUCTIONS,
    min_diameter: float | None = None,
    frame_duration: timedelta | None = None,
    alignment: str = "start",
    connectivity: int = 1,
    dimensions: Mapping[str, str] | None = None,
    return_labels: bool = False,
) -> gpd.GeoDataFrame | tuple[gpd.GeoDataFrame, Any]:
    """
    Finds the features of a mask at each time: the contiguous regions of its
    cells that are `True` (such as the cells of a nature run whose cloud-top
    temperature is below a threshold), each as the polygon of its cells, which
    exists for the frame of its time (from `start` to `end`), as
    `tatc.analysis.collect_feature_observations` takes them.

    Cells are contiguous if they share an edge (`connectivity=1`) or also a
    corner (`connectivity=2`), and across the antimeridian if the grid goes
    around the globe. The mask may be on a regular grid (latitude and
    longitude dimensions) or a curvilinear grid (two-dimensional latitude and
    longitude coordinates; see `tatc.analysis.sample_field`), whose cells'
    corners are midway between their centers.

    Gives each feature's `frame` (the index of its time), `label` (its label
    in the frame, as in the labels returned with `return_labels`), `time`,
    `start` and `end` (the frame of its time: from the time for a duration,
    with `alignment="start"`, or centered on it, `"center"`, as for means
    over intervals centered on their times), number of cells (`n_cells`),
    `area_m2` (square meters, of its cells, as `tatc.analysis.sample_field`
    weighs them), `diameter_m` (meters, of a circle of its area), and centroid (`lon` and
    `lat`, weighted by area), with the reductions of the values of
    `variables` over its cells (ignoring cells without values), each named
    `{reduction}_{variable}`: `min`, `max`, `mean` (weighted by area), and
    `integral` (the sum of values times areas, as of a precipitation rate
    over the feature). Features smaller than `min_diameter` are left out
    (and unlabeled).

    Args:
        mask (xarray.DataArray): The mask (`True` in the cells of features,
            `False` or missing elsewhere), with a time dimension.
        variables (Mapping[str, xarray.DataArray] | None): Variables over the
            same times and grid as the mask, by name, to reduce over each
            feature's cells.
        reductions (Sequence[str]): The reductions of the variables: any of
            `"min"`, `"max"`, `"mean"`, and `"integral"`.
        min_diameter (float | None): The least diameter (meters, of a circle
            of a feature's area) of features, or `None` for any.
        frame_duration (datetime.timedelta | None): The duration of each
            frame, or `None` for the median step between the mask's times.
        alignment (str): `"start"` for frames that start at their times, or
            `"center"` for frames centered on them.
        connectivity (int): 1 for cells that share an edge, 2 for cells that
            share an edge or a corner.
        dimensions (Mapping[str, str] | None): The names of the mask's `time`
            dimension and of its latitudes (`lat`) and longitudes (`lon`),
            if not the defaults (see `tatc.analysis.sample_field`).
        return_labels (bool): `True` to also return the label of each cell
            at each time (zero outside features).

    Returns:
        geopandas.GeoDataFrame | tuple[geopandas.GeoDataFrame, xarray.DataArray]:
            the features (and, with `return_labels`, the labels of the cells)
    """
    # pylint: disable=import-outside-toplevel,too-many-locals
    from scipy import ndimage

    if alignment not in ("start", "center"):
        raise ValueError(f"alignment must be 'start' or 'center', not {alignment!r}")
    if connectivity not in (1, 2):
        raise ValueError(f"connectivity must be 1 or 2, not {connectivity!r}")
    unknown = [reduction for reduction in reductions if reduction not in _REDUCTIONS]
    if unknown:
        raise ValueError(f"reductions must be among {list(_REDUCTIONS)}, not {unknown}")
    dimensions = dict(dimensions or {})
    time_dim = _get_time_dimension(mask, dimensions.get("time"))
    grid = FieldGrid.from_field(
        mask, dimensions.get("lat"), dimensions.get("lon"), time_dim
    )
    if grid.kind == "unstructured":
        raise ValueError(
            "features need a regular or curvilinear grid, not unstructured"
        )
    other = [dim for dim in mask.dims if dim != time_dim and dim not in grid.dims]
    if other:
        raise ValueError(
            f"the mask must have no dimensions besides time and its grid, not {other}"
        )
    # a regular grid in the order of its coordinates
    indexers: dict[str, np.ndarray] = {}
    around = False
    if grid.kind == "regular":
        indexers = {
            grid.dims[0]: np.argsort(grid.latitude, kind="stable"),
            grid.dims[1]: np.argsort(grid.longitude, kind="stable"),
        }
        grid = FieldGrid(
            "regular",
            grid.dims,
            grid.latitude[indexers[grid.dims[0]]],
            grid.longitude[indexers[grid.dims[1]]],
        )
        step = np.median(np.diff(grid.longitude)) if len(grid.longitude) > 1 else 0
        around = len(grid.longitude) * step >= 360 - 1e-6
    ordered = mask.isel(indexers).transpose(time_dim, *grid.dims)
    values = {
        name: np.asarray(
            variable.isel(indexers).transpose(time_dim, *grid.dims).to_numpy(),
            dtype=float,
        )
        for name, variable in (variables or {}).items()
    }
    for name, variable in values.items():
        if variable.shape != ordered.shape:
            raise ValueError(
                f"variable {name!r} must have the mask's shape {ordered.shape}, "
                + f"not {variable.shape}"
            )
    times = pd.DatetimeIndex(pd.to_datetime(ordered[time_dim].to_numpy(), utc=True))
    duration = _get_frame_duration(times, frame_duration)
    areas = grid.cell_areas(mask.isel(indexers)).reshape(grid.shape)
    corner_longitude, corner_latitude = _get_corners(grid)
    cell_latitude, cell_longitude = grid.cells()
    cell_vectors = unit_vectors(cell_longitude, cell_latitude)
    structure = ndimage.generate_binary_structure(2, connectivity)
    data = np.nan_to_num(np.asarray(ordered.to_numpy(), dtype=float)) != 0
    all_labels = np.zeros(data.shape, dtype=np.int32)
    frames = []
    for frame, time in enumerate(times):
        labels, count = ndimage.label(data[frame], structure)
        if count == 0:
            continue
        if around:
            labels = _merge_around(labels, count, connectivity)
            count = int(labels.max())
        flat = labels.ravel()
        cells = np.bincount(flat, minlength=count + 1)[1:]
        area = np.bincount(flat, weights=areas.ravel(), minlength=count + 1)[1:]
        diameter = np.sqrt(4 * area / np.pi)
        keep = np.ones(count, dtype=bool)
        if min_diameter is not None:
            keep = diameter >= min_diameter
        # the area-weighted centroid, from unit vectors
        weighted = np.column_stack(
            [
                np.bincount(
                    flat,
                    weights=cell_vectors[:, axis] * areas.ravel(),
                    minlength=count + 1,
                )[1:]
                for axis in range(3)
            ]
        )
        lon, lat = _to_degrees(weighted)
        polygons = _get_polygons(labels, keep, corner_longitude, corner_latitude, lon)
        reduced = _reduce(
            labels, count, areas, {k: v[frame] for k, v in values.items()}, reductions
        )
        kept = np.flatnonzero(keep)
        start = time if alignment == "start" else time - duration / 2
        table = pd.DataFrame(
            {
                "frame": frame,
                "label": kept + 1,
                "time": time,
                "start": start,
                "end": start + duration,
                "n_cells": cells[kept],
                "area_m2": area[kept],
                "diameter_m": diameter[kept],
                "lon": lon[kept],
                "lat": lat[kept],
                **{column: reduction[kept] for column, reduction in reduced.items()},
            }
        )
        frames.append(
            gpd.GeoDataFrame(
                table, geometry=[polygons[i] for i in kept], crs="EPSG:4326"
            )
        )
        mapping = np.concatenate([[0], np.where(keep, np.arange(1, count + 1), 0)])
        all_labels[frame] = mapping[labels]
    if frames:
        features = pd.concat(frames, ignore_index=True)
        features = gpd.GeoDataFrame(features, geometry="geometry", crs="EPSG:4326")
    else:
        features = _empty_features(values, reductions)
    if not return_labels:
        return features
    labels = ordered.copy(data=all_labels).rename("label")
    if indexers:
        # in the order of the mask's coordinates
        inverse = {
            dim: np.argsort(order, kind="stable") for dim, order in indexers.items()
        }
        labels = labels.isel(inverse).transpose(*mask.dims)
    return features, labels


def _empty_features(
    values: Mapping[str, Any], reductions: Sequence[str]
) -> gpd.GeoDataFrame:
    """Gets a data frame of no features, with the columns of features."""
    columns = {
        "frame": pd.Series(dtype=int),
        "label": pd.Series(dtype=int),
        "time": pd.Series(dtype="datetime64[ns, UTC]"),
        "start": pd.Series(dtype="datetime64[ns, UTC]"),
        "end": pd.Series(dtype="datetime64[ns, UTC]"),
        "n_cells": pd.Series(dtype=int),
        "area_m2": pd.Series(dtype=float),
        "diameter_m": pd.Series(dtype=float),
        "lon": pd.Series(dtype=float),
        "lat": pd.Series(dtype=float),
    }
    for name in values:
        for reduction in reductions:
            columns[f"{reduction}_{name}"] = pd.Series(dtype=float)
    return gpd.GeoDataFrame(columns, geometry=gpd.GeoSeries([], crs="EPSG:4326"))

"""
Private helpers for the horizontal grids of fields (such as the variables of
a nature run): regular grids of latitude and longitude dimensions,
curvilinear grids of two-dimensional latitude and longitude coordinates (as
of a regional model on a projected grid), and unstructured grids of cells
with latitude and longitude coordinates (as of a model on a mesh).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Literal

import numpy as np
import numpy.typing as npt

from .. import constants

LATITUDE_NAMES = ("lat", "latitude", "XLAT", "latCell", "nav_lat")
"""The usual names of the latitudes of a field's grid."""

LONGITUDE_NAMES = ("lon", "longitude", "XLONG", "lonCell", "nav_lon")
"""The usual names of the longitudes of a field's grid."""


_AREA_NAMES = ("cell_area", "area", "areaCell")
"""The names of coordinates of the areas of a field's cells."""

_REACH_TOLERANCE = 1.1
"""The factor of the reach of a cell (about half its diagonal) beyond which a
place is outside a grid that is not regular, as cells vary in size (as the
widths in longitude of a grid's rows narrow toward the poles)."""


def unit_vectors(
    longitude: npt.ArrayLike, latitude: npt.ArrayLike
) -> npt.NDArray[np.float64]:
    """
    Gets the Earth-centered unit vectors of places on a sphere.

    Args:
        longitude (numpy.typing.ArrayLike): The longitudes (degrees).
        latitude (numpy.typing.ArrayLike): The latitudes (degrees).

    Returns:
        numpy.typing.NDArray: the unit vectors, of shape (n, 3)
    """
    lon = np.radians(np.ravel(np.asarray(longitude, dtype=float)))
    lat = np.radians(np.ravel(np.asarray(latitude, dtype=float)))
    return np.column_stack(
        [np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)]
    )


def _find_name(
    names: Sequence[str], candidates: Sequence[str], given: str | None
) -> str | None:
    """Gets the given name, or the first candidate among names."""
    if given is not None:
        return given if given in names else None
    return next((name for name in candidates if name in names), None)


def _in_degrees(coordinate: Any) -> npt.NDArray[np.float64]:
    """Gets the values of a coordinate in degrees (from radians, if its units
    say so, as of the cells of an MPAS mesh)."""
    values = np.asarray(coordinate.to_numpy(), dtype=float)
    if "radian" in str(coordinate.attrs.get("units", "")).lower():
        return np.degrees(values)
    return values


class FieldGrid:
    """
    The horizontal grid of a field: `regular` (latitude and longitude
    dimensions), `curvilinear` (two-dimensional latitude and longitude
    coordinates over two dimensions), or `unstructured` (latitude and
    longitude coordinates over one dimension of cells).

    Args:
        kind (str): The kind of grid.
        dims (tuple[str, ...]): The horizontal dimensions: latitude and
            longitude (regular), the two dimensions of the coordinates
            (curvilinear), or that of the cells (unstructured).
        latitude (numpy.typing.NDArray): The latitudes (degrees) of the
            cells: of the latitude dimension (regular), or of each cell
            (others).
        longitude (numpy.typing.NDArray): The longitudes (degrees) of the
            cells, as the latitudes.
    """

    def __init__(
        self,
        kind: Literal["regular", "curvilinear", "unstructured"],
        dims: tuple[str, ...],
        latitude: npt.NDArray[np.float64],
        longitude: npt.NDArray[np.float64],
    ):
        self.kind = kind
        self.dims = dims
        self.latitude = latitude
        self.longitude = longitude
        self._tree: Any = None
        self._reach: npt.NDArray[np.float64] | None = None
        self._valid: npt.NDArray[np.intp] | None = None

    @classmethod
    def from_field(
        cls,
        field: Any,
        lat: str | None = None,
        lon: str | None = None,
        time: str | None = None,
    ) -> FieldGrid:
        """
        Gets the horizontal grid of a field: regular, if it has latitude and
        longitude dimensions, or else curvilinear or unstructured, by the
        dimensions of its latitude and longitude coordinates (of the first
        time, if they vary in time, as those of a WRF run).

        Args:
            field (xarray.DataArray): The field.
            lat (str | None): The name of the latitude dimension or
                coordinate, or `None` for the first of the usual names (see
                `LATITUDE_NAMES`).
            lon (str | None): The name of the longitude dimension or
                coordinate, or `None` for the first of the usual names (see
                `LONGITUDE_NAMES`).
            time (str | None): The name of the time dimension, if any.

        Returns:
            FieldGrid: the grid
        """
        dims = list(field.dims)
        lat_dim = _find_name(dims, LATITUDE_NAMES, lat)
        lon_dim = _find_name(dims, LONGITUDE_NAMES, lon)
        if lat_dim is not None and lon_dim is not None:
            return cls(
                "regular",
                (lat_dim, lon_dim),
                _in_degrees(field[lat_dim]),
                _in_degrees(field[lon_dim]),
            )
        coords = list(field.coords)
        lat_name = _find_name(coords, LATITUDE_NAMES, lat)
        lon_name = _find_name(coords, LONGITUDE_NAMES, lon)
        if lat_name is None or lon_name is None:
            raise ValueError(
                "the field has no latitude and longitude dimensions or "
                + f"coordinates (named any of {list(LATITUDE_NAMES)} and "
                + f"{list(LONGITUDE_NAMES)}, or as given)"
            )
        latitude, longitude = field.coords[lat_name], field.coords[lon_name]
        if time is not None and time in latitude.dims:
            latitude = latitude.isel({time: 0})
        if time is not None and time in longitude.dims:
            longitude = longitude.isel({time: 0})
        if latitude.dims != longitude.dims:
            raise ValueError(
                f"the field's latitudes {latitude.dims} and longitudes "
                + f"{longitude.dims} must have the same dimensions"
            )
        if latitude.ndim not in (1, 2):
            raise ValueError(
                "the field's latitude and longitude coordinates must have one "
                + f"(unstructured) or two (curvilinear) dimensions, not {latitude.dims}"
            )
        return cls(
            "curvilinear" if latitude.ndim == 2 else "unstructured",
            tuple(str(dim) for dim in latitude.dims),
            _in_degrees(latitude),
            _in_degrees(longitude),
        )

    @property
    def shape(self) -> tuple[int, ...]:
        """tuple[int, ...]: The shape of the horizontal dimensions."""
        if self.kind == "regular":
            return (len(self.latitude), len(self.longitude))
        return self.latitude.shape

    def _cells(self) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        """Gets the latitudes and longitudes of every cell, flattened in the
        order of the horizontal dimensions."""
        if self.kind == "regular":
            lon, lat = np.meshgrid(self.longitude, self.latitude)
            return lat.ravel(), lon.ravel()
        return self.latitude.ravel(), self.longitude.ravel()

    def _build(self) -> None:
        """
        Builds the tree of the cells with finite coordinates and the reach of
        each: about half its diagonal, from the distances to its nearest
        distinct neighbor and to its third (in the other direction of a grid
        of rows and columns, or the next of a mesh), as half their hypotenuse.
        """
        # pylint: disable=import-outside-toplevel
        from scipy.spatial import cKDTree

        latitude, longitude = self._cells()
        valid = np.flatnonzero(np.isfinite(latitude) & np.isfinite(longitude))
        if len(valid) == 0:
            raise ValueError("the field's grid has no cells with coordinates")
        vectors = unit_vectors(longitude[valid], latitude[valid])
        tree = cKDTree(vectors)
        reach = np.full(len(valid), np.inf)
        if len(valid) > 1:
            # distinct neighbors (cells may repeat, as at a pole)
            distance, _ = tree.query(vectors, k=min(len(valid), 7))
            distance = np.sort(np.where(distance > 0, distance, np.inf), axis=1)
            first = distance[:, 0]
            third = distance[:, min(2, distance.shape[1] - 1)]
            third = np.where(np.isfinite(third), third, first)
            reach = np.hypot(first, third) / 2
        self._tree, self._reach, self._valid = tree, reach, valid

    @property
    def tree(self) -> Any:
        """scipy.spatial.cKDTree: The tree of the unit vectors of the cells
        with coordinates (see `valid`), built when first used."""
        if self._tree is None:
            self._build()
        return self._tree

    @property
    def reach(self) -> npt.NDArray[np.float64]:
        """numpy.typing.NDArray: The chord distance (on the unit sphere) from
        each cell with coordinates to its corners, about half its diagonal
        (see `_build`)."""
        if self._reach is None:
            self._build()
        return self._reach  # type: ignore[return-value]

    @property
    def valid(self) -> npt.NDArray[np.intp]:
        """numpy.typing.NDArray: The flat indices of the cells with
        coordinates, in the order of the tree."""
        if self._valid is None:
            self._build()
        return self._valid  # type: ignore[return-value]

    def nearest(
        self, longitude: npt.ArrayLike, latitude: npt.ArrayLike
    ) -> tuple[npt.NDArray[np.intp], npt.NDArray[np.bool_]]:
        """
        Finds the nearest cell to each place, and whether the place is
        outside the grid. Of a regular grid, the nearest cell is that of the
        nearest latitude and longitude, and a place is outside by more than
        half a spacing beyond its outer latitudes or (unless they go around
        the globe) longitudes. Of other grids, the nearest cell is found on the
        sphere, and a place is outside if farther from it than about half the
        diagonal of the cell (its `reach`, with a tolerance for cells that
        vary in size, `_REACH_TOLERANCE`), as beyond the edge of a regional
        grid by more than about half a cell.

        Args:
            longitude (numpy.typing.ArrayLike): The longitudes (degrees).
            latitude (numpy.typing.ArrayLike): The latitudes (degrees).

        Returns:
            tuple[numpy.typing.NDArray, numpy.typing.NDArray]: the flat
                indices of the nearest cells (in the order of the horizontal
                dimensions), and whether each place is outside
        """
        longitude = np.ravel(np.asarray(longitude, dtype=float))
        latitude = np.ravel(np.asarray(latitude, dtype=float))
        if len(longitude) == 0:
            return np.zeros(0, dtype=np.intp), np.zeros(0, dtype=bool)
        if self.kind == "regular":
            rows, row_outside = _nearest_coordinate(self.latitude, latitude, False)
            columns, column_outside = _nearest_coordinate(
                self.longitude, self.to_field_longitudes(longitude), True
            )
            return rows * len(self.longitude) + columns, row_outside | column_outside
        distance, index = self.tree.query(unit_vectors(longitude, latitude))
        outside = distance > _REACH_TOLERANCE * self.reach[index]
        return self.valid[index], outside

    def to_field_longitudes(self, longitude: npt.ArrayLike) -> npt.NDArray[np.float64]:
        """
        Converts longitudes to the convention of the grid: [0, 360) if its
        longitudes exceed 180 degrees, else [-180, 180).

        Args:
            longitude (numpy.typing.ArrayLike): The longitudes (degrees).

        Returns:
            numpy.typing.NDArray: the longitudes (degrees)
        """
        longitude = np.asarray(longitude, dtype=float)
        if np.nanmax(self.longitude) > 180:
            return np.mod(longitude, 360)
        return (longitude + 180) % 360 - 180

    def unravel(self, index: npt.ArrayLike) -> tuple[npt.NDArray[np.intp], ...]:
        """
        Gets the indices along each horizontal dimension of flat indices.

        Args:
            index (numpy.typing.ArrayLike): The flat indices.

        Returns:
            tuple[numpy.typing.NDArray, ...]: the index along each dimension
        """
        return tuple(np.unravel_index(np.asarray(index, dtype=np.intp), self.shape))

    def positions(
        self, index: npt.ArrayLike
    ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        """
        Gets the longitudes and latitudes (degrees) of cells by their flat
        indices (`NaN` for negative indices, of no cell).

        Args:
            index (numpy.typing.ArrayLike): The flat indices.

        Returns:
            tuple[numpy.typing.NDArray, numpy.typing.NDArray]: the longitudes
                and latitudes
        """
        index = np.asarray(index, dtype=np.intp)
        missing = index < 0
        index = np.where(missing, 0, index)
        if self.kind == "regular":
            rows, columns = np.divmod(index, len(self.longitude))
            longitude, latitude = self.longitude[columns], self.latitude[rows]
        else:
            longitude = self.longitude.ravel()[index]
            latitude = self.latitude.ravel()[index]
        return (
            np.where(missing, np.nan, longitude),
            np.where(missing, np.nan, latitude),
        )

    def cell_areas(self, field: Any = None) -> npt.NDArray[np.float64]:
        """
        Gets the area (square meters) of each cell, flattened in the order of
        the horizontal dimensions: that of a coordinate of the field over the
        grid's dimensions named `cell_area`, `area`, or `areaCell` (as of an
        MPAS mesh), if given; or else, of a regular grid, the area between
        the midpoints of neighboring latitudes and longitudes; of a
        curvilinear grid, the area spanned by the steps to neighboring cells
        along each dimension; and of an unstructured grid, a square of the
        cell's diagonal (twice its reach squared).

        Args:
            field (xarray.DataArray | None): The field, with any coordinate of
                the cells' areas.

        Returns:
            numpy.typing.NDArray: the areas (square meters)
        """
        if field is not None:
            for name in _AREA_NAMES:
                if name in field.coords and set(field.coords[name].dims) == set(
                    self.dims
                ):
                    return np.ravel(
                        field.coords[name].transpose(*self.dims).to_numpy()
                    ).astype(float)
        radius = constants.EARTH_MEAN_RADIUS
        if self.kind == "regular":
            row_step = _coordinate_steps(self.latitude)
            column_step = _coordinate_steps(self.longitude)
            areas = (
                radius**2
                * np.radians(row_step)[:, None]
                * np.radians(column_step)[None, :]
                * np.cos(np.radians(self.latitude))[:, None]
            )
            return np.abs(areas).ravel()
        if self.kind == "curvilinear":
            vectors = unit_vectors(self.longitude, self.latitude).reshape(
                self.shape + (3,)
            )
            steps = [
                np.gradient(vectors, axis=axis) if self.shape[axis] > 1 else None
                for axis in (0, 1)
            ]
            if steps[0] is not None and steps[1] is not None:
                return (
                    radius**2 * np.linalg.norm(np.cross(steps[0], steps[1]), axis=-1)
                ).ravel()
        areas = np.full(int(np.prod(self.shape)), np.nan)
        areas[self.valid] = 2 * (self.reach * radius) ** 2
        return areas

    def candidates(
        self, longitude: npt.ArrayLike, latitude: npt.ArrayLike, reach: npt.ArrayLike
    ) -> npt.NDArray[np.intp]:
        """
        Finds the cells near each place: those within a distance (`reach`,
        meters) of it, and others of a window around it (of a regular grid, the
        rows and columns that span the distance), by their flat indices,
        with `-1` for none (padding each place's cells to the same number).

        Args:
            longitude (numpy.typing.ArrayLike): The longitudes (degrees).
            latitude (numpy.typing.ArrayLike): The latitudes (degrees).
            reach (numpy.typing.ArrayLike): The distance (meters) from each
                place.

        Returns:
            numpy.typing.NDArray: the flat indices of the cells, with shape
                (places, cells)
        """
        longitude = np.ravel(np.asarray(longitude, dtype=float))
        latitude = np.ravel(np.asarray(latitude, dtype=float))
        reach = np.broadcast_to(np.asarray(reach, dtype=float), longitude.shape)
        if len(longitude) == 0:
            return np.zeros((0, 0), dtype=np.intp)
        if self.kind == "regular":
            return self._window_candidates(longitude, latitude, reach)
        return self._tree_candidates(longitude, latitude, reach)

    def _window_candidates(
        self,
        longitude: npt.NDArray[np.float64],
        latitude: npt.NDArray[np.float64],
        reach: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.intp]:
        """
        Finds the cells of a regular grid in the rows and columns around each
        place that span a distance from it (see `candidates`).
        """
        angle = np.degrees(reach / constants.EARTH_MEAN_RADIUS)
        row_order = np.argsort(self.latitude, kind="stable")
        column_order = np.argsort(self.longitude, kind="stable")
        rows_sorted = self.latitude[row_order]
        columns_sorted = self.longitude[column_order]
        row_count, column_count = len(rows_sorted), len(columns_sorted)
        row_step = _median_step(rows_sorted)
        column_step = _median_step(columns_sorted)
        around = column_count * column_step >= 360 - 1e-6
        row = np.searchsorted(rows_sorted, latitude)
        row = _nearer(rows_sorted, latitude, row)
        column = _nearer(
            columns_sorted,
            self.to_field_longitudes(longitude),
            np.searchsorted(columns_sorted, self.to_field_longitudes(longitude)),
            around,
        )
        row_width = np.ceil(angle / row_step).astype(int)
        # columns narrow toward the poles, at the farthest latitude reached
        farthest = np.minimum(np.abs(latitude) + angle, 89.9)
        column_width = np.ceil(
            angle / (column_step * np.cos(np.radians(farthest)))
        ).astype(int)
        column_width = np.minimum(
            column_width, (column_count - 1) // 2 if around else column_count
        )
        rows_offset = np.arange(-row_width.max(), row_width.max() + 1)
        columns_offset = np.arange(-column_width.max(), column_width.max() + 1)
        rows = row[:, None] + rows_offset[None, :]
        row_valid = (
            (rows >= 0)
            & (rows < row_count)
            & (np.abs(rows_offset)[None, :] <= row_width[:, None])
        )
        columns = column[:, None] + columns_offset[None, :]
        column_valid = np.abs(columns_offset)[None, :] <= column_width[:, None]
        if around:
            columns = np.mod(columns, column_count)
        else:
            column_valid &= (columns >= 0) & (columns < column_count)
        rows = row_order[np.clip(rows, 0, row_count - 1)]
        columns = column_order[np.clip(columns, 0, column_count - 1)]
        index = rows[:, :, None] * column_count + columns[:, None, :]
        valid = row_valid[:, :, None] & column_valid[:, None, :]
        return np.where(valid, index, -1).reshape(len(longitude), -1)

    def _tree_candidates(
        self,
        longitude: npt.NDArray[np.float64],
        latitude: npt.NDArray[np.float64],
        reach: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.intp]:
        """
        Finds the cells of a grid that is not regular within a distance of
        each place (see `candidates`), as the nearest cells (by a number
        estimated from the cells' size, doubled while any place has that many
        within its distance).
        """
        vectors = unit_vectors(longitude, latitude)
        angle = reach / constants.EARTH_MEAN_RADIUS
        chord = 2 * np.sin(np.minimum(angle, np.pi) / 2)
        cell = np.median(2 * self.reach[np.isfinite(self.reach)] ** 2)
        total = len(self.valid)
        count = int(min(total, np.ceil(1.5 * np.pi * np.max(angle) ** 2 / cell) + 4))
        while True:
            distance, index = self.tree.query(
                vectors, k=count, distance_upper_bound=float(np.max(chord))
            )
            distance = np.reshape(distance, (len(vectors), -1))
            index = np.reshape(index, (len(vectors), -1))
            if count >= total or not np.any(np.isfinite(distance[:, -1])):
                break
            count = min(total, 2 * count)
        found = np.isfinite(distance) & (distance <= chord[:, None])
        return np.where(found, self.valid[np.minimum(index, total - 1)], -1)


def _coordinate_steps(values: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Gets the spacing of a coordinate at each value: the distance between
    the midpoints to its neighbors (in their sorted order)."""
    if len(values) < 2:
        return np.ones(len(values))
    order = np.argsort(values, kind="stable")
    steps = np.empty(len(values))
    steps[order] = np.gradient(values[order])
    return steps


def _median_step(values: npt.NDArray[np.float64]) -> float:
    """Gets the median spacing of sorted coordinate values."""
    if len(values) < 2:
        return np.inf
    return float(np.median(np.diff(values)))


def _nearer(
    ordered: npt.NDArray[np.float64],
    values: npt.NDArray[np.float64],
    after: npt.NDArray[np.intp],
    around: bool = False,
) -> npt.NDArray[np.intp]:
    """
    Gets the index of the nearest of sorted coordinate values to each value,
    from the index after it (see `numpy.searchsorted`), around the globe in
    longitude, if the values go around it.
    """
    count = len(ordered)
    if count == 1:
        return np.zeros(len(values), dtype=np.intp)
    if around:
        below = np.mod(after - 1, count)
        above = np.mod(after, count)
        gap_below = np.mod(values - ordered[below], 360)
        gap_above = np.mod(ordered[above] - values, 360)
        return np.where(gap_below < gap_above, below, above)
    below = np.clip(after - 1, 0, count - 1)
    above = np.clip(after, 0, count - 1)
    return np.where(
        np.abs(values - ordered[below]) <= np.abs(ordered[above] - values),
        below,
        above,
    )


def _nearest_coordinate(
    coordinate: npt.NDArray[np.float64],
    values: npt.NDArray[np.float64],
    longitude: bool,
) -> tuple[npt.NDArray[np.intp], npt.NDArray[np.bool_]]:
    """
    Gets the index of the nearest value of a coordinate of a regular grid to
    each value, and whether each is outside the coordinate's extent: half a
    spacing beyond its outer values (unless, of longitudes, they go around
    the globe).
    """
    order = np.argsort(coordinate, kind="stable")
    ordered = coordinate[order]
    step = _median_step(ordered)
    around = longitude and len(ordered) * step >= 360 - 1e-6
    nearest = _nearer(ordered, values, np.searchsorted(ordered, values), around)
    if around or len(ordered) < 2:
        return order[nearest], np.zeros(len(values), dtype=bool)
    low = ordered[0] - (ordered[1] - ordered[0]) / 2
    high = ordered[-1] + (ordered[-1] - ordered[-2]) / 2
    return order[nearest], (values < low) | (values > high)

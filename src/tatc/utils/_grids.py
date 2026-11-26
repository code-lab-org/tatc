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

LATITUDE_NAMES = ("lat", "latitude", "XLAT", "latCell", "nav_lat")
"""The usual names of the latitudes of a field's grid."""

LONGITUDE_NAMES = ("lon", "longitude", "XLONG", "lonCell", "nav_lon")
"""The usual names of the longitudes of a field's grid."""


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
        outside the grid: farther from its nearest cell than about half the
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
        vectors = unit_vectors(longitude, latitude)
        if len(vectors) == 0:
            return np.zeros(0, dtype=np.intp), np.zeros(0, dtype=bool)
        distance, index = self.tree.query(vectors)
        outside = distance > _REACH_TOLERANCE * self.reach[index]
        return self.valid[index], outside

    def unravel(self, index: npt.ArrayLike) -> tuple[npt.NDArray[np.intp], ...]:
        """
        Gets the indices along each horizontal dimension of flat indices.

        Args:
            index (numpy.typing.ArrayLike): The flat indices.

        Returns:
            tuple[numpy.typing.NDArray, ...]: the index along each dimension
        """
        return tuple(np.unravel_index(np.asarray(index, dtype=np.intp), self.shape))

"""
Geometry utility functions.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from typing import overload

import geopandas as gpd
import numpy as np
from pyproj import Geod
from shapely import make_valid
from shapely.geometry import (
    GeometryCollection,
    LineString,
    MultiPolygon,
    Point,
    Polygon,
)
from shapely.ops import split

# WGS 84 ellipsoid geodesic solver, shared across calls
_WGS84_GEOD = Geod(ellps="WGS84")


def geodesic_distance(
    longitude_1: float, latitude_1: float, longitude_2: float, latitude_2: float
) -> float:
    """
    Computes the geodesic surface distance between two longitude/latitude
    points on the WGS 84 ellipsoid: the length of the shortest path
    between them that stays on the ellipsoid surface. Unlike a
    longitude/latitude-based (planar) distance, this accounts for the
    Earth's oblateness and the convergence of meridians toward the poles,
    so it remains accurate at any latitude or longitude separation
    (including antipodal-ish or antimeridian-spanning point pairs).

    Args:
        longitude_1 (float): Longitude (degrees) of the first point.
        latitude_1 (float): Latitude (degrees) of the first point.
        longitude_2 (float): Longitude (degrees) of the second point.
        latitude_2 (float): Latitude (degrees) of the second point.

    Returns:
        float: The geodesic distance (meters) between the two points.
    """
    _, _, distance = _WGS84_GEOD.inv(longitude_1, latitude_1, longitude_2, latitude_2)
    return distance


@overload
def project_polygon_to_elevation(polygon: Polygon, elevation: float) -> Polygon: ...


@overload
def project_polygon_to_elevation(
    polygon: MultiPolygon, elevation: float
) -> MultiPolygon: ...


def project_polygon_to_elevation(
    polygon: Polygon | MultiPolygon, elevation: float
) -> Polygon | MultiPolygon:
    """
    Assigns a fixed z-coordinate (elevation) to every coordinate of a
    polygon or multipolygon, including exterior and interior (hole) rings.
    Any existing z-coordinate is overwritten, not offset.

    Args:
        polygon (shapely.geometry.Polygon | shapely.geometry.MultiPolygon): The polygon to project.
        elevation (float): The elevation (meters) above the WGS 84 geoid to
            assign to every coordinate.

    Returns:
        shapely.geometry.Polygon | shapely.geometry.MultiPolygon: The projected
        polygon, matching the input type.
    """
    if isinstance(polygon, Polygon):
        return Polygon(
            [(p[0], p[1], elevation) for p in polygon.exterior.coords],
            [[(p[0], p[1], elevation) for p in i.coords] for i in polygon.interiors],
        )
    return MultiPolygon(
        [project_polygon_to_elevation(g, elevation) for g in polygon.geoms]
    )


def _flatten_polygons(pgons: list[Polygon | MultiPolygon]) -> list[Polygon]:
    """
    Flattens a list of Polygon/MultiPolygon geometries into a flat list of
    Polygon, unpacking any MultiPolygon into its constituent polygons.

    Args:
       pgons (list[shapely.geometry.Polygon | shapely.geometry.MultiPolygon]):
           The geometries to flatten.

    Returns:
       list[shapely.geometry.Polygon]: The flattened list of polygons.
    """
    return [g for p in pgons for g in (p.geoms if isinstance(p, MultiPolygon) else [p])]


def _wrap_polygon_over_pole(
    polygon: Polygon | MultiPolygon, pole: int
) -> Polygon | MultiPolygon:
    """
    Wraps polygon coordinates over a pole (`pole` = 1 for the North pole,
    -1 for the South pole). Due to buffering and projection, sometimes
    latitudes exceed the pole (i.e. exceed 90 * pole degrees). This method
    wraps them to the correct latitude between -90 and 90 degrees and
    adjusts the longitude by 180 degrees. Only coordinates exceeding the
    pole are shifted; other coordinates are left unchanged.

    This method requires a polygon exceeding the pole to be confined to one
    side of the prime meridian: this is guaranteed by `_split_polygon_over_pole`,
    which is the only caller, since it splits off any piece straddling the
    prime meridian before wrapping. A polygon violating this precondition
    would produce a self-intersecting (invalid) ring.

    Note: this method only changes coordinates: it does not create a MultiPolygon.

    Args:
       polygon (shapely.geometry.Polygon | shapely.geometry.MultiPolygon): The polygon to wrap.
       pole (int): 1 for the North pole, -1 for the South pole.

    Returns:
       shapely.geometry.Polygon | shapely.geometry.MultiPolygon: The wrapped polygon.
    """
    if isinstance(polygon, Polygon):
        if all(c[1] * pole <= 90 for c in polygon.exterior.coords):
            # no wrapping necessary
            return polygon
        # map latitudes beyond the pole back between -90 and 90, adjusting longitude by 180 degrees
        lat_shift = 180 if all(c[0] <= 0 for c in polygon.exterior.coords) else -180
        return Polygon(
            [
                [
                    c[0] + lat_shift if c[1] * pole >= 90 else c[0],
                    pole * 180 - c[1] if c[1] * pole >= 90 else c[1],
                ]
                for c in polygon.exterior.coords
            ],
            [
                [
                    [
                        c[0] + lat_shift if c[1] * pole >= 90 else c[0],
                        pole * 180 - c[1] if c[1] * pole >= 90 else c[1],
                    ]
                    for c in i.coords
                ]
                for i in polygon.interiors
            ],
        )
    # recursive call for each polygon
    return MultiPolygon(
        _flatten_polygons([_wrap_polygon_over_pole(p, pole) for p in polygon.geoms])
    )


def _split_polygon_over_pole(
    polygon: Polygon | MultiPolygon, pole: int
) -> Polygon | MultiPolygon:
    """
    Splits a polygon that encompasses a pole (`pole` = 1 for the North
    pole, -1 for the South pole; i.e. exceeds 90 * pole degrees latitude)
    into a valid MultiPolygon on the standard (-180, -90, 180, 90) plane.
    The polygon is first split along the pole latitude; any resulting
    piece that also straddles the prime meridian (0 degrees longitude)
    while beyond the pole is split again there, since such a piece cannot
    be wrapped to one side in a single step. Each piece exceeding the pole
    latitude is then wrapped to its correct latitude/longitude. A polygon
    that does not exceed the pole latitude is returned unchanged.

    Args:
       polygon (shapely.geometry.Polygon | shapely.geometry.MultiPolygon): The polygon to split.
       pole (int): 1 for the North pole, -1 for the South pole.

    Returns:
       shapely.geometry.Polygon | shapely.geometry.MultiPolygon: The split polygon.
    """
    if isinstance(polygon, Polygon):
        if all(c[1] * pole <= 90 for c in polygon.exterior.coords):
            # no splitting necessary
            return polygon
        # split polygon along the pole
        pole_line = LineString([(-360, 90 * pole), (360, 90 * pole)])
        parts = split(polygon, pole_line)
        # check and split part over prime meridian if necessary
        meridian_line = LineString([(0, 90 * pole), (0, 180 * pole)])
        for part in parts.geoms:
            if part.crosses(meridian_line):
                parts = GeometryCollection(
                    [g for g in parts.geoms if g != part]
                    + list(split(part, meridian_line).geoms)
                )
        # convert to a multi polygon
        if isinstance(parts, GeometryCollection):
            parts = _convert_collection_to_polygon(parts)
        # return polygon with components wrapped over the pole
        return _wrap_polygon_over_pole(parts, pole)
    # recursive call for each polygon
    return MultiPolygon(
        _flatten_polygons([_split_polygon_over_pole(p, pole) for p in polygon.geoms])
    )


@overload
def _wrap_polygon_over_antimeridian(polygon: Polygon) -> Polygon: ...


@overload
def _wrap_polygon_over_antimeridian(polygon: MultiPolygon) -> MultiPolygon: ...


def _wrap_polygon_over_antimeridian(
    polygon: Polygon | MultiPolygon,
) -> Polygon | MultiPolygon:
    """
    Wraps polygon coordinates over the antimeridian. Due to buffering and projection,
    sometimes longitudes exceed 180 degrees. This method wraps them to
    the correct longitude between -180 and 180 degrees by adding or
    subtracting 360 degrees to every coordinate.

    This method requires all coordinates to be at or beyond the same side
    of the antimeridian (all at or below -180 degrees, or all at or above
    180 degrees): this is guaranteed by `_split_polygon_antimeridian`, which
    is the only caller, since it splits along a single meridian line before
    wrapping, confining each resulting piece to one side.

    Note: this method only changes coordinates: it does not create a MultiPolygon.

    Args:
       polygon (shapely.geometry.Polygon | shapely.geometry.MultiPolygon): The polygon to wrap.

    Returns:
       shapely.geometry.Polygon | shapely.geometry.MultiPolygon: The wrapped polygon.
    """
    if isinstance(polygon, Polygon):
        if all(c[0] >= -180 and c[0] <= 180 for c in polygon.exterior.coords):
            # no wrapping necessary
            return polygon
        if all(c[0] <= -180 for c in polygon.exterior.coords):
            # map longitudes from (-540, -180] to (-180, 180]
            return Polygon(
                [[c[0] + 360, c[1]] for c in polygon.exterior.coords],
                [[[c[0] + 360, c[1]] for c in i.coords] for i in polygon.interiors],
            )
        # map longitudes from [180, 540) to [-180, 180)
        return Polygon(
            [[c[0] - 360, c[1]] for c in polygon.exterior.coords],
            [[[c[0] - 360, c[1]] for c in i.coords] for i in polygon.interiors],
        )
    # recursive call for each polygon
    return MultiPolygon(
        _flatten_polygons([_wrap_polygon_over_antimeridian(p) for p in polygon.geoms])
    )


def _convert_collection_to_polygon(
    collection: GeometryCollection,
) -> Polygon | MultiPolygon:
    """
    Converts a GeometryCollection to a Polygon or MultiPolygon. Quick clipping
    can create dirty results with points or lines on boundaries. This method
    drops any points or lines from a GeometryCollection, flattens any
    MultiPolygon members into their constituent polygons, and returns a
    single Polygon if only one remains, or a MultiPolygon otherwise
    (empty if no polygons remain).

    Args:
       collection (shapely.geometry.GeometryCollection): The geometry collection to convert.

    Returns:
       shapely.geometry.Polygon | shapely.geometry.MultiPolygon: The converted polygon.
    """
    pgons = [p for p in collection.geoms if isinstance(p, Polygon)] + [
        p for g in collection.geoms if isinstance(g, MultiPolygon) for p in g.geoms
    ]
    if len(pgons) == 1:
        return pgons[0]
    return MultiPolygon(pgons)


def _split_polygon_antimeridian(
    polygon: Polygon | MultiPolygon,
) -> Polygon | MultiPolygon:
    """
    Splits a polygon that crosses the antimeridian (180 degrees longitude)
    into a valid MultiPolygon on the standard (-180, 180) longitude range.
    A crossing is detected when adjacent exterior vertices jump by 180
    degrees of longitude or more.

    Two cases are handled differently:

    - If the polygon's vertex longitudes wrap all the way around the globe
      (e.g. a polar cap that does not itself exceed +/-90 degrees
      latitude), it is reconstructed with a flattened edge at the pole and
      split along the prime meridian instead of the antimeridian. Note:
      the raw result of this case may be reported as invalid (the two
      pieces touch along the shared prime-meridian cut edge); the public
      `split_polygon` function repairs this via `shapely.make_valid`.
    - Otherwise, coordinates are "unrolled" past +/-180 degrees according to
      the cumulative crossing direction, split along the antimeridian, and
      wrapped back with `_wrap_polygon_over_antimeridian`.

    A polygon with no detected crossing is returned unchanged. Note: this
    function only supports polygons that span LESS than 360 degrees longitude.

    Args:
       polygon (shapely.geometry.Polygon | shapely.geometry.MultiPolygon): The polygon to split.

    Returns:
       shapely.geometry.Polygon | shapely.geometry.MultiPolygon: The split polygon.
    """
    if isinstance(polygon, Polygon):
        lon = np.array([c[0] for c in polygon.exterior.coords])
        # check if any longitudes cross the anti-meridian
        # (adjacent coordinate longitude differs by more than 180 degrees)
        if all(np.abs(np.diff(lon)) < 180):
            return polygon
        # check if this polygon contains a pole
        if Polygon(zip(np.cos(np.radians(lon)), np.sin(np.radians(lon)))).contains(
            Point(0, 0)
        ):
            # extract (lon, lat) only, discarding any z-dimension, and sort by longitude
            coords = [(c[0], c[1]) for c in polygon.exterior.coords[0:-1]]
            coords.sort(key=lambda r: r[0])
            # determine if contains north or south pole based on sign of mean latitude
            n_s = 1 if np.array(coords)[:, 1].mean() > 0 else -1
            # interpolate latitude at antimeridian
            lat = np.interp(
                180, [coords[-1][0], coords[0][0] + 180], [coords[-1][1], coords[0][1]]
            )
            # reconstruct polygon (ccw) with added coords on antimeridian
            pgon = Polygon(
                [(-180, 90 * n_s), (-180, lat)]
                + coords
                + [(180, lat), (180, 90 * n_s), (-180, 90 * n_s)],
                [
                    [(c[0], c[1]) for c in interior.coords]
                    for interior in polygon.interiors
                ],
            )
            # return polygon split down prime meridian to improve handling
            parts = split(pgon, LineString([(0, -180), (0, 180)]))
            # convert to multi polygon
            if isinstance(parts, GeometryCollection):
                parts = _convert_collection_to_polygon(parts)
            return parts
        # find anti-meridian crossings and calculate shift direction
        # coords from W -> E (shift < 0) will add 360 degrees to E component
        # coords from E -> W (shift > 0) will subtract 360 degrees from W component
        shift = np.insert(np.cumsum(np.around(np.diff(lon) / 360)), 0, 0)
        pgon = Polygon(
            [
                (c[0] - 360 * shift[i], c[1])
                for i, c in enumerate(polygon.exterior.coords)
            ],
            [
                [
                    (
                        ic[0]
                        - 360 * np.interp(ic[0], np.sort(lon), shift[np.argsort(lon)]),
                        ic[1],
                    )
                    for ic in i.coords
                ]
                for i in polygon.interiors
            ],
        )
        # split along the anti-meridian (-180 for shift > 0; 180 for shift < 0)
        shift_dir = -180 if shift.max() >= 1 else 180
        parts = split(pgon, LineString([(shift_dir, -180), (shift_dir, 180)]))
        # convert to multi polygon
        if isinstance(parts, GeometryCollection):
            parts = _convert_collection_to_polygon(parts)
        # return polygon with components wrapped over anti-meridian
        return _wrap_polygon_over_antimeridian(parts)
    if isinstance(polygon, MultiPolygon):
        # recursive call for each polygon
        return MultiPolygon(
            _flatten_polygons([_split_polygon_antimeridian(p) for p in polygon.geoms])
        )
    raise ValueError("Unknown geometry: " + str(type(polygon)))


def split_polygon(
    polygon: Polygon | MultiPolygon,
) -> Polygon | MultiPolygon:
    """
    Splits a Polygon into a MultiPolygon if it crosses the anti-meridian
    (180 degrees longitude), exceeds the north pole (90 degrees latitude), or
    exceeds the south pole (-90 degrees latitude). Note: this function
    only supports polygons that span LESS than 360 degrees longitude.
    Operates on (longitude, latitude) only: any z-dimension on the input
    is discarded. Use `project_polygon_to_elevation` to add elevation back
    after splitting.

    Args:
        polygon (shapely.geometry.Polygon | shapely.geometry.MultiPolygon): The polygon to split.

    Returns:
        shapely.geometry.Polygon | shapely.geometry.MultiPolygon: The split polygon.
    """
    polygon = _split_polygon_over_pole(
        _split_polygon_over_pole(_split_polygon_antimeridian(polygon), pole=-1),
        pole=1,
    )
    # invalid polygons can arise from narrow sensor geometries in polar regions
    if not polygon.is_valid:
        # try to fix geometry
        polygon = make_valid(polygon)  # type: ignore
        if isinstance(polygon, GeometryCollection):
            polygon = _convert_collection_to_polygon(polygon)
    return polygon


def get_planar_bounds(
    mask: Polygon | MultiPolygon | None,
) -> tuple[float, float, float, float]:
    """
    Generates a tuple of bounds for a polygon mask.

    Known limitation: this method assumes `mask` lies on the planar
    (non-antimeridian-crossing) longitude domain. A mask crossing the
    antimeridian with longitude below -180 (e.g. bounds spanning
    -190 to -180, representing the same region as 170 to 180) only has
    its max_longitude corrected to 180; min_longitude is left unadjusted,
    so the resulting bounds do not coherently describe such a mask.
    Fully supporting antimeridian-crossing masks can be achieved by
    pre-processing the geometry (e.g. via `split_polygon`).

    Args:
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | None):
            Geometric shape using WGS84 (EPSG:4326)
            geodetic coordinates in a Polygon or MultiPolygon.

    Returns:
        tuple[float, float, float, float]: min longitude (degrees),
            min latitude (degrees), max longitude (degrees), max latitude (degrees)
    """
    if isinstance(mask, (Polygon, MultiPolygon)):
        if not mask.is_valid:
            raise ValueError("Mask is not a valid Polygon or MultiPolygon.")
        total_bounds = mask.bounds
    else:
        total_bounds = [-180, -90, 180, 90]
    min_longitude = total_bounds[0]
    min_latitude = total_bounds[1]
    max_longitude = 180 if total_bounds[2] == -180 else total_bounds[2]
    max_latitude = total_bounds[3]
    return (min_longitude, min_latitude, max_longitude, max_latitude)


def normalize_geometry(
    geometry: Polygon | MultiPolygon | gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    """
    Normalize geometry to a GeoDataFrame with antimeridian wrapping.

    Args:
        geometry (shapely.geometry.Polygon | shapely.geometry.MultiPolygon |
            geopandas.GeoDataFrame): The geometry to normalize.

    Returns:
        geopandas.GeoDataFrame: The normalized geometry.
    """
    if isinstance(geometry, (Polygon, MultiPolygon)):
        if not geometry.is_valid:
            raise ValueError("Geometry is not a valid Polygon or MultiPolygon.")
        geometry = gpd.GeoDataFrame(geometry=gpd.GeoSeries([geometry]), crs="EPSG:4326")
    elif isinstance(geometry, gpd.GeoSeries):
        geometry = gpd.GeoDataFrame(geometry=geometry, crs="EPSG:4326")
    if isinstance(geometry, gpd.GeoDataFrame):
        geometry["geometry"] = geometry.apply(
            lambda r: split_polygon(r.geometry),
            axis=1,
        )
    return geometry

"""
Geometry utility functions.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import hashlib
from typing import Any, overload

import geopandas as gpd
import numpy as np
import numpy.typing as npt
import shapely
from pyproj import Geod
from shapely import make_valid
from shapely.geometry import (
    GeometryCollection,
    LineString,
    MultiPolygon,
    Point,
    Polygon,
)
from shapely.geometry.base import BaseGeometry
from shapely.ops import split

from .ellipsoid import _get_surface_directions

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


def geodesic_destination(
    longitude: float, latitude: float, azimuth: float, distance: float
) -> tuple[float, float]:
    """
    Computes the longitude/latitude reached by traveling a specified
    geodesic distance from a starting point at a specified initial
    azimuth, on the WGS 84 ellipsoid.

    Args:
        longitude (float): Longitude (degrees) of the starting point.
        latitude (float): Latitude (degrees) of the starting point.
        azimuth (float): Initial azimuth (degrees, clockwise from north)
            of travel.
        distance (float): Geodesic distance (meters) to travel.

    Returns:
        tuple[float, float]: The destination (longitude, latitude), in degrees.
    """
    destination_longitude, destination_latitude, _ = _WGS84_GEOD.fwd(
        longitude, latitude, azimuth, distance
    )
    return destination_longitude, destination_latitude


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
    return shapely.force_3d(shapely.force_2d(polygon), elevation)


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
        exterior = shapely.get_coordinates(polygon.exterior)
        if np.all(exterior[:, 1] * pole <= 90):
            # no wrapping necessary
            return polygon
        # map latitudes beyond the pole back between -90 and 90, adjusting longitude by 180 degrees
        lat_shift = 180 if np.all(exterior[:, 0] <= 0) else -180

        def wrap(coords: npt.NDArray) -> npt.NDArray:
            beyond = coords[:, 1] * pole >= 90
            return np.column_stack(
                [
                    np.where(beyond, coords[:, 0] + lat_shift, coords[:, 0]),
                    np.where(beyond, pole * 180 - coords[:, 1], coords[:, 1]),
                ]
            )

        return shapely.transform(polygon, wrap)
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
        if np.all(shapely.get_coordinates(polygon.exterior)[:, 1] * pole <= 90):
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
        longitude = shapely.get_coordinates(polygon.exterior)[:, 0]
        if np.all((longitude >= -180) & (longitude <= 180)):
            # no wrapping necessary
            return polygon
        # map longitudes from (-540, -180] to (-180, 180], or from [180, 540)
        # to [-180, 180)
        offset = 360 if np.all(longitude <= -180) else -360
        return shapely.transform(polygon, lambda coords: coords + [offset, 0])
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
    A crossing is either a jump, where adjacent exterior vertices differ by
    more than 180 degrees of longitude (the shorter way around the globe
    crosses the antimeridian), or a coordinate beyond +/-180 degrees (the
    exterior continues past the antimeridian). A seam edge between vertices
    both on the antimeridian (e.g. from -180 to 180 degrees along the edge of
    a polygon spanning all longitudes) is not a jump, unless the exterior
    encircles a pole.

    Two cases are handled differently:

    - If the exterior encircles a pole (its longitude winds a full 360
      degrees around the globe, e.g. a polar cap that does not itself
      exceed +/-90 degrees latitude), it is reconstructed with a flattened
      edge at the pole and split along the prime meridian instead of the
      antimeridian. Note: the raw result of this case may be reported as
      invalid (the two pieces touch along the shared prime-meridian cut
      edge); the public `split_polygon` function repairs this via
      `shapely.make_valid`.
    - Otherwise, coordinates are "unrolled" past +/-180 degrees according to
      the cumulative jump direction, split along the antimeridian that the
      unrolled coordinates cross, and wrapped back with
      `_wrap_polygon_over_antimeridian`.

    A polygon with no crossing is returned unchanged. Note: this function
    only supports polygons that span LESS than 360 degrees longitude, other
    than polar caps and polygons spanning exactly -180 to 180 degrees.

    Args:
       polygon (shapely.geometry.Polygon | shapely.geometry.MultiPolygon): The polygon to split.

    Returns:
       shapely.geometry.Polygon | shapely.geometry.MultiPolygon: The split polygon.
    """
    if isinstance(polygon, Polygon):
        exterior = shapely.get_coordinates(polygon.exterior)
        lon = exterior[:, 0]
        diff = np.diff(lon)
        # jumps: adjacent coordinate longitudes differ by more than 180 degrees,
        # distinguishing seam edges between vertices both on the anti-meridian
        jump = np.abs(diff) > 180
        seam = (np.abs(lon[:-1]) == 180) & (np.abs(lon[1:]) == 180)
        # net longitude winding of the exterior (+/-360 degrees around a pole),
        # taking each jump (including seam edges) the shorter way around
        winding = np.sum(diff - 360 * np.where(jump, np.round(diff / 360), 0))
        # check if this polygon contains a pole
        if np.abs(winding) > 180:
            # extract (lon, lat) only, discarding any z-dimension, and sort by longitude
            coords = exterior[:-1][np.argsort(exterior[:-1, 0], kind="stable")]
            # determine if contains north or south pole based on sign of mean latitude
            n_s = 1 if coords[:, 1].mean() > 0 else -1
            # interpolate latitude at antimeridian
            lat = np.interp(
                180, [coords[-1, 0], coords[0, 0] + 180], [coords[-1, 1], coords[0, 1]]
            )
            # reconstruct polygon (ccw) with added coords on antimeridian and a
            # flattened edge along the pole, spanning -180 to 180 degrees
            return Polygon(
                np.concatenate(
                    [
                        [[-180, 90 * n_s], [-180, lat]],
                        coords,
                        [[180, lat], [180, 90 * n_s], [-180, 90 * n_s]],
                    ]
                ),
                [shapely.get_coordinates(interior) for interior in polygon.interiors],
            )
        # find anti-meridian crossings and calculate shift direction
        # coords from W -> E (shift < 0) will add 360 degrees to E component
        # coords from E -> W (shift > 0) will subtract 360 degrees from W component
        # (seam edges are not crossings outside of a polar cap)
        jump &= ~seam
        shift = np.insert(np.cumsum(np.where(jump, np.round(diff / 360), 0)), 0, 0)
        unrolled = lon - 360 * shift
        if np.all(np.abs(unrolled) <= 180):
            # no jumps and no coordinates beyond the anti-meridian
            return polygon
        pgon = Polygon(
            np.column_stack([unrolled, exterior[:, 1]]),
            [
                np.column_stack(
                    [
                        ic[:, 0]
                        - 360
                        * np.interp(ic[:, 0], np.sort(lon), shift[np.argsort(lon)]),
                        ic[:, 1],
                    ]
                )
                for ic in (shapely.get_coordinates(i) for i in polygon.interiors)
            ],
        )
        # split along the anti-meridian that the unrolled coordinates cross
        shift_dir = -180 if unrolled.min() < -180 else 180
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


def _is_within_planar_domain(polygon: Polygon | MultiPolygon) -> bool:
    """
    Checks, vectorized across coordinates, whether a polygon lies within the
    standard (-180, -90, 180, 90) longitude and latitude domain without
    crossing the anti-meridian (no adjacent coordinates, even across rings,
    differ by more than 180 degrees of longitude), in which case it needs no
    splitting (see `split_polygon`).

    Args:
        polygon (shapely.geometry.Polygon | shapely.geometry.MultiPolygon): The polygon.

    Returns:
        bool: True, if the polygon needs no splitting
    """
    if not isinstance(polygon, (Polygon, MultiPolygon)):
        return False
    coordinates = shapely.get_coordinates(polygon)
    return bool(
        np.all(np.abs(coordinates[:, 0]) <= 180)
        and np.all(np.abs(coordinates[:, 1]) <= 90)
        and np.all(np.abs(np.diff(coordinates[:, 0])) <= 180)
    )


def split_polygon(
    polygon: Polygon | MultiPolygon,
) -> Polygon | MultiPolygon:
    """
    Splits a Polygon into a MultiPolygon if it crosses the anti-meridian
    (180 degrees longitude), exceeds the north pole (90 degrees latitude), or
    exceeds the south pole (-90 degrees latitude). The anti-meridian may be
    crossed either by a jump between adjacent vertices (e.g. from 170 to
    -170 degrees) or by coordinates beyond +/-180 degrees (e.g. from 170 to
    190 degrees). Note: this function only supports polygons that span LESS
    than 360 degrees longitude, other than polar caps (encircling a pole)
    and polygons spanning exactly -180 to 180 degrees. Operates on (longitude, latitude) only: any z-dimension on the input
    is discarded. Use `project_polygon_to_elevation` to add elevation back
    after splitting.

    Args:
        polygon (shapely.geometry.Polygon | shapely.geometry.MultiPolygon): The polygon to split.

    Returns:
        shapely.geometry.Polygon | shapely.geometry.MultiPolygon: The split polygon.
    """
    if not _is_within_planar_domain(polygon):
        if isinstance(polygon, MultiPolygon) and any(
            len(part.interiors) > 0 for part in polygon.geoms
        ):
            # split each part (with its holes) separately
            polygon = MultiPolygon(
                _flatten_polygons([split_polygon(part) for part in polygon.geoms])
            )
        elif isinstance(polygon, Polygon) and len(polygon.interiors) > 0:
            # split the exterior and each hole as polygons of their own (so
            # that holes across the anti-meridian or around a pole are split
            # like exteriors) and subtract the holes from the exterior
            polygon = split_polygon(Polygon(polygon.exterior)).difference(
                shapely.union_all(
                    [split_polygon(Polygon(ring)) for ring in polygon.interiors]
                )
            )
            if isinstance(polygon, GeometryCollection):
                polygon = _convert_collection_to_polygon(polygon)
        else:
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
        if isinstance(polygon, MultiPolygon):
            # drop degenerate (zero-area) parts left by the repair
            parts = [p for p in polygon.geoms if p.area > 1e-20]
            polygon = parts[0] if len(parts) == 1 else MultiPolygon(parts)
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


def hash_geometry(geometry: BaseGeometry) -> str:
    """
    Computes a compact hash of a geometry, to identify it (for example, the
    region of interest of an observation) without storing it: the first 16
    hexadecimal digits (64 bits) of the BLAKE2b digest of its normalized
    well-known binary (little-endian, with z coordinates, if any). Geometries
    with the same normalized coordinates have the same hash, regardless of
    the order of their vertices, rings, or parts.

    Args:
        geometry (shapely.geometry.base.BaseGeometry): The geometry.

    Returns:
        str: The hash, as 16 hexadecimal digits.
    """
    return hashlib.blake2b(
        shapely.to_wkb(shapely.normalize(geometry), byte_order=1), digest_size=8
    ).hexdigest()


def _get_boundary_arcs(
    geometry: Polygon | MultiPolygon, elevation: float = 0, max_segment: float = 1
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """
    Gets the boundary of a (split, see `split_polygon`) geometry as great
    circle arcs between geocentric unit vectors (see
    `tatc.utils.ellipsoid._get_surface_directions`). Edges, straight in longitude and latitude,
    are first divided into segments of at most `max_segment` degrees, over
    which a great circle arc departs from them by about a thousandth of a
    degree at most. Degenerate edges (for example, along a pole) and seams
    along the anti-meridian between parts of a split geometry, which are
    not boundaries, are omitted.

    Args:
        geometry (shapely.geometry.Polygon | shapely.geometry.MultiPolygon): The geometry.
        elevation (float): The elevation (meters) above the WGS 84 ellipsoid.
        max_segment (float): The maximum length (degrees) of a segment.

    Returns:
        tuple[numpy.typing.NDArray[numpy.float64], numpy.typing.NDArray[numpy.float64]]:
            the unit vectors (shape (3, E)) at the start and end of each arc
    """
    densified = shapely.segmentize(geometry, max_segment)
    polygons = densified.geoms if isinstance(densified, MultiPolygon) else [densified]
    rings = [
        np.asarray(ring.coords)[:, :2]
        for polygon in polygons
        for ring in [polygon.exterior, *polygon.interiors]
    ]
    if len(rings) == 0:
        return np.empty((3, 0)), np.empty((3, 0))
    start = np.concatenate([ring[:-1] for ring in rings])
    end = np.concatenate([ring[1:] for ring in rings])
    # seams: edges along the anti-meridian where the geometry continues
    # across it (on the opposite side of the anti-meridian)
    on_meridian = (np.abs(start[:, 0]) == 180) & (end[:, 0] == start[:, 0])
    seam = np.zeros(len(start), dtype=bool)
    seam[on_meridian] = shapely.intersects_xy(
        geometry,
        -start[on_meridian, 0],
        (start[on_meridian, 1] + end[on_meridian, 1]) / 2,
    )
    start_u = _get_surface_directions(start[:, 0], start[:, 1], elevation)
    end_u = _get_surface_directions(end[:, 0], end[:, 1], elevation)
    degenerate = np.linalg.norm(np.cross(start_u.T, end_u.T), axis=1) < 1e-12
    keep = ~seam & ~degenerate
    return start_u[:, keep], end_u[:, keep]


def _get_nearest_arc_points(
    directions: npt.NDArray[np.float64],
    arcs: tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]],
    chunk_size: int = 1024,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """
    Computes the minimum angular distance from each of a set of directions
    to a set of great circle arcs (see `_get_boundary_arcs`), each shorter
    than a half circle, and the nearest point on the arcs.

    Args:
        directions (numpy.typing.NDArray[numpy.float64]): The unit vectors (shape (3, N)).
        arcs (tuple[numpy.typing.NDArray[numpy.float64], numpy.typing.NDArray[numpy.float64]]):
            The unit vectors (shape (3, E)) at the start and end of each arc.
        chunk_size (int): The number of directions processed at once, to limit memory.

    Returns:
        tuple[numpy.typing.NDArray[numpy.float64], numpy.typing.NDArray[numpy.float64]]:
            the angular distances (radians, shape (N,)) and the unit vectors
            toward the nearest points (shape (3, N); the directions
            themselves if there are no arcs)
    """
    start, end = arcs
    distance = np.full(directions.shape[1], np.pi)
    nearest = np.array(directions, dtype=float, copy=True)
    if start.shape[1] == 0:
        return distance, nearest
    # unit normal to each arc's great circle and, within its plane, the
    # normals to the arc's ends (pointing into the arc)
    normal = np.cross(start.T, end.T).T
    normal /= np.linalg.norm(normal, axis=0)
    after_start = np.cross(normal.T, start.T).T
    before_end = np.cross(end.T, normal.T).T
    for i in range(0, directions.shape[1], chunk_size):
        u = directions[:, i : i + chunk_size]
        # sine of the distance to each great circle (shape (E, M))
        sine = normal.T @ u
        # the direction's projection onto a great circle lies within its
        # arc if on the inner side of both of the arc's ends
        within = (after_start.T @ u >= 0) & (before_end.T @ u >= 0)
        to_circle = np.arcsin(np.clip(np.abs(sine), 0, 1))
        cos_start, cos_end = start.T @ u, end.T @ u
        to_ends = np.arccos(np.clip(np.maximum(cos_start, cos_end), -1, 1))
        candidates = np.where(within, to_circle, to_ends)
        arc = np.argmin(candidates, axis=0)
        column = np.arange(u.shape[1])
        distance[i : i + chunk_size] = candidates[arc, column]
        # nearest point: the projection onto the great circle within the
        # arc, or otherwise the nearer end of the arc
        projection = u - normal[:, arc] * sine[arc, column]
        projection /= np.linalg.norm(projection, axis=0)
        end_point = np.where(
            cos_start[arc, column] >= cos_end[arc, column], start[:, arc], end[:, arc]
        )
        nearest[:, i : i + chunk_size] = np.where(
            within[arc, column], projection, end_point
        )
    return distance, nearest


def _get_angular_distance_to_arcs(
    directions: npt.NDArray[np.float64],
    arcs: tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]],
    chunk_size: int = 1024,
) -> npt.NDArray[np.float64]:
    """
    Computes the minimum angular distance from each of a set of directions
    to a set of great circle arcs (see `_get_nearest_arc_points`).

    Args:
        directions (numpy.typing.NDArray[numpy.float64]): The unit vectors (shape (3, N)).
        arcs (tuple[numpy.typing.NDArray[numpy.float64], numpy.typing.NDArray[numpy.float64]]):
            The unit vectors (shape (3, E)) at the start and end of each arc.
        chunk_size (int): The number of directions processed at once, to limit memory.

    Returns:
        numpy.typing.NDArray[numpy.float64]: the angular distances (radians, shape (N,))
    """
    return _get_nearest_arc_points(directions, arcs, chunk_size)[0]


def _get_point_coordinates(point: Any) -> tuple[float, float, float]:
    """
    Gets the geodetic coordinates of a point: either a TAT-C `Point` (or a
    subclass, such as `GroundStation`) or a shapely `Point`, whose x, y, and
    (optional) z coordinates are its longitude (degrees), latitude
    (degrees), and elevation (meters) in the WGS 84 coordinate system.

    Args:
        point (tatc.schemas.Point | shapely.geometry.Point): The point.

    Returns:
        tuple[float, float, float]: the longitude (degrees), latitude
            (degrees), and elevation (meters)
    """
    if isinstance(point, Point):
        if point.is_empty:
            raise ValueError("Point is empty.")
        longitude, latitude = point.x, point.y
        elevation = point.z if point.has_z else 0.0
        if not -90 <= latitude <= 90:
            raise ValueError(f"Point latitude {latitude} is not within [-90, 90].")
        return float(longitude), float(latitude), float(elevation)
    return float(point.longitude), float(point.latitude), float(point.elevation)

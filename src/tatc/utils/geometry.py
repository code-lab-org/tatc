"""
Geometry utility functions.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""
from __future__ import annotations

import geopandas as gpd
import numpy as np
from shapely import make_valid
from shapely.geometry import (
    GeometryCollection,
    LineString,
    MultiPolygon,
    Point,
    Polygon,
)
from shapely.ops import split


def project_polygon_to_elevation(
    polygon: Polygon | MultiPolygon, elevation: float
) -> Polygon | MultiPolygon:
    """
    Projects a polygon to a specified elevation (z-coordinate).

    Args:
        polygon (shapely.geometry.Polygon | shapely.geometry.MultiPolygon): The polygon to project.
        elevation (float): The elevation (meters) above the WGS 84 geoid.

    Returns:
        shapely.geometry.Polygon | shapely.geometry.MultiPolygon: The projected polygon.
    """
    if isinstance(polygon, Polygon):
        return Polygon(
            [(p[0], p[1], elevation) for p in polygon.exterior.coords],
            [[(p[0], p[1], elevation) for p in i.coords] for i in polygon.interiors],
        )
    return MultiPolygon(
        [project_polygon_to_elevation(g, elevation) for g in polygon.geoms]
    )


def _wrap_polygon_over_north_pole(
    polygon: Polygon | MultiPolygon,
) -> Polygon | MultiPolygon:
    """
    Wraps polygon coordinates over the North pole. Due to buffering and projection,
    sometimes latitudes exceed 90 degrees. This method wraps them to the correct
    latitude between -90 and 90 degrees and adjusts the longitude by 180 degrees.
    This method requires a polygon above 90 degrees latitude to be only on one
    side of the prime meridian.

    Note: this method only changes coordinates: it does not create a MultiPolygon.

    Args:
       polygon (shapely.geometry.Polygon | shapely.geometry.MultiPolygon): The polygon to wrap.

    Returns:
       shapely.geometry.Polygon | shapely.geometry.MultiPolygon: The wrapped polygon.
    """
    if isinstance(polygon, Polygon):
        if all(c[1] <= 90 for c in polygon.exterior.coords):
            # no wrapping necessary
            return polygon
        # map latitudes from [90, 180) to [90, -90), adjusting longitude by 180 degrees
        lat_shift = 180 if all(c[0] <= 0 for c in polygon.exterior.coords) else -180
        pgon = Polygon(
            [
                [
                    c[0] + lat_shift if c[1] >= 90 else c[0],
                    180 - c[1] if c[1] >= 90 else c[1],
                ]
                for c in polygon.exterior.coords
            ],
            [
                [
                    [
                        c[0] + lat_shift if c[1] >= 90 else c[0],
                        180 - c[1] if c[1] >= 90 else c[1],
                    ]
                    for c in i.coords
                ]
                for i in polygon.interiors
            ],
        )
        # give up and return original polygon if invalid
        if not pgon.is_valid:
            return polygon
        return pgon
    if isinstance(polygon, MultiPolygon):
        # recursive call for each polygon
        polygons = [_wrap_polygon_over_north_pole(p) for p in polygon.geoms]
        return MultiPolygon(
            [
                g
                for p in polygons
                for g in (p.geoms if isinstance(p, MultiPolygon) else [p])
            ]
        )
    raise ValueError("Unknown geometry: " + str(type(polygon)))


def _split_polygon_north_pole(
    polygon: Polygon | MultiPolygon,
) -> Polygon | MultiPolygon:
    """
    Splits a Polygon into a MultiPolygon if it crosses north pole.

    Args:
       polygon (shapely.geometry.Polygon | shapely.geometry.MultiPolygon): The polygon to split.

    Returns:
       shapely.geometry.Polygon | shapely.geometry.MultiPolygon: The split polygon.
    """
    if isinstance(polygon, Polygon):
        if all(c[1] <= 90 for c in polygon.exterior.coords):
            # no splitting necessary
            return polygon
        # split polygon along north pole
        parts = split(polygon, LineString([(-360, 90), (360, 90)]))
        # check and split part over prime meridian if necessary
        for part in parts.geoms:
            if part.crosses(LineString([(0, 90), (0, 180)])):
                parts = GeometryCollection(
                    [g for g in parts.geoms if g != part]
                    + [g for g in split(part, LineString([(0, 90), (0, 180)])).geoms]
                )
        # convert to a multi polygon
        if isinstance(parts, GeometryCollection):
            parts = _convert_collection_to_polygon(parts)
        # return polygon with components wrapped over north pole
        return _wrap_polygon_over_north_pole(parts)
    if isinstance(polygon, MultiPolygon):
        # recursive call for each polygon
        pgons = [_split_polygon_north_pole(p) for p in polygon.geoms]
        return MultiPolygon(
            [
                g
                for p in pgons
                for g in (p.geoms if isinstance(p, MultiPolygon) else [p])
            ]
        )
    raise ValueError("Unknown geometry: " + str(type(polygon)))


def _wrap_polygon_over_south_pole(
    polygon: Polygon | MultiPolygon,
) -> Polygon | MultiPolygon:
    """
    Wraps polygon coordinates over the South pole. Due to buffering and projection,
    sometimes latitudes exceed -90 degrees. This method wraps them to the correct
    latitude between -90 and 90 degrees and adjusts the longitude by 180 degrees.
    This method requires a polygon above 90 degrees latitude to be only on one
    side of the prime meridian.

    Note: this method only changes coordinates: it does not create a MultiPolygon.

    Args:
       polygon (shapely.geometry.Polygon | shapely.geometry.MultiPolygon): The polygon to wrap.

    Returns:
       shapely.geometry.Polygon | shapely.geometry.MultiPolygon: The wrapped polygon.
    """
    if isinstance(polygon, Polygon):
        if all(c[1] >= -90 for c in polygon.exterior.coords):
            # no splitting necessary
            return polygon
        # map latitudes from [-90, -180) to [-90, 90), adjusting longitude by 180 degrees
        lat_shift = 180 if all(c[0] <= 0 for c in polygon.exterior.coords) else -180
        pgon = Polygon(
            [
                [
                    c[0] + lat_shift if c[1] <= -90 else c[0],
                    -180 - c[1] if c[1] <= -90 else c[1],
                ]
                for c in polygon.exterior.coords
            ],
            [
                [
                    [
                        (c[0] + lat_shift if c[1] <= -90 else c[0],),
                        -180 - c[1] if c[1] <= -90 else c[1],
                    ]
                    for c in i.coords
                ]
                for i in polygon.interiors
            ],
        )
        # give up and return original polygon if invalid
        if not pgon.is_valid:
            return polygon
        return pgon
    if isinstance(polygon, MultiPolygon):
        # recursive call for each polygon
        polygons = [_wrap_polygon_over_south_pole(p) for p in polygon.geoms]
        return MultiPolygon(
            [
                g
                for p in polygons
                for g in (p.geoms if isinstance(p, MultiPolygon) else [p])
            ]
        )
    raise ValueError("Unknown geometry: " + str(type(polygon)))


def _split_polygon_south_pole(
    polygon: Polygon | MultiPolygon,
) -> Polygon | MultiPolygon:
    """
    Splits a Polygon into a MultiPolygon if it crosses south pole.

    Args:
       polygon (shapely.geometry.Polygon | shapely.geometry.MultiPolygon): The polygon to split.

    Returns:
       shapely.geometry.Polygon | shapely.geometry.MultiPolygon: The split polygon.
    """
    if isinstance(polygon, Polygon):
        lat = np.array([c[1] for c in polygon.exterior.coords])
        if np.all(lat >= -90):
            return polygon
        # split polygon along south pole
        parts = split(polygon, LineString([(-360, -90), (360, -90)]))
        # check and split part over prime meridian if necessary
        for part in parts.geoms:
            if part.crosses(LineString([(0, -90), (0, -180)])):
                parts = GeometryCollection(
                    [g for g in parts.geoms if g != part]
                    + [g for g in split(part, LineString([(0, -90), (0, -180)])).geoms]
                )
        # convert to a multi polygon
        if isinstance(parts, GeometryCollection):
            parts = _convert_collection_to_polygon(parts)
        # return polygon with components wrapped over south pole
        return _wrap_polygon_over_south_pole(parts)
    if isinstance(polygon, MultiPolygon):
        # recursive call for each polygon
        pgons = [_split_polygon_south_pole(p) for p in polygon.geoms]
        return MultiPolygon(
            [
                g
                for p in pgons
                for g in (p.geoms if isinstance(p, MultiPolygon) else [p])
            ]
        )
    raise ValueError("Unknown geometry: " + str(type(polygon)))


def _wrap_polygon_over_antimeridian(
    polygon: Polygon | MultiPolygon,
) -> Polygon | MultiPolygon:
    """
    Wraps polygon coordinates over the antimeridian. Due to buffering and projection,
    sometimes longitudes exceed 180 degrees. This method wraps them to
    the correct longitude between -180 and 180 degrees.

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
            pgon = Polygon(
                [[c[0] + 360, c[1]] for c in polygon.exterior.coords],
                [[[c[0] + 360, c[1]] for c in i.coords] for i in polygon.interiors],
            )
        if all(c[0] >= 180 for c in polygon.exterior.coords):
            # map longitudes from [180, 540) to [-180, 180)
            pgon = Polygon(
                [[c[0] - 360, c[1]] for c in polygon.exterior.coords],
                [[[c[0] - 360, c[1]] for c in i.coords] for i in polygon.interiors],
            )
        # give up and return original polygon if invalid
        if not pgon.is_valid:
            return polygon
        return pgon
    if isinstance(polygon, MultiPolygon):
        # recursive call for each polygon
        pgons = [_wrap_polygon_over_antimeridian(p) for p in polygon.geoms]
        return MultiPolygon(
            [
                g
                for p in pgons
                for g in (p.geoms if isinstance(p, MultiPolygon) else [p])
            ]
        )
    raise ValueError("Unknown geometry: " + str(type(polygon)))


def _convert_collection_to_polygon(
    collection: GeometryCollection,
) -> Polygon | MultiPolygon:
    """
    Converts a GeometryCollection to a Polygon or MultiPolygon. Quick clipping
    can create dirty results with points or lines on boundaries. This method
    drops and lines or points from a GeometryCollection to return only the
    Polygon or MultiPolygon geometry.

    Args:
       polygon (shapely.geometry.Polygon | shapely.geometry.MultiPolygon): The polygon to convert.

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
    Splits a Polygon into a MultiPolygon if it crosses the anti-meridian after
    wrapping its coordinates using `wrap_coordinates_antimeridian`. Note: this
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
            # extract and sort coords by longitude
            coords = polygon.exterior.coords[0:-1]
            coords.sort(key=lambda r: r[0])
            # determine if contains north or south pole based on sign of mean latitude
            n_s = 1 if np.array(coords)[:, 1].mean() > 0 else -1
            # interpolate latitude at antimeridian
            lat = np.interp(
                180, [coords[-1][0], coords[0][0] + 180], [coords[-1][1], coords[0][1]]
            )
            # reconstruct polygon (ccw) with added coords on antimeridian
            # TODO potential problem if provided polygon has z-dimension
            pgon = Polygon(
                [(-180, 90 * n_s), (-180, lat)]
                + coords
                + [(180, lat), (180, 90 * n_s), (-180, 90 * n_s)],
                polygon.interiors,
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
        pgons = [_split_polygon_antimeridian(p) for p in polygon.geoms]
        return MultiPolygon(
            [
                g
                for p in pgons
                for g in (p.geoms if isinstance(p, MultiPolygon) else [p])
            ]
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

    Args:
        polygon (shapely.geometry.Polygon | shapely.geometry.MultiPolygon): The polygon to split.

    Returns:
        shapely.geometry.Polygon | shapely.geometry.MultiPolygon: The split polygon.
    """
    polygon = _split_polygon_north_pole(
        _split_polygon_south_pole(_split_polygon_antimeridian(polygon))
    )
    # invalid polygons can arise from narrow sensor geometries in polar regions
    if not polygon.is_valid:
        # try to fix geometry
        polygon = make_valid(polygon)
        if isinstance(polygon, GeometryCollection):
            polygon = _convert_collection_to_polygon(polygon)
    return polygon


def normalize_geometry(
    geometry: Polygon | MultiPolygon | gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    """
    Normalize geometry to a GeoDataFrame with antimeridian wrapping.

    Args:
        geometry (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame): The geometry to normalize.

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

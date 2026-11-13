"""
Methods to collect the observed ground footprints of airborne instruments
along a flight trajectory.

The methods take a generic trajectory (a data frame of aircraft states) and
do not depend on how it was produced (a flight plan, a flight recording, or a
simulation).

@author: Naeun Kim <hkim439@asu.edu>
"""

from __future__ import annotations

from typing import Literal

import geopandas as gpd
import numpy as np
import numpy.typing as npt
import pandas as pd
import shapely
from shapely.geometry import Polygon

from ..utils.ellipsoid import (
    compute_ellipsoid_intersection,
    geodetic_to_rectangular,
    rectangular_to_geodetic,
)
from ..utils.geometry import geodesic_destination

# trajectory columns that every footprint method requires
_REQUIRED_COLUMNS = (
    "time",
    "latitude",
    "longitude",
    "altitude",
    "track",
    "ground_speed",
)

# beam and pitch angles (degrees) beyond which the flat-surface formulas are
# not valid (the footprint is flagged as not valid)
_MAX_FLAT_BEAM_ANGLE = 80.0

# samples further apart than this multiple of the observation interval are
# not continuous (the footprint is extrapolated along the track instead)
_CONTINUITY_FACTOR = 1.5


def _get_empty_airborne_track() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for airborne ground track results.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "time": pd.Series([], dtype="datetime64[ns, utc]"),
        "aircraft": pd.Series([], dtype="str"),
        "instrument": pd.Series([], dtype="str"),
        "roll": pd.Series([], dtype="float"),
        "pitch": pd.Series([], dtype="float"),
        "valid_obs": pd.Series([], dtype="bool"),
        "geometry": pd.Series([], dtype="object"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def _offset_position(
    longitude: npt.ArrayLike,
    latitude: npt.ArrayLike,
    east: npt.ArrayLike,
    north: npt.ArrayLike,
) -> tuple[npt.NDArray, npt.NDArray]:
    """
    Moves positions by east and north offsets along WGS 84 geodesics.

    Args:
        longitude (numpy.typing.ArrayLike): The longitudes (degrees).
        latitude (numpy.typing.ArrayLike): The latitudes (degrees).
        east (numpy.typing.ArrayLike): The eastward offsets (meters).
        north (numpy.typing.ArrayLike): The northward offsets (meters).

    Returns:
        tuple[numpy.typing.NDArray, numpy.typing.NDArray]: the longitudes and
            latitudes (degrees) of the moved positions
    """
    moved_longitude, moved_latitude = geodesic_destination(
        np.asarray(longitude, dtype=float),
        np.asarray(latitude, dtype=float),
        np.degrees(np.arctan2(east, north)),
        np.hypot(east, north),
    )
    return np.asarray(moved_longitude), np.asarray(moved_latitude)


def _get_flat_scan_offsets(
    altitude: npt.NDArray,
    roll: npt.NDArray,
    pitch: npt.NDArray,
    half_angle: float,
) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray, npt.NDArray]:
    """
    Gets the ground offsets (meters) of the two scan edges of a cross-track
    scanner over a flat surface, relative to the aircraft's heading.

    The body frame has its x axis forward, y axis to the right wing, and z
    axis down. A beam at scan angle `a` (positive to the right) reaches the
    ground ahead by `altitude * tan(pitch)` and to the right by
    `altitude * tan(a - roll) / cos(pitch)`, so a positive roll (right wing
    down) tilts the scan to the left, the outside of a right turn.

    Args:
        altitude (numpy.typing.NDArray): The heights (meters) above the surface.
        roll (numpy.typing.NDArray): The roll angles (degrees, right wing down
            positive).
        pitch (numpy.typing.NDArray): The pitch angles (degrees, nose up
            positive).
        half_angle (float): The cross-track half angle (degrees) of the scan.

    Returns:
        tuple[
            numpy.typing.NDArray,
            numpy.typing.NDArray,
            numpy.typing.NDArray,
            numpy.typing.NDArray,
        ]:
            the forward offsets, the rightward offsets of the left and right
            edges (meters), and whether the flat-surface formulas are valid
    """
    phi, theta, alpha = np.radians(roll), np.radians(pitch), np.radians(half_angle)
    limit = np.radians(_MAX_FLAT_BEAM_ANGLE)
    beam_left, beam_right = -alpha - phi, alpha - phi
    valid = (
        (np.abs(beam_left) < limit)
        & (np.abs(beam_right) < limit)
        & (np.abs(theta) < limit)
        & (altitude > 0)
    )
    forward = altitude * np.tan(np.clip(theta, -limit, limit))
    left = altitude * np.tan(np.clip(beam_left, -limit, limit)) / np.cos(theta)
    right = altitude * np.tan(np.clip(beam_right, -limit, limit)) / np.cos(theta)
    return forward, left, right, valid


def _get_flat_edge_points(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    longitude: npt.NDArray,
    latitude: npt.NDArray,
    heading: npt.NDArray,
    forward: npt.NDArray,
    left: npt.NDArray,
    right: npt.NDArray,
) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray, npt.NDArray]:
    """
    Gets the longitudes and latitudes (degrees) of the left and right scan
    edges from their offsets relative to the heading (see
    `_get_flat_scan_offsets`).

    Returns:
        tuple[
            numpy.typing.NDArray,
            numpy.typing.NDArray,
            numpy.typing.NDArray,
            numpy.typing.NDArray,
        ]:
            the longitudes and latitudes of the left edges, then of the right edges
    """
    psi = np.radians(heading)
    # unit vectors (east, north) ahead of and to the right of the heading
    ahead_east, ahead_north = np.sin(psi), np.cos(psi)
    right_east, right_north = np.cos(psi), -np.sin(psi)
    left_lon, left_lat = _offset_position(
        longitude,
        latitude,
        forward * ahead_east + left * right_east,
        forward * ahead_north + left * right_north,
    )
    right_lon, right_lat = _offset_position(
        longitude,
        latitude,
        forward * ahead_east + right * right_east,
        forward * ahead_north + right * right_north,
    )
    return left_lon, left_lat, right_lon, right_lat


def _get_body_to_ned(
    roll: npt.NDArray, pitch: npt.NDArray, yaw: npt.NDArray
) -> npt.NDArray:
    """
    Gets the rotation matrices from the aircraft body frame (x forward, y
    right wing, z down) to the local north-east-down frame, for yaw (heading),
    pitch (nose up positive), and roll (right wing down positive) applied in
    that order about the body axes (the aerospace 3-2-1 sequence).

    Args:
        roll (numpy.typing.NDArray): The roll angles (degrees).
        pitch (numpy.typing.NDArray): The pitch angles (degrees).
        yaw (numpy.typing.NDArray): The yaw angles (degrees), the heading
            clockwise from north.

    Returns:
        numpy.typing.NDArray: the matrices (shape (N, 3, 3))
    """
    r, p, y = np.radians(roll), np.radians(pitch), np.radians(yaw)
    cr, sr = np.cos(r), np.sin(r)
    cp, sp = np.cos(p), np.sin(p)
    cy, sy = np.cos(y), np.sin(y)
    matrix = np.empty((len(r), 3, 3))
    matrix[:, 0, 0] = cp * cy
    matrix[:, 0, 1] = sr * sp * cy - cr * sy
    matrix[:, 0, 2] = cr * sp * cy + sr * sy
    matrix[:, 1, 0] = cp * sy
    matrix[:, 1, 1] = sr * sp * sy + cr * cy
    matrix[:, 1, 2] = cr * sp * sy - sr * cy
    matrix[:, 2, 0] = -sp
    matrix[:, 2, 1] = sr * cp
    matrix[:, 2, 2] = cr * cp
    return matrix


def _get_ned_to_earth_fixed(
    longitude: npt.NDArray, latitude: npt.NDArray
) -> npt.NDArray:
    """
    Gets the rotation matrices from the local north-east-down frame to the
    Earth-fixed frame, whose columns are the local north, east, and down
    axes (for geodetic latitudes and longitudes).

    Args:
        longitude (numpy.typing.NDArray): The longitudes (degrees).
        latitude (numpy.typing.NDArray): The latitudes (degrees).

    Returns:
        numpy.typing.NDArray: the matrices (shape (N, 3, 3))
    """
    lon, lat = np.radians(longitude), np.radians(latitude)
    north = np.stack(
        [-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)], -1
    )
    east = np.stack([-np.sin(lon), np.cos(lon), np.zeros_like(lon)], -1)
    down = np.stack(
        [-np.cos(lat) * np.cos(lon), -np.cos(lat) * np.sin(lon), -np.sin(lat)], -1
    )
    return np.stack([north, east, down], -1)


def compute_airborne_edge_points(  # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
    latitude: npt.ArrayLike,
    longitude: npt.ArrayLike,
    altitude: npt.ArrayLike,
    roll: npt.ArrayLike,
    pitch: npt.ArrayLike,
    heading: npt.ArrayLike,
    cross_track_half_angle_deg: float,
    elevation: float = 0,
) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray]:
    """
    Compute the points on the ground where the two edges of a cross-track
    scan of an airborne instrument meet the WGS 84 ellipsoid (extended by an
    elevation). The instrument looks down the aircraft's body z axis and
    scans across its wings, so its edges are the rays at scan angles of minus
    (left) and plus (right) the half angle. The rays are rotated from the
    body frame by the aircraft's attitude (heading, then pitch, then roll)
    and intersected with the ellipsoid in 3D, so large roll and pitch angles
    and the Earth's curvature are modeled. The intersection uses the
    ellipsoid utilities of TAT-C (`tatc.utils.ellipsoid`) and does not call
    SPICE; the points agree with the SPICE surface intersection to within
    nanometers.

    Args:
        latitude (numpy.typing.ArrayLike): The aircraft geodetic latitudes (degrees).
        longitude (numpy.typing.ArrayLike): The aircraft longitudes (degrees).
        altitude (numpy.typing.ArrayLike): The aircraft altitudes (meters)
            above the WGS 84 ellipsoid.
        roll (numpy.typing.ArrayLike): The roll angles (degrees, right wing
            down positive).
        pitch (numpy.typing.ArrayLike): The pitch angles (degrees, nose up
            positive).
        heading (numpy.typing.ArrayLike): The headings (degrees clockwise
            from north) of the aircraft's nose.
        cross_track_half_angle_deg (float): The cross-track half angle
            (degrees) of the scan, from nadir to one edge.
        elevation (float): The elevation (meters) of the surface above the
            WGS 84 ellipsoid, in the same reference as the altitudes. The
            surface at an elevation is approximated by the ellipsoid with
            both semi-axes extended by the elevation (see
            `tatc.utils.ellipsoid.compute_ellipsoid_intersection`); its
            geodetic height differs from the elevation by less than 1 cm up
            to 5 km and about 3 cm at 20 km (largest near 45 degrees
            latitude, exact at the equator and the poles).

    Returns:
        tuple[numpy.typing.NDArray, numpy.typing.NDArray, numpy.typing.NDArray]:
            the latitudes and longitudes (degrees) of the edge points and
            whether each edge ray meets the surface (shape (N, 2), left edge
            first; rays that miss, or from an aircraft at or below the
            surface, are `NaN` and not valid)
    """
    lat, lon, alt, roll, pitch, heading = (
        np.atleast_1d(np.asarray(a, dtype=float))
        for a in (latitude, longitude, altitude, roll, pitch, heading)
    )
    scan = np.radians([-cross_track_half_angle_deg, cross_track_half_angle_deg])
    # the edge rays in the body frame (shape (2, 3))
    rays = np.stack([np.zeros(2), np.sin(scan), np.cos(scan)], -1)
    body_to_earth = np.einsum(
        "nij,njk->nik",
        _get_ned_to_earth_fixed(lon, lat),
        _get_body_to_ned(roll, pitch, heading),
    )
    directions = np.einsum("nij,kj->nki", body_to_earth, rays)
    points, found = compute_ellipsoid_intersection(
        np.repeat(geodetic_to_rectangular(lon, lat, alt), 2, axis=1),
        directions.reshape(-1, 3).T,
        elevation=elevation,
    )
    # an aircraft at or below the surface cannot observe it (a ray from a
    # point on the surface has an unstable first intersection)
    found &= np.repeat(alt > elevation, 2)
    edge_lon, edge_lat, _ = rectangular_to_geodetic(points)
    shape = (len(lat), 2)
    return (
        np.where(found, edge_lat, np.nan).reshape(shape),
        np.where(found, edge_lon, np.nan).reshape(shape),
        found.reshape(shape),
    )


def _prepare_trajectory(trajectory: pd.DataFrame) -> pd.DataFrame:
    """
    Prepares a trajectory: a copy sorted by (UTC) time with a continuous index.
    """
    prepared = trajectory.copy()
    prepared["time"] = pd.to_datetime(prepared["time"], utc=True)
    return prepared.sort_values("time").reset_index(drop=True)


def _get_next_sample_index(
    trajectory: pd.DataFrame, observation_interval_s: float
) -> tuple[npt.NDArray, npt.NDArray]:
    """
    Gets the index of the next sample of each sample (the last sample is its
    own next sample) and whether the next sample is continuous (no more
    than 1.5 observation intervals later).
    """
    n = len(trajectory)
    seconds = (
        (trajectory["time"] - trajectory["time"].iloc[0]).dt.total_seconds().to_numpy()
    )
    following = np.minimum(np.arange(n) + 1, n - 1)
    continuous = (np.arange(n) < n - 1) & (
        (seconds[following] - seconds) <= _CONTINUITY_FACTOR * observation_interval_s
    )
    return following, continuous


def _get_attitude(
    trajectory: pd.DataFrame, use_attitude: bool
) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray]:
    """
    Gets the roll, pitch, and heading (degrees) of each sample: zero roll and
    pitch if the attitude is not used or not given, and the track as the
    heading if the heading is not given (in still air an aircraft points
    along its track; with wind its heading differs from its track).
    """
    n = len(trajectory)
    roll, pitch = (
        (
            trajectory[name].to_numpy(dtype=float)
            if use_attitude and name in trajectory
            else np.zeros(n)
        )
        for name in ("roll", "pitch")
    )
    heading = (
        trajectory["heading"].to_numpy(dtype=float)
        if "heading" in trajectory
        else trajectory["track"].to_numpy(dtype=float)
    )
    return roll, pitch, heading


def _get_hull_footprints(
    corners: npt.NDArray, include: npt.NDArray
) -> tuple[npt.NDArray, npt.NDArray]:
    """
    Gets the convex hulls of the corner points of each sample as polygons.
    The hull (rather than the polygon through the corners in order) avoids a
    twisted footprint where the swath edges cross.

    Args:
        corners (numpy.typing.NDArray): The corner longitudes and latitudes
            of each sample (shape (N, 4, 2)).
        include (numpy.typing.NDArray): Whether to build each sample's hull
            (samples with non-finite corners are never built).

    Returns:
        tuple[numpy.typing.NDArray, numpy.typing.NDArray]: the polygons (an
            empty polygon where not built) and whether each has an area
    """
    geometries = np.empty(len(corners), dtype=object)
    geometries[:] = [Polygon() for _ in range(len(corners))]
    has_area = np.zeros(len(corners), dtype=bool)
    rows = np.flatnonzero(include & np.isfinite(corners).all(axis=(1, 2)))
    if len(rows) > 0:
        hulls = shapely.convex_hull(
            shapely.multipoints(
                corners[rows].reshape(-1, 2),
                indices=np.repeat(np.arange(len(rows)), corners.shape[1]),
            )
        )
        # a collinear or repeated set of corners has a point or line hull
        is_polygon = shapely.get_type_id(hulls) == shapely.GeometryType.POLYGON
        geometries[rows[is_polygon]] = hulls[is_polygon]
        has_area[rows[is_polygon]] = True
    return geometries, has_area


def _collect_level_footprints(  # pylint: disable=too-many-locals
    trajectory: pd.DataFrame,
    half_angle: float,
    observation_interval_s: float,
    elevation: float = 0,
) -> tuple[npt.NDArray, npt.NDArray]:
    """
    Collects level-flight footprints: for each sample, the rectangle that is
    centered on the aircraft, oriented along its track, as long as the
    distance flown during an observation interval (at least 1 meter) and as
    wide as the swath of a nadir-pointing scan over a flat level at an
    elevation above the flat ground. A sample at or below that level
    (altitude no greater than the elevation) is not valid and has an empty
    footprint.
    """
    longitude = trajectory["longitude"].to_numpy(dtype=float)
    latitude = trajectory["latitude"].to_numpy(dtype=float)
    track = np.radians(trajectory["track"].to_numpy(dtype=float))
    along_half = (
        np.maximum(
            trajectory["ground_speed"].to_numpy(dtype=float) * observation_interval_s,
            1.0,
        )
        / 2
    )
    # the height above the observed level (the aircraft altitude and the
    # elevation are both measured from the flat ground)
    height = trajectory["altitude"].to_numpy(dtype=float) - elevation
    valid = height > 0
    cross_half = height * np.tan(np.radians(half_angle))
    corners = []
    for cross_sign, along_sign in ((1, 1), (-1, 1), (-1, -1), (1, -1)):
        cross, along = cross_sign * cross_half, along_sign * along_half
        corners.append(
            np.stack(
                _offset_position(
                    longitude,
                    latitude,
                    along * np.sin(track) + cross * np.cos(track),
                    along * np.cos(track) - cross * np.sin(track),
                ),
                -1,
            )
        )
    footprints = shapely.polygons(np.stack(corners, axis=1))
    footprints[~valid] = Polygon()
    return footprints, valid


def _collect_flat_strips(  # pylint: disable=too-many-locals
    trajectory: pd.DataFrame,
    half_angle: float,
    observation_interval_s: float,
    attitude: tuple[npt.NDArray, npt.NDArray, npt.NDArray],
    elevation: float = 0,
) -> tuple[npt.NDArray, npt.NDArray]:
    """
    Collects flat-surface footprints with attitude: for each sample, the
    strip between the scan edges at the sample and at the next sample (or,
    where the next sample is not continuous, at the position reached after
    one observation interval along the track with the same attitude), on a
    flat level at an elevation above the flat ground. The strips follow each
    other without gaps on the outside of a turn. A sample at or below that
    level (altitude no greater than the elevation) is not valid and has an
    empty footprint.
    """
    roll, pitch, heading = attitude
    longitude = trajectory["longitude"].to_numpy(dtype=float)
    latitude = trajectory["latitude"].to_numpy(dtype=float)
    track = np.radians(trajectory["track"].to_numpy(dtype=float))
    distance = trajectory["ground_speed"].to_numpy(dtype=float) * observation_interval_s
    # the height above the observed level (the aircraft altitude and the
    # elevation are both measured from the flat ground)
    height = trajectory["altitude"].to_numpy(dtype=float) - elevation
    forward, left, right, valid = _get_flat_scan_offsets(
        height, roll, pitch, half_angle
    )
    edges = _get_flat_edge_points(longitude, latitude, heading, forward, left, right)
    # the edges after one observation interval, as extrapolated along the track
    ahead_lon, ahead_lat = _offset_position(
        longitude, latitude, distance * np.sin(track), distance * np.cos(track)
    )
    ahead = _get_flat_edge_points(ahead_lon, ahead_lat, heading, forward, left, right)
    following, continuous = _get_next_sample_index(trajectory, observation_interval_s)
    next_edges = [
        np.where(continuous, edge[following], extrapolated)
        for edge, extrapolated in zip(edges, ahead)
    ]
    corners = np.stack(
        [
            np.stack([edges[0], edges[1]], -1),
            np.stack([edges[2], edges[3]], -1),
            np.stack([next_edges[2], next_edges[3]], -1),
            np.stack([next_edges[0], next_edges[1]], -1),
        ],
        axis=1,
    )
    above = height > 0
    geometries, has_area = _get_hull_footprints(corners, above)
    # without an area (as when stationary), use the triangle of the edges
    # and the aircraft position
    for i in np.flatnonzero(above & ~has_area & np.isfinite(corners).all(axis=(1, 2))):
        geometries[i] = Polygon(
            [corners[i, 0], corners[i, 1], (longitude[i], latitude[i])]
        )
    return geometries, valid & has_area


def _collect_ellipsoid_strips(  # pylint: disable=too-many-locals
    trajectory: pd.DataFrame,
    half_angle: float,
    observation_interval_s: float,
    attitude: tuple[npt.NDArray, npt.NDArray, npt.NDArray],
    elevation: float = 0,
) -> tuple[npt.NDArray, npt.NDArray]:
    """
    Collects ellipsoid footprints (see `compute_airborne_edge_points`): for
    each sample, the strip between the scan edges at the sample and at the
    next sample (or, where the next sample is not continuous, at the
    position reached after one observation interval along the track with the
    same attitude), where the scan rays meet the surface at an elevation
    above the WGS 84 ellipsoid. The aircraft is not moved; only the surface
    is raised. A sample is not valid if any edge ray misses that surface or
    the aircraft is at or below it.
    """
    roll, pitch, heading = attitude
    longitude = trajectory["longitude"].to_numpy(dtype=float)
    latitude = trajectory["latitude"].to_numpy(dtype=float)
    altitude = trajectory["altitude"].to_numpy(dtype=float)
    lat, lon, found = compute_airborne_edge_points(
        latitude, longitude, altitude, roll, pitch, heading, half_angle, elevation
    )
    # the edges after one observation interval, as extrapolated along the track
    ahead_lon, ahead_lat = geodesic_destination(
        longitude,
        latitude,
        trajectory["track"].to_numpy(dtype=float),
        trajectory["ground_speed"].to_numpy(dtype=float) * observation_interval_s,
    )
    ahead_lat, ahead_lon, ahead_found = compute_airborne_edge_points(
        ahead_lat, ahead_lon, altitude, roll, pitch, heading, half_angle, elevation
    )
    following, continuous = _get_next_sample_index(trajectory, observation_interval_s)
    next_lat = np.where(continuous[:, None], lat[following], ahead_lat)
    next_lon = np.where(continuous[:, None], lon[following], ahead_lon)
    valid = found.all(axis=1) & np.where(
        continuous, found[following].all(axis=1), ahead_found.all(axis=1)
    )
    corners = np.stack(
        [
            np.stack([lon[:, 0], lat[:, 0]], -1),
            np.stack([lon[:, 1], lat[:, 1]], -1),
            np.stack([next_lon[:, 1], next_lat[:, 1]], -1),
            np.stack([next_lon[:, 0], next_lat[:, 0]], -1),
        ],
        axis=1,
    )
    geometries, has_area = _get_hull_footprints(corners, valid)
    return geometries, valid & has_area


def collect_airborne_ground_track(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    trajectory: pd.DataFrame,
    cross_track_half_angle_deg: float,
    observation_interval_s: float,
    aircraft_name: str = "aircraft",
    instrument_name: str = "instrument",
    use_attitude: bool = False,
    geometry: Literal["flat", "ellipsoid", "spice"] = "flat",
    elevation: float = 0,
) -> gpd.GeoDataFrame:
    """
    Collect the observed ground footprints of a cross-track scanning
    instrument that looks down from an aircraft, along its flight trajectory.
    The altitude is the height above the ground (flat) or above the WGS 84
    ellipsoid (ellipsoid), and the footprints are not clipped to a region.

    By default the footprints are on the ground (an elevation of zero). With
    an elevation, they are where the scan meets a level at that height
    instead, for example the altitude of a layer of a storm; the aircraft
    stays at its altitude and only the observed level is raised. In the
    `"flat"` geometry the elevation, like the altitude, is a height above the
    flat ground; in the `"ellipsoid"` geometry it is a height above the WGS
    84 ellipsoid, in the same reference as the altitude (see
    `compute_airborne_edge_points` for how that surface is approximated).
    The footprint is the geometric overlap of the instrument's field of view
    with that level: whether the instrument detects anything there (its
    sensitivity, or attenuation along the beam) is not modeled.

    The geometry sets how a footprint is computed:

    - `"flat"` without attitude: the rectangle centered on the aircraft,
      oriented along its track, as long as the distance flown during an
      observation interval and as wide as the swath of a nadir-pointing
      scan over a flat surface.
    - `"flat"` with attitude (`use_attitude`): the strip between the scan
      edges at a sample and at the next sample, from formulas for a flat
      surface in which the roll and pitch angles shift the swath. Samples
      with beams (or a pitch) beyond 80 degrees from nadir are not valid.
    - `"ellipsoid"`: the same strip with the scan edges found by rotating
      the scan rays with the aircraft attitude and intersecting them with
      the WGS 84 ellipsoid in 3D (see `compute_airborne_edge_points`), which
      models large roll and pitch angles and the Earth's curvature. Samples
      with an edge ray that misses the Earth are not valid and have an
      empty footprint.
    - `"spice"`: an alias of `"ellipsoid"` for backward compatibility with
      the earlier prototype of this method, which intersected the rays with
      the ellipsoid by calling SPICE. The same computation as `"ellipsoid"`
      is used, with the ellipsoid utilities of TAT-C
      (`tatc.utils.ellipsoid`) instead of SPICE (they agree to within
      nanometers).

    In every geometry, a sample at or below the observed level (an altitude
    no greater than the elevation; zero or less by default) is not valid and
    has an empty footprint, as the instrument looks down. Where the next sample
    is more than 1.5 observation intervals later (as between flight legs),
    the strip is extrapolated along the track for one observation interval
    instead of spanning the gap.

    Args:
        trajectory (pandas.DataFrame): The aircraft states, one row per
            sample, with columns `time` (timezone-aware or UTC), `latitude`
            and `longitude` (degrees), `altitude` (meters), `track` (degrees
            clockwise from north, the direction of motion over the ground)
            and `ground_speed` (meters per second); and optionally `roll`
            (degrees, right wing down positive), `pitch` (degrees, nose up
            positive), `heading` (degrees clockwise from north, the
            direction of the aircraft's nose, by default the track) and
            `segment_name` (copied to the results with attitude).
        cross_track_half_angle_deg (float): The cross-track half angle
            (degrees) of the scan, from nadir to one edge.
        observation_interval_s (float): The time (seconds) between
            observations.
        aircraft_name (str): The name of the aircraft.
        instrument_name (str): The name of the instrument.
        use_attitude (bool): `True` to apply the roll and pitch angles of
            the trajectory (zero if not given) to the footprints, otherwise
            the aircraft is level. Always `False` for the `"flat"` rectangle.
        geometry (str): The footprint geometry, `"flat"` or `"ellipsoid"`
            (`"spice"` is an alias of `"ellipsoid"`).
        elevation (float): The height (meters) of the observed level: above
            the flat ground for `"flat"`, above the WGS 84 ellipsoid (in the
            reference of the altitude) for `"ellipsoid"`. Zero by default
            (the ground).

    Returns:
        geopandas.GeoDataFrame: The data frame of collected airborne ground
            track results, with columns `time`, `aircraft`, `instrument`,
            `roll`, `pitch` (degrees; the input angles for the level
            rectangle, otherwise those applied), `heading` (degrees, except
            for the level rectangle), `valid_obs` and `geometry`.
    """
    if trajectory.empty:
        return _get_empty_airborne_track()
    missing = set(_REQUIRED_COLUMNS).difference(trajectory.columns)
    if missing:
        raise ValueError(f"Trajectory is missing required columns: {sorted(missing)}")
    if observation_interval_s <= 0:
        raise ValueError("observation_interval_s must be positive.")
    if not 0 < cross_track_half_angle_deg < 90:
        raise ValueError("cross_track_half_angle_deg must be between 0 and 90 degrees.")
    if geometry not in ("flat", "ellipsoid", "spice"):
        raise ValueError(
            "geometry must be 'flat' or 'ellipsoid' ('spice' is an alias of 'ellipsoid')."
        )
    if geometry == "spice":
        geometry = "ellipsoid"
    if not np.isfinite(elevation):
        raise ValueError("elevation must be a finite number of meters.")

    trajectory = _prepare_trajectory(trajectory)
    results = pd.DataFrame(
        {
            "time": trajectory["time"],
            "aircraft": aircraft_name,
            "instrument": instrument_name,
        }
    )
    if geometry == "flat" and not use_attitude:
        footprints, valid = _collect_level_footprints(
            trajectory, cross_track_half_angle_deg, observation_interval_s, elevation
        )
        for name in ("roll", "pitch"):
            results[name] = (
                trajectory[name].to_numpy(dtype=float)
                if name in trajectory
                else np.zeros(len(trajectory))
            )
        results["valid_obs"] = valid
    else:
        attitude = _get_attitude(trajectory, use_attitude)
        collect = (
            _collect_flat_strips if geometry == "flat" else _collect_ellipsoid_strips
        )
        footprints, valid = collect(
            trajectory,
            cross_track_half_angle_deg,
            observation_interval_s,
            attitude,
            elevation,
        )
        results["roll"], results["pitch"], results["heading"] = attitude
        results["valid_obs"] = valid
        if "segment_name" in trajectory:
            results["segment_name"] = trajectory["segment_name"].to_numpy()
    return gpd.GeoDataFrame(results, geometry=list(footprints), crs="EPSG:4326")

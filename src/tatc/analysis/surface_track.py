"""
Methods to collect the static tracks of instruments on fixed surface platforms.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import warnings
from datetime import datetime, timezone

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import MultiPolygon, Polygon
from skyfield.positionlib import Geocentric

from ..schemas import Instrument, SurfacePlatform
from .check import (
    _check_platforms,
    _check_times,
    _combine_results,
    _dissolve_footprints,
    _get_instrument_indices,
    _is_single,
)
from .ground_track import (
    _assemble_pixels,
    _assemble_track,
    _clip_track,
    _get_angle_columns,
    _get_empty_track,
    _get_mask_geometry,
)

_STATIC_TIME = datetime(2020, 1, 1, tzinfo=timezone.utc)
"""
Time of the track of a fixed platform analyzed without times, whose views
and footprints do not depend on time (see `SurfacePlatform.get_track`).
"""


def _get_empty_surface_track(
    angles: list[str] | tuple[str, ...] = (), timed: bool = False
) -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for surface track results.

    Args:
        angles (list[str] | tuple[str, ...]): The angle columns requested.
        timed (bool): `True` to include a `time` column.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    return _get_empty_track(
        "surface_platform", angles, timed=timed, extra={"elevation": "float"}
    )


def _get_target_elevation(platform: SurfacePlatform, elevation: float | None) -> float:
    """
    Gets the elevation at which to project a platform's footprints: the
    given elevation, or the platform's ground elevation if None, which must
    differ from the elevation of its instruments (below them, looking down,
    or above them, looking up).

    Args:
        platform (SurfacePlatform): The platform.
        elevation (float | None): The target elevation (meters) above the
            WGS 84 ellipsoid, or None for the platform's ground elevation.

    Returns:
        float: the target elevation (meters).
    """
    target = platform.elevation if elevation is None else elevation
    if target == platform.get_altitude():
        raise ValueError(
            f"target elevation ({target} m) must differ from the elevation of "
            f"the instruments of platform {platform.name!r} "
            f"({platform.get_altitude()} m)"
        )
    return target


def _clip_to_region(
    footprint: Polygon | MultiPolygon, region: Polygon | MultiPolygon, elevation: float
) -> Polygon | MultiPolygon:
    """
    Clips a footprint to a region (such as the region visible over a terrain
    mask), keeping only its polygonal parts.

    Args:
        footprint (shapely.geometry.Polygon | shapely.geometry.MultiPolygon): The footprint.
        region (shapely.geometry.Polygon | shapely.geometry.MultiPolygon): The region.
        elevation (float): The elevation (meters) of the footprint.

    Returns:
        shapely.geometry.Polygon | shapely.geometry.MultiPolygon: The clipped footprint.
    """
    clipped = shapely.intersection(shapely.force_2d(footprint), region)
    parts = [
        part
        for part in shapely.get_parts(shapely.get_parts(clipped))
        if isinstance(part, Polygon) and not part.is_empty
    ]
    if len(parts) == 0:
        return Polygon()
    return shapely.force_3d(
        parts[0] if len(parts) == 1 else MultiPolygon(parts), elevation
    )


def _get_surface_platform_altaz(
    track: Geocentric,
    instrument: Instrument,
    elevation: float,
    ground_elevation: float,
    empty: bool,
) -> dict[str, float]:
    """
    Gets the altitude/azimuth angles (degrees) of a platform's instrument
    seen from the center of its view at an elevation (see
    `collect_surface_track`).

    Args:
        track (skyfield.positionlib.Geocentric): The platform's static track.
        instrument (Instrument): The instrument.
        elevation (float): The target elevation (meters).
        ground_elevation (float): The elevation (meters) of the ground below
            the platform.
        empty (bool): `True` if the footprint is empty, without angles.

    Returns:
        dict[str, float]: the angles (`surface_platform_alt` and `surface_platform_az`)
    """
    if empty:
        return {"surface_platform_alt": np.nan, "surface_platform_az": np.nan}
    center = instrument.compute_footprint_center(
        track, elevation=elevation, ground_elevation=ground_elevation
    )
    alt, az, _ = (track - center.at(track.t)).altaz()
    return {
        "surface_platform_alt": float(np.ravel(alt.degrees)[0]),
        "surface_platform_az": float(np.ravel(az.degrees)[0]),
    }


def collect_surface_track(
    platforms: SurfacePlatform | list[SurfacePlatform],
    times: list[datetime] | None = None,
    instrument_index: int | None = 0,
    elevation: float | None = None,
    mask: Polygon | MultiPolygon | gpd.GeoDataFrame | gpd.GeoSeries | None = None,
    surface_platform_altaz: bool = False,
    solar_altaz: bool = False,
) -> gpd.GeoDataFrame:
    """
    Collect the static footprint of each instrument on fixed surface
    platforms (e.g. towers or tethered balloons), projected to a target
    elevation by computing the exact ray/WGS 84 ellipsoid intersection, as
    for a satellite's ground track (see `collect_ground_track`). Rays that
    miss a target below the instruments (looking above the horizon) are
    projected to the horizon. A target above the instruments (as for clouds
    observed by an instrument looking up, e.g. with a roll angle of 180
    degrees) is observed only above the horizon: views are clipped to the
    lines of sight tangent to the ground, which blocks the rest (a view
    entirely below the horizon has an empty footprint), and to the
    platform's terrain mask, if any (see
    `SurfacePlatform.compute_visible_region`). An empty footprint (as of a
    view entirely blocked by the ground or terrain) is not a valid
    observation (`valid_obs`), as for satellites and aircraft.

    A platform does not move in the Earth-fixed frame, so its footprints do
    not change with time. Without `times`, this returns one row per platform
    and instrument, without a time, and an instrument's time-varying
    requirements (such as a sunlit target, see
    `Instrument.is_valid_observation`) are not evaluated. With `times`, it
    returns one row per time, platform, and instrument, as
    `collect_ground_track` does for satellites, with the validity of each
    observation evaluated at its time (at the footprint center) and, if
    `solar_altaz`, the solar angles.

    Args:
        platforms (SurfacePlatform | list[SurfacePlatform]): The platform(s).
        times (list[datetime.datetime] | None): The times to sample, or
                `None` for a static footprint of each instrument.
        instrument_index (int | None): The index of the observing instrument
                on each platform, or `None` for every instrument of each
                platform.
        elevation (float | None): The elevation (meters) above the WGS 84
                ellipsoid of the observed target for which to compute
                footprints, below (looking down) or above (looking up) each
                platform's instruments. Defaults to each platform's ground
                elevation.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame | geopandas.GeoSeries | None):
                An optional mask, always interpreted in WGS84 (lon/lat)
                coordinates, to constrain results.
        surface_platform_altaz (bool): `True` to include the altitude/azimuth
                angles (`surface_platform_alt` and `surface_platform_az`) of
                the instruments seen from the center of their view at the
                target elevation (as `sat_altaz` does for a satellite),
                `NaN` for an empty footprint. Below the instruments, the
                altitude is the elevation angle of the view; above them, it
                is negative.
        solar_altaz (bool): `True` to include the solar altitude/azimuth
                angles (`solar_alt` and `solar_az`) and the glint angle
                (`glint_angle`, see `tatc.utils.compute_glint_angle`) at the
                footprint center, `NaN` for an empty footprint. Requires
                `times`.

    Returns:
        geopandas.GeoDataFrame: The data frame of collected surface track
            results, with columns `surface_platform`, `instrument`,
            `elevation`, `valid_obs`, and `geometry` (and the angles
            requested): without `times`, in the order of the platforms and
            their instruments; with `times`, also `time`, for a single
            platform and instrument (an integer `instrument_index`) sorted
            by time, otherwise those of each platform and instrument,
            concatenated and sorted by time.
    """
    return _collect_surface(
        platforms,
        times,
        instrument_index,
        elevation,
        mask,
        surface_platform_altaz,
        solar_altaz,
    )


def _collect_surface(
    platforms: SurfacePlatform | list[SurfacePlatform],
    times: list[datetime] | None,
    instrument_index: int | None,
    elevation: float | None,
    mask: Polygon | MultiPolygon | gpd.GeoDataFrame | gpd.GeoSeries | None,
    surface_platform_altaz: bool,
    solar_altaz: bool,
    pixels: bool = False,
) -> gpd.GeoDataFrame:
    """
    Collect the footprints (see `collect_surface_track`) or, if `pixels`,
    the pixels (see `tatc.analysis.footprints.collect_footprints`) of
    instruments on fixed surface platforms.
    """
    single = _is_single(platforms, instrument_index=instrument_index)
    platforms = _check_platforms(platforms, SurfacePlatform, "platforms")
    _check_times(times)
    if solar_altaz and times is None:
        raise ValueError(
            "solar_altaz requires times: the solar angles of a fixed "
            "platform's footprint change with time"
        )
    if mask is not None:
        mask = _get_mask_geometry(mask)
    platform_altaz = "surface_platform" if surface_platform_altaz else None
    angles = _get_angle_columns(platform_altaz, solar_altaz)
    timed = times is not None
    if timed and len(times) == 0:
        return _get_empty_surface_track(angles, timed)
    frames = []
    for platform in platforms:
        target = _get_target_elevation(platform, elevation)
        # the static track, from which the footprints are projected
        static = platform.get_track()
        # the track at each time (or at one time, without times)
        track = platform.get_track(times if timed else [_STATIC_TIME])
        # the region visible over the terrain mask, toward a target above
        region = (
            platform.compute_visible_region(target)
            if platform.terrain_mask is not None and target > platform.get_altitude()
            else None
        )
        for index in _get_instrument_indices(
            platform, instrument_index, "surface platform"
        ):
            instrument = platform.get_instrument(index)
            if pixels:
                frames.append(
                    _assemble_pixels(
                        "surface_platform",
                        platform.name,
                        instrument,
                        pd.to_datetime(times, utc=True) if timed else None,
                        track,
                        target,
                        mask,
                        platform_altaz,
                        solar_altaz,
                        ground_elevation=platform.elevation,
                        extra={"elevation": target},
                        region=region,
                    )
                )
                continue
            footprint = instrument.compute_footprint(
                static, elevation=target, ground_elevation=platform.elevation
            )[0]
            if region is not None:
                footprint = _clip_to_region(footprint, region, target)
            if timed:
                frames.append(
                    _assemble_track(
                        "surface_platform",
                        platform.name,
                        instrument,
                        pd.to_datetime(times, utc=True),
                        track,
                        [footprint] * len(times),
                        target,
                        mask,
                        platform_altaz,
                        solar_altaz,
                        ground_elevation=platform.elevation,
                        extra={"elevation": target},
                    )
                )
                continue
            record = {
                "surface_platform": platform.name,
                "instrument": instrument.name,
                "elevation": target,
                "valid_obs": not footprint.is_empty,
                "geometry": footprint,
            }
            if surface_platform_altaz:
                record.update(
                    _get_surface_platform_altaz(
                        static,
                        instrument,
                        target,
                        platform.elevation,
                        footprint.is_empty,
                    )
                )
            frames.append(
                _clip_track(gpd.GeoDataFrame([record], crs="EPSG:4326"), mask)
            )
    if timed:
        return _combine_results(
            frames, single, "time", lambda: _get_empty_surface_track(angles, timed)
        )
    frames = [frame for frame in frames if not frame.empty]
    if len(frames) == 0:
        return _get_empty_surface_track(angles)
    return pd.concat(frames, ignore_index=True)


def compute_surface_track(
    platforms: SurfacePlatform | list[SurfacePlatform],
    times: list[datetime] | None = None,
    instrument_index: int | None = 0,
    elevation: float | None = None,
    mask: Polygon | MultiPolygon | gpd.GeoDataFrame | gpd.GeoSeries | None = None,
) -> gpd.GeoDataFrame:
    """
    Compute the aggregated (dissolved) coverage of instruments on fixed
    surface platforms into a single composite geometry (e.g. a network
    coverage map): unlike `collect_surface_track` (one footprint polygon per
    platform and instrument, and time), this dissolves the footprints of
    valid observations into one geometry (see
    `tatc.analysis.footprints.dissolve_footprints`).

    Args:
        platforms (SurfacePlatform | list[SurfacePlatform]): The platform(s).
        times (list[datetime.datetime] | None): The times to sample (see
                `collect_surface_track`), or `None` for static footprints
                (whose time-varying requirements are not evaluated).
        instrument_index (int | None): The index of the observing instrument
                on each platform, or `None` for every instrument of each
                platform.
        elevation (float | None): The elevation (meters) above the WGS 84
                ellipsoid of the observed target for which to compute
                footprints, below (looking down) or above (looking up) each
                platform's instruments. Defaults to each platform's ground
                elevation.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame | geopandas.GeoSeries | None):
                An optional mask, always interpreted in WGS84 (lon/lat)
                coordinates, to constrain results.

    Returns:
        geopandas.GeoDataFrame: The data frame of the aggregated surface
            track: one row with the observed area, or none if no observation
            is valid.
    """
    return _dissolve_footprints(
        collect_surface_track(
            platforms,
            times,
            instrument_index=instrument_index,
            elevation=elevation,
            mask=mask,
        )
    )


def _warn_renamed(old: str, new: str) -> None:
    """
    Warns that a function was renamed.

    Args:
        old (str): The deprecated name.
        new (str): The new name.
    """
    warnings.warn(
        f"{old} is deprecated and will be removed in a future release: use {new} "
        "(whose results name each platform in a surface_platform column)",
        DeprecationWarning,
        stacklevel=3,
    )


def collect_platform_track(
    platforms: SurfacePlatform | list[SurfacePlatform],
    instrument_index: int | None = 0,
    elevation: float | None = None,
    mask: Polygon | MultiPolygon | gpd.GeoDataFrame | gpd.GeoSeries | None = None,
    **kwargs,
) -> gpd.GeoDataFrame:
    """
    Deprecated: use `collect_surface_track`, of which this is an alias (with
    the former order of arguments, without times). Its results name each
    platform in a `surface_platform` column (formerly `platform`).
    """
    _warn_renamed("collect_platform_track", "collect_surface_track")
    return collect_surface_track(
        platforms,
        instrument_index=instrument_index,
        elevation=elevation,
        mask=mask,
        **kwargs,
    )


def compute_platform_track(
    platforms: SurfacePlatform | list[SurfacePlatform],
    instrument_index: int | None = 0,
    elevation: float | None = None,
    mask: Polygon | MultiPolygon | gpd.GeoDataFrame | gpd.GeoSeries | None = None,
) -> gpd.GeoDataFrame:
    """
    Deprecated: use `compute_surface_track`, of which this is an alias (with
    the former order of arguments, without times).
    """
    _warn_renamed("compute_platform_track", "compute_surface_track")
    return compute_surface_track(
        platforms, instrument_index=instrument_index, elevation=elevation, mask=mask
    )

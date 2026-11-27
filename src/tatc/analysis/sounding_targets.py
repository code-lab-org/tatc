"""
Methods to collect the soundings of limb sounders and radio occultation
receivers as observations of targets.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timedelta
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from shapely import geometry as geo

from ..constants import EARTH_MEAN_RADIUS
from ..schemas import LimbSounder, OccultationReceiver, Satellite
from ..utils.geometry import _get_point_coordinates, hash_geometry
from .architecture import _split_targets
from .check import (
    _check_instrument_index,
    _check_satellites,
    _check_time_window,
)
from .limb_sampling import collect_limb_soundings
from .ro_sampling import collect_occultations
from .sampling import _get_empty_coverage_frame


def _collect_soundings(
    satellite: Satellite,
    start: datetime,
    end: datetime,
    instrument_index: int | None,
    transmitters: Satellite | list[Satellite] | None,
    time_step: timedelta,
) -> pd.DataFrame:
    """
    Collects the soundings of a satellite's sounders (see
    `collect_sounding_observations`): the representative point (`position`)
    and time of each, with its `satellite`, `instrument`, `start`, and `end`.
    """
    # pylint: disable=too-many-arguments,too-many-positional-arguments
    if instrument_index is None:
        kinds = {type(i) for i in satellite.instruments}
    else:
        kinds = {
            type(
                satellite.instruments[
                    _check_instrument_index(satellite, instrument_index)
                ]
            )
        }
        if not kinds & {LimbSounder, OccultationReceiver}:
            instrument = satellite.instruments[instrument_index]
            raise ValueError(
                f"instrument_index {instrument_index} of {satellite.name!r} is the "
                f"{type(instrument).__name__} {instrument.name!r}, not a sounder"
            )
    frames = []
    if LimbSounder in kinds:
        frames.append(collect_limb_soundings(satellite, start, end, instrument_index))
    if OccultationReceiver in kinds:
        if transmitters is None:
            raise ValueError(
                f"{satellite.name!r} carries an OccultationReceiver, whose "
                "occultations need transmitters (see tatc.catalog.get_transmitters)"
            )
        frames.append(
            collect_occultations(
                satellite,
                transmitters,
                start,
                end,
                instrument_index,
                time_step=time_step,
            ).rename(columns={"receiver": "satellite"})
        )
    columns = ["satellite", "instrument", "position", "start", "end", "time"]
    frames = [frame[columns] for frame in frames if not frame.empty]
    if not frames:
        return pd.DataFrame(columns=columns)
    return pd.concat(frames, ignore_index=True)


def _get_distances(
    longitude: np.ndarray, latitude: np.ndarray, point: tuple[float, float]
) -> np.ndarray:
    """
    Gets the great-circle distances (meters, on a sphere of the Earth's mean
    radius) from a point (longitude and latitude) to positions.
    """
    lon0, lat0 = np.radians(point)
    lon, lat = np.radians(longitude), np.radians(latitude)
    haversine = (
        np.sin((lat - lat0) / 2) ** 2
        + np.cos(lat0) * np.cos(lat) * np.sin((lon - lon0) / 2) ** 2
    )
    return 2 * EARTH_MEAN_RADIUS * np.arcsin(np.sqrt(np.clip(haversine, 0, 1)))


def collect_sounding_observations(
    targets: Any,
    platform: Any,
    start: datetime,
    end: datetime,
    *,
    radius: float,
    transmitters: Satellite | list[Satellite] | None = None,
    instrument_index: int | None = None,
    time_step: timedelta = timedelta(seconds=10),
) -> gpd.GeoDataFrame:
    """
    Collects the soundings of the limb sounders and radio occultation receivers
    that a satellite carries (see `collect_limb_soundings` and
    `collect_occultations`) as observations of points and regions, with the
    columns of the observations of `collect_observations`: a point is observed
    by each sounding whose representative point (its `position`, at the
    sounder's sample elevation) lies within a `radius` of it (a great-circle
    distance, on a sphere of the Earth's mean radius), and a region by each
    sounding whose representative point lies within it. Each observation's
    `epoch` is the sounding's time, its `start` and `end` the sounding's, and
    its `geometry` the point (with its elevation as its z coordinate) or, for
    a region, the sounding's representative point (the part of the region
    observed), with the target's `target_hash`. The angles of the satellite
    seen from a target (`platform_alt` and `platform_az`) are `NaN`, as a
    sounding does not view its targets.

    It is an analysis of one platform for a trade study's observations (see
    `tatc.tradespace.Observations`), called as `analysis(targets, platform,
    start, end)` with its options given by `functools.partial` (for example,
    `partial(collect_sounding_observations, radius=250e3, transmitters=gnss)`),
    so that a study's metrics of revisits and coverage count soundings. A
    platform other than a satellite, or a satellite without sounders, has no
    such observations. The observations of other instruments are collected by
    `collect_platform_observations` (whose `instrument_index` of `None` skips
    sounders), which an analysis may combine with these.

    Args:
        targets (Any): The point(s) and region(s) of interest (see
            `collect_platform_observations`).
        platform (Any): The observing platform.
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.
        radius (float): The distance (meters) from a point within which a
            sounding observes it.
        transmitters (Satellite | list[Satellite] | None): The transmitters of
            radio occultations (see `tatc.catalog.get_transmitters`), needed
            for a satellite that carries an occultation receiver.
        instrument_index (int | None): The index of the sounder on the
            satellite, or `None` for every sounder it carries.
        time_step (datetime.timedelta): The time step with which each
            occultation is sampled.

    Returns:
        geopandas.GeoDataFrame: The observations, sorted by start time.
    """
    # pylint: disable=too-many-arguments,too-many-locals
    if radius <= 0:
        raise ValueError("radius must be positive")
    _check_time_window(start, end)
    empty = _get_empty_coverage_frame(True)
    if not isinstance(platform, Satellite):
        return empty
    satellite = _check_satellites(platform)[0]
    points, regions = _split_targets(targets)
    soundings = _collect_soundings(
        satellite, start, end, instrument_index, transmitters, time_step
    )
    if soundings.empty:
        return empty
    positions = soundings["position"].to_numpy()
    longitude, latitude = shapely.get_x(positions), shapely.get_y(positions)
    rows = []
    for point in points:
        target = geo.Point(*_get_point_coordinates(point))
        near = _get_distances(longitude, latitude, (target.x, target.y)) <= radius
        rows.append((hash_geometry(target), [target] * int(near.sum()), near))
    for region in regions:
        inside = shapely.contains_xy(region, longitude, latitude)
        rows.append((hash_geometry(region), list(positions[inside]), inside))
    frames = [
        pd.DataFrame(
            {
                "target_hash": target_hash,
                "geometry": geometries,
                "satellite": soundings["satellite"].to_numpy()[selected],
                "instrument": soundings["instrument"].to_numpy()[selected],
                "start": soundings["start"].to_numpy()[selected],
                "epoch": soundings["time"].to_numpy()[selected],
                "end": soundings["end"].to_numpy()[selected],
                "platform_alt": np.nan,
                "platform_az": np.nan,
            }
        )
        for target_hash, geometries, selected in rows
        if selected.any()
    ]
    if not frames:
        return empty
    observations = pd.concat(frames, ignore_index=True)
    for column in ("start", "epoch", "end"):
        observations[column] = pd.to_datetime(observations[column], utc=True)
    return gpd.GeoDataFrame(
        observations.sort_values("start", kind="stable", ignore_index=True),
        geometry="geometry",
        crs="EPSG:4326",
    )

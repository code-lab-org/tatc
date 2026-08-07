"""
Methods to perform coverage analysis.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timedelta, timezone

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import geometry as geo
from skyfield.api import wgs84

from ..constants import EARTH_MEAN_RADIUS, de421, timescale
from ..schemas import Point, PointedInstrument, Satellite
from ..utils.observation import (
    compute_max_access_time,
    compute_min_elevation_angle,
)
from ..utils.orbital import compute_apoapsis_radius
from ..utils.projection import compute_footprint


def _get_visible_interval_series(
    point: Point,
    satellite: Satellite,
    min_elevation_angle: float,
    max_altitude: float,
    start: datetime,
    end: datetime,
) -> pd.Series:
    """
    Get the series of visible intervals based on altitude angle constraints.

    Args:
        point (Point): Point to observe.
        satellite (Satellite): Satellite doing the observation.
        min_elevation_angle (float): Minimum elevation angle (degrees) for valid observation.
        max_altitude (float): A conservative upper-bound satellite altitude
                (meters, e.g. the orbit's apogee altitude), used only to
                compute a generously large `max_access_time` bound for
                matching rise events to their corresponding set events.
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.

    Returns:
        pandas.Series: Series of observation intervals.
    """
    # compute the maximum access time to filter bad data
    max_access_time = timedelta(
        seconds=compute_max_access_time(max_altitude, min_elevation_angle)
    )
    # find the set of observation events
    times, events = satellite.orbit.to_gp_orbit().get_observation_events(
        point, start, end, min_elevation_angle
    )

    # build the observation periods
    obs_periods = []
    if len(events) == 0:
        # no rise, culminate, or set event was captured in [start, end]. This
        # means the elevation angle never crossed min_elevation_angle and had
        # no interior local maximum in this window -- which happens both
        # when the point is never visible, and when [start, end] falls
        # entirely within a longer visible pass (no rise/set inside the
        # window, and the window is too narrow, or off-center, to contain
        # the pass's culmination). Disambiguate by sampling the true
        # elevation angle at the window's midpoint.
        mid = start + (end - start) / 2
        topos = wgs84.latlon(point.latitude, point.longitude, point.elevation)
        orbit_track = satellite.orbit.to_gp_orbit().get_orbit_track(mid)
        elevation_angle = (
            orbit_track - topos.at(timescale.from_datetime(mid))
        ).altaz()[0].degrees
        if elevation_angle > min_elevation_angle:
            # continuously visible for the entire window
            obs_periods += [
                pd.Interval(
                    left=pd.Timestamp(start.astimezone(tz=timezone.utc)),
                    right=pd.Timestamp(end.astimezone(tz=timezone.utc)),
                )
            ]
    elif np.all(events == 1):
        # if all events are type 1 (culminate), create a period from start to end
        obs_periods += [
            pd.Interval(
                left=pd.Timestamp(start.astimezone(tz=timezone.utc)),
                right=pd.Timestamp(end.astimezone(tz=timezone.utc)),
            )
        ]
    else:
        # otherwise, match rise/set events
        rises = times[events == 0]
        sets = times[events == 2]
        if (
            len(sets) > 0
            and (len(rises) == 0 or sets[0].utc_datetime() < rises[0].utc_datetime())
            and start < sets[0].utc_datetime()
        ):
            # if first event is a set, create a period from the start
            obs_periods += [
                pd.Interval(
                    left=pd.Timestamp(start.astimezone(tz=timezone.utc)),
                    right=pd.Timestamp(sets[0].utc_datetime()),
                )
            ]
        # create an observation period to match with each rise event if
        # there is a following set event within twice the maximum access time
        obs_periods += [
            pd.Interval(
                left=pd.Timestamp(rise.utc_datetime()),
                right=pd.Timestamp(
                    sets[
                        np.logical_and(
                            rise.utc_datetime() < sets.utc_datetime(),
                            sets.utc_datetime()
                            < rise.utc_datetime() + 2 * max_access_time,
                        )
                    ][0].utc_datetime()
                ),
            )
            for rise in rises
            if np.any(
                np.logical_and(
                    rise.utc_datetime() < sets.utc_datetime(),
                    sets.utc_datetime() < rise.utc_datetime() + 2 * max_access_time,
                )
            )
        ]
        if (
            len(rises) > 0
            and (len(sets) == 0 or rises[-1].utc_datetime() > sets[-1].utc_datetime())
            and rises[-1].utc_datetime() < end
        ):
            # if last event is a rise, create a period to the end
            obs_periods += [
                pd.Interval(
                    left=pd.Timestamp(rises[-1].utc_datetime()),
                    right=pd.Timestamp(end.astimezone(tz=timezone.utc)),
                )
            ]
    return pd.Series(obs_periods, dtype="interval")


def _get_empty_coverage_frame(omit_solar: bool) -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for coverage analysis results.

    Args:
        omit_solar (bool): `True`, to omit solar angles to improve performance.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "point_id": pd.Series([], dtype="int"),
        "geometry": pd.Series([], dtype="object"),
        "satellite": pd.Series([], dtype="str"),
        "instrument": pd.Series([], dtype="str"),
        "start": pd.Series([], dtype="datetime64[ns, utc]"),
        "epoch": pd.Series([], dtype="datetime64[ns, utc]"),
        "end": pd.Series([], dtype="datetime64[ns, utc]"),
        "sat_alt": pd.Series(dtype="float"),
        "sat_az": pd.Series(dtype="float"),
    }
    if not omit_solar:
        columns = {
            **columns,
            "sat_sunlit": pd.Series(dtype="bool"),
            "solar_alt": pd.Series(dtype="float"),
            "solar_az": pd.Series(dtype="float"),
            "solar_time": pd.Series(dtype="float"),
        }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def collect_observations(
    point: Point,
    satellite: Satellite,
    start: datetime,
    end: datetime,
    instrument_index: int = 0,
    omit_solar: bool = True,
) -> gpd.GeoDataFrame:
    """
    Collect single satellite observations of a geodetic point of interest.

    Args:
        point (Point): The ground point of interest.
        satellite (Satellite): The observing satellite.
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.
        instrument_index (int): The index of the observing instrument in satellite.
        omit_solar (bool): `True`, to omit solar angles to improve performance.

    Returns:
        geopandas.GeoDataFrame: The data frame with recorded observations.
    """
    instrument = satellite.instruments[instrument_index]
    # use the apogee altitude as a conservative upper bound for computing access times
    max_altitude = (
        compute_apoapsis_radius(
            satellite.orbit.get_semimajor_axis(), satellite.orbit.get_eccentricity()
        )
        - EARTH_MEAN_RADIUS
    )
    # compute the minimum altitude angle required for observation
    min_elevation_angle = compute_min_elevation_angle(
        max_altitude,
        instrument.field_of_regard,
    )
    records = [
        {
            "point_id": point.id,
            "geometry": geo.Point(point.longitude, point.latitude, point.elevation),
            "satellite": satellite.name,
            "instrument": instrument.name,
            "start": (
                period.left
                if not instrument.access_time_fixed
                else period.mid - instrument.min_access_time / 2
            ),
            "end": (
                period.right
                if not instrument.access_time_fixed
                else period.mid + instrument.min_access_time / 2
            ),
            "epoch": period.mid,
        }
        for period in _get_visible_interval_series(
            point, satellite, min_elevation_angle, max_altitude, start, end
        )
        # instrument validity (illumination, footprint containment) below is
        # only checked at each coarse period's midpoint, as an approximation
        # of the whole interval; a more general approach would refine the
        # exact observation period boundaries with Skyfield's find_discrete
        # using the instrument's own validity condition, but that is out of
        # scope for now
        if (
            instrument.min_access_time <= period.right - period.left
            and instrument.is_valid_observation(
                (
                    orbit_track := satellite.orbit.to_gp_orbit().get_orbit_track(
                        period.mid
                    )
                ),
                wgs84.latlon(point.latitude, point.longitude, point.elevation),
            ).all()
            and (
                not isinstance(instrument, PointedInstrument)
                or compute_footprint(
                    orbit_track=orbit_track,
                    cross_track_field_of_view=instrument.cross_track_field_of_view,
                    along_track_field_of_view=instrument.along_track_field_of_view,
                    roll_angle=instrument.roll_angle,
                    pitch_angle=instrument.pitch_angle,
                    is_rectangular=instrument.is_rectangular,
                    elevation=point.elevation,
                )[0].contains(
                    geo.Point(point.longitude, point.latitude)
                )
            )
        )
    ]

    # build the dataframe
    if len(records) > 0:
        gdf = gpd.GeoDataFrame(records, crs="EPSG:4326")
        topos = wgs84.latlon(point.latitude, point.longitude, point.elevation)
        ts = timescale.from_datetimes(gdf.epoch)
        orbit_track = satellite.orbit.to_gp_orbit().get_orbit_track(gdf.epoch.tolist())
        # append satellite altitude/azimuth columns
        sat_altaz = (orbit_track - topos.at(ts)).altaz()
        gdf["sat_alt"] = sat_altaz[0].degrees  # type: ignore
        gdf["sat_az"] = sat_altaz[1].degrees  # type: ignore
        if not omit_solar:
            # append satellite sunlit column
            gdf["sat_sunlit"] = orbit_track.is_sunlit(de421)
            # append solar altitude/azimuth columns
            sun_altaz = (
                (de421["earth"] + topos).at(ts).observe(de421["sun"]).apparent().altaz()
            )
            gdf["solar_alt"] = sun_altaz[0].degrees
            gdf["solar_az"] = sun_altaz[1].degrees
            # append local solar time column
            gdf["solar_time"] = (de421["earth"] + topos).at(ts).observe(
                de421["sun"]
            ).apparent().hadec()[0].hours + 12
    else:
        gdf = _get_empty_coverage_frame(omit_solar)
    return gdf


def collect_multi_observations(
    point: Point,
    satellites: Satellite | list[Satellite],
    start: datetime,
    end: datetime,
    omit_solar: bool = True,
) -> gpd.GeoDataFrame:
    """
    Collect multiple satellite observations of a geodetic point of interest:
    calls `collect_observations` for every instrument on every satellite in
    `satellites`, and concatenates the results into one data frame.

    Args:
        point (Point): The ground point of interest.
        satellites (Satellite | list[Satellite]): The observing satellite(s),
                each contributing an observation per instrument it carries.
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.
        omit_solar (bool): `True`, to omit solar angles to improve performance.

    Returns:
        geopandas.GeoDataFrame: The data frame with all recorded observations.
    """
    gdfs = [
        collect_observations(point, satellite, start, end, instrument_index, omit_solar)
        for satellite in (satellites if isinstance(satellites, list) else [satellites])
        for instrument_index in range(len(satellite.instruments))
    ]
    if len(gdfs) == 0:
        # an empty `satellites` list leaves nothing to concatenate
        return _get_empty_coverage_frame(omit_solar)
    # concatenate into one data frame, sort by start time, and re-index
    return pd.concat(gdfs).sort_values("start").reset_index(drop=True)


def _get_empty_aggregate_frame() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for aggregated coverage analysis results.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "point_id": pd.Series([], dtype="int"),
        "geometry": pd.Series([], dtype="object"),
        "satellite": pd.Series([], dtype="str"),
        "instrument": pd.Series([], dtype="str"),
        "start": pd.Series([], dtype="datetime64[ns, utc]"),
        "epoch": pd.Series([], dtype="datetime64[ns, utc]"),
        "end": pd.Series([], dtype="datetime64[ns, utc]"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def aggregate_observations(observations: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Aggregate constellation observations. Interleaves observations by multiple
    satellites to compute aggregate performance metrics including access
    (observation duration) and revisit (duration between observations).
    Overlapping (including fully nested) observations for the same point,
    possibly from different satellites/instruments, are merged into a single
    continuous coverage period; `satellite`/`instrument` record every
    contributor to that period, comma-separated. `epoch` is reassigned to
    the midpoint of the merged period's `start`/`end` (a representative
    instant), not the mean of the constituent observations' own epochs.
    Per-observation columns that lose their meaning once merged across
    satellites and over a potentially much longer period -- e.g. `sat_alt`,
    `sat_az`, `sat_sunlit`, `solar_alt`, `solar_az`, `solar_time` -- are
    intentionally dropped, even if present on `observations`.

    Args:
        observations (geopandas.GeoDataFrame): The collected observations.

    Returns:
        geopandas.GeoDataFrame: The data frame with aggregated observations.
    """
    if observations.empty:
        return _get_empty_aggregate_frame()
    gdfs = []
    # split into constituent data frames based on point_id
    for _, gdf in observations.groupby("point_id"):
        # sort the values by start datetime
        gdf = gdf.sort_values("start")
        # assign the observation group number based on overlapping start/end times
        gdf["obs"] = (gdf["start"] > gdf["end"].shift().cummax()).cumsum()
        # perform the aggregation to group overlapping observations
        gdf = gdf.dissolve(
            "obs",
            aggfunc={
                "point_id": "first",
                "satellite": ", ".join,  # type: ignore
                "instrument": ", ".join,  # type: ignore
                "start": "min",
                "end": "max",
            },
        )
        # reassign epoch to the midpoint of the merged period, as a single
        # representative instant, rather than the mean of the constituent
        # observations' own (pre-merge) epochs
        gdf["epoch"] = gdf["start"] + (gdf["end"] - gdf["start"]) / 2
        # compute access and revisit metrics
        gdf["access"] = gdf["end"] - gdf["start"]
        gdf["revisit"] = gdf["start"] - gdf["end"].shift()
        # append to the list of data frames
        gdfs.append(gdf)
    # return a concatenated data frame and re-index
    return pd.concat(gdfs).reset_index(drop=True)


def _get_empty_reduce_frame() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for reduced coverage analysis results.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "point_id": pd.Series([], dtype="int"),
        "geometry": pd.Series([], dtype="object"),
        "access": pd.Series([], dtype="timedelta64[ns]"),
        "revisit": pd.Series([], dtype="timedelta64[ns]"),
        "samples": pd.Series([], dtype="int"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def reduce_observations(aggregated_observations: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Reduce constellation observations: for each unique point_id in
    `aggregated_observations`, computes the mean access period, the mean
    revisit period, and the total number of samples (aggregated periods)
    over the analysis period. The first sample's revisit is undefined (no
    prior observation to measure a gap from) and is excluded from the mean
    rather than counted as zero, which would otherwise bias the mean
    downward; a point with only one sample accordingly has an undefined
    (NaT) mean revisit.

    Args:
        aggregated_observations (geopandas.GeoDataFrame): The aggregated observations.

    Returns:
        geopandas.GeoDataFrame: The data frame with reduced observations.
    """
    if aggregated_observations.empty:
        return _get_empty_reduce_frame()
    # operate on a copy of the data frame
    gdf = aggregated_observations.copy()
    # convert access and revisit to numeric values before aggregation
    gdf["access"] = gdf["access"].dt.total_seconds()
    gdf["revisit"] = gdf["revisit"].dt.total_seconds()
    # assign each record to one observation
    gdf["samples"] = 1
    # perform the aggregation operation
    gdf = gdf.dissolve(
        "point_id",
        aggfunc={
            "access": "mean",
            "revisit": "mean",
            "samples": "sum",
        },
    ).reset_index()
    # convert access and revisit from numeric values after aggregation
    gdf["access"] = pd.to_timedelta(gdf["access"], unit="s")
    gdf["revisit"] = pd.to_timedelta(gdf["revisit"], unit="s")
    return gdf


def grid_observations(
    reduced_observations: gpd.GeoDataFrame, cells: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """
    Grid reduced observations to cells: for every cell, sums the number of
    samples across every point it contains, and combines those points'
    access/revisit statistics into a single representative value per cell.
    Both access (a per-event duration) and revisit (a time-between-events
    duration, i.e. the reciprocal of a sampling rate) use a sample-weighted
    mean -- arithmetic for access, harmonic for revisit, since revisit
    needs to be averaged as a rate to stay a representative statistic for
    a typical point in the cell: unlike the harmonic mean, summing
    reciprocal rates directly (without normalizing by sample count) would
    make a cell's reported revisit shrink simply because more (possibly
    near-identical) points happen to fall inside it, which is a property
    of the input point density, not of the underlying coverage geometry.

    Args:
        reduced_observations (geopandas.GeoDataFrame): The reduced observations.
        cells (geopandas.GeoDataFrame): The cell specification.

    Returns:
        geopandas.GeoDataFrame: The data frame with gridded observations.
    """
    if reduced_observations.empty:
        gdf = cells.copy()
        gdf["samples"] = 0
        gdf["access"] = None
        gdf["revisit"] = None
        return gdf
    # operate on a copy of the data frame
    gdf = reduced_observations.copy()
    # convert access and revisit to numeric values before aggregation
    gdf["access"] = gdf["access"].dt.total_seconds()
    gdf["revisit"] = gdf["revisit"].dt.total_seconds()
    # pre-transform so the means below reduce to plain sums: groupby().agg()
    # with a dict of {column: function} only ever hands a custom callable
    # its own column's Series, never a sibling column like "samples" needed
    # to compute a weighted statistic within the callable. The weighted
    # harmonic mean of revisit is sum(samples) / sum(samples/revisit).
    gdf["access_x_samples"] = gdf["access"] * gdf["samples"]
    gdf["samples_over_revisit"] = gdf["samples"] / gdf["revisit"]
    gdf = (
        cells.sjoin(gdf, how="inner", predicate="contains")
        .dissolve(
            by="cell_id",
            aggfunc={
                "samples": "sum",
                "access_x_samples": "sum",
                "samples_over_revisit": "sum",
            },
        )
        .reset_index()
    )
    # finish the aggregation: sample-weighted arithmetic mean for access,
    # sample-weighted harmonic mean for revisit
    gdf["access"] = gdf["access_x_samples"] / gdf["samples"]
    gdf["revisit"] = gdf["samples"] / gdf["samples_over_revisit"]
    gdf = gdf.drop(columns=["access_x_samples", "samples_over_revisit"])
    # convert access and revisit from numeric values after aggregation
    gdf["access"] = pd.to_timedelta(gdf["access"], unit="s")
    gdf["revisit"] = pd.to_timedelta(gdf["revisit"], unit="s")
    return gdf

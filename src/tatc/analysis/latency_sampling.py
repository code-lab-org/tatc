"""
Methods to perform latency analysis.

@author: Isaac Feldman
@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timedelta
from typing import Literal

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import geometry as geo

from ..schemas import GroundStation, Satellite
from ..schemas.architecture import Relay
from ..utils.computation import _run_together
from .check import (
    _check_satellites,
    _check_time_window,
    _get_observation_times,
    _get_platform_columns,
)
from .coverage_metrics import _grid_means, _reduce_by_target
from .point_sampling import _compute_access_periods
from .space_sampling import collect_space_observations


def _get_empty_downlinks_frame() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for downlink results.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "ground_station": pd.Series([], dtype="str"),
        "geometry": pd.Series([], dtype="object"),
        "satellite": pd.Series([], dtype="str"),
        "start": pd.Series([], dtype="datetime64[ns, utc]"),
        "epoch": pd.Series([], dtype="datetime64[ns, utc]"),
        "end": pd.Series([], dtype="datetime64[ns, utc]"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def collect_downlinks(
    stations: GroundStation | list[GroundStation],
    satellites: Satellite | list[Satellite],
    start: datetime,
    end: datetime,
) -> gpd.GeoDataFrame:
    """
    Collect satellite downlink opportunities to ground station(s) of
    interest. The access periods of every station and satellite are
    computed together (see `tatc.utils.computation._run_together`).

    Args:
        stations (GroundStation | list[GroundStation]): The ground stations.
        satellites (Satellite | list[Satellite]): The observing satellite(s).
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.

    Returns:
        geopandas.GeoDataFrame: The data frame of collected downlink results,
            each naming its ground station in a `ground_station` column,
            sorted by start time.
    """
    satellites = _check_satellites(satellites)
    _check_time_window(start, end)
    stations = [stations] if isinstance(stations, GroundStation) else stations
    pairs = [(station, satellite) for satellite in satellites for station in stations]
    access_periods = _run_together(
        [
            _compute_access_periods(
                station, satellite, start, end, station.min_elevation_angle
            )
            for station, satellite in pairs
        ]
    )
    # collect the records of ground station overpasses
    records = [
        {
            "ground_station": station.name,
            "geometry": geo.Point(
                station.longitude, station.latitude, station.elevation
            ),
            "satellite": satellite.name,
            "start": period.left,
            "end": period.right,
            "epoch": period.mid,
        }
        for (station, satellite), periods in zip(pairs, access_periods)
        for period in periods
        if (station.min_access_time <= period.right - period.left)
    ]
    # build the dataframe
    if len(records) > 0:
        gdf = gpd.GeoDataFrame(records, crs="EPSG:4326")
        # times in nanoseconds (as in the empty frame and the observations)
        for column in ["start", "end", "epoch"]:
            gdf[column] = gdf[column].dt.as_unit("ns")
        gdf = gdf.sort_values("start").reset_index(drop=True)
    else:
        gdf = _get_empty_downlinks_frame()
    return gdf


_RELAY_SEARCH = timedelta(seconds=60)
"""The shortest contact with a relay that the search guarantees to detect,
if the relay's `min_access_time` is shorter (see `collect_relay_contacts`)."""


def _get_empty_relay_contacts_frame() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for relay contacts (see `collect_relay_contacts`).

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "ground_station": pd.Series([], dtype="str"),
        "relay": pd.Series([], dtype="str"),
        "relay_satellite": pd.Series([], dtype="str"),
        "geometry": pd.Series([], dtype="object"),
        "satellite": pd.Series([], dtype="str"),
        "start": pd.Series([], dtype="datetime64[ns, utc]"),
        "epoch": pd.Series([], dtype="datetime64[ns, utc]"),
        "end": pd.Series([], dtype="datetime64[ns, utc]"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def collect_relay_contacts(
    relays: Relay | list[Relay],
    satellites: Satellite | list[Satellite],
    start: datetime,
    end: datetime,
) -> gpd.GeoDataFrame:
    """
    Collect the contacts of satellites with data relays (see
    `tatc.schemas.Relay`), which act as downlinks (see `compute_latencies`):
    the periods when a satellite has a line of sight to a relay satellite
    that grazes the Earth no lower than the relay's `min_grazing_altitude`
    (and within its `max_range`, see
    `tatc.analysis.collect_space_observations`), and, if the relay has
    ground stations, when the relay satellite is also in view of one of
    them (its access periods, see `collect_downlinks`), lasting at least the
    relay's `min_access_time`. The search guarantees to detect a line of
    sight that lasts at least the relay's `min_access_time`, or 60 seconds
    if that is shorter.

    Args:
        relays (Relay | list[Relay]): The relays.
        satellites (Satellite | list[Satellite]): The satellites.
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.

    Returns:
        geopandas.GeoDataFrame: The contacts, each naming its `relay`, the
            `relay_satellite`, and the `ground_station` that receives its
            data (missing for a relay without ground stations), with the
            geometry of the lines of sight at the start and end of the line
            of sight (see `collect_space_observations`), sorted by start time.
    """
    satellites = _check_satellites(satellites)
    _check_time_window(start, end)
    relays = [relays] if isinstance(relays, Relay) else list(relays)
    frames = []
    for relay in relays:
        sight = collect_space_observations(
            satellites,
            start,
            end,
            to_satellites=relay.satellites,
            min_range=None,
            max_range=relay.max_range,
            min_grazing_altitude=relay.min_grazing_altitude,
            min_duration=max(relay.min_access_time, _RELAY_SEARCH),
        )
        if sight.empty:
            continue
        sight = pd.DataFrame(
            {
                "relay_satellite": sight["to_satellite"].to_numpy(),
                "geometry": sight.geometry.to_numpy(),
                "satellite": sight["from_satellite"].to_numpy(),
                "start": sight["start"].dt.as_unit("ns"),
                "end": sight["end"].dt.as_unit("ns"),
            }
        )
        if relay.ground_stations:
            # the parts of each line of sight when the relay satellite is in
            # view of a ground station
            access = collect_downlinks(
                relay.ground_stations, relay.satellites, start, end
            )
            access = pd.DataFrame(
                {
                    "relay_satellite": access["satellite"].to_numpy(),
                    "ground_station": access["ground_station"].to_numpy(),
                    "access_start": access["start"].dt.as_unit("ns"),
                    "access_end": access["end"].dt.as_unit("ns"),
                }
            )
            sight = sight.merge(access, on="relay_satellite")
            sight["start"] = sight[["start", "access_start"]].max(axis=1)
            sight["end"] = sight[["end", "access_end"]].min(axis=1)
            sight = sight.drop(columns=["access_start", "access_end"])
        else:
            sight["ground_station"] = None
        sight = sight[
            sight["end"] - sight["start"] >= max(relay.min_access_time, timedelta(0))
        ]
        sight = sight[sight["end"] > sight["start"]]
        if sight.empty:
            continue
        sight["relay"] = relay.name
        frames.append(sight)
    if not frames:
        return _get_empty_relay_contacts_frame()
    contacts = pd.concat(frames, ignore_index=True)
    contacts["epoch"] = contacts["start"] + (contacts["end"] - contacts["start"]) / 2
    contacts = gpd.GeoDataFrame(
        contacts[list(_get_empty_relay_contacts_frame().columns)],
        geometry="geometry",
        crs="EPSG:4326",
    )
    for column in ["start", "epoch", "end"]:
        contacts[column] = contacts[column].dt.as_unit("ns")
    return contacts.sort_values("start", kind="stable").reset_index(drop=True)


def _get_empty_latency_frame(
    platform_columns: list[str] | None = None, angles: bool = True
) -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for latency results.

    Args:
        platform_columns (list[str] | None): The columns that name the
            observing platforms (see `_PLATFORM_COLUMNS`), by default
            `satellite`.
        angles (bool): `True` to include the platform's angles
            (`platform_alt` and `platform_az`), as of observations that
            have them.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "target_hash": pd.Series([], dtype="str"),
        "geometry": pd.Series([], dtype="object"),
        **{
            column: pd.Series([], dtype="str")
            for column in platform_columns or ["satellite"]
        },
        "instrument": pd.Series([], dtype="str"),
        **(
            {
                "platform_alt": pd.Series([], dtype="float"),
                "platform_az": pd.Series([], dtype="float"),
            }
            if angles
            else {}
        ),
        "ground_station": pd.Series([], dtype="str"),
        "downlinked": pd.Series([], dtype="datetime64[ns, utc]"),
        "latency": pd.Series([], dtype="timedelta64[ns]"),
        "observed": pd.Series([], dtype="datetime64[ns, utc]"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def compute_latencies(
    observations: gpd.GeoDataFrame,
    downlinks: gpd.GeoDataFrame,
    during_contact: Literal["next", "end", "immediate"] = "end",
) -> gpd.GeoDataFrame:
    """
    Collect latencies between an observation and the first downlink opportunity.

    An observation that ends before a downlink starts is downlinked at that
    downlink's epoch (midpoint). The `during_contact` option sets how an
    observation that ends while a downlink is in progress is downlinked:
    `"end"` (default) downlinks it at the end of the downlink in progress
    (stored data is played back after the data recorded before the contact);
    `"next"` waits for the next downlink to start (stored data is only
    played back from the start of a contact); `"immediate"` downlinks it as it is
    observed (real-time downlink), at the later of the observation epoch and
    the start of the downlink in progress. Each observation is assigned the
    earliest of these downlink times.

    The observations are those of each satellite (not aggregated, see
    `aggregate_observations`, whose periods may merge several satellites'):
    of points (see `collect_observations` or `collect_footprint_observations`),
    of regions (see `collect_region_observations`), or of features (see
    `collect_feature_observations`), each at its time (its `epoch`, or else
    its `time`), ending at its `end` (or its time, if none). Only satellites
    downlink (see `collect_downlinks`): the observations of other platforms
    (aircraft, fixed surface platforms, and radar stations, named in the
    column of their kind) are kept without a downlink, as an observation
    without a matching downlink is (`NaT` latency).

    Args:
        observations (geopandas.GeoDataFrame): The data frame of observations to downlink.
        downlinks (geopandas.GeoDataFrame): The data frame of downlink opportunities.
        during_contact (str): Downlink of observations that end during a
            downlink opportunity: `"end"` (default), `"next"`, or `"immediate"`.

    Returns:
        geopandas.GeoDataFrame: The data frame of collected latency results, sorted by the 'observed'
        column in ascending order. It includes the following columns, of
        those that the observations have (with every observation, `NaT`
        latency without a matching downlink, as without any downlinks):

        - 'target_hash' (str): Hash of the observation point (see
          `tatc.utils.geometry.hash_geometry`).
        - 'geometry' (geometry): Geometry representing the observation point.
        - 'satellite' (object): Name or identifier of the satellite (and
          the other platform columns of the observations, if any).
        - 'instrument' (object): Name or identifier of the instrument.
        - 'platform_alt' (float64): Altitude of the satellite at the time of observation.
        - 'platform_az' (float64): Azimuth of the satellite at the time of observation.
        - 'ground_station' (object): Name or identifier of the ground station for downlink.
        - 'downlinked' (datetime64[ns, UTC]): Timestamp when the observation data was downlinked.
        - 'latency' (timedelta64[ns]): Latency between observation and downlink.
        - 'observed' (datetime64[ns, UTC]): Timestamp when the observation was made.
    """
    if during_contact not in ("next", "end", "immediate"):
        raise ValueError(
            f"during_contact must be 'next', 'end', or 'immediate', not {during_contact!r}"
        )
    if "revisit" in observations.columns:
        raise ValueError(
            "observations must be those of each satellite, not aggregated "
            "(see aggregate_observations): compute their latencies before "
            "aggregating them"
        )
    # the columns of the result, of those the observations have
    empty = _get_empty_latency_frame(
        _get_platform_columns(observations),
        # (by default, those of satellites' observations of points)
        angles="platform_alt" in observations.columns or len(observations.columns) == 0,
    )
    if observations.empty:
        return empty
    obs = observations.copy()
    # the time of each observation, and its end (or its time, if none)
    obs["epoch"] = _get_observation_times(observations).dt.as_unit("ns").to_numpy()
    if "end" not in obs.columns:
        obs["end"] = obs["epoch"]
    found = _find_downlinks(obs, downlinks, during_contact)
    obs = obs.assign(**{column: found[column] for column in found.columns})
    obs = obs.sort_values(by="end").reset_index(drop=True)

    # compute latency
    obs["latency"] = obs["downlinked"] - obs["epoch"]
    obs = obs.rename(columns={"epoch": "observed"})

    # select relevant columns, of those the observations have (and the relay
    # of each downlink, if the downlinks name relays)
    columns = list(empty.columns)
    if "relay" in obs.columns:
        columns.insert(columns.index("ground_station") + 1, "relay")
    obs = obs[[column for column in columns if column in obs]]

    # ensure the result is a GeoDataFrame with the observations' CRS
    obs = gpd.GeoDataFrame(
        obs,
        geometry="geometry",
        crs=(
            observations.crs
            if isinstance(observations, gpd.GeoDataFrame) and observations.crs
            else "EPSG:4326"
        ),
    )

    # sort observations by observed time
    obs = obs.sort_values(by="observed").reset_index(drop=True)
    return obs


def _find_downlinks(
    observations: pd.DataFrame,
    downlinks: pd.DataFrame,
    during_contact: Literal["next", "end", "immediate"] = "end",
) -> pd.DataFrame:
    """
    Finds the downlink of each observation of satellites (see
    `compute_latencies`): its `ground_station` (and `relay`, if the downlinks
    name relays, see `collect_relay_contacts`) and the time it is
    `downlinked` (missing without one, as for the observations of other
    platforms).

    Args:
        observations (pandas.DataFrame): The observations, each with its
            time (`epoch`) and `end`.
        downlinks (pandas.DataFrame): The downlink opportunities.
        during_contact (str): Downlink of observations that end during a
            downlink opportunity (see `compute_latencies`).

    Returns:
        pandas.DataFrame: the downlink of each observation, with the
            observations' index
    """
    named = ["ground_station"] + (["relay"] if "relay" in downlinks.columns else [])
    if downlinks.empty:
        downlinks = _get_empty_downlinks_frame().assign(
            **{column: pd.Series([], dtype="str") for column in named}
        )
    obs = pd.DataFrame(
        {
            # the satellite of each observation (missing for other platforms),
            # as objects, so that it is compared with the downlinks' whatever
            # the dtypes
            "_satellite": pd.Series(
                (
                    observations["satellite"].to_numpy()
                    if "satellite" in observations.columns
                    else None
                ),
                index=range(len(observations)),
                dtype=object,
            ),
            "epoch": observations["epoch"].array,
            # merge keys must share a unit (frames may come from other
            # sources, such as files)
            "end": observations["end"].dt.as_unit("ns").array,
            "_row": np.arange(len(observations)),
        }
    )
    obs = obs.sort_values(by="end").reset_index(drop=True)
    contacts = downlinks[["satellite", *named, "start", "epoch", "end"]].rename(
        columns={
            "satellite": "_satellite",
            "start": "downlink_start",
            "epoch": "downlinked",
            "end": "downlink_end",
        }
    )
    for column in ["_satellite", *named]:
        contacts[column] = contacts[column].astype(object)
    for column in ["downlink_start", "downlinked", "downlink_end"]:
        contacts[column] = contacts[column].dt.as_unit("ns")
    # pair each observation with the first downlink that starts after it ends
    obs = pd.merge_asof(
        obs,
        contacts.sort_values(by="downlink_start"),
        by="_satellite",
        left_on="end",
        right_on="downlink_start",
        direction="forward",
    )
    if during_contact != "next":
        # find the first downlink that ends after the observation ends, which
        # is in progress if it started before the observation ends
        current = pd.merge_asof(
            obs[["_satellite", "end"]],
            contacts.sort_values(by="downlink_end"),
            by="_satellite",
            left_on="end",
            right_on="downlink_end",
            direction="forward",
        )
        in_progress = current["downlink_start"] <= current["end"]
        if during_contact == "end":
            downlinked = current["downlink_end"]
        else:
            downlinked = current["downlink_start"].where(
                current["downlink_start"] > obs["epoch"], obs["epoch"]
            )
        earlier = in_progress & (
            obs["downlinked"].isna() | (downlinked < obs["downlinked"])
        )
        for column in named:
            obs.loc[earlier, column] = current.loc[earlier, column]
        obs.loc[earlier, "downlinked"] = downlinked[earlier]
    obs = obs.sort_values(by="_row")
    return pd.DataFrame(
        {column: obs[column].array for column in [*named, "downlinked"]},
        index=observations.index,
    )


def reduce_latencies(
    latency_observations: gpd.GeoDataFrame,
    targets: gpd.GeoDataFrame | gpd.GeoSeries | None = None,
) -> gpd.GeoDataFrame:
    """
    Reduce observation latencies: for each unique target (`target_hash`,
    see `tatc.analysis.coverage_metrics._get_target_keys`) in
    `latency_observations`, computes the mean latency and the total number
    of samples (observation/downlink pairs). An observation with no
    matching downlink has an undefined (NaT) latency (see
    `compute_latencies`), which is excluded from the mean rather than
    counted as zero, while still counting toward `samples`. With `targets`,
    the result has a row for every target, in their order, as for
    `tatc.analysis.coverage_metrics.reduce_observations`: a target not
    observed has no samples (zero) and an undefined (NaT) latency.

    Args:
        latency_observations (geopandas.GeoDataFrame): The latency observations.
        targets (geopandas.GeoDataFrame | geopandas.GeoSeries | None): The
            targets, including those not observed, with the other columns of
            a data frame of targets, or `None` for those observed.

    Returns:
        geopandas.GeoDataFrame: The data frame with reduced latencies.
    """
    return _reduce_by_target(latency_observations, targets, ["latency"])


def grid_latencies(
    reduced_latencies: gpd.GeoDataFrame, cells: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """
    Grid reduced latencies to cells: for every cell, sums the number of
    samples across every point it contains, and combines those points'
    latency into a single sample-weighted arithmetic mean per cell. Every
    cell is gridded: a cell that contains no observed point has no samples
    (zero) and an undefined (NaT) latency.

    Args:
        reduced_latencies (geopandas.GeoDataFrame): The reduced latencies.
        cells (geopandas.GeoDataFrame): The cell specification, with a
            `cell_id` column.

    Returns:
        geopandas.GeoDataFrame: The cells, with the `samples` and `latency`
            of each.
    """
    return _grid_means(reduced_latencies, cells, ["latency"])

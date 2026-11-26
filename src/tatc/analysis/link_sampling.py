"""
Methods to sample the geometry of the contacts of satellites with ground
stations and relays on a grid of time, and the rates of their links (see
`tatc.utils.compute_link_rate`), for data flow analyses (see
`tatc.analysis.compute_data_flow`).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from datetime import datetime, timedelta

import numpy as np
import pandas as pd
from skyfield.api import wgs84

from ..schemas import GroundStation, Relay, Satellite
from ..schemas.communication import Receiver, Transmitter
from ..utils.link_budget import compute_link_rate
from ..utils.time import _index_orbit_track, _to_time_from_offsets

_SAMPLE_COLUMNS = {
    "contact": "int64",
    "satellite": "object",
    "ground_station": "object",
    "relay": "object",
    "relay_satellite": "object",
    "cell": "int64",
    "time": "datetime64[ns, UTC]",
    "fraction": "float64",
    "distance": "float64",
    "elevation": "float64",
}


@dataclass(frozen=True)
class _Terminal:
    """
    A transmitter of a satellite for data flow: a `Transmitter`, or a fixed
    `rate` (from its resources' downlink rate) to ground stations in any band.
    """

    name: str
    link: str
    transmitter: Transmitter | None = None
    rate: float | None = None


def _get_terminals(satellite: Satellite) -> list[_Terminal]:
    """
    Gets the transmitters of a satellite for data flow: those of its
    communications, or else a downlink at its resources' downlink rate (if
    given).
    """
    if satellite.communications is not None and satellite.communications.transmitters:
        return [
            _Terminal(name=t.get_name(), link=t.link, transmitter=t)
            for t in satellite.communications.transmitters
        ]
    if satellite.resources is not None and satellite.resources.downlink_rate:
        return [
            _Terminal(
                name="downlink", link="downlink", rate=satellite.resources.downlink_rate
            )
        ]
    return []


def get_generation_rate(satellite: Satellite) -> float | None:
    """
    Gets the rate (bits per second) at which a satellite generates data while
    active: its communications' generation `rate`, or else the total data
    rate of its instruments (see `SpaceSystem.get_payload_resources`).

    Args:
        satellite (Satellite): The satellite.

    Returns:
        float | None: The rate, or `None` if not given.
    """
    if (
        satellite.communications is not None
        and satellite.communications.generation.rate is not None
    ):
        return satellite.communications.generation.rate
    return satellite.get_payload_resources().data_rate


def _get_receivers(endpoint: GroundStation | Relay) -> list[Receiver]:
    """Gets the receivers of a ground station or relay (none if not given)."""
    if endpoint.communications is None:
        return []
    return list(endpoint.communications.receivers)


def _get_rate(
    terminal: _Terminal,
    receivers: list[Receiver],
    distance: np.ndarray,
    elevation: np.ndarray | None,
) -> np.ndarray | None:
    """
    Gets the rates of a terminal's link to an endpoint with receivers (any band
    without receivers; the best receiver of its band otherwise), or `None` if
    the link is not compatible (no receiver in its band) or has no known rate.
    """
    if terminal.transmitter is None:
        # a fixed rate in any band, at most the fastest receiver's maximum
        caps = [r.max_data_rate for r in receivers]
        rate = terminal.rate
        if receivers and all(cap is not None for cap in caps):
            rate = min(rate, max(caps))
        return np.full(distance.shape, rate)
    transmitter = terminal.transmitter
    if not receivers:
        candidates: list[Receiver | None] = [None]
    else:
        candidates = [r for r in receivers if r.band == transmitter.band]
        if not candidates:
            return None
    rates = [
        np.asarray(compute_link_rate(transmitter, r, distance, elevation), dtype=float)
        for r in candidates
    ]
    rate = np.max(rates, axis=0)
    if not np.all(np.isfinite(rate)):
        # a modem without a maximum and without a budget: no known rate
        return None
    return rate


def _to_utc(time: datetime) -> pd.Timestamp:
    """Gets a time as a timestamp in UTC (taken as UTC without a time zone)."""
    stamp = pd.Timestamp(time)
    return stamp.tz_convert("UTC") if stamp.tzinfo else stamp.tz_localize("UTC")


def _get_cells(
    contacts: pd.DataFrame, start: datetime, step: timedelta
) -> pd.DataFrame:
    """
    Gets the cells of the grid of time (`start + k * step`) that each contact
    overlaps: the contact's row (`contact`), the cell's index (`cell`), the
    fraction of the cell that the contact overlaps (`fraction`), and the
    offset (seconds) from the start of the middle of the overlap (`offset`).
    """
    step_s = step.total_seconds()
    reference = _to_utc(start)
    first = (
        (pd.to_datetime(contacts["start"], utc=True) - reference).dt.total_seconds()
    ).to_numpy(dtype=float)
    last = (
        (pd.to_datetime(contacts["end"], utc=True) - reference).dt.total_seconds()
    ).to_numpy(dtype=float)
    k0 = np.floor(first / step_s).astype(np.int64)
    k1 = np.maximum(np.ceil(last / step_s).astype(np.int64), k0)
    counts = k1 - k0
    contact = np.repeat(np.arange(len(contacts)), counts)
    within = np.arange(counts.sum()) - np.repeat(np.cumsum(counts) - counts, counts)
    cell = np.repeat(k0, counts) + within
    low = np.maximum(np.repeat(first, counts), cell * step_s)
    high = np.minimum(np.repeat(last, counts), (cell + 1) * step_s)
    keep = high > low
    return pd.DataFrame(
        {
            "contact": contact[keep],
            "cell": cell[keep],
            "fraction": ((high - low) / step_s)[keep],
            "offset": ((low + high) / 2)[keep],
        }
    )


def _empty_samples() -> pd.DataFrame:
    return pd.DataFrame(
        {
            column: pd.Series([], dtype=dtype)
            for column, dtype in _SAMPLE_COLUMNS.items()
        }
    )


def sample_link_geometry(
    contacts: pd.DataFrame,
    satellites: Sequence[Satellite],
    ground_stations: Sequence[GroundStation] = (),
    relays: Sequence[Relay] = (),
    *,
    start: datetime,
    step: timedelta = timedelta(seconds=30),
) -> pd.DataFrame:
    """
    Samples the geometry of contacts (see `collect_downlinks` and
    `collect_relay_contacts`) on a grid of time from `start` with `step`: for
    each cell of the grid that a contact overlaps, the fraction it overlaps
    and, at the middle of the overlap, the distance between the satellite and
    the ground station (with the satellite's elevation angle seen from it) or
    the relay satellite.

    Args:
        contacts (pandas.DataFrame): The contacts, each naming its `satellite`
            and its `ground_station` or its `relay` and `relay_satellite`,
            with its `start` and `end`.
        satellites (Sequence[Satellite]): The satellites.
        ground_stations (Sequence[GroundStation]): The ground stations.
        relays (Sequence[Relay]): The relays.
        start (datetime.datetime): The start of the grid of time.
        step (datetime.timedelta): The step of the grid of time.

    Returns:
        pandas.DataFrame: One row per contact and cell: the row of the contact
            (`contact`, its position in `contacts`), its `satellite`,
            `ground_station`, `relay`, and `relay_satellite` (missing if not
            given or for a contact with a ground station), the `cell` of the
            grid, the `time` of the sample, the `fraction` of the cell, the
            `distance` (meters), and the `elevation` (degrees, missing for a
            relay).

    Raises:
        ValueError: If a contact names a member that is not given.
    """
    if contacts.empty:
        return _empty_samples()
    contacts = contacts.reset_index(drop=True)
    cells = _get_cells(contacts, start, step)
    if cells.empty:
        return _empty_samples()
    columns = {}
    for column in ["satellite", "ground_station", "relay", "relay_satellite"]:
        values = (
            contacts[column].astype(object).to_numpy()
            if column in contacts.columns
            else np.full(len(contacts), None, dtype=object)
        )
        values = np.where(pd.isna(values), None, values)
        columns[column] = values[cells["contact"].to_numpy()]
    samples = pd.DataFrame(
        {
            "contact": cells["contact"].to_numpy(),
            **columns,
            "cell": cells["cell"].to_numpy(),
            "fraction": cells["fraction"].to_numpy(),
        }
    )
    offsets = cells["offset"].to_numpy()
    samples["time"] = _to_utc(start) + pd.to_timedelta(offsets, unit="s")
    distance = np.full(len(samples), np.nan)
    elevation = np.full(len(samples), np.nan)
    satellites_by_name = {s.name: s for s in satellites}
    stations = {s.name: s for s in ground_stations}
    relays_by_name = {r.name: r for r in relays}
    is_relay = samples["relay"].notna().to_numpy()
    for name, rows in samples.groupby("satellite", sort=False).indices.items():
        if name not in satellites_by_name:
            raise ValueError(f"contacts name a satellite {name!r} that is not given")
        t = _to_time_from_offsets(start, offsets[rows])
        track = satellites_by_name[name].get_track(t)
        position = track.position.m
        # contacts with ground stations
        ground = rows[~is_relay[rows]]
        for station_name in pd.unique(samples["ground_station"].to_numpy()[ground]):
            if station_name not in stations:
                raise ValueError(
                    f"contacts name a ground station {station_name!r} that is "
                    + "not given"
                )
            station = stations[station_name]
            select = samples["ground_station"].to_numpy()[rows] == station_name
            select &= ~is_relay[rows]
            topos = wgs84.latlon(station.latitude, station.longitude, station.elevation)
            index = np.flatnonzero(select)
            sub = _index_orbit_track(track, index)
            alt, _, dist = (sub - topos.at(sub.t)).altaz()
            elevation[rows[index]] = alt.degrees
            distance[rows[index]] = dist.m
        # contacts with relays
        relayed = rows[is_relay[rows]]
        if len(relayed):
            pairs = samples.iloc[relayed][["relay", "relay_satellite"]]
            for (relay_name, relay_satellite), group in pairs.groupby(
                ["relay", "relay_satellite"], sort=False
            ).indices.items():
                if relay_name not in relays_by_name:
                    raise ValueError(
                        f"contacts name a relay {relay_name!r} that is not given"
                    )
                relay = relays_by_name[relay_name]
                relay_satellites = {s.name: s for s in relay.satellites}
                if relay_satellite not in relay_satellites:
                    raise ValueError(
                        f"contacts name a satellite {relay_satellite!r} that "
                        + f"relay {relay_name!r} does not have"
                    )
                index = np.flatnonzero(np.isin(rows, relayed[group]))
                sub_t = _to_time_from_offsets(start, offsets[rows[index]])
                other = relay_satellites[relay_satellite].get_track(sub_t).position.m
                distance[rows[index]] = np.linalg.norm(
                    position[:, index] - other, axis=0
                )
    samples["distance"] = distance
    samples["elevation"] = elevation
    return samples[list(_SAMPLE_COLUMNS)].astype(
        {column: dtype for column, dtype in _SAMPLE_COLUMNS.items() if column != "time"}
    )


def compute_link_rates(
    samples: pd.DataFrame,
    satellites: Sequence[Satellite],
    ground_stations: Sequence[GroundStation] = (),
    relays: Sequence[Relay] = (),
) -> pd.DataFrame:
    """
    Computes the rate of each transmitter of a satellite at each sample of
    its contacts (see `sample_link_geometry`), from a link budget where the
    transmitter and the receiver of its band give what one needs (see
    `tatc.utils.compute_link_rate`). A satellite's transmitters are those of
    its `communications`, or else a downlink at its resources'
    `downlink_rate` in any band. A transmitter's `downlink` reaches ground
    stations and its `crosslink` relays, whose receivers (if given) must
    include its band.

    Args:
        samples (pandas.DataFrame): The samples (see `sample_link_geometry`).
        satellites (Sequence[Satellite]): The satellites.
        ground_stations (Sequence[GroundStation]): The ground stations.
        relays (Sequence[Relay]): The relays.

    Returns:
        pandas.DataFrame: One row per sample and transmitter whose rate is
            positive: the columns of the samples, the `transmitter` (its
            name), and the `rate` (bits per second).
    """
    empty = samples.iloc[:0].assign(
        transmitter=pd.Series([], dtype=object), rate=pd.Series([], dtype=float)
    )
    if samples.empty:
        return empty
    satellites_by_name = {s.name: s for s in satellites}
    endpoints: dict[tuple[str, str], GroundStation | Relay] = {
        **{("ground_station", s.name): s for s in ground_stations},
        **{("relay", r.name): r for r in relays},
    }
    is_relay = samples["relay"].notna().to_numpy()
    keys = np.where(
        is_relay,
        samples["relay"].astype(object).to_numpy(),
        samples["ground_station"].astype(object).to_numpy(),
    )
    frames = []
    groups = (
        pd.DataFrame(
            {
                "satellite": samples["satellite"].to_numpy(),
                "relay": is_relay,
                "key": keys,
            }
        )
        .groupby(["satellite", "relay", "key"], sort=False)
        .indices
    )
    for (satellite_name, relayed, key), rows in groups.items():
        terminals = _get_terminals(satellites_by_name[satellite_name])
        endpoint = endpoints[("relay" if relayed else "ground_station", key)]
        receivers = _get_receivers(endpoint)
        sub = samples.iloc[rows]
        distance = sub["distance"].to_numpy()
        elevation = None if relayed else sub["elevation"].to_numpy()
        for terminal in terminals:
            if terminal.link != ("crosslink" if relayed else "downlink"):
                continue
            rate = _get_rate(terminal, receivers, distance, elevation)
            if rate is None:
                continue
            positive = rate > 0
            if np.any(positive):
                frames.append(
                    sub[positive].assign(transmitter=terminal.name, rate=rate[positive])
                )
    if not frames:
        return empty
    return pd.concat(frames, ignore_index=True)

"""
Methods to analyze the flow of data from satellites to the ground: the data
each satellite generates and stores, the links scheduled through ground
stations and relays with a limited number of antennas, and the time each
observation's data reaches the ground (its latency) or whether storage
overflows and it is lost. A coarse model for concept studies, on a grid of
time whose step sets the accuracy of times.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timedelta

import geopandas as gpd
import numpy as np
import pandas as pd
from skyfield.functions import angle_between

from .._warnings import warn
from ..constants import de421
from ..schemas import GroundStation, Relay, Satellite
from ..schemas.communication import DataFlow
from ..utils.time import _to_time_from_offsets
from .check import _get_observation_times, _get_platform_columns
from .latency_sampling import _find_downlinks, _get_empty_latency_frame
from .link_sampling import (
    _get_terminals,
    _to_utc,
    compute_link_rates,
    get_generation_rate,
    sample_link_geometry,
)

_ENDPOINT_COLUMNS = ["ground_station", "relay", "relay_satellite"]


def _get_number_cells(start: datetime, end: datetime, step: timedelta) -> int:
    """Gets the number of cells of the grid of time from start to end."""
    return max(0, math.ceil((end - start) / step))


def _get_cell_lengths(start: datetime, end: datetime, step: timedelta) -> np.ndarray:
    """Gets the length (seconds) of each cell, the last one ending at the end."""
    lengths = np.full(_get_number_cells(start, end, step), step.total_seconds())
    if len(lengths):
        lengths[-1] = (end - start).total_seconds() - step.total_seconds() * (
            len(lengths) - 1
        )
    return lengths


def compute_generation(
    satellite: Satellite,
    start: datetime,
    end: datetime,
    step: timedelta = timedelta(seconds=30),
) -> np.ndarray:
    """
    Computes the volume of data a satellite generates in each cell of a grid
    of time from `start` to `end` with `step`: its rate of generation (see
    `link_sampling.get_generation_rate`) times its duty cycle over the length
    of the cell, while it is active (its activity sampled at the middle of
    each cell; see `schemas.DataGeneration`).

    Args:
        satellite (Satellite): The satellite.
        start (datetime.datetime): The start of the grid.
        end (datetime.datetime): The end of the grid.
        step (datetime.timedelta): The step of the grid.

    Returns:
        numpy.ndarray: The volume (bits) of each cell (zeros without a rate).
    """
    lengths = _get_cell_lengths(start, end, step)
    rate = get_generation_rate(satellite)
    if rate is None or len(lengths) == 0:
        return np.zeros_like(lengths)
    generation = (
        satellite.communications.generation
        if satellite.communications is not None
        else None
    )
    volume = rate * lengths * (generation.duty_cycle if generation else 1)
    if generation is not None and generation.activity != "always":
        middles = step.total_seconds() * np.arange(len(lengths)) + lengths / 2
        track = satellite.get_track(_to_time_from_offsets(start, middles))
        if generation.activity == "sunlit":
            active = np.asarray(track.is_sunlit(de421), dtype=bool)
        else:
            # the point below is in daylight if the sun is within 90 degrees of
            # the satellite's direction from the Earth's center
            sun = de421["earth"].at(track.t).observe(de421["sun"]).position.m
            active = angle_between(track.position.m, sun) < np.pi / 2
        volume = np.where(active, volume, 0.0)
    return volume


def _get_endpoint_key(row: tuple) -> tuple:
    """Gets the key of the endpoint of a link: a ground station or relay satellite."""
    ground_station, relay, relay_satellite = row
    if relay is not None:
        return ("relay", relay, relay_satellite)
    return ("ground_station", ground_station, None)


def _empty_schedule() -> pd.DataFrame:
    """Gets an empty schedule of links (see `_simulate`)."""
    return pd.DataFrame(
        {
            **{
                column: pd.Series([], dtype=object)
                for column in ["satellite", "transmitter", *_ENDPOINT_COLUMNS]
            },
            "contact": pd.Series([], dtype="int64"),
            "cell": pd.Series([], dtype="int64"),
            "fraction": pd.Series([], dtype=float),
            "rate": pd.Series([], dtype=float),
            "capacity": pd.Series([], dtype=float),
            "volume": pd.Series([], dtype=float),
        }
    )


@dataclass(frozen=True)
class _Storage:
    """
    The storage of a satellite on the grid of time: the volume generated,
    lost, and delivered in each cell, and the backlog at the end of each.
    """

    generated: np.ndarray
    lost: np.ndarray
    delivered: np.ndarray
    backlog: np.ndarray
    initial_backlog: float


@dataclass(frozen=True)
class _Endpoints:
    """The ground stations and relay satellites that receive links, by key."""

    index: dict[tuple, int]
    antennas: list[float]
    min_access_times: list[float]
    feeder_rates: list[float]

    @classmethod
    def of(
        cls, ground_stations: Sequence[GroundStation], relays: Sequence[Relay]
    ) -> _Endpoints:
        """Gets the endpoints of ground stations and relays."""
        endpoints = cls({}, [], [], [])
        for station in ground_stations:
            endpoints.add(
                ("ground_station", station.name, None),
                station.communications.antennas if station.communications else np.inf,
                station.min_access_time,
                np.inf,
            )
        for relay in relays:
            communications = relay.communications
            for relay_satellite in relay.satellites:
                endpoints.add(
                    ("relay", relay.name, relay_satellite.name),
                    communications.antennas if communications else np.inf,
                    relay.min_access_time,
                    (
                        communications.feeder_rate
                        if communications and communications.feeder_rate
                        else np.inf
                    ),
                )
        return endpoints

    def add(
        self, key: tuple, antennas: float, min_access_time: timedelta, feeder: float
    ) -> None:
        """Adds an endpoint."""
        self.index[key] = len(self.antennas)
        self.antennas.append(antennas)
        self.min_access_times.append(min_access_time.total_seconds())
        self.feeder_rates.append(feeder)


class _Links:
    """
    The links that satellites can take in each cell of the grid of time (see
    `link_sampling.compute_link_rates`), sorted by cell and then by
    decreasing rate (a satellite's faster links first), from which links are
    scheduled cell by cell.
    """

    def __init__(
        self,
        rates: pd.DataFrame,
        satellites: Sequence[Satellite],
        endpoints: _Endpoints,
        n_cells: int,
        step: timedelta,
    ):
        frame = rates[(rates["cell"] >= 0) & (rates["cell"] < n_cells)]
        sat_index = {s.name: i for i, s in enumerate(satellites)}
        frame = frame[frame["satellite"].isin(sat_index)].reset_index(drop=True)
        values = frame[_ENDPOINT_COLUMNS].astype(object)
        values = values.where(values.notna(), None)
        keys = [_get_endpoint_key(row) for row in values.itertuples(index=False)]
        known = np.array([key in endpoints.index for key in keys], dtype=bool)
        frame = frame[known].reset_index(drop=True)
        endpoint = np.array(
            [endpoints.index[key] for key, ok in zip(keys, known) if ok], dtype=int
        )
        satellite = frame["satellite"].map(sat_index).to_numpy(dtype=int)
        transmitter = pd.factorize(
            pd.Series(list(zip(satellite, frame["transmitter"])), dtype=object)
        )[0]
        seconds = frame["fraction"].to_numpy(dtype=float) * step.total_seconds()
        # the time left in each contact of a transmitter from each cell
        remaining = (
            pd.Series(seconds[::-1])
            .groupby([frame["contact"].to_numpy()[::-1], transmitter[::-1]], sort=False)
            .cumsum()
            .to_numpy()[::-1]
        )
        order = np.lexsort(
            (-frame["rate"].to_numpy(dtype=float), frame["cell"].to_numpy())
        )
        self.frame = frame
        self.order = order
        self.bounds = np.searchsorted(
            frame["cell"].to_numpy()[order], np.arange(n_cells + 1)
        ).tolist()
        self.satellite = satellite[order].tolist()
        self.transmitter = transmitter[order].tolist()
        self.endpoint = endpoint[order].tolist()
        self.contact = frame["contact"].to_numpy()[order].tolist()
        self.rate = frame["rate"].to_numpy(dtype=float)[order].tolist()
        self.seconds = seconds[order].tolist()
        self.remaining = remaining[order].tolist()
        self.endpoints = endpoints

    def schedule(
        self,
        cell: int,
        wants: np.ndarray,
        waiting: np.ndarray,
        active: dict[int, int],
    ) -> list[tuple[int, float]]:
        """
        Schedules the links of a cell: links in progress (`active`, the contact
        of each transmitter) continue while their contacts last, then
        satellites that want links take free antennas, the longest `waiting`
        first (updated for those that request one), each its fastest link
        first, if its contact lasts the endpoint's minimum access time.

        Returns:
            list[tuple[int, float]]: The links taken (their rows, in sorted
                order) and their rates (bits per second).
        """
        first, last = self.bounds[cell], self.bounds[cell + 1]
        if last <= first:
            active.clear()
            return []
        endpoints = self.endpoints
        occupied: dict[int, int] = {}
        fed: dict[int, float] = {}
        busy: set[int] = set()
        grants: list[tuple[int, float]] = []

        def take(row: int) -> bool:
            endpoint = self.endpoint[row]
            rate = min(
                self.rate[row], endpoints.feeder_rates[endpoint] - fed.get(endpoint, 0)
            )
            if occupied.get(endpoint, 0) >= endpoints.antennas[endpoint] or rate <= 0:
                return False
            occupied[endpoint] = occupied.get(endpoint, 0) + 1
            fed[endpoint] = fed.get(endpoint, 0) + rate
            busy.add(self.transmitter[row])
            grants.append((row, rate))
            return True

        available = {
            (self.contact[r], self.transmitter[r]): r for r in range(first, last)
        }
        for transmitter, contact in list(active.items()):
            row = available.get((contact, transmitter))
            if row is None or not take(row):
                del active[transmitter]
        requests = [
            r
            for r in range(first, last)
            if self.transmitter[r] not in busy
            and wants[self.satellite[r]]
            and self.remaining[r] >= endpoints.min_access_times[self.endpoint[r]]
        ]
        for r in requests:
            waiting[self.satellite[r]] = min(waiting[self.satellite[r]], cell)
        requests.sort(
            key=lambda r: (waiting[self.satellite[r]], self.satellite[r], -self.rate[r])
        )
        for row in requests:
            if self.transmitter[row] not in busy and take(row):
                active[self.transmitter[row]] = self.contact[row]
        return grants


def _simulate(
    rates: pd.DataFrame,
    satellites: Sequence[Satellite],
    endpoints: _Endpoints,
    generation: np.ndarray,
    model: DataFlow,
) -> tuple[pd.DataFrame, dict[str, _Storage]]:
    """
    Simulates the links and storage of satellites cell by cell of the grid of
    time (see `compute_data_flow`).

    Args:
        rates (pandas.DataFrame): The rates of links in each cell (see
            `link_sampling.compute_link_rates`).
        satellites (Sequence[Satellite]): The satellites.
        endpoints (_Endpoints): The ground stations and relay satellites.
        generation (numpy.ndarray): The volume each satellite generates in each
            cell (shape: satellites by cells).
        model (DataFlow): The options of the analysis.

    Returns:
        tuple[pandas.DataFrame, dict[str, _Storage]]: The cells of the links
            taken, and the storage of each satellite.
    """
    n_sats, n_cells = generation.shape
    limits = np.array(
        [
            s.resources.storage if s.resources and s.resources.storage else np.inf
            for s in satellites
        ]
    )
    backlog = np.minimum(np.full(n_sats, float(model.initial_backlog)), limits)
    initial = backlog.copy()
    lost = np.zeros((n_sats, n_cells))
    delivered = np.zeros((n_sats, n_cells))
    stored = np.zeros((n_sats, n_cells))
    links = _Links(rates, satellites, endpoints, n_cells, model.step)

    active: dict[int, int] = {}  # the contact of each transmitter's link
    waiting = np.full(n_sats, np.inf)  # the cell since which each has waited
    taken: list[tuple[int, float, float]] = []
    for k in range(n_cells):
        # satellites want links while they store or generate more than the
        # minimum backlog (or always, holding links to the ends of contacts)
        wants = (
            backlog + generation[:, k] > model.min_backlog
            if model.release == "empty"
            else np.ones(n_sats, dtype=bool)
        )
        grants = links.schedule(k, wants, waiting, active)
        capacity = np.zeros(n_sats)
        for row, rate in grants:
            capacity[links.satellite[row]] += rate * links.seconds[row]
        waiting[(capacity > 0) | ~wants] = np.inf
        # first in, first out storage, losing data generated while it is full
        level = np.maximum(0.0, backlog + generation[:, k] - capacity)
        lost[:, k] = np.maximum(0.0, level - limits)
        level = np.minimum(level, limits)
        delivered[:, k] = np.maximum(
            0.0, backlog + generation[:, k] - lost[:, k] - level
        )
        stored[:, k] = backlog = level
        for row, rate in grants:
            s = links.satellite[row]
            # the volume delivered through each link, by its share of capacity
            share = rate * links.seconds[row] / capacity[s] if capacity[s] > 0 else 0
            taken.append((links.order[row], rate, share * delivered[s, k]))
            if model.release == "empty" and backlog[s] <= 0:
                # released once the satellite's storage empties
                active.pop(links.transmitter[row], None)
    storages = {
        s.name: _Storage(
            generated=generation[i],
            lost=lost[i],
            delivered=delivered[i],
            backlog=stored[i],
            initial_backlog=float(initial[i]),
        )
        for i, s in enumerate(satellites)
    }
    if not taken:
        return _empty_schedule(), storages
    rows, taken_rates, volumes = zip(*taken)
    schedule = links.frame.iloc[list(rows)][
        ["satellite", "transmitter", *_ENDPOINT_COLUMNS, "contact", "cell", "fraction"]
    ].reset_index(drop=True)
    schedule["rate"] = taken_rates
    schedule["capacity"] = (
        schedule["rate"] * schedule["fraction"] * model.step.total_seconds()
    )
    schedule["volume"] = volumes
    return schedule, storages


@dataclass(frozen=True)
class DataFlowResult:
    """
    The result of a data flow analysis (see `compute_data_flow`).

    Attributes:
        schedule (pandas.DataFrame): The links taken, as runs of consecutive
            cells: `satellite`, `transmitter`, `ground_station`, `relay`,
            `relay_satellite`, `start`, `end`, the `volume` (bits) delivered,
            and the `mean_rate` (bits per second) of the link.
        storage (pandas.DataFrame): The storage of each satellite at the end
            of each cell of the grid of time: `satellite`, `time`, `backlog`
            (bits), and the cumulative volumes `generated`, `delivered`, and
            `lost` (bits).
        summary (pandas.DataFrame): For each satellite: the volumes
            `generated`, `delivered`, and `lost` (bits), the `loss_fraction`,
            the `peak_backlog` and `final_backlog` (bits), and the `storage`
            (bits, missing without a limit).
        utilization (pandas.DataFrame): For each ground station and relay
            satellite that links: its `antennas` (missing without a limit),
            the `busy_fraction` of the antennas' time taken, and the `volume`
            (bits) delivered through it.
    """

    schedule: pd.DataFrame
    storage: pd.DataFrame
    summary: pd.DataFrame
    utilization: pd.DataFrame
    start: datetime
    step: timedelta
    _storages: dict[str, _Storage]
    _endpoints: dict[str, tuple[np.ndarray, np.ndarray]]

    @property
    def satellites(self) -> list[str]:
        """The names of the satellites whose data flow is analyzed."""
        return list(self._storages)

    def find(self, observations: pd.DataFrame) -> pd.DataFrame:
        """
        Finds the time each observation of an analyzed satellite reaches the
        ground (`downlinked`), through which `ground_station` (and `relay`),
        and whether it is `lost` (generated while the storage was full). Its
        data is the volume generated up to its `end`, delivered first in,
        first out.

        Args:
            observations (pandas.DataFrame): The observations, each naming its
                `satellite`, with its `end`.

        Returns:
            pandas.DataFrame: The `ground_station`, `relay`, `downlinked`, and
                `lost` of each observation, with the observations' index
                (missing for the observations of other satellites, and the
                downlink of those not delivered by the end of the analysis).
        """
        n = len(observations)
        ground_station = np.full(n, None, dtype=object)
        relay = np.full(n, None, dtype=object)
        downlinked = np.full(n, np.nan)
        lost = np.zeros(n, dtype=bool)
        if n and "satellite" in observations.columns and self._storages:
            step_s = self.step.total_seconds()
            ends = (
                (pd.to_datetime(observations["end"], utc=True) - _to_utc(self.start))
                .dt.total_seconds()
                .to_numpy(dtype=float)
            )
            names = observations["satellite"].astype(object).to_numpy()
            for name, storage in self._storages.items():
                rows = np.flatnonzero(names == name)
                if len(rows) and len(storage.generated):
                    found = _find_deliveries(storage, ends[rows], step_s)
                    downlinked[rows], lost[rows], cell = found
                    endpoints = self._endpoints[name]
                    delivered = ~np.isnan(downlinked[rows])
                    ground_station[rows[delivered]] = endpoints[0][cell[delivered]]
                    relay[rows[delivered]] = endpoints[1][cell[delivered]]
        times = _to_utc(self.start) + pd.to_timedelta(downlinked, unit="s")
        return pd.DataFrame(
            {
                "ground_station": ground_station,
                "relay": relay,
                "downlinked": pd.Series(times, dtype="datetime64[ns, UTC]").array,
                "lost": lost,
            },
            index=observations.index,
        )


def _find_deliveries(
    storage: _Storage, ends: np.ndarray, step_s: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Finds the delivery of the data generated up to times (seconds from the
    start of the grid of time) through storage: the time (seconds, `NaN` if
    lost or not delivered), whether it is lost, and the cell of its delivery.
    Volumes are linear within each cell.
    """
    n_cells = len(storage.generated)
    # the cell in which the data up to each time was generated
    cell = np.clip(np.ceil(ends / step_s).astype(int) - 1, 0, n_cells - 1)
    within = np.clip(ends / step_s - cell, 0, 1)
    # data generated after the storage fills within a cell is lost (tail drop)
    kept = storage.generated[cell] - storage.lost[cell]
    produced = storage.generated[cell] * within
    lost = (storage.lost[cell] > 0) & (produced > kept)
    accepted = np.r_[0, np.cumsum(storage.generated - storage.lost)]
    delivered = np.r_[0, np.cumsum(storage.delivered)]
    position = storage.initial_backlog + accepted[cell] + np.minimum(produced, kept)
    # the first time the cumulative volume delivered reaches the position
    index = np.clip(np.searchsorted(delivered, position, side="left"), 1, n_cells)
    low, high = delivered[index - 1], delivered[index]
    share = np.clip((position - low) / np.where(high > low, high - low, 1), 0, 1)
    time = np.maximum((index - 1 + share) * step_s, ends)
    done = ~lost & (position <= delivered[-1] * (1 + 1e-12))
    return np.where(done, time, np.nan), lost, index - 1


def _get_runs(schedule: pd.DataFrame, start: datetime, step: timedelta) -> pd.DataFrame:
    """
    Gets the runs of a schedule of links (see `_simulate`): consecutive
    cells of a contact and transmitter, with their `start`, `end`, `volume`,
    and `mean_rate`.
    """
    columns = ["satellite", "transmitter", *_ENDPOINT_COLUMNS]
    if schedule.empty:
        return pd.DataFrame(
            {
                **{column: pd.Series([], dtype=object) for column in columns},
                "start": pd.Series([], dtype="datetime64[ns, UTC]"),
                "end": pd.Series([], dtype="datetime64[ns, UTC]"),
                "volume": pd.Series([], dtype=float),
                "mean_rate": pd.Series([], dtype=float),
            }
        )
    frame = schedule.sort_values(["contact", "transmitter", "cell"], kind="stable")
    new_run = (
        (frame["contact"].diff() != 0)
        | (frame["transmitter"] != frame["transmitter"].shift())
        | (frame["cell"].diff() != 1)
    )
    frame = frame.assign(
        _run=new_run.cumsum(), _seconds=frame["fraction"] * step.total_seconds()
    )
    runs = frame.groupby("_run", sort=False).agg(
        **{column: (column, "first") for column in columns},
        _first=("cell", "first"),
        _last=("cell", "last"),
        volume=("volume", "sum"),
        _capacity=("capacity", "sum"),
        _seconds=("_seconds", "sum"),
    )
    reference = _to_utc(start)
    runs["start"] = reference + runs["_first"].to_numpy() * step
    runs["end"] = reference + (runs["_last"].to_numpy() + 1) * step
    runs["mean_rate"] = runs["_capacity"] / runs["_seconds"]
    runs = runs.sort_values(["start", "satellite"], kind="stable").reset_index(
        drop=True
    )
    return runs[[*columns, "start", "end", "volume", "mean_rate"]]


def _get_endpoint_by_cell(
    schedule: pd.DataFrame, n_cells: int
) -> tuple[np.ndarray, np.ndarray]:
    """
    Gets the ground station and relay of the fastest link of a satellite in
    each cell (missing in cells without one).
    """
    ground_station = np.full(n_cells, None, dtype=object)
    relay = np.full(n_cells, None, dtype=object)
    if len(schedule):
        best = schedule.sort_values("rate", ascending=False, kind="stable")
        best = best.drop_duplicates("cell")
        cells = best["cell"].to_numpy()
        ground_station[cells] = best["ground_station"].to_numpy()
        relay[cells] = best["relay"].to_numpy()
    return ground_station, relay


def _get_utilization(
    schedule: pd.DataFrame,
    ground_stations: Sequence[GroundStation],
    relays: Sequence[Relay],
    duration: float,
    step: timedelta,
) -> pd.DataFrame:
    """
    Gets the utilization of each ground station and relay satellite: its
    antennas, the fraction of their time taken, and the volume received.
    """
    seconds = schedule["fraction"].to_numpy() * step.total_seconds()
    is_relay = schedule["relay"].notna().to_numpy()
    members = []
    for station in ground_stations:
        taken = ~is_relay & (schedule["ground_station"] == station.name).to_numpy()
        antennas = station.communications.antennas if station.communications else None
        members.append((station.name, None, None, antennas, taken))
    for relay in relays:
        antennas = relay.communications.antennas if relay.communications else None
        for relay_satellite in relay.satellites:
            taken = (schedule["relay"] == relay.name).to_numpy() & (
                schedule["relay_satellite"] == relay_satellite.name
            ).to_numpy()
            members.append((None, relay.name, relay_satellite.name, antennas, taken))
    return pd.DataFrame(
        [
            {
                "ground_station": station,
                "relay": relay,
                "relay_satellite": relay_satellite,
                "antennas": antennas,
                "busy_fraction": (
                    seconds[taken].sum() / (antennas * duration)
                    if antennas and duration > 0
                    else np.nan
                ),
                "volume": float(schedule["volume"].to_numpy()[taken].sum()),
            }
            for station, relay, relay_satellite, antennas, taken in members
        ],
        columns=[
            "ground_station",
            "relay",
            "relay_satellite",
            "antennas",
            "busy_fraction",
            "volume",
        ],
    )


def is_data_flow_modeled(satellite: Satellite) -> bool:
    """
    Checks whether the data flow of a satellite can be modeled: whether it
    generates data at a known rate (see `link_sampling.get_generation_rate`)
    and has a transmitter (of its communications, or at its resources'
    `downlink_rate`).

    Args:
        satellite (Satellite): The satellite.

    Returns:
        bool: `True` if its data flow can be modeled
    """
    return get_generation_rate(satellite) is not None and bool(
        _get_terminals(satellite)
    )


def analyze_data_flow(
    contacts: pd.DataFrame,
    satellites: Sequence[Satellite],
    ground_stations: Sequence[GroundStation] = (),
    relays: Sequence[Relay] = (),
    *,
    start: datetime,
    end: datetime,
    model: DataFlow = DataFlow(),
    samples: pd.DataFrame | None = None,
    periods: Mapping[str, tuple[datetime, datetime]] | None = None,
) -> DataFlowResult:
    """
    Analyzes the flow of data from satellites through their contacts with
    ground stations and relays (see `compute_data_flow`), without the
    observations.

    Args:
        contacts (pandas.DataFrame): The contacts with ground stations and
            relays (see `collect_downlinks` and `collect_relay_contacts`).
        satellites (Sequence[Satellite]): The satellites, in order of priority
            among contacts that start together.
        ground_stations (Sequence[GroundStation]): The ground stations.
        relays (Sequence[Relay]): The relays.
        start (datetime.datetime): The start of the analysis.
        end (datetime.datetime): The end of the analysis.
        model (DataFlow): The options of the analysis.
        samples (pandas.DataFrame | None): The samples of the contacts'
            geometry on the model's grid of time (see `sample_link_geometry`),
            if already computed.
        periods (Mapping[str, tuple[datetime, datetime]] | None): The period
            in which each satellite is active (and generates data), by name,
            if not the whole analysis.

    Returns:
        DataFlowResult: The result.

    Raises:
        ValueError: If the data flow of a satellite cannot be modeled and the
            model is `strict`.
    """
    modeled = [s for s in satellites if is_data_flow_modeled(s)]
    others = [s.name for s in satellites if not is_data_flow_modeled(s)]
    if others:
        message = (
            f"the data flow of satellites {others} cannot be modeled (without a "
            + "rate of generation or a transmitter): their latencies are found "
            + "by timing alone, and they take no antennas"
        )
        if model.strict:
            raise ValueError(message)
        warn(message)
    names = {s.name for s in modeled}
    step = model.step
    lengths = _get_cell_lengths(start, end, step)
    n_cells = len(lengths)
    if samples is None:
        if len(contacts):
            contacts = contacts[contacts["satellite"].isin(names)]
        samples = sample_link_geometry(
            contacts, modeled, ground_stations, relays, start=start, step=step
        )
    elif len(samples):
        samples = samples[samples["satellite"].isin(names)]
    rates = compute_link_rates(samples, modeled, ground_stations, relays)
    generation = np.zeros((len(modeled), n_cells))
    for i, satellite in enumerate(modeled):
        generation[i] = compute_generation(satellite, start, end, step)
        if periods is not None and satellite.name in periods:
            # generated only while active (in cells whose middles are in its period)
            first, last = (
                (_to_utc(time) - _to_utc(start)).total_seconds()
                for time in periods[satellite.name]
            )
            middles = step.total_seconds() * np.arange(n_cells) + lengths / 2
            generation[i, (middles < first) | (middles >= last)] = 0
    schedule, storages = _simulate(
        rates, modeled, _Endpoints.of(ground_stations, relays), generation, model
    )
    endpoints = {}
    summaries = []
    series = []
    edges = _to_utc(start) + pd.to_timedelta(np.cumsum(lengths), unit="s")
    for satellite in modeled:
        storage = storages[satellite.name]
        endpoints[satellite.name] = _get_endpoint_by_cell(
            schedule[schedule["satellite"] == satellite.name], n_cells
        )
        generated = float(storage.generated.sum())
        summaries.append(
            {
                "satellite": satellite.name,
                "generated": generated,
                "delivered": float(storage.delivered.sum()),
                "lost": float(storage.lost.sum()),
                "loss_fraction": (
                    float(storage.lost.sum()) / generated if generated > 0 else 0.0
                ),
                "peak_backlog": float(
                    max(storage.initial_backlog, storage.backlog.max(initial=0))
                ),
                "final_backlog": float(
                    storage.backlog[-1] if n_cells else storage.initial_backlog
                ),
                "storage": (
                    satellite.resources.storage if satellite.resources else None
                ),
            }
        )
        series.append(
            pd.DataFrame(
                {
                    "satellite": satellite.name,
                    "time": edges,
                    "backlog": storage.backlog,
                    "generated": np.cumsum(storage.generated),
                    "delivered": np.cumsum(storage.delivered),
                    "lost": np.cumsum(storage.lost),
                }
            )
        )
    storage_columns = ["satellite", "time", "backlog", "generated", "delivered", "lost"]
    return DataFlowResult(
        schedule=_get_runs(schedule, start, step),
        storage=(
            pd.concat(series, ignore_index=True)
            if series
            else pd.DataFrame(columns=storage_columns)
        ),
        summary=pd.DataFrame(
            summaries,
            columns=[
                "satellite",
                "generated",
                "delivered",
                "lost",
                "loss_fraction",
                "peak_backlog",
                "final_backlog",
                "storage",
            ],
        ),
        utilization=_get_utilization(
            schedule, ground_stations, relays, float(lengths.sum()), step
        ),
        start=start,
        step=step,
        _storages=storages,
        _endpoints=endpoints,
    )


def find_data_flow_downlinks(
    observations: pd.DataFrame, contacts: pd.DataFrame, result: DataFlowResult
) -> pd.DataFrame:
    """
    Finds the downlink of each observation: through the data flow of its
    satellite, if analyzed (see `DataFlowResult.find`), or else by timing
    alone (with `during_contact="end"`, see `compute_latencies`).

    Args:
        observations (pandas.DataFrame): The observations, each with its time
            (`epoch`) and `end`.
        contacts (pandas.DataFrame): The contacts with ground stations and
            relays.
        result (DataFlowResult): The data flow.

    Returns:
        pandas.DataFrame: The `ground_station`, `relay`, `downlinked`, and
            `lost` of each observation, with the observations' index.
    """
    found = result.find(observations)
    if "satellite" not in observations.columns or len(contacts) == 0:
        return found
    timed = ~observations["satellite"].isin(result.satellites).to_numpy()
    if np.any(timed):
        contacts = contacts.copy()
        for column in ["relay", "relay_satellite"]:
            if column not in contacts.columns:
                contacts[column] = None
        alone = _find_downlinks(observations[timed], contacts, "end")
        for column in ["ground_station", "relay", "downlinked"]:
            found.loc[timed, column] = alone[column].to_numpy()
    return found


def compute_data_flow(
    observations: gpd.GeoDataFrame,
    contacts: pd.DataFrame,
    satellites: Sequence[Satellite],
    ground_stations: Sequence[GroundStation] = (),
    relays: Sequence[Relay] = (),
    *,
    start: datetime,
    end: datetime,
    model: DataFlow = DataFlow(),
) -> tuple[gpd.GeoDataFrame, DataFlowResult]:
    """
    Computes the latencies of observations through the flow of data from
    satellites to the ground, rather than by timing alone (see
    `compute_latencies`). Each satellite generates data (see
    `compute_generation`), stores it, and sends it first in, first out over
    links through its contacts with ground stations and relays (see
    `collect_downlinks` and `collect_relay_contacts`), at rates from link
    budgets (see `link_sampling.compute_link_rates`). Data generated while its
    storage (its resources' `storage`) is full is lost. An observation
    reaches the ground when the data generated up to its end is delivered.

    Links are scheduled cell by cell of the grid of time. Each transmitter
    sends over one link at a time, and each ground station and relay
    satellite receives over at most its `antennas` at a time (any number
    without communications); a relay satellite forwards at most its
    `feeder_rate`, shared by its links, and stores nothing. A link continues
    while its contact lasts, until the satellite's storage empties (with the
    model's `release="empty"`, after which it requests a link again once more
    than `min_backlog` is stored) or to the end of the contact (with
    `release="end"`). Satellites that request links take free antennas in
    the order they began to wait (then in the order given), each its fastest
    link first, if its contact lasts at least the endpoint's
    `min_access_time`.

    The analysis is on a grid of time from `start` with the model's `step`,
    which sets the accuracy of times, starting with the model's
    `initial_backlog` stored on each satellite (an analysis of at least a day
    is less sensitive to it). The latencies of satellites whose data flow
    cannot be modeled (see `is_data_flow_modeled`) are found by timing alone,
    unless the model is `strict`.

    Args:
        observations (geopandas.GeoDataFrame): The observations of each
            satellite (see `compute_latencies`).
        contacts (pandas.DataFrame): The contacts with ground stations and
            relays.
        satellites (Sequence[Satellite]): The satellites, in order of priority
            among contacts that start together.
        ground_stations (Sequence[GroundStation]): The ground stations.
        relays (Sequence[Relay]): The relays.
        start (datetime.datetime): The start of the analysis.
        end (datetime.datetime): The end of the analysis.
        model (DataFlow): The options of the analysis.

    Returns:
        tuple[geopandas.GeoDataFrame, DataFlowResult]: The latencies (with the
            columns of `compute_latencies`, the `relay` of each downlink, and
            whether each observation is `lost`), and the data flow.
    """
    result = analyze_data_flow(
        contacts,
        satellites,
        ground_stations,
        relays,
        start=start,
        end=end,
        model=model,
    )
    empty = _get_empty_latency_frame(
        _get_platform_columns(observations),
        angles="platform_alt" in observations.columns or len(observations.columns) == 0,
    )
    columns = list(empty.columns)
    columns.insert(columns.index("ground_station") + 1, "relay")
    columns.append("lost")
    if observations.empty:
        empty = empty.assign(
            relay=pd.Series([], dtype=object), lost=pd.Series([], dtype=bool)
        )
        return empty[columns], result
    obs = observations.copy()
    obs["epoch"] = _get_observation_times(observations).dt.as_unit("ns").to_numpy()
    if "end" not in obs.columns:
        obs["end"] = obs["epoch"]
    found = find_data_flow_downlinks(obs, contacts, result)
    obs = obs.assign(**{column: found[column] for column in found.columns})
    obs["latency"] = obs["downlinked"] - obs["epoch"]
    obs = obs.rename(columns={"epoch": "observed"})
    obs = gpd.GeoDataFrame(
        obs[[column for column in columns if column in obs]],
        geometry="geometry",
        crs=(
            observations.crs
            if isinstance(observations, gpd.GeoDataFrame) and observations.crs
            else "EPSG:4326"
        ),
    )
    return obs.sort_values(by="observed").reset_index(drop=True), result

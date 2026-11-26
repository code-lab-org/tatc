"""
Evaluations of mission architectures (see `tatc.schemas.Architecture`): the
observations of targets by platforms of every kind, the contacts of
satellites with ground stations and relays, the latency of each observation,
and metrics of the observations, with the value of each member measured by
leaving it out or adding it.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import dataclasses
import numbers
from collections.abc import Callable, Iterable, Mapping
from datetime import datetime, timedelta
from typing import Any, Literal, Protocol

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import geometry as geo

from ..schemas import (
    Aircraft,
    AllPlatforms,
    GroundStation,
    Point,
    RadarStation,
    Satellite,
)
from .._warnings import warn
from ..schemas.architecture import ActiveWindow, Architecture, Member, Relay, _field_of
from ..schemas.communication import DataFlow
from .check import (
    _PLATFORM_KINDS,
    _check_time_window,
    _get_observation_times,
    _get_platform_kind,
)
from .data_flow import DataFlowResult, analyze_data_flow, find_data_flow_downlinks
from .footprints import (
    _get_empty_footprint_observations,
    collect_footprint_observations,
    collect_footprints,
)
from .latency_sampling import (
    _find_downlinks,
    _get_empty_relay_contacts_frame,
    collect_downlinks,
    collect_relay_contacts,
)
from .point_sampling import collect_observations
from .region_sampling import collect_region_observations
from .sampling import _get_empty_coverage_frame

_KINDS = {
    "satellites": "satellite",
    "aircraft": "aircraft",
    "surface_platforms": "surface_platform",
    "radar_stations": "radar_station",
    "ground_stations": "ground_station",
    "relays": "relay",
}
"""The kind of the members of each field of an architecture, as named in the
results (the column of each kind of platform, see `_PLATFORM_KINDS`)."""

DATA_PATH_COLUMNS = ("ground_station", "relay", "downlinked", "latency")
"""The columns of the downlink of each observation of an architecture with
ground stations or relays (see `evaluate_architecture`)."""

DATA_FLOW_COLUMNS = ("lost",)
"""The columns of the data flow of each observation of an architecture
evaluated with a data flow model, besides its downlink (see
`evaluate_architecture`)."""


def _kind_of(member: Member) -> str:
    """Gets the kind of a member, as named in results (such as `satellite`)."""
    return _KINDS[_field_of(member)]


def _split_targets(targets: Any) -> tuple[list, list]:
    """
    Splits targets into points (TAT-C or shapely points) and regions
    (shapely polygons and multipolygons), as `collect_observations` and
    `collect_region_observations` take them.
    """
    if isinstance(targets, gpd.GeoDataFrame):
        targets = targets.geometry
    if isinstance(targets, gpd.GeoSeries):
        targets = list(targets)
    targets = list(targets) if isinstance(targets, (list, tuple)) else [targets]
    points, regions = [], []
    for target in targets:
        if isinstance(target, (Point, geo.Point)):
            points.append(target)
        elif isinstance(target, (geo.Polygon, geo.MultiPolygon)):
            regions.append(target)
        else:
            raise TypeError(
                "target must be a Point, shapely Point, Polygon, or MultiPolygon, "
                f"not a {type(target).__name__}"
            )
    return points, regions


def _concat(
    frames: list[pd.DataFrame], empty: Callable[[], gpd.GeoDataFrame]
) -> gpd.GeoDataFrame:
    """
    Concatenates the frames that are not empty, sorted (stably) by start, or
    gets an empty frame if all are.
    """
    frames = [frame for frame in frames if len(frame) > 0]
    if not frames:
        return empty()
    if len(frames) == 1:
        frame = frames[0]
    else:
        frame = pd.concat(frames, ignore_index=True)
    frame = frame.sort_values("start", kind="stable").reset_index(drop=True)
    if not isinstance(frame, gpd.GeoDataFrame):
        frame = gpd.GeoDataFrame(frame, geometry="geometry", crs="EPSG:4326")
    return frame


def _sample_times(
    platform: AllPlatforms,
    start: datetime,
    end: datetime,
    time_step: timedelta | None,
) -> list[datetime] | None:
    """
    Gets the times at which to sample the footprints of a platform other
    than a satellite over a period (see `collect_platform_observations`), or
    `None` for all those of an aircraft's trajectory.
    """
    if time_step is None:
        if isinstance(platform, Aircraft):
            return None
        return [start, end]
    times = list(pd.date_range(start, end, freq=time_step).to_pydatetime())
    if isinstance(platform, Aircraft):
        contained = np.atleast_1d(platform.trajectory.contains(times))
        times = [time for time, inside in zip(times, contained) if inside]
    return times


def collect_platform_observations(
    targets: Any,
    platform: AllPlatforms,
    start: datetime,
    end: datetime,
    *,
    time_step: timedelta | None = None,
    instrument_index: int | None = 0,
    elevation: float | Mapping[str, float | None] | None = None,
) -> gpd.GeoDataFrame:
    """
    Collect the observations of points and regions by one observing platform
    of any kind over a period, with the columns of each kind's observations:
    of a satellite, from its orbit (see `collect_observations` for points and
    `collect_region_observations` for regions); of an aircraft, a fixed
    surface platform, or a radar station, from its footprints (see
    `collect_footprints` and `collect_footprint_observations`, with its
    angles seen from points): those at each `time_step` from the start of
    the period (within an aircraft's trajectory), or by default the times of
    an aircraft's trajectory within the period, or the start and end of the
    period for a fixed platform, whose footprint does not change, so that it
    observes a target in view throughout the period once.

    It is the analysis of each platform of an architecture (see
    `evaluate_architecture`) and of a trade study's observations (see
    `tatc.tradespace.Observations`).

    Args:
        targets (Any): The point(s) and region(s) of interest (see
            `collect_footprint_observations`).
        platform (AllPlatforms): The observing platform.
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.
        time_step (datetime.timedelta | None): The time between the samples
            of the footprints of a platform other than a satellite, if any.
        instrument_index (int | None): The index of the observing instrument
            on the platform, or `None` for every instrument (a radar station
            has none).
        elevation (float | Mapping[str, float | None] | None): The elevation
            (meters) above the WGS 84 ellipsoid of the footprints of
            platforms other than satellites, or that of each kind (see
            `collect_footprints`): required for a radar station.

    Returns:
        geopandas.GeoDataFrame: The observations, sorted by start time.
    """
    # pylint: disable=too-many-locals,too-many-arguments
    kind = _get_platform_kind(platform)
    if kind is None:
        names = ", ".join(kind.type.__name__ for kind in _PLATFORM_KINDS)
        raise TypeError(f"platform must be a {names}, not a {type(platform).__name__}")
    _check_time_window(start, end)
    if isinstance(platform, Satellite):
        points, regions = _split_targets(targets)
        frames = []
        if points:
            frames.append(
                collect_observations(points, platform, start, end, instrument_index)
            )
        if regions:
            frames.append(
                collect_region_observations(
                    regions, platform, start, end, instrument_index
                )
            )
        return _concat(
            frames,
            lambda: _get_empty_coverage_frame(True),
        )
    if isinstance(platform, RadarStation):
        instrument_index = 0
    times = _sample_times(platform, start, end, time_step)
    empty = _get_empty_footprint_observations([kind.column], True)
    if times is not None and len(times) == 0:
        return empty
    footprints = collect_footprints(
        platform,
        times,
        instrument_index=instrument_index,
        elevation=elevation,
        time_step=time_step,
    )
    if times is None and not footprints.empty:
        # the times of an aircraft's trajectory within the period
        time = pd.to_datetime(footprints["time"], utc=True)
        footprints = footprints[
            (time >= pd.Timestamp(start)) & (time <= pd.Timestamp(end))
        ]
    if footprints.empty:
        return empty
    return collect_footprint_observations(targets, footprints, platforms=platform)


class _Analyses(Protocol):
    """
    The analyses of the members of an architecture (see
    `ArchitectureEvaluation`): the observations of a platform, and the
    contacts of a satellite with ground stations and relays, each over a
    period.
    """

    def observe(
        self, platform: AllPlatforms, start: datetime, end: datetime
    ) -> pd.DataFrame:
        """Gets the observations of the targets by a platform over a period."""
        raise NotImplementedError

    def contact(
        self,
        satellite: Satellite,
        ground_stations: list[GroundStation],
        relays: list[Relay],
        start: datetime,
        end: datetime,
    ) -> pd.DataFrame:
        """Gets the contacts of a satellite with ground stations and relays."""
        raise NotImplementedError


def _collect_contacts(
    satellite: Satellite,
    ground_stations: list[GroundStation],
    relays: list[Relay],
    start: datetime,
    end: datetime,
) -> pd.DataFrame:
    """
    Collects the contacts of a satellite with ground stations (see
    `collect_downlinks`) and relays (see `collect_relay_contacts`) over a
    period, with the columns of relay contacts (no relay for a ground
    station's).
    """
    frames = []
    if ground_stations:
        downlinks = collect_downlinks(ground_stations, satellite, start, end)
        if len(downlinks) > 0:
            frames.append(downlinks.assign(relay=None, relay_satellite=None))
    if relays:
        frames.append(collect_relay_contacts(relays, satellite, start, end))
    columns = list(_get_empty_relay_contacts_frame().columns)
    contacts = _concat(frames, _get_empty_relay_contacts_frame)
    return contacts[columns]


@dataclasses.dataclass(frozen=True)
class _Direct:
    """
    The analyses of the members of an architecture, without a cache (see
    `evaluate_architecture`).
    """

    targets: Any
    time_step: timedelta | None
    instrument_index: int | None
    elevation: float | Mapping[str, float | None] | None

    def observe(
        self, platform: AllPlatforms, start: datetime, end: datetime
    ) -> pd.DataFrame:
        """Gets the observations of the targets by a platform over a period."""
        return collect_platform_observations(
            self.targets,
            platform,
            start,
            end,
            time_step=self.time_step,
            instrument_index=self.instrument_index,
            elevation=self.elevation,
        )

    def contact(
        self,
        satellite: Satellite,
        ground_stations: list[GroundStation],
        relays: list[Relay],
        start: datetime,
        end: datetime,
    ) -> pd.DataFrame:
        """Gets the contacts of a satellite with ground stations and relays."""
        # pylint: disable=too-many-arguments,too-many-positional-arguments
        return _collect_contacts(satellite, ground_stations, relays, start, end)


def _clip_contacts(
    contacts: pd.DataFrame,
    window: ActiveWindow,
    min_access_time: timedelta,
) -> pd.DataFrame:
    """
    Clips contacts to the window in which their ground station or relay is
    active, keeping those that last at least its minimum access time, with
    each epoch at the midpoint of the part kept.
    """
    if window.start is None and window.end is None:
        return contacts
    start, end = contacts["start"], contacts["end"]
    if window.start is not None:
        start = start.clip(lower=pd.Timestamp(window.start))
    if window.end is not None:
        end = end.clip(upper=pd.Timestamp(window.end))
    duration = end - start
    kept = (duration > timedelta(0)) & (duration >= min_access_time)
    clipped = contacts.assign(start=start, end=end, epoch=start + duration / 2)
    return clipped[kept.to_numpy()]


@dataclasses.dataclass(frozen=True)
class _Parts:
    """
    The analyses of the members of an architecture over a period, from which
    the evaluations of it and of other architectures of its members (and
    others) are assembled (see `ArchitectureEvaluation`).

    Attributes:
        analyses (_Analyses): The analyses of members.
        start (datetime): The start of the period.
        end (datetime): The end of the period.
        observed (dict): The observations of each platform (by name and
            content: a platform that changes is analyzed again).
        contacts (dict): The contacts of each satellite (by name and content)
            with each ground station and relay (by name and content), over
            the part of the period in which the satellite is active.
    """

    analyses: _Analyses
    start: datetime
    end: datetime
    observed: dict = dataclasses.field(default_factory=dict)
    contacts: dict = dataclasses.field(default_factory=dict)

    def observations_of(
        self, architecture: Architecture, platform: AllPlatforms
    ) -> pd.DataFrame | None:
        """
        Gets the observations of a platform of an architecture over the part
        of the period in which it is active (analyzed if needed), or `None`
        if it is not active in the period.
        """
        period = architecture.get_active_period(platform.name, self.start, self.end)
        if period is None:
            return None
        key = (platform.name, _geometric_content(platform), period)
        if key not in self.observed:
            self.observed[key] = self.analyses.observe(platform, *period)
        return self.observed[key]

    def contacts_of(
        self, architecture: Architecture, satellite: Satellite
    ) -> list[pd.DataFrame]:
        """
        Gets the contacts of a satellite of an architecture with its ground
        stations and relays (each analyzed if needed, all those of the
        satellite that are not together), clipped to the windows in which
        the ground stations and relays are active.
        """
        period = architecture.get_active_period(satellite.name, self.start, self.end)
        if period is None:
            return []
        endpoints = [*architecture.ground_stations, *architecture.relays]
        keys = {
            id(endpoint): (
                satellite.name,
                _geometric_content(satellite),
                period,
                endpoint.name,
                _geometric_content(endpoint),
            )
            for endpoint in endpoints
        }
        missing = [e for e in endpoints if keys[id(e)] not in self.contacts]
        if missing:
            stations = [e for e in missing if isinstance(e, GroundStation)]
            relays = [e for e in missing if isinstance(e, Relay)]
            contacts = self.analyses.contact(satellite, stations, relays, *period)
            for endpoint in missing:
                if isinstance(endpoint, Relay):
                    mine = contacts["relay"].eq(endpoint.name)
                else:
                    mine = contacts["relay"].isna() & contacts["ground_station"].eq(
                        endpoint.name
                    )
                self.contacts[keys[id(endpoint)]] = contacts[mine.to_numpy()]
        return [
            _clip_contacts(
                self.contacts[keys[id(endpoint)]],
                architecture.get_window(endpoint.name),
                endpoint.min_access_time,
            )
            for endpoint in endpoints
        ]


def _content(model: Any) -> str:
    """Gets a key of a model's content (its JSON)."""
    return model.model_dump_json()


def _geometric_content(model: Any) -> str:
    """
    Gets a key of the content of a member that sets its observations or
    contacts (its JSON without its communications and resources, which do not
    change them, nor those of a relay's satellites and ground stations).
    """
    exclude: dict[str, Any] = {
        name: True
        for name in ("communications", "resources")
        if name in type(model).model_fields
    }
    if isinstance(model, Relay):
        exclude["satellites"] = {"__all__": {"communications", "resources"}}
        exclude["ground_stations"] = {"__all__": {"communications"}}
    return model.model_dump_json(exclude=exclude)


def _assemble(
    architecture: Architecture,
    parts: _Parts,
    during_contact: Literal["next", "end", "immediate"],
    data_path: bool,
    data_flow: DataFlow | None = None,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, DataFlowResult | None]:
    """
    Assembles the observations of an architecture's platforms and the
    downlinks of its satellites from the analyses of its members (see
    `_Parts`), with the downlink of each observation if it has a
    `data_path`: by timing alone, or through the flow of data of its
    satellites with a `data_flow` model (see `compute_data_flow`), whose
    result is also given.
    """
    # pylint: disable=too-many-locals
    frames = []
    for platform in architecture.platforms:
        observed = parts.observations_of(architecture, platform)
        if observed is not None:
            frames.append(observed)
    columns = [
        kind.column
        for kind in _PLATFORM_KINDS
        if any(_get_platform_kind(p) is kind for p in architecture.platforms)
    ]
    observations = _concat(
        frames,
        lambda: _get_empty_footprint_observations(columns or ["satellite"], True),
    )
    downlinks = _concat(
        [
            contacts
            for satellite in architecture.satellites
            for contacts in parts.contacts_of(architecture, satellite)
        ],
        _get_empty_relay_contacts_frame,
    )
    if data_path:
        observed = observations.assign(
            epoch=_get_observation_times(observations).dt.as_unit("ns").to_numpy()
        )
        if "end" not in observed.columns:
            observed["end"] = observed["epoch"]
        if data_flow is None:
            found = _find_downlinks(observed, downlinks, during_contact)
        else:
            periods = {
                satellite.name: period
                for satellite in architecture.satellites
                if (
                    period := architecture.get_active_period(
                        satellite.name, parts.start, parts.end
                    )
                )
                is not None
            }
            flow = analyze_data_flow(
                downlinks,
                [s for s in architecture.satellites if s.name in periods],
                architecture.ground_stations,
                architecture.relays,
                start=parts.start,
                end=parts.end,
                model=data_flow,
                periods=periods,
            )
            found = find_data_flow_downlinks(observed, downlinks, flow)
        observations = observations.assign(
            ground_station=found["ground_station"],
            relay=found["relay"],
            downlinked=found["downlinked"],
            latency=found["downlinked"] - observed["epoch"],
        )
        if data_flow is not None:
            observations = observations.assign(lost=found["lost"].to_numpy())
            return observations, downlinks, flow
    return observations, downlinks, None


def _members_of(values: Iterable[Any]) -> list[Member]:
    """Flattens members given alone or in lists."""
    flat = []
    for value in values:
        if isinstance(value, (list, tuple)):
            flat.extend(value)
        else:
            flat.append(value)
    return flat


def _contribution(value: Any, other: Any) -> Any:
    """
    Gets the difference of two values of a metric (`value - other`), if both
    are numbers (or times), or else `None`.
    """
    if isinstance(value, bool) or isinstance(other, bool):
        return None
    if value is None or other is None:
        return None
    if isinstance(value, (numbers.Number, np.number, timedelta, pd.Timedelta)):
        try:
            return value - other
        except TypeError:
            return None
    return None


@dataclasses.dataclass(frozen=True, eq=False)
class ArchitectureEvaluation:  # pylint: disable=too-many-instance-attributes
    """
    The evaluation of an architecture over a period (see
    `evaluate_architecture`): the observations of the targets by its
    platforms, the downlinks of its satellites, and the metrics of the
    observations. The value of each member is measured by evaluating the
    architecture without it (`leave_one_out`) and of candidates by
    evaluating it with each (`add_one`), from the same analyses of its
    members (each member is analyzed once).

    Attributes:
        architecture (Architecture): The architecture.
        start (datetime): The start of the period.
        end (datetime): The end of the period.
        observations (geopandas.GeoDataFrame): The observations of every
            platform (each over the part of the period in which it is
            active), sorted by start time, with the platform columns of each
            kind (see `collect_platform_observations`) and, if the
            architecture has ground stations or relays, the downlink of
            each (see `DATA_PATH_COLUMNS`): its `ground_station`, `relay`,
            the time it is `downlinked`, and its `latency` (missing without
            one, as for the observations of platforms other than satellites,
            see `tatc.analysis.compute_latencies`), and with a data flow
            model, whether each is `lost` (see `DATA_FLOW_COLUMNS`).
        downlinks (geopandas.GeoDataFrame): The contacts of the satellites
            with the ground stations (see `collect_downlinks`) and the relays
            (see `collect_relay_contacts`), with the columns of relay
            contacts (no relay for a ground station's), each over the part
            of the period in which both are active.
        metrics (dict | None): The metrics of the observations, if any.
        flow (DataFlowResult | None): The flow of data of the satellites to
            the ground stations and relays, with a data flow model (see
            `compute_data_flow`).
    """

    architecture: Architecture
    start: datetime
    end: datetime
    observations: gpd.GeoDataFrame
    downlinks: gpd.GeoDataFrame
    metrics: dict | None
    _: dataclasses.KW_ONLY
    metric: Callable[[pd.DataFrame], Mapping[str, Any]] | None = dataclasses.field(
        default=None, repr=False
    )
    during_contact: Literal["next", "end", "immediate"] = "end"
    data_path: bool = False
    data_flow: DataFlow | None = None
    flow: DataFlowResult | None = dataclasses.field(default=None, repr=False)
    _parts: _Parts | None = dataclasses.field(default=None, repr=False)

    def _of(self, architecture: Architecture) -> ArchitectureEvaluation:
        """
        Evaluates another architecture over the period, with the analyses of
        this one's members (and of others as needed), with the columns of the
        downlinks of observations if this evaluation has them (or the other
        architecture has ground stations or relays).
        """
        data_path = self.data_path or bool(
            architecture.ground_stations or architecture.relays
        )
        observations, downlinks, flow = _assemble(
            architecture, self._parts, self.during_contact, data_path, self.data_flow
        )
        return ArchitectureEvaluation(
            architecture,
            self.start,
            self.end,
            observations,
            downlinks,
            None if self.metric is None else dict(self.metric(observations)),
            metric=self.metric,
            during_contact=self.during_contact,
            data_path=data_path,
            data_flow=self.data_flow,
            flow=flow,
            _parts=self._parts,
        )

    def without_members(self, *names: str | Iterable[str]) -> ArchitectureEvaluation:
        """
        Evaluates the architecture without members (see
        `Architecture.without_members`), from the analyses of this
        evaluation.

        Args:
            *names (str | Iterable[str]): The names of the members (or lists
                of them).

        Returns:
            ArchitectureEvaluation: the evaluation
        """
        return self._of(self.architecture.without_members(*names))

    def with_members(
        self,
        *members: Member | Iterable[Member],
        windows: Mapping[str, ActiveWindow] | None = None,
    ) -> ArchitectureEvaluation:
        """
        Evaluates the architecture with other members (see
        `Architecture.with_members`), analyzing only the new members (and
        the contacts of the satellites with new ground stations and relays).

        Args:
            *members (Member | Iterable[Member]): The members (or lists of
                them).
            windows (Mapping[str, ActiveWindow] | None): The windows of the
                new members, if any.

        Returns:
            ArchitectureEvaluation: the evaluation
        """
        return self._of(self.architecture.with_members(*members, windows=windows))

    def _check_metrics(self) -> None:
        if self.metric is None:
            raise ValueError(
                "the value of members is measured by metrics: evaluate the "
                + "architecture with metrics"
            )

    def _compare(
        self, evaluations: dict[str, tuple[str, ArchitectureEvaluation]], sign: int
    ) -> pd.DataFrame:
        """
        Tabulates the metrics of other evaluations (by label, with a kind)
        and their contributions: the differences from this one's (or from
        theirs, with a negative sign).
        """
        rows = []
        for label, (kind, evaluation) in evaluations.items():
            row: dict[str, Any] = {"member": label, "kind": kind}
            row.update(evaluation.metrics)
            for name, value in evaluation.metrics.items():
                base = self.metrics.get(name)
                row[f"{name} contribution"] = (
                    _contribution(value, base)
                    if sign > 0
                    else _contribution(base, value)
                )
            rows.append(row)
        names = list(self.metrics)
        columns = ["member", "kind", *names, *(f"{n} contribution" for n in names)]
        frame = pd.DataFrame(rows, columns=columns if not rows else None)
        if rows:
            ordered = [c for c in columns if c in frame.columns]
            frame = frame[ordered + [c for c in frame.columns if c not in ordered]]
        return frame.set_index("member")

    def leave_one_out(
        self, members: Iterable[str] | Mapping[str, str | Iterable[str]] | None = None
    ) -> pd.DataFrame:
        """
        Measures the value of members of the architecture by leaving each
        out: the metrics of the architecture without it and its
        contribution to each numeric metric (the metric with it, this
        evaluation's, minus that without it), such as of a satellite, a
        ground station, a relay, or a group of members (such as a radar
        network, or a constellation's satellites). The architecture without
        each is assembled from the analyses of this evaluation.

        Args:
            members (Iterable[str] | Mapping[str, str | Iterable[str]] | None):
                The names of the members to leave out, or groups of members
                (each by a label) to leave out together, or `None` for each
                member.

        Returns:
            pandas.DataFrame: a row of each member (or group, of kind
                `group`) by name (or label), `member`, with its `kind`, the
                metrics without it, and the contribution to each metric
                (`<metric> contribution`, missing for metrics that are not
                numbers or times)

        Raises:
            ValueError: if the evaluation has no metrics.
            KeyError: if a name is not of a member.
        """
        self._check_metrics()
        if members is None:
            members = [member.name for member in self.architecture.members]
        if isinstance(members, Mapping):
            groups = {
                str(label): [names] if isinstance(names, str) else list(names)
                for label, names in members.items()
            }
        else:
            groups = {name: [name] for name in members}
        evaluations = {}
        for label, names in groups.items():
            kinds = {_kind_of(self.architecture.get_member(name)) for name in names}
            kind = kinds.pop() if len(kinds) == 1 and len(names) == 1 else "group"
            evaluations[label] = (kind, self.without_members(names))
        return self._compare(evaluations, -1)

    def add_one(
        self, candidates: Iterable[Member] | Mapping[str, Member | Iterable[Member]]
    ) -> pd.DataFrame:
        """
        Measures the value of candidates to add to the architecture by adding
        each: the metrics of the architecture with it and its contribution to
        each numeric metric (the metric with it minus this evaluation's),
        such as of a satellite, a ground station, a relay, or a group of
        members (such as a constellation's satellites). Only the candidates
        (and the contacts of satellites with new ground stations and relays)
        are analyzed, each once.

        Args:
            candidates (Iterable[Member] | Mapping[str, Member | Iterable[Member]]):
                The members to add, each by its name, or groups of members
                (each by a label) to add together. A candidate active only in
                a window is given in a group by a mapping of the window
                (`{"label": (member, ActiveWindow(...))}`).

        Returns:
            pandas.DataFrame: a row of each candidate (or group, of kind
                `group`) by name (or label), `member`, with its `kind`, the
                metrics with it, and the contribution to each metric
                (`<metric> contribution`, missing for metrics that are not
                numbers or times)

        Raises:
            ValueError: if the evaluation has no metrics.
        """
        self._check_metrics()
        if isinstance(candidates, Mapping):
            groups = dict(candidates)
        else:
            groups = {member.name: member for member in candidates}
        evaluations = {}
        for label, value in groups.items():
            windows = None
            if (
                isinstance(value, tuple)
                and len(value) == 2
                and isinstance(value[1], ActiveWindow)
            ):
                value, window = value
                windows = {member.name: window for member in _members_of([value])}
            members = _members_of([value])
            kinds = {_kind_of(member) for member in members}
            kind = kinds.pop() if len(kinds) == 1 and len(members) == 1 else "group"
            evaluations[str(label)] = (
                kind,
                self.with_members(members, windows=windows),
            )
        return self._compare(evaluations, 1)


def _evaluate(
    architecture: Architecture,
    parts: _Parts,
    metrics: Callable[[pd.DataFrame], Mapping[str, Any]] | None,
    during_contact: Literal["next", "end", "immediate"],
    data_flow: DataFlow | None = None,
) -> ArchitectureEvaluation:
    """
    Evaluates an architecture from the analyses of its members (see
    `_Parts`), such as cached analyses (see `tatc.tradespace.Observations`).
    """
    if during_contact not in ("next", "end", "immediate"):
        raise ValueError(
            f"during_contact must be 'next', 'end', or 'immediate', not {during_contact!r}"
        )
    if data_flow is not None and during_contact != "end":
        warn(
            "during_contact is ignored with a data flow model, through which "
            + "data is delivered as links allow"
        )
    seed = ArchitectureEvaluation(
        architecture,
        parts.start,
        parts.end,
        _get_empty_footprint_observations(["satellite"], True),
        _get_empty_relay_contacts_frame(),
        None,
        metric=metrics,
        during_contact=during_contact,
        data_path=False,
        data_flow=data_flow,
        _parts=parts,
    )
    return seed._of(architecture)  # pylint: disable=protected-access


def evaluate_architecture(
    architecture: Architecture,
    targets: Any,
    start: datetime,
    end: datetime,
    metrics: Callable[[pd.DataFrame], Mapping[str, Any]] | None = None,
    *,
    time_step: timedelta | None = None,
    instrument_index: int | None = 0,
    elevation: float | Mapping[str, float | None] | None = None,
    during_contact: Literal["next", "end", "immediate"] = "end",
    data_flow: DataFlow | None = None,
) -> ArchitectureEvaluation:
    """
    Evaluates an architecture over a period: the observations of targets by
    each of its platforms over the part of the period in which it is active
    (see `collect_platform_observations`), the contacts of its satellites
    with its ground stations and relays while both are active (see
    `collect_downlinks` and `collect_relay_contacts`), the downlink of each
    observation (see `compute_latencies`), and the metrics of the
    observations. The value of each member is then measured by leaving it
    out, and of candidates by adding each (see
    `ArchitectureEvaluation.leave_one_out` and
    `ArchitectureEvaluation.add_one`), without analyzing the architecture's
    members again.

    Args:
        architecture (Architecture): The architecture.
        targets (Any): The point(s) and region(s) of interest (see
            `collect_footprint_observations`).
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.
        metrics (Callable | None): A function of the observations (which may
            be empty) that gives their metrics (a mapping), if any.
        time_step (datetime.timedelta | None): The time between the samples
            of the footprints of platforms other than satellites, if any
            (see `collect_platform_observations`).
        instrument_index (int | None): The index of the observing instrument
            on each platform, or `None` for every instrument.
        elevation (float | Mapping[str, float | None] | None): The elevation
            (meters) above the WGS 84 ellipsoid of the footprints of
            platforms other than satellites, or that of each kind (see
            `collect_footprints`): required for radar stations.
        during_contact (str): Downlink of observations that end during a
            contact: `"end"` (default), `"next"`, or `"immediate"` (see
            `compute_latencies`), by timing alone.
        data_flow (DataFlow | None): A model of the flow of data, through
            which to find the downlink of each observation of satellites
            (see `compute_data_flow`), or `None` to find it by timing alone.
            Satellites that generate data and ground stations and relays that
            receive it share links through their communications.

    Returns:
        ArchitectureEvaluation: the evaluation
    """
    # pylint: disable=too-many-arguments
    _check_time_window(start, end)
    if not isinstance(architecture, Architecture):
        raise TypeError(
            f"architecture must be an Architecture, not a {type(architecture).__name__}"
        )
    analyses = _Direct(targets, time_step, instrument_index, elevation)
    return _evaluate(
        architecture, _Parts(analyses, start, end), metrics, during_contact, data_flow
    )

"""
Methods to task agile satellites: to schedule images of requested targets
(each a place, with the time from which it may be imaged, an optional time by
which it must be, and a priority) by satellites that turn (slew) their bodies
to point a narrow instrument at each target in turn, within their fields of
regard, taking time to slew between targets (see `tatc.schemas.Agility`).

The opportunities of each request are the periods in which a satellite may
point at its target (`collect_imaging_opportunities`). A schedule assigns
requests to opportunities so that each satellite has time to slew between its
images (`schedule_tasks`), and gives the images as observations (with the
columns of `collect_observations`, so that latencies and responses follow, see
`compute_latencies` and `compute_event_responses`) and the attitude of each
satellite over time (a `tatc.schemas.AttitudeSchedule`).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from typing import Any

import geopandas as gpd
import numpy as np
import numpy.typing as npt
import pandas as pd
from skyfield.api import wgs84

from ..constants import timescale
from ..schemas.agility import Agility
from ..schemas.attitude import Attitude, AttitudeSchedule, TargetAttitude
from ..schemas.instrument import Instrument
from ..schemas.space import Satellite
from ..schemas.space.base_constellation import BaseConstellation
from ..schemas.surface.point import Point
from ..utils.geometry import _hash_geometries
from ..utils.projection import compute_view_angles
from .point_sampling import collect_observations

_ENVELOPE_FIELDS = (
    "field_of_regard",
    "req_self_sunlit",
    "req_target_sunlit",
    "min_target_solar_elevation",
    "max_target_solar_elevation",
    "min_target_glint_angle",
    "max_target_glint_angle",
    "nadir_reference",
)
"""The fields of an instrument that bound the opportunities to point it at
a target: its field of regard (the cone about nadir within which the body may
point it) and the conditions of its observations."""

OPPORTUNITY_COLUMNS = (
    "request",
    "target_hash",
    "geometry",
    "satellite",
    "start",
    "end",
)
"""tuple[str, ...]: The columns of imaging opportunities (see
`collect_imaging_opportunities`)."""


def _get_members(satellites) -> list[Satellite]:
    """
    Gets the satellites of satellites and constellations (their members).
    """
    if isinstance(satellites, (Satellite, BaseConstellation)):
        satellites = [satellites]
    members = []
    for satellite in satellites:
        if isinstance(satellite, BaseConstellation):
            members.extend(satellite.generate_members())
        else:
            members.append(satellite)
    return members


def _get_window(
    requests: gpd.GeoDataFrame, start: datetime, end: datetime
) -> tuple[pd.Series, pd.Series]:
    """
    Gets the time from which each request may be imaged (its `onset`, or the
    start) and by which it must be (its `expires` time, or the end), within
    the period (timezone-aware, UTC).
    """
    first = pd.Series(pd.Timestamp(start), index=requests.index).dt.tz_convert("UTC")
    last = pd.Series(pd.Timestamp(end), index=requests.index).dt.tz_convert("UTC")
    if "onset" in requests.columns:
        onset = pd.to_datetime(requests["onset"], utc=True)
        first = onset.where(onset > first, first).fillna(first)
    if "expires" in requests.columns:
        expires = pd.to_datetime(requests["expires"], utc=True)
        last = expires.where(expires < last, last).fillna(last)
    return first, last


def _get_envelope(instrument: Instrument) -> Instrument:
    """
    Gets an instrument that observes a target whenever the body may point the
    instrument at it: within its field of regard, with its conditions.
    """
    values = {field: getattr(instrument, field) for field in _ENVELOPE_FIELDS}
    return Instrument(name=instrument.name, **values)


def collect_imaging_opportunities(
    satellites,
    requests: gpd.GeoDataFrame,
    start: datetime,
    end: datetime,
    *,
    instrument_index: int = 0,
) -> gpd.GeoDataFrame:
    """
    Collects the opportunities to image requested targets: the periods in
    which a satellite may point its instrument at a request's target, within
    the instrument's field of regard (a cone about nadir, the reach of the
    body's slews) and its conditions (such as a sunlit target), from the time
    from which the request may be imaged (its `onset`) to the time by which
    it must be (its `expires` time), within the period.

    Args:
        satellites (Satellite | BaseConstellation | list): The satellites (or
            constellations).
        requests (geopandas.GeoDataFrame): The requests: the target of each
            (a point), and optionally its `onset` and `expires` times (`NaT`
            for none; times without time zones are UTC).
        start (datetime.datetime): The start of the period.
        end (datetime.datetime): The end of the period.
        instrument_index (int): The index of the instrument of each satellite
            that images the targets.

    Returns:
        geopandas.GeoDataFrame: a row per opportunity (see
            `OPPORTUNITY_COLUMNS`), with the label of its request (the index
            of the requests), in order of start
    """
    first, last = _get_window(requests, start, end)
    hashes = pd.Series(_hash_geometries(requests.geometry), index=requests.index)
    places = requests.geometry[~hashes.duplicated().to_numpy()]
    frames = []
    for satellite in _get_members(satellites):
        instrument = satellite.instruments[instrument_index]
        envelope = satellite.model_copy(
            update={"instruments": [_get_envelope(instrument)], "attitude": None}
        )
        access = collect_observations(places, envelope, start, end)
        if len(access) == 0:
            continue
        matched = (
            pd.DataFrame({"request": requests.index, "target_hash": hashes.to_numpy()})
            .merge(
                pd.DataFrame(
                    {
                        "target_hash": access["target_hash"].to_numpy(),
                        "window_start": pd.to_datetime(access["start"], utc=True).array,
                        "window_end": pd.to_datetime(access["end"], utc=True).array,
                    }
                ),
                on="target_hash",
            )
            .assign(satellite=satellite.name)
        )
        matched["start"] = np.maximum(
            matched["window_start"].array, first.loc[matched["request"]].array
        )
        matched["end"] = np.minimum(
            matched["window_end"].array, last.loc[matched["request"]].array
        )
        frames.append(matched[matched["start"] < matched["end"]])
    if not frames:
        return gpd.GeoDataFrame(
            {
                column: (
                    pd.Series(dtype="datetime64[ns, UTC]")
                    if column in ("start", "end")
                    else pd.Series(dtype=object)
                )
                for column in OPPORTUNITY_COLUMNS
                if column != "geometry"
            },
            geometry=gpd.GeoSeries([], crs="EPSG:4326"),
        )[list(OPPORTUNITY_COLUMNS)]
    found = pd.concat(frames, ignore_index=True)
    found["geometry"] = requests.geometry.loc[found["request"]].to_numpy()
    found = found.sort_values(["start", "satellite"], kind="stable").reset_index(
        drop=True
    )
    return gpd.GeoDataFrame(
        found[list(OPPORTUNITY_COLUMNS)], geometry="geometry", crs="EPSG:4326"
    )


def _get_direction(roll: npt.ArrayLike, pitch: npt.ArrayLike) -> npt.NDArray:
    """
    Gets the unit vector of a view's center in the nadir frame (along-track,
    cross-track, nadir) of its roll and pitch angles (degrees, see
    `tatc.utils.compute_view_angles`): rolled about the along-track axis, then
    pitched about the rolled cross-track axis.
    """
    r, p = np.radians(roll), np.radians(pitch)
    return np.stack([np.sin(p), -np.sin(r) * np.cos(p), np.cos(r) * np.cos(p)], axis=-1)


NADIR = np.array([0.0, 0.0, 1.0])
"""The direction of nadir in the nadir frame."""


def _get_angle(first: npt.NDArray, second: npt.NDArray) -> float:
    """
    Gets the angle (degrees) between two unit vectors.
    """
    return float(np.degrees(np.arccos(np.clip(np.dot(first, second), -1, 1))))


class _Pointing:
    """
    The pointing of a satellite's instrument at a target over a period:
    the roll and pitch angles of the body that point its view at the target,
    sampled over the period and interpolated between.
    """

    def __init__(
        self,
        satellite: Satellite,
        instrument: Instrument,
        target,
        start: pd.Timestamp,
        end: pd.Timestamp,
        step: timedelta,
    ):
        count = max(2, int(np.ceil((end - start) / step)) + 1)
        self.times = pd.date_range(start, end, periods=count)
        self.start = start
        track = satellite.get_track(
            timescale.from_datetimes(self.times.to_pydatetime())
        )
        position = wgs84.latlon(target.y, target.x, target.z if target.has_z else 0)
        roll, pitch = compute_view_angles(
            track, position, nadir_reference=instrument.nadir_reference
        )
        self.seconds = np.asarray((self.times - start).total_seconds(), dtype=float)
        self.roll = np.atleast_1d(roll)
        self.pitch = np.atleast_1d(pitch)

    def at(self, time: pd.Timestamp) -> tuple[float, float]:
        """
        Gets the roll and pitch angles (degrees) at a time.
        """
        second = (time - self.start).total_seconds()
        return (
            float(np.interp(second, self.seconds, self.roll)),
            float(np.interp(second, self.seconds, self.pitch)),
        )

    def direction(self, time: pd.Timestamp) -> npt.NDArray:
        """
        Gets the direction of the view's center in the nadir frame at a time.
        """
        return _get_direction(*self.at(time))


_ANGLE_TOLERANCE = 0.01
"""The tolerance (degrees) of the off-nadir angle of an image beyond the
field of regard, for the interpolation of its pointing."""

_ITERATIONS = 20
"""The most iterations of the time of an image after a slew."""


@dataclass
class _Image:
    """
    An image scheduled on a satellite's timeline: its request, period, and
    the pointing of its view (at its start and end).
    """

    request: Any
    satellite: str
    opportunity: int
    start: pd.Timestamp
    end: pd.Timestamp
    first: npt.NDArray
    last: npt.NDArray


@dataclass
class TaskingPlan:
    """
    A schedule of images of requested targets by agile satellites (see
    `schedule_tasks`).
    """

    observations: gpd.GeoDataFrame
    """geopandas.GeoDataFrame: The images, as observations (with the columns
    of `collect_observations`: `target_hash`, `geometry`, `satellite`,
    `instrument`, `start`, `epoch`, and `end` of each image, and the
    elevation and azimuth of the satellite seen from its target at its
    epoch, `platform_alt` and `platform_az`), with the label of each one's
    `request`, the `roll`, `pitch`, and `off_nadir` angles (degrees) of the
    body that point its view at its target at its epoch, and the angle
    (degrees) and time of the slew to it from the previous image (or from
    nadir), `slew_angle` and `slew_time`, in order of start."""
    unscheduled: gpd.GeoDataFrame
    """geopandas.GeoDataFrame: The requests without images."""
    attitudes: dict[str, AttitudeSchedule]
    """dict[str, AttitudeSchedule]: The attitude of each satellite over time
    (by name): pointing at the target of each image during it (see
    `tatc.schemas.TargetAttitude`), and at nadir, without observing, between
    images (as while slewing)."""
    members: list[Satellite] = field(repr=False, default_factory=list)
    """list[Satellite]: The satellites (members of constellations)."""

    def satellites(self) -> list[Satellite]:
        """
        Gets the satellites with the attitudes of the plan, from which to
        analyze what they observe (such as with `collect_observations`).

        Returns:
            list[Satellite]: the satellites
        """
        return [
            satellite.model_copy(update={"attitude": self.attitudes[satellite.name]})
            for satellite in self.members
        ]


def _get_agility(agility: Agility | Mapping[str, Agility], name: str) -> Agility:
    """
    Gets the agility of a satellite (by name), of one agility or a mapping of
    satellites to agilities.
    """
    if isinstance(agility, Agility):
        return agility
    if name not in agility:
        raise ValueError(f"no agility of satellite {name!r}")
    return agility[name]


def _get_order(requests: gpd.GeoDataFrame) -> list:
    """
    Gets the labels of requests in the order in which to schedule them: by
    priority (highest first), then by expiry (earliest first, then those
    that do not expire), then by onset (earliest first).
    """
    never = pd.Timestamp.max.tz_localize("UTC")
    keys = pd.DataFrame(
        {
            "priority": (
                -pd.to_numeric(requests["priority"]).to_numpy()
                if "priority" in requests.columns
                else np.zeros(len(requests))
            ),
            "expires": (
                pd.to_datetime(requests["expires"], utc=True).fillna(never).array
                if "expires" in requests.columns
                else [never] * len(requests)
            ),
            "onset": (
                pd.to_datetime(requests["onset"], utc=True).fillna(never).array
                if "onset" in requests.columns
                else [never] * len(requests)
            ),
            "position": np.arange(len(requests)),
        },
        index=requests.index,
    )
    return list(keys.sort_values(["priority", "expires", "onset", "position"]).index)


def _fit(
    timeline: list[_Image],
    opportunity: Any,
    pointing: _Pointing,
    agility: Agility,
    dwell: timedelta,
    start: pd.Timestamp,
    reach: float,
) -> tuple[int, pd.Timestamp] | None:
    """
    Finds the earliest time at which an image of an opportunity fits a
    satellite's timeline (its images, in order): after the slew to it from
    the previous image (or from nadir, at the start of the period), within
    the opportunity and the reach of the body (the field of regard), and
    leaving time to slew to the next image. Returns the position at which to
    insert it in the timeline and its start, or `None` if it does not fit.
    """
    for position in range(len(timeline) + 1):
        previous = timeline[position - 1] if position > 0 else None
        following = timeline[position] if position < len(timeline) else None
        free = start if previous is None else previous.end
        if following is not None and following.start <= opportunity.start:
            continue
        if free >= opportunity.end:
            break
        origin = NADIR if previous is None else previous.last
        time = max(opportunity.start, free).ceil("us")
        for _ in range(_ITERATIONS):
            slew = agility.get_slew_time(_get_angle(origin, pointing.direction(time)))
            later = max(opportunity.start, free + pd.Timedelta(seconds=slew)).ceil("us")
            if later <= time:
                break
            time = later
        finish = time + pd.Timedelta(dwell)
        if finish > opportunity.end:
            continue
        if any(
            np.degrees(np.arccos(np.clip(pointing.direction(t)[2], -1, 1)))
            > reach + _ANGLE_TOLERANCE
            for t in (time, finish)
        ):
            continue
        if following is not None:
            onward = agility.get_slew_time(
                _get_angle(pointing.direction(finish), following.first)
            )
            if finish + pd.Timedelta(seconds=onward) > following.start:
                continue
        return position, time
    return None


def schedule_tasks(
    satellites,
    requests: gpd.GeoDataFrame,
    start: datetime,
    end: datetime,
    agility: Agility | Mapping[str, Agility],
    *,
    instrument_index: int = 0,
    dwell: timedelta = timedelta(0),
    pointing_step: timedelta = timedelta(seconds=1),
) -> TaskingPlan:
    """
    Schedules images of requested targets by agile satellites, each of which
    slews its body to point its instrument at one target at a time, within
    the instrument's field of regard (a cone about nadir), and images it for
    a dwell time (see `collect_imaging_opportunities`). Requests are
    scheduled in order of priority (highest first), then of expiry (earliest
    first), then of onset, each at the earliest time of its opportunities
    (on any satellite) at which its image fits the satellite's timeline: after
    the slew from the previous image (from nadir, at the start of the
    period), and leaving time to slew to the next (see
    `tatc.schemas.Agility`). A request without such a time is not scheduled.
    This greedy insertion favors priorities over the number of images, and
    does not revisit earlier choices.

    The angle of a slew is that between the directions of the view (from
    the satellite, relative to its nadir and direction of motion) at the end
    of one image and the start of the next: the turn of the body relative to
    the orbit's frame, which itself turns with the orbit. Each request is
    imaged once; the momentum, power, and data of the satellites do not
    limit their images.

    Args:
        satellites (Satellite | BaseConstellation | list): The satellites (or
            constellations), each with its instrument.
        requests (geopandas.GeoDataFrame): The requests: the target of each
            (a point), and optionally its `onset` and `expires` times (`NaT`
            for none; times without time zones are UTC) and `priority`
            (higher first; by default, 1).
        start (datetime.datetime): The start of the period.
        end (datetime.datetime): The end of the period.
        agility (Agility | Mapping[str, Agility]): The agility of the
            satellites, or of each (by name).
        instrument_index (int): The index of the instrument of each satellite
            that images the targets.
        dwell (datetime.timedelta): The time to image a target.
        pointing_step (datetime.timedelta): The step at which to sample the
            pointing at each target over its opportunities.

    Returns:
        TaskingPlan: the plan
    """
    members = _get_members(satellites)
    by_name = {satellite.name: satellite for satellite in members}
    opportunities = collect_imaging_opportunities(
        members, requests, start, end, instrument_index=instrument_index
    )
    period_start = pd.Timestamp(start).tz_convert("UTC")
    timelines: dict[str, list[_Image]] = {name: [] for name in by_name}
    pointings: dict[int, _Pointing] = {}

    def pointing_of(index: int) -> _Pointing:
        if index not in pointings:
            row = opportunities.loc[index]
            satellite = by_name[row["satellite"]]
            pointings[index] = _Pointing(
                satellite,
                satellite.instruments[instrument_index],
                row["geometry"],
                row["start"],
                row["end"],
                pointing_step,
            )
        return pointings[index]

    by_request = opportunities.groupby("request", sort=False).groups
    scheduled = []
    for request in _get_order(requests):
        for index in by_request.get(request, []):
            opportunity = opportunities.loc[index]
            name = opportunity["satellite"]
            satellite = by_name[name]
            pointing = pointing_of(index)
            fitted = _fit(
                timelines[name],
                opportunity,
                pointing,
                _get_agility(agility, name),
                dwell,
                period_start,
                satellite.instruments[instrument_index].field_of_regard / 2,
            )
            if fitted is None:
                continue
            position, time = fitted
            finish = time + pd.Timedelta(dwell)
            timelines[name].insert(
                position,
                _Image(
                    request,
                    name,
                    index,
                    time,
                    finish,
                    pointing.direction(time),
                    pointing.direction(finish),
                ),
            )
            scheduled.append(request)
            break
    observations = _get_observations(
        timelines, opportunities, pointings, by_name, agility, instrument_index
    )
    attitudes = {
        name: _get_attitude_schedule(
            timelines[name],
            opportunities,
            period_start,
            by_name[name].instruments[instrument_index].field_of_regard / 2,
        )
        for name in by_name
    }
    unscheduled = requests.loc[~requests.index.isin(scheduled)]
    return TaskingPlan(observations, unscheduled, attitudes, members)


def _get_observations(
    timelines: dict[str, list[_Image]],
    opportunities: gpd.GeoDataFrame,
    pointings: dict[int, _Pointing],
    by_name: dict[str, Satellite],
    agility: Agility | Mapping[str, Agility],
    instrument_index: int,
) -> gpd.GeoDataFrame:
    """
    Gets the images of the timelines of satellites as observations (see
    `TaskingPlan.observations`).
    """
    rows = []
    for name, timeline in timelines.items():
        satellite = by_name[name]
        origin = NADIR
        for image in timeline:
            opportunity = opportunities.loc[image.opportunity]
            epoch = (image.start + (image.end - image.start) / 2).floor("us")
            roll, pitch = pointings[image.opportunity].at(epoch)
            angle = _get_angle(origin, image.first)
            target = opportunity["geometry"]
            track = satellite.get_track(timescale.from_datetime(epoch.to_pydatetime()))
            topos = wgs84.latlon(target.y, target.x, target.z if target.has_z else 0)
            altitude, azimuth, _ = (track - topos.at(track.t)).altaz()
            rows.append(
                {
                    "target_hash": opportunity["target_hash"],
                    "geometry": target,
                    "satellite": name,
                    "instrument": satellite.instruments[instrument_index].name,
                    "start": image.start,
                    "epoch": epoch,
                    "end": image.end,
                    "platform_alt": float(altitude.degrees),
                    "platform_az": float(azimuth.degrees),
                    "request": image.request,
                    "roll": roll,
                    "pitch": pitch,
                    "off_nadir": _get_angle(NADIR, _get_direction(roll, pitch)),
                    "slew_angle": angle,
                    "slew_time": pd.Timedelta(
                        seconds=_get_agility(agility, name).get_slew_time(angle)
                    ),
                }
            )
            origin = image.last
    columns = [
        "target_hash",
        "geometry",
        "satellite",
        "instrument",
        "start",
        "epoch",
        "end",
        "platform_alt",
        "platform_az",
        "request",
        "roll",
        "pitch",
        "off_nadir",
        "slew_angle",
        "slew_time",
    ]
    if not rows:
        frame = pd.DataFrame({column: [] for column in columns})
        for column in ("start", "epoch", "end"):
            frame[column] = pd.Series(dtype="datetime64[ns, UTC]")
        frame["slew_time"] = pd.Series(dtype="timedelta64[ns]")
        frame["geometry"] = gpd.GeoSeries([])
    else:
        frame = pd.DataFrame(rows, columns=columns)
        frame = frame.sort_values(["start", "satellite"], kind="stable")
    return gpd.GeoDataFrame(
        frame.reset_index(drop=True), geometry="geometry", crs="EPSG:4326"
    )


def _get_attitude_schedule(
    timeline: list[_Image],
    opportunities: gpd.GeoDataFrame,
    start: pd.Timestamp,
    reach: float,
) -> AttitudeSchedule:
    """
    Gets the attitude of a satellite over time from its images: pointing at
    the target of each image during it, and at nadir, without observing,
    between them.
    """
    starts, attitudes, observing = [start.to_pydatetime()], [Attitude()], [False]
    for image in timeline:
        target = opportunities.loc[image.opportunity, "geometry"]
        point = Point(
            latitude=target.y,
            longitude=target.x,
            elevation=target.z if target.has_z else 0,
        )
        starts.append(image.start.floor("us").to_pydatetime())
        attitudes.append(
            TargetAttitude(
                targets=[point],
                max_off_nadir_angle=min(reach + _ANGLE_TOLERANCE, 89.99),
            )
        )
        observing.append(True)
        starts.append(image.end.floor("us").to_pydatetime())
        attitudes.append(Attitude())
        observing.append(False)
    return AttitudeSchedule.from_attitudes(starts, attitudes, observing)

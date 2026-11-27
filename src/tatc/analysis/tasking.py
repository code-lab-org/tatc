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

from datetime import datetime, timedelta

import geopandas as gpd
import numpy as np
import numpy.typing as npt
import pandas as pd
from skyfield.api import wgs84

from ..constants import timescale
from ..schemas.instrument import Instrument
from ..schemas.space import Satellite
from ..schemas.space.base_constellation import BaseConstellation
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
    values = {
        field: getattr(instrument, field)
        for field in _ENVELOPE_FIELDS
        if field in Instrument.model_fields and hasattr(instrument, field)
    }
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

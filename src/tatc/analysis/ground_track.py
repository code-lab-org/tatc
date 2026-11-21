"""
Methods to collect the ground tracks of instruments on satellites.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timedelta

import geopandas as gpd
import numpy as np
import numpy.typing as npt
import pandas as pd
from shapely.geometry import MultiPolygon, Polygon
from skyfield.positionlib import Geocentric

from ..schemas import AllInstruments, Satellite
from ..utils.time import _index_orbit_track
from ._tracks import (
    Mask,
    _assemble,
    _check_pixel_instrument,
    _collect_tracks,
    _compute_step_footprints,
    _FootprintRequest,
    _get_platform_altaz,
)
from .check import (
    SATELLITE,
    _check_distinct_names,
    _check_satellites,
    _dissolve_footprints,
    _get_instrument_indices,
)
from .region_sampling import compute_region_access_periods


def _cull_orbit_track(
    satellite: Satellite,
    instrument: AllInstruments,
    times: list[datetime],
    half_spans: npt.NDArray,
    mask: Polygon | MultiPolygon,
    elevation: float,
) -> tuple[Geocentric, list[Polygon | MultiPolygon], list[datetime | None]] | None:
    """
    Propagates the orbit track only at the times when the instrument's
    footprint intersects a mask. The times are first culled to those whose
    step (the time, or its span for a swept instrument; see
    `_get_half_spans`) overlaps a period when the instrument's field of regard
    may observe any part of the mask (see
    `tatc.analysis.region_sampling.compute_region_access_periods`), so that
    the orbit is not propagated far from the mask, and then to those whose
    footprint intersects the mask.

    Args:
        satellite (Satellite): The observing satellite.
        instrument (AllInstruments): The observing instrument.
        times (list[datetime.datetime]): The times to sample.
        half_spans (numpy.typing.NDArray): The half span (seconds) of each
            time's step.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon):
                The mask, in WGS84 (lon/lat) coordinates (see `_get_mask_geometry`).
        elevation (float): The elevation (meters) above the WGS 84 ellipsoid
                at which to project footprints.

    Returns:
        tuple[skyfield.positionlib.Geocentric, list[shapely.geometry.Polygon | shapely.geometry.MultiPolygon], list[datetime.datetime | None]] | None:
            the culled orbit track, its footprints, and the times at which to
            evaluate their validity if not their own (see
            `_compute_step_footprints`), or None if no footprint intersects
            the mask
    """
    periods = compute_region_access_periods(
        mask,
        satellite,
        min(times),
        max(times),
        instrument.field_of_regard,
        elevation,
    )
    if periods.empty:
        return None
    # select the times whose steps overlap a period (the periods are sorted
    # and disjoint)
    sample = pd.DatetimeIndex(pd.to_datetime(times, utc=True))
    span = pd.to_timedelta(np.asarray(half_spans, dtype=float), unit="s")
    left = pd.DatetimeIndex([period.left for period in periods])
    right = pd.DatetimeIndex([period.right for period in periods])
    index = np.searchsorted(left, sample + span, side="right") - 1
    selected = (index >= 0) & (sample - span <= right[np.maximum(index, 0)])
    if not any(selected):
        return None
    chosen = np.flatnonzero(selected)
    orbit_track = satellite.get_track([times[i] for i in chosen])
    # cull orbit track to observable footprint
    footprint, represent = _compute_step_footprints(
        satellite,
        instrument,
        [times[i] for i in chosen],
        np.asarray(half_spans)[chosen],
        elevation,
        orbit_track,
    )
    keep = np.flatnonzero([mask.intersects(f) for f in footprint])
    if len(keep) == 0:
        return None
    # keep the Earth orientation quantities cached on the times
    return (
        _index_orbit_track(orbit_track, keep),
        [footprint[i] for i in keep],
        [represent[i] for i in keep],
    )


def collect_ground_track(
    satellites: Satellite | list[Satellite],
    times: list[datetime],
    instrument_index: int | None = 0,
    elevation: float | None = None,
    mask: Mask | None = None,
    platform_altaz: bool = False,
    solar_altaz: bool = False,
    time_step: timedelta | None = None,
    *,
    sat_altaz: bool | None = None,
) -> gpd.GeoDataFrame:
    """
    Collect the instrument's viewable ground footprint at each requested
    time, projected to a specified elevation by computing the exact
    ray/WGS 84 ellipsoid intersection (see
    `tatc.utils.projection.compute_footprint`). The footprint of a swept
    instrument (see `PointedInstrument.is_swept`) is the strip its
    cross-track line sweeps over the time step centered on each time, trimmed
    to the times at which it observes (see `tatc.schemas.Attitude.blind_intervals`).

    Args:
        satellites (Satellite | list[Satellite]): The observing satellite(s).
        times (list[datetime.datetime]): The times to sample.
        instrument_index (int | None): The index of the observing instrument
                on each satellite, or `None` for every instrument of each
                satellite.
        elevation (float | None): The elevation (meters) above the WGS 84
                ellipsoid at which to compute footprints, or `None` for the
                ellipsoid (zero).
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame | geopandas.GeoSeries | None):
                An optional mask, always interpreted in WGS84 (lon/lat)
                coordinates, to constrain results. Providing a mask limits
                the propagated orbit to the periods when the instrument may
                observe the mask (see `compute_region_access_periods`),
                improving performance for small areas.
        platform_altaz (bool): `True` to include the altitude/azimuth angles
                of the satellite (`platform_alt` and `platform_az`) seen
                from the footprint center.
        solar_altaz (bool): `True` to include the solar altitude/azimuth
                angles (`solar_alt` and `solar_az`) and the glint angle
                (`glint_angle`, see `tatc.utils.compute_glint_angle`) at the
                footprint center. The angles are `NaN` for an empty
                footprint.
        time_step (datetime.timedelta | None): The time step that each time
                represents, over which a swept instrument's footprint sweeps,
                or `None` to take half the smaller gap to the neighboring
                times on either side of each time (so that a gap between
                times does not stretch the footprints beside it, and a single
                time is a view at an instant).
        sat_altaz (bool | None): Deprecated: use `platform_altaz` (whose
                angles were named `sat_alt` and `sat_az`).

    Returns:
        geopandas.GeoDataFrame: The data frame of collected ground track results:
            for a single satellite and instrument (an integer
            `instrument_index`), sorted by time; otherwise, those of each
            satellite and instrument, concatenated and sorted by time.
    """
    return _collect_tracks(
        satellites,
        SATELLITE,
        _FootprintRequest.create(
            times,
            instrument_index,
            elevation,
            mask,
            _get_platform_altaz(platform_altaz, sat_altaz),
            solar_altaz,
            time_step,
        ),
        _collect_satellite,
    )


def _collect_satellite(
    satellite: Satellite, instrument_index: int, request: _FootprintRequest
) -> gpd.GeoDataFrame | None:
    """
    Collect a satellite instrument's footprints (see `collect_ground_track`)
    or, if the request is of pixels, its pixels (see
    `collect_ground_pixels`), or `None` if there are none within the mask.

    Args:
        satellite (Satellite): The satellite.
        instrument_index (int): The index of the instrument.
        request (_FootprintRequest): The request.

    Returns:
        geopandas.GeoDataFrame | None: the footprints or pixels
    """
    instrument = satellite.get_instrument(instrument_index)
    times: list[datetime] = request.times  # type: ignore
    elevation = request.get_elevation(satellite)
    if request.pixels:
        # before propagating the orbit
        _check_pixel_instrument(instrument)
    # pixel arrays are views at an instant (steps of no span)
    half_spans = np.zeros(len(times)) if request.pixels else request.half_spans
    if request.mask is not None and len(times) > 1:
        # propagate orbit only where the footprint intersects the mask,
        # reusing the footprints computed to cull it
        culled = _cull_orbit_track(
            satellite, instrument, times, half_spans, request.mask, elevation
        )
        if culled is None:
            return None
        orbit_track, geometries, represent = culled
    else:
        orbit_track = satellite.get_track(request.skyfield_times)
        geometries, represent = (
            ([], [])
            if request.pixels
            # compute footprints (exact ray/WGS 84 ellipsoid intersection)
            else _compute_step_footprints(
                satellite, instrument, times, half_spans, elevation, orbit_track
            )
        )
    return _assemble(
        satellite,
        instrument,
        orbit_track.t.utc_datetime(),  # type: ignore
        orbit_track,
        elevation,
        request,
        geometries,
        represent,
    )


def _get_empty_orbit_coverage() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for the coverage of each orbit (see
    `compute_ground_track`).

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "satellite": pd.Series([], dtype="str"),
        "instrument": pd.Series([], dtype="str"),
        "orbit_id": pd.Series([], dtype="int"),
        "start": pd.Series([], dtype="datetime64[ns, utc]"),
        "end": pd.Series([], dtype="datetime64[ns, utc]"),
        "geometry": pd.Series([], dtype="object"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def _compute_orbit_coverage(
    satellite: Satellite, track: gpd.GeoDataFrame, first: datetime
) -> gpd.GeoDataFrame:
    """
    Dissolves the valid footprints of a satellite's instrument within each
    orbit (see `compute_ground_track`): the number of periods, of the element
    that propagates the orbit at each time (see `get_element_index`), since
    the first time.

    Args:
        satellite (Satellite): The satellite.
        track (geopandas.GeoDataFrame): The valid footprints of an instrument.
        first (datetime.datetime): The first time.

    Returns:
        geopandas.GeoDataFrame: the coverage of each orbit
    """
    orbit = satellite.orbit.to_gp_orbit()
    orbit_id = [
        int(
            (time - first)
            // orbit.elements[orbit.get_element_index(time)].get_orbit_period()
        )
        for time in track["time"]
    ]
    track = track.assign(orbit_id=orbit_id)
    times = track.groupby("orbit_id")["time"]
    dissolved = track[["orbit_id", "geometry"]].dissolve(by="orbit_id")
    return gpd.GeoDataFrame(
        {
            "satellite": satellite.name,
            "instrument": track["instrument"].iloc[0],
            "orbit_id": dissolved.index.to_numpy(),
            "start": times.min().to_numpy(),
            "end": times.max().to_numpy(),
        },
        geometry=dissolved.geometry.to_numpy(),
        crs=track.crs,
    )


def compute_ground_track(
    satellites: Satellite | list[Satellite],
    times: list[datetime],
    instrument_index: int | None = 0,
    elevation: float | None = None,
    mask: Mask | None = None,
    dissolve_orbits: bool = True,
    time_step: timedelta | None = None,
) -> gpd.GeoDataFrame:
    """
    Compute the aggregated (dissolved) coverage of instruments on satellites:
    unlike `collect_ground_track` (one footprint per time), this dissolves
    the footprints of valid observations into one geometry (for example,
    the area observed over the times) or, if not `dissolve_orbits`, into one
    geometry for each orbit of each satellite's instrument.

    Args:
        satellites (Satellite | list[Satellite]): The observing satellite(s).
        times (list[datetime.datetime]): The times to sample.
        instrument_index (int | None): The index of the observing instrument
                on each satellite, or `None` for every instrument of each
                satellite.
        elevation (float | None): The elevation (meters) above the WGS 84
                ellipsoid at which to compute footprints, or `None` for the
                ellipsoid (zero).
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame | geopandas.GeoSeries | None):
                An optional mask, always interpreted in WGS84 (lon/lat)
                coordinates, to constrain results.
        dissolve_orbits (bool): `True` to dissolve the coverage of every orbit,
                satellite, and instrument into one geometry; `False` for
                one geometry for each orbit of each satellite's instrument.
        time_step (datetime.timedelta | None): The time step that each time
                represents, over which a swept instrument's footprint sweeps
                (see `collect_ground_track`), or `None` to infer it from the
                neighboring times.

    Returns:
        geopandas.GeoDataFrame: The data frame of the aggregated ground
            track: if `dissolve_orbits`, one row with the observed area (or
            none if no observation is valid); otherwise, one row for each
            orbit of each satellite's instrument with a valid observation,
            with columns `satellite`, `instrument`, `orbit_id` (the number of
            orbital periods since the first time), `start` and `end` (the
            first and last times of its valid footprints), and `geometry`,
            sorted by `start`.
    """
    request = _FootprintRequest.create(
        times, instrument_index, elevation, mask, time_step=time_step
    )
    if dissolve_orbits:
        return _dissolve_footprints(
            _collect_tracks(satellites, SATELLITE, request, _collect_satellite)
        )
    checked = _check_satellites(satellites)
    _check_distinct_names(checked, SATELLITE, SATELLITE.argument)
    orbits = []
    for satellite in checked if len(times) > 0 else []:
        for index in _get_instrument_indices(satellite, instrument_index):
            track = _collect_satellite(satellite, index, request)
            if track is not None:
                track = track[track["valid_obs"]]
                if not track.empty:
                    orbits.append(_compute_orbit_coverage(satellite, track, min(times)))
    if len(orbits) == 0:
        return _get_empty_orbit_coverage()
    return (
        pd.concat(orbits, ignore_index=True)
        .sort_values("start", kind="stable")
        .reset_index(drop=True)
    )


def collect_ground_pixels(
    satellites: Satellite | list[Satellite],
    times: list[datetime],
    instrument_index: int | None = 0,
    elevation: float | None = None,
    mask: Mask | None = None,
    platform_altaz: bool = False,
    solar_altaz: bool = False,
    *,
    sat_altaz: bool | None = None,
) -> gpd.GeoDataFrame:
    """
    Collect the instrument's individual ground sample point (pixel)
    locations at each requested time, rather than an aggregate observable
    geometry (see `collect_ground_track`). Only supported for a rectangular
    `PointedInstrument`, whose `cross_track_pixels` x `along_track_pixels`
    grid defines the pixel array. For the pixels of aircraft and fixed
    surface platforms, see `tatc.analysis.footprints.collect_footprints`.

    Args:
        satellites (Satellite | list[Satellite]): The observing satellite(s).
        times (list[datetime.datetime]): The times to sample.
        instrument_index (int | None): The index of the observing instrument
                on each satellite, or `None` for every instrument of each
                satellite.
        elevation (float | None): The elevation (meters) above the WGS 84
                ellipsoid at which to compute pixels, or `None` for the
                ellipsoid (zero).
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame | geopandas.GeoSeries | None):
                An optional mask, always interpreted in WGS84 (lon/lat)
                coordinates, to constrain results. Providing a mask limits
                the propagated orbit to the periods when the instrument may
                observe the mask (see `compute_region_access_periods`),
                improving performance for small areas.
        platform_altaz (bool): `True` to include the altitude/azimuth angles
                of the satellite (`platform_alt` and `platform_az`) seen
                from each pixel.
        solar_altaz (bool): `True` to include the solar altitude/azimuth
                angles (`solar_alt` and `solar_az`) and the glint angle
                (`glint_angle`, see `tatc.utils.compute_glint_angle`) at
                each pixel.
        sat_altaz (bool | None): Deprecated: use `platform_altaz` (whose
                angles were named `sat_alt` and `sat_az`).

    Returns:
        geopandas.GeoDataFrame: The data frame of collected ground pixels results:
            for a single satellite and instrument (an integer
            `instrument_index`), sorted by time; otherwise, those of each
            satellite and instrument, concatenated and sorted by time.
    """
    return _collect_tracks(
        satellites,
        SATELLITE,
        _FootprintRequest.create(
            times,
            instrument_index,
            elevation,
            mask,
            _get_platform_altaz(platform_altaz, sat_altaz),
            solar_altaz,
            pixels=True,
        ),
        _collect_satellite,
    )

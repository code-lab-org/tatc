"""
Methods to compute the sizes of the samples of observations: the extent on
the ground of an instrument's pixel at each observed target, and the
spacing between pixels, which bound the detail an observation can resolve.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import shapely
from skyfield.api import wgs84

from ..schemas import ConicalInstrument, PointedInstrument, Satellite
from ..utils.projection import VelocityFrame, compute_view_angles
from ..utils.time import _to_time
from .check import _check_satellites

_COLUMNS = [
    "sat_range",
    "scan_angle",
    "sample_along",
    "sample_cross",
    "spacing_along",
    "spacing_cross",
]


def compute_sample_sizes(
    observations: gpd.GeoDataFrame,
    satellites: Satellite | list[Satellite],
) -> gpd.GeoDataFrame:
    """
    Compute the sizes of the samples of observations: for each observation,
    the extent on the ground of the observing instrument's pixel at the
    target, along track and across track (`sample_along` and
    `sample_cross`, meters), and the spacing between the centers of
    neighboring pixels (`spacing_along` and `spacing_cross`, meters), which
    is larger than the pixel if there are gaps between pixels (a negative
    oversampling). An observation cannot resolve detail finer than the
    larger of the two. The sizes are those of a `PointedInstrument`'s
    pixels (see `PointedInstrument.compute_sample_size` and
    `PointedInstrument.compute_sample_spacing`) or of a
    `ConicalInstrument`'s beam and samples along its scanned arc (see
    `ConicalInstrument.compute_sample_size` and
    `ConicalInstrument.compute_sample_spacing`, `NaN` unless its
    `beam_width` and `scan_sample_spacing` are given), and `NaN` for other
    instruments, whose samples are not defined.

    Also computes the satellite's slant range from the target (`sat_range`,
    meters) and the target's scan angle (`scan_angle`, degrees): its angle
    across track from nadir, positive to the left of the direction of motion
    (see `tatc.utils.compute_view_angles`), measured from the body's axes
    if the satellite has an attitude, which gives its position in the swath.

    The observations are those of points (see `collect_observations`), at
    their `epoch`, or of features (see `collect_feature_observations`), at
    their `time`, whose target is the centroid of the viewed part.

    Args:
        observations (geopandas.GeoDataFrame): The observations, with the
                `satellite` and `instrument` that observe each.
        satellites (Satellite | list[Satellite]): The observing satellites,
                including those named by the observations.

    Returns:
        geopandas.GeoDataFrame: the observations, with the sample sizes,
            spacings, slant ranges, and scan angles.
    """
    satellites = _check_satellites(satellites)
    by_name = {satellite.name: satellite for satellite in satellites}
    result = observations.copy()
    for column in _COLUMNS:
        result[column] = np.nan
    if len(observations) == 0:
        return result
    time_column = "epoch" if "epoch" in observations.columns else "time"
    # the targets: observed points, or the centroids of observed regions
    geometry = observations.geometry.to_numpy()
    is_point = shapely.get_type_id(geometry) == shapely.GeometryType.POINT
    centroid = shapely.centroid(geometry)
    longitude, latitude = shapely.get_x(centroid), shapely.get_y(centroid)
    elevation = np.where(is_point, np.nan_to_num(shapely.get_z(geometry)), 0)
    position = np.arange(len(observations))
    for (satellite_name, instrument_name), group in observations.groupby(
        ["satellite", "instrument"], sort=False
    ):
        satellite = by_name.get(satellite_name)
        if satellite is None:
            raise ValueError(f"no satellite named {satellite_name!r}")
        instrument_index = next(
            (
                i
                for i, instrument in enumerate(satellite.instruments)
                if instrument.name == instrument_name
            ),
            None,
        )
        if instrument_index is None:
            raise ValueError(
                f"satellite {satellite_name!r} has no instrument named {instrument_name!r}"
            )
        # the instrument, with its views rotated from the body's axes
        instrument = satellite.get_instrument(instrument_index)
        index = position[observations.index.get_indexer(group.index)]
        orbit_track = satellite.orbit.to_gp_orbit().get_orbit_track_at_time(
            _to_time(list(group[time_column]))
        )
        target = wgs84.latlon(latitude[index], longitude[index], elevation[index])
        distance = (orbit_track - target.at(orbit_track.t)).altaz()[2]
        result.iloc[index, result.columns.get_loc("sat_range")] = np.atleast_1d(
            distance.m
        )
        velocity_frame = (
            instrument.velocity_frame
            if isinstance(instrument, (PointedInstrument, ConicalInstrument))
            else VelocityFrame.EARTH_FIXED
        )
        roll, _ = compute_view_angles(
            orbit_track,
            target,
            velocity_frame,
            instrument.nadir_reference,
            attitude=instrument.platform_attitude,
        )
        result.iloc[index, result.columns.get_loc("scan_angle")] = np.atleast_1d(roll)
        if isinstance(instrument, (PointedInstrument, ConicalInstrument)):
            for prefix, method in (
                ("sample", instrument.compute_sample_size),
                ("spacing", instrument.compute_sample_spacing),
            ):
                along, cross = method(orbit_track, target)
                result.iloc[index, result.columns.get_loc(f"{prefix}_along")] = (
                    np.atleast_1d(along)
                )
                result.iloc[index, result.columns.get_loc(f"{prefix}_cross")] = (
                    np.atleast_1d(cross)
                )
    return result

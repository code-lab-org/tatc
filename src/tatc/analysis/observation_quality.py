"""
Methods to assess the quality of observations: the resolution of their
samples, and the noise of their measurements relative to that of an
instrument's nominal performance (see `tatc.schemas.InstrumentPerformance`),
from which to simulate observation errors in an observing system simulation
experiment (OSSE).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Mapping

import geopandas as gpd
import numpy as np
import pandas as pd

from ..schemas import AllInstruments, AllPlatforms
from .check import _get_instruments, _get_platform_lookup, _get_platforms
from .sample_sizes import _get_sized_instrument, compute_sample_sizes

_QUALITY_COLUMNS = [
    "resolution_along",
    "resolution_cross",
    "ground_sample_distance",
    "noise_scale",
    "signal_to_noise_ratio",
]
"""The columns of the quality of observations (besides their sample sizes)."""


def compute_observation_quality(
    observations: gpd.GeoDataFrame,
    platforms: AllPlatforms | list[AllPlatforms],
    instrument: AllInstruments | Mapping[str, AllInstruments] | None = None,
) -> gpd.GeoDataFrame:
    """
    Computes the quality of observations: the sizes of their samples (see
    `compute_sample_sizes`: the extent and spacing of each pixel or beam
    along and across track, its azimuth, the slant range, and the scan
    angle), the resolution they bound, and the scaling of the observing
    instrument's nominal performance (see `tatc.schemas.InstrumentPerformance`)
    to the slant range of each observation.

    The resolution along and across track (`resolution_along` and
    `resolution_cross`, meters) is the larger of the size of a sample and the
    spacing between samples, the finest detail an observation can resolve;
    where the samples of an instrument are not modeled (`NaN` sizes, as of an
    instrument without pixels), it is the nominal ground sample distance.

    The performance of the instrument (if given) is scaled by its power laws
    to the slant range (`platform_range`) of each observation, which equals
    the altitude at nadir and grows off nadir (as the range of a lidar's or
    radar's echo does):

    - `ground_sample_distance` (meters): the nominal ground sample distance at
      the slant range (as at nadir at that altitude, without the lengthening
      of an oblique view, which the sample sizes include);
    - `noise_scale`: the factor by which noise relative to the signal grows
      from the reference altitude (`(platform_range / reference_altitude) **
      signal_to_noise_exponent`, 1 for a passive instrument), with which to
      scale an observation error stated at the reference altitude (see
      `simulate_observation_errors`);
    - `signal_to_noise_ratio`: the nominal signal-to-noise ratio divided by
      the noise scale.

    Each is `NaN` for observations without the performance it needs: of
    instruments without performance (or without a ground sample distance or
    signal-to-noise ratio), of radar stations (which have no instruments),
    and of observations without an `instrument` column (such as limb
    soundings), unless `instrument` gives one.

    Args:
        observations (geopandas.GeoDataFrame): The observations, with the
                platform and `instrument` that observe each (see
                `compute_sample_sizes`).
        platforms (AllPlatforms | list[AllPlatforms]): The observing
                platforms, including those named by the observations.
        instrument (AllInstruments | Mapping[str, AllInstruments] | None):
                The instrument whose samples and performance to assess in
                place of the observing instrument, or one for each observing
                instrument's name, or `None` for the observing instruments.

    Returns:
        geopandas.GeoDataFrame: the observations, with the sample sizes (see
            `compute_sample_sizes`), resolutions, and scaled performance.
    """
    result = compute_sample_sizes(observations, platforms, instrument)
    for column in _QUALITY_COLUMNS:
        result[column] = np.nan
    if len(result) == 0:
        return result
    lookup = _get_platform_lookup(platforms)
    distance = result["platform_range"].to_numpy(dtype=float)
    position = np.arange(len(result))
    columns = {column: result.columns.get_loc(column) for column in _QUALITY_COLUMNS}
    for (key, name), index in (
        pd.DataFrame(
            {
                "key": _get_platforms(observations),
                "name": _get_instruments(observations),
            }
        )
        .groupby(["key", "name"], sort=False, dropna=False)
        .indices.items()
    ):
        index = position[index]
        sized = _get_sized_instrument(lookup[key], name, instrument)
        performance = getattr(sized, "performance", None)
        if performance is None:
            continue
        reach = distance[index]
        result.iloc[index, columns["noise_scale"]] = np.atleast_1d(
            performance.get_noise_scale(reach)
        )
        if performance.ground_sample_distance is not None:
            result.iloc[index, columns["ground_sample_distance"]] = np.atleast_1d(
                performance.get_ground_sample_distance(reach)
            )
        if performance.signal_to_noise_ratio is not None:
            result.iloc[index, columns["signal_to_noise_ratio"]] = np.atleast_1d(
                performance.get_signal_to_noise_ratio(reach)
            )
    nominal = result["ground_sample_distance"].to_numpy(dtype=float)
    for axis in ("along", "cross"):
        resolution = np.fmax(
            result[f"sample_{axis}"].to_numpy(dtype=float),
            result[f"spacing_{axis}"].to_numpy(dtype=float),
        )
        result[f"resolution_{axis}"] = np.where(
            np.isnan(resolution), nominal, resolution
        )
    return result

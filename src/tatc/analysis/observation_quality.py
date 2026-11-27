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
    and of observations without an `instrument` column (such as those of
    `collect_limb_observations`), unless `instrument` gives one. The
    observations of a limb sounder (see `collect_limb_soundings`) have its
    performance at the slant range to the representative tangent point.

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


def _get_values(
    frame: pd.DataFrame, value: float | str | None, default: float
) -> np.ndarray:
    """
    Gets a value for each row of a frame: a number for every row, the values
    of a column (by name), or a default (for `None`).
    """
    if value is None:
        return np.full(len(frame), default)
    if isinstance(value, str):
        if value not in frame.columns:
            raise ValueError(f"the observations have no column {value!r}")
        return frame[value].to_numpy(dtype=float)
    return np.full(len(frame), float(value))


def simulate_observation_errors(
    sampled: gpd.GeoDataFrame,
    name: str,
    standard_deviation: float | str,
    noise_scale: bool | str = True,
    bias: float | str | None = None,
    seed: int | np.random.Generator | None = None,
) -> gpd.GeoDataFrame:
    """
    Simulates the errors of observations of sampled values (such as a field
    of a nature run sampled at observations, see `sample_field`): random
    errors drawn from a normal distribution with a standard deviation (an
    instrument's error at its reference altitude, as a noise equivalent
    temperature difference, which the user gives), scaled by each
    observation's noise scale (see `compute_observation_quality`), plus an
    optional bias. Gives the standard deviation of each observation's error
    (`{name}_error_std`) and its observed value (`{name}_observed`: the
    sampled value plus the bias and the error), as the observations of an
    observing system simulation experiment (OSSE) and their error, which a
    data assimilation system takes. The errors are independent; errors
    correlated between observations (such as those of a scan) are not
    simulated.

    Args:
        sampled (geopandas.GeoDataFrame): The observations, with the sampled
            values (`name`).
        name (str): The column of the sampled values.
        standard_deviation (float | str): The standard deviation of the
            errors at the reference altitude of the instrument's performance,
            for every observation or by a column (as of errors that depend on
            the scene).
        noise_scale (bool | str): `True` (default) to scale the standard
            deviation by each observation's `noise_scale` (see
            `compute_observation_quality`; observations without one, `NaN`,
            are not scaled), `False` not to, or the name of a column of the
            scale.
        bias (float | str | None): The bias of the observations, for every
            observation or by a column, or `None` for none.
        seed (int | numpy.random.Generator | None): The seed of the random
            errors (or a generator of them), for reproducible errors.

    Returns:
        geopandas.GeoDataFrame: the observations, with the standard deviation
            of their errors and their observed values
    """
    if name not in sampled.columns:
        raise ValueError(f"the observations have no column {name!r}")
    if noise_scale is True:
        if "noise_scale" not in sampled.columns:
            raise ValueError(
                "the observations have no noise_scale: see "
                + "tatc.analysis.compute_observation_quality (or noise_scale=False)"
            )
        scale = np.nan_to_num(sampled["noise_scale"].to_numpy(dtype=float), nan=1.0)
    elif noise_scale is False:
        scale = np.ones(len(sampled))
    else:
        scale = _get_values(sampled, noise_scale, 1.0)
    error_std = _get_values(sampled, standard_deviation, 0.0) * scale
    if np.any(error_std < 0):
        raise ValueError("the standard deviation of errors must not be negative")
    generator = (
        seed if isinstance(seed, np.random.Generator) else np.random.default_rng(seed)
    )
    error = generator.standard_normal(len(sampled)) * error_std
    result = sampled.copy()
    result[f"{name}_error_std"] = error_std
    result[f"{name}_observed"] = (
        sampled[name].to_numpy(dtype=float) + _get_values(sampled, bias, 0.0) + error
    )
    return result

"""
Methods to assess the statistical sufficiency of observations: how well the
samples of a target (such as a point) estimate a statistic of a varying
field, and how their times are distributed.

Observations sample a field (such as precipitation or temperature) at
discrete times, from which a statistic (such as the mean over a month) is
estimated. Its error depends on how many samples there are, and on how they
are spread in time relative to how quickly the field decorrelates: samples
closer together than the decorrelation time repeat each other. The error
can be estimated for a model of the field (`estimate_sampling_error`) or
measured against a simulation of it (a nature run, see `sample_field` and
`compare_sampled_means`). The distribution of the samples in time
(`compute_revisit_statistics`) and in local time of day
(`compute_local_time_coverage`, for fields with daily cycles) explains it.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from datetime import datetime, timedelta
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd

from ..utils.geometry import _hash_geometries


def _get_points(geometries: gpd.GeoSeries) -> tuple[np.ndarray, np.ndarray]:
    """
    Gets the longitudes and latitudes (degrees) of geometries: points, or the
    representative points of other geometries.
    """
    points = geometries.representative_point()
    return points.x.to_numpy(), points.y.to_numpy()


def _get_seconds(times: pd.Series | pd.Index, origin: datetime) -> np.ndarray:
    """
    Gets times as seconds after an origin.
    """
    times = pd.DatetimeIndex(pd.to_datetime(times, utc=True))
    return np.asarray((times - pd.Timestamp(origin)).total_seconds(), dtype=float)


def _get_targets(
    observations: pd.DataFrame, targets: gpd.GeoDataFrame | gpd.GeoSeries | None
) -> gpd.GeoDataFrame:
    """
    Gets the targets: those given (identified by the hash of their geometry),
    or those observed.
    """
    if targets is None:
        first = ~observations["target_hash"].duplicated().to_numpy()
        result = observations.loc[first, ["target_hash", "geometry"]]
        result = result.reset_index(drop=True)
        return result if result.crs is not None else result.set_crs("EPSG:4326")
    geometry = targets if isinstance(targets, gpd.GeoSeries) else targets.geometry
    return gpd.GeoDataFrame(
        {"target_hash": _hash_geometries(geometry)},
        geometry=geometry.to_numpy(),
        crs="EPSG:4326",
    )


def _exponential_integral(t: np.ndarray, length: float, scale: float) -> np.ndarray:
    """
    Integrates the exponential correlation exp(-|t - s| / scale) over s in
    [0, length], for each time t.
    """
    inside = (t >= 0) & (t <= length)
    before = np.exp(-np.abs(t) / scale) - np.exp(-np.abs(length - t) / scale)
    within = 2 - np.exp(-t / scale) - np.exp(-(length - t) / scale)
    return scale * np.where(inside, within, np.abs(before))


def sampling_error_variance(
    times: Sequence[float] | np.ndarray,
    length: float,
    decorrelation_time: float,
) -> float:
    """
    Computes the variance (as a share of the variance of the field) of the
    error of estimating the mean of a field over a period from its samples
    at given times: the mean of the samples, compared with the mean over the
    period. The field is stationary, with an exponential correlation in time
    (`exp(-|dt| / decorrelation_time)`, a first-order Markov process), for
    which the variance is exact::

        var(e) / var(x) = (1/n^2) sum_ij r(t_i - t_j)
                          - (2/(n L)) sum_i integral_0^L r(t_i - s) ds
                          + (1/L^2) integral_0^L integral_0^L r(s - s') ds ds'

    Samples farther apart than the decorrelation time are nearly independent,
    and the variance approaches 1/n; samples closer together repeat each
    other, and the variance falls less.

    Args:
        times (Sequence[float]): The times of the samples (seconds after the
            start of the period).
        length (float): The length of the period (seconds).
        decorrelation_time (float): The decorrelation time (seconds).

    Returns:
        float: the variance of the error, as a share of the field's variance
            (`NaN` without samples)
    """
    return _error_terms(times, length, decorrelation_time)[1]


def _error_terms(
    times: Sequence[float] | np.ndarray, length: float, decorrelation_time: float
) -> tuple[float, float]:
    """
    Computes the mean correlation between pairs of samples (the variance of
    their mean, as a share of the field's) and the variance of the sampling
    error (see `sampling_error_variance`).
    """
    t = np.asarray(times, dtype=float)
    if len(t) == 0:
        return np.nan, np.nan
    if length <= 0 or decorrelation_time <= 0:
        raise ValueError("the length and decorrelation time must be positive")
    scale = decorrelation_time
    pairs = np.exp(-np.abs(t[:, None] - t[None, :]) / scale).mean()
    cross = _exponential_integral(t, length, scale).sum() / (len(t) * length)
    within = 2 * (scale * length - scale**2 * (1 - np.exp(-length / scale))) / length**2
    return float(pairs), float(max(pairs - 2 * cross + within, 0.0))


def _error_terms_by_group(
    seconds: np.ndarray,
    codes: np.ndarray,
    groups: int,
    length: float,
    scales: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Computes the error terms (see `_error_terms`) of the samples of each
    group (such as a target) at once: the samples' times (seconds), the code
    of the group of each (0 to `groups` - 1), and the decorrelation time of
    each group. Gives the number of samples of each group, the mean
    correlation between pairs of its samples, and the variance of its
    sampling error (`NaN` without samples).

    The sum of the correlations between pairs of samples follows, for the
    samples of a group in order of time, from the recursion `S_j = exp(-(t_j
    - t_(j-1)) / scale) * (S_(j-1) + 1)` (the sum of the correlations of
    sample `j` with those before it), stepped for every group at once.
    """
    if length <= 0 or np.any(scales <= 0):
        raise ValueError("the length and decorrelation time must be positive")
    order = np.lexsort((seconds, codes))
    t, c = seconds[order], codes[order]
    counts = np.bincount(c, minlength=groups)
    first = np.concatenate([[0], np.cumsum(counts)[:-1]])
    position = np.arange(len(t)) - first[c]
    width = int(counts.max()) if len(counts) > 0 and len(t) > 0 else 0
    table = np.zeros((groups, max(width, 1)))
    table[c, position] = t
    recent = np.zeros(groups)  # S of the latest sample of each group
    total = np.zeros(groups)  # sum of S over the samples of each group
    for k in range(1, width):
        active = counts > k
        gap = table[active, k] - table[active, k - 1]
        recent[active] = np.exp(-gap / scales[active]) * (recent[active] + 1)
        total[active] += recent[active]
    with np.errstate(divide="ignore", invalid="ignore"):
        pairs = (counts + 2 * total) / counts.astype(float) ** 2
        cross = np.bincount(
            c, weights=_exponential_integral(t, length, scales[c]), minlength=groups
        ) / (counts * length)
        within = (
            2
            * (scales * length - scales**2 * (1 - np.exp(-length / scales)))
            / length**2
        )
        variance = np.maximum(pairs - 2 * cross + within, 0.0)
    empty = counts == 0
    pairs[empty] = np.nan
    variance[empty] = np.nan
    return counts, pairs, variance


def estimate_sampling_error(
    observations: gpd.GeoDataFrame,
    start: datetime,
    end: datetime,
    decorrelation_time: timedelta | Mapping[str, timedelta] | pd.Series,
    standard_deviation: float | Mapping[str, float] | pd.Series = 1.0,
    targets: gpd.GeoDataFrame | gpd.GeoSeries | None = None,
    time_column: str = "epoch",
) -> gpd.GeoDataFrame:
    """
    Estimates the sampling error of each target: the error of estimating the
    mean of a field over a period from the target's samples (observations)
    for a model of the field, stationary with a standard deviation and an
    exponential correlation in time with a decorrelation time (see
    `sampling_error_variance`). The decorrelation time of a field is the time
    over which its values become nearly independent: hours for convective
    precipitation, days for synoptic temperature.

    Gives for each target the number of `samples`, the standard deviation of
    the error (`sampling_error`, in the units of the field's standard
    deviation), the error relative to the field's standard deviation
    (`relative_error`), and the effective number of independent samples
    (`effective_samples`): the number of independent samples whose mean
    would vary as much as the mean of the samples (the inverse of the mean
    correlation between pairs of samples), which equals the number of
    samples if they are farther apart than the decorrelation time, and is
    less if they repeat each other. Targets without samples have no error
    estimate (`NaN`) and no effective samples.

    Args:
        observations (geopandas.GeoDataFrame): The observations (see
            `collect_observations`), with the `target_hash` and the time
            (`time_column`) of each sample.
        start (datetime.datetime): The start of the period.
        end (datetime.datetime): The end of the period.
        decorrelation_time (datetime.timedelta | Mapping | pandas.Series): The
            decorrelation time of the field, or of each target (by
            `target_hash`), as for a field that decorrelates faster in some
            places (such as over land) than others.
        standard_deviation (float | Mapping | pandas.Series): The standard
            deviation of the field, or of each target (by `target_hash`).
        targets (geopandas.GeoDataFrame | geopandas.GeoSeries | None): The
            targets, including those not observed, or `None` for those
            observed.
        time_column (str): The column of the times of the samples.

    Returns:
        geopandas.GeoDataFrame: a row per target
    """
    length = (end - start).total_seconds()
    result = _get_targets(observations, targets)
    seconds = _get_seconds(observations[time_column], start)
    codes, keys = pd.factorize(observations["target_hash"].to_numpy())
    # the decorrelation time of each observed target
    if isinstance(decorrelation_time, (Mapping, pd.Series)):
        scales = {
            key: pd.Timedelta(value).total_seconds()
            for key, value in dict(decorrelation_time).items()
        }
        per_group = np.array([scales[key] for key in keys], dtype=float)
    else:
        per_group = np.full(len(keys), decorrelation_time.total_seconds())
    counts, pairs, variance = _error_terms_by_group(
        seconds, codes, len(keys), length, per_group
    )
    # the group of each target (-1 for targets without samples)
    index = pd.Index(keys).get_indexer(result["target_hash"].to_numpy())
    found = index >= 0
    picked = np.where(found, index, 0)
    samples = np.where(found, counts[picked], 0).astype(int)
    relative_error = np.sqrt(np.where(found, variance[picked], np.nan))
    with np.errstate(divide="ignore"):
        effective = np.where(found, 1 / pairs[picked], 0.0)
    if isinstance(standard_deviation, (Mapping, pd.Series)):
        sigma = result["target_hash"].map(standard_deviation).to_numpy(dtype=float)
    else:
        sigma = float(standard_deviation)
    return result.assign(
        samples=samples,
        relative_error=relative_error,
        sampling_error=relative_error * sigma,
        effective_samples=effective,
    )


def compute_revisit_statistics(
    observations: gpd.GeoDataFrame,
    start: datetime,
    end: datetime,
    quantiles: Sequence[float] = (0.5, 0.9),
    targets: gpd.GeoDataFrame | gpd.GeoSeries | None = None,
    time_column: str = "epoch",
) -> gpd.GeoDataFrame:
    """
    Computes statistics of the revisits of each target over a period: the
    number of `samples`, the mean and quantiles of the times between
    consecutive samples (`revisit`, and `revisit p50`, `revisit p90`, ...),
    and the longest time without a sample (`max_gap`), including the times
    from the start of the period to the first sample and from the last
    sample to the end, which a mean revisit hides. Targets without samples
    have a gap of the whole period.

    Args:
        observations (geopandas.GeoDataFrame): The observations (see
            `collect_observations`), with the `target_hash` and the time
            (`time_column`) of each sample.
        start (datetime.datetime): The start of the period.
        end (datetime.datetime): The end of the period.
        quantiles (Sequence[float]): The quantiles of the revisits.
        targets (geopandas.GeoDataFrame | geopandas.GeoSeries | None): The
            targets, including those not observed, or `None` for those
            observed.
        time_column (str): The column of the times of the samples.

    Returns:
        geopandas.GeoDataFrame: a row per target
    """
    length = (end - start).total_seconds()
    result = _get_targets(observations, targets)
    seconds = pd.Series(
        _get_seconds(observations[time_column], start), index=observations.index
    )
    rows = {}
    for key, group in seconds.groupby(observations["target_hash"].to_numpy()):
        t = np.sort(group.to_numpy())
        gaps = np.diff(t)
        edges = np.diff(np.concatenate([[0.0], t, [length]]))
        row = {
            "samples": len(t),
            "revisit": gaps.mean() if len(gaps) else np.nan,
            "max_gap": edges.max(),
        }
        for q in quantiles:
            row[f"revisit p{round(100 * q)}"] = (
                np.quantile(gaps, q) if len(gaps) else np.nan
            )
        rows[key] = row
    table = pd.DataFrame.from_dict(rows, orient="index")
    columns = ["samples", "revisit"] + [f"revisit p{round(100 * q)}" for q in quantiles]
    for column in columns + ["max_gap"]:
        values = (
            result["target_hash"].map(table[column])
            if column in table
            else pd.Series(np.nan, index=result.index)
        )
        result[column] = values
    result["samples"] = result["samples"].fillna(0).astype(int)
    result["max_gap"] = result["max_gap"].fillna(length)
    for column in columns[1:] + ["max_gap"]:
        result[column] = pd.to_timedelta(result[column], unit="s")
    return result


def compute_local_time_coverage(
    observations: gpd.GeoDataFrame,
    bins: int = 8,
    targets: gpd.GeoDataFrame | gpd.GeoSeries | None = None,
    time_column: str = "epoch",
) -> gpd.GeoDataFrame:
    """
    Computes the distribution of the samples of each target in local (mean
    solar) time of day, which a field with a daily cycle (such as convective
    precipitation, or land surface temperature) needs to be sampled across:
    the number of samples in each of `bins` equal intervals of the day (`local
    00:00`, ...), the number of intervals sampled (`local_bins_sampled`), and
    their share (`local_time_coverage`). The local time of a sample is the
    universal time plus the target's longitude (15 degrees per hour).

    Args:
        observations (geopandas.GeoDataFrame): The observations (see
            `collect_observations`), with the `target_hash` and the time
            (`time_column`) of each sample.
        bins (int): The number of intervals of the day.
        targets (geopandas.GeoDataFrame | geopandas.GeoSeries | None): The
            targets, including those not observed, or `None` for those
            observed.
        time_column (str): The column of the times of the samples.

    Returns:
        geopandas.GeoDataFrame: a row per target
    """
    if bins < 1:
        raise ValueError("bins must be at least 1")
    result = _get_targets(observations, targets)
    longitude, _ = _get_points(observations.geometry)
    times = pd.to_datetime(observations[time_column], utc=True)
    hours = (
        times.dt.hour + times.dt.minute / 60 + times.dt.second / 3600
    ).to_numpy() + longitude / 15
    index = np.floor((hours % 24) / (24 / bins)).astype(int) % bins
    counts = (
        pd.crosstab(observations["target_hash"].to_numpy(), index)
        .reindex(columns=range(bins), fill_value=0)
        .reindex(result["target_hash"], fill_value=0)
    )
    width = 24 / bins
    for k in range(bins):
        hour = k * width
        label = f"local {int(hour):02d}:{round(60 * (hour % 1)):02d}"
        result[label] = counts[k].to_numpy()
    sampled = (counts.to_numpy() > 0).sum(axis=1)
    result["local_bins_sampled"] = sampled
    result["local_time_coverage"] = sampled / bins
    return result


def _get_dimension(field: Any, names: Sequence[str], given: str | None) -> str:
    """Gets the name of a dimension of a field."""
    if given is not None:
        return given
    for name in names:
        if name in field.dims:
            return name
    raise ValueError(f"the field has no dimension named any of {list(names)}")


def _field_times(field: Any, dimension: str, times: pd.Series) -> np.ndarray:
    """
    Converts times to those of the time coordinate of a field: naive (UTC) if
    the field's are.
    """
    times = pd.to_datetime(times, utc=True)
    if getattr(field[dimension].dtype, "tz", None) is None:
        times = times.dt.tz_convert("UTC").dt.tz_localize(None)
    return times.to_numpy()


def _field_longitudes(field: Any, dimension: str, longitude: np.ndarray) -> np.ndarray:
    """
    Converts longitudes to the convention of a field: [0, 360) if its
    longitudes exceed 180 degrees, else [-180, 180).
    """
    if float(field[dimension].max()) > 180:
        return np.mod(longitude, 360)
    return (longitude + 180) % 360 - 180


def sample_field(
    observations: gpd.GeoDataFrame,
    field: Any,
    name: str | None = None,
    method: str = "nearest",
    time_column: str = "epoch",
    dimensions: Mapping[str, str] | None = None,
) -> gpd.GeoDataFrame:
    """
    Samples a field (such as a variable of a nature run) at the time and place
    of each observation: at its point (or the representative point of its
    region), by the nearest value or by linear interpolation in time and
    space. Observations outside the period of the field have no value
    (`NaN`). Requires xarray.

    Args:
        observations (geopandas.GeoDataFrame): The observations, with the time
            (`time_column`) of each sample.
        field (xarray.DataArray): The field, with dimensions of time,
            latitude, and longitude (by default named `time`, `lat` or
            `latitude`, and `lon` or `longitude`), in degrees, with longitudes
            in [-180, 180) or [0, 360).
        name (str | None): The column of the sampled values, or `None` for the
            field's name (or `"value"`).
        method (str): `"nearest"` or `"linear"`.
        time_column (str): The column of the times of the samples.
        dimensions (Mapping[str, str] | None): The names of the field's
            dimensions of `time`, `lat`, and `lon`, if not the defaults.

    Returns:
        geopandas.GeoDataFrame: the observations, with the sampled values
    """
    # pylint: disable=import-outside-toplevel
    import xarray as xr

    if method not in ("nearest", "linear"):
        raise ValueError(f"method must be 'nearest' or 'linear', not {method!r}")
    dimensions = dict(dimensions or {})
    time_dim = _get_dimension(field, ["time"], dimensions.get("time"))
    lat_dim = _get_dimension(field, ["lat", "latitude"], dimensions.get("lat"))
    lon_dim = _get_dimension(field, ["lon", "longitude"], dimensions.get("lon"))
    name = name or field.name or "value"
    result = observations.copy()
    if len(observations) == 0:
        result[name] = pd.Series(dtype=float)
        return result
    longitude, latitude = _get_points(observations.geometry)
    where = {
        time_dim: xr.DataArray(
            _field_times(field, time_dim, observations[time_column]), dims="sample"
        ),
        lat_dim: xr.DataArray(latitude, dims="sample"),
        lon_dim: xr.DataArray(
            _field_longitudes(field, lon_dim, longitude), dims="sample"
        ),
    }
    if method == "nearest":
        values = field.sel(where, method="nearest").to_numpy().astype(float)
        times = where[time_dim].to_numpy()
        coordinate = field[time_dim].to_numpy()
        outside = (times < coordinate.min()) | (times > coordinate.max())
        values[outside] = np.nan
    else:
        values = field.interp(where).to_numpy().astype(float)
    result[name] = values
    return result


def compare_sampled_means(
    sampled: gpd.GeoDataFrame,
    field: Any,
    start: datetime,
    end: datetime,
    name: str | None = None,
    targets: gpd.GeoDataFrame | gpd.GeoSeries | None = None,
    method: str = "nearest",
    dimensions: Mapping[str, str] | None = None,
) -> gpd.GeoDataFrame:
    """
    Compares the mean of each target's samples of a field (see
    `sample_field`) with the true mean of the field at the target over a
    period (the mean of its values at the target's place, nearest or
    interpolated, at every time of the field in the period): the number of
    `samples`, the `sampled_mean`, the `true_mean`, and the `error` (sampled
    less true), which measures the sampling error against a simulation of
    the field. Targets without samples have no sampled mean (`NaN`). Requires
    xarray.

    Args:
        sampled (geopandas.GeoDataFrame): The observations with sampled
            values.
        field (xarray.DataArray): The field (see `sample_field`).
        start (datetime.datetime): The start of the period.
        end (datetime.datetime): The end of the period.
        name (str | None): The column of the sampled values, or `None` for the
            field's name (or `"value"`).
        targets (geopandas.GeoDataFrame | geopandas.GeoSeries | None): The
            targets, including those not observed, or `None` for those
            observed.
        method (str): `"nearest"` or `"linear"`, for the field at the targets.
        dimensions (Mapping[str, str] | None): The names of the field's
            dimensions of `time`, `lat`, and `lon`, if not the defaults.

    Returns:
        geopandas.GeoDataFrame: a row per target
    """
    # pylint: disable=import-outside-toplevel
    import xarray as xr

    name = name or field.name or "value"
    dimensions = dict(dimensions or {})
    time_dim = _get_dimension(field, ["time"], dimensions.get("time"))
    lat_dim = _get_dimension(field, ["lat", "latitude"], dimensions.get("lat"))
    lon_dim = _get_dimension(field, ["lon", "longitude"], dimensions.get("lon"))
    result = _get_targets(sampled, targets)
    window = _field_times(field, time_dim, pd.Series([start, end]))
    period = field.sel({time_dim: slice(window[0], window[1])})
    longitude, latitude = _get_points(result.geometry)
    where = {
        lat_dim: xr.DataArray(latitude, dims="target"),
        lon_dim: xr.DataArray(
            _field_longitudes(field, lon_dim, longitude), dims="target"
        ),
    }
    at_targets = (
        period.sel(where, method="nearest")
        if method == "nearest"
        else period.interp(where)
    )
    result["true_mean"] = at_targets.mean(dim=time_dim).to_numpy().astype(float)
    groups = sampled.groupby("target_hash")[name]
    result["samples"] = result["target_hash"].map(groups.count()).fillna(0).astype(int)
    result["sampled_mean"] = result["target_hash"].map(groups.mean()).astype(float)
    result["error"] = result["sampled_mean"] - result["true_mean"]
    return result

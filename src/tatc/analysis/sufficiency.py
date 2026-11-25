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

from .check import _get_observation_times
from .coverage_metrics import _get_naive_utc, _get_targets, _group_overlapping


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


def _get_samples(
    observations: pd.DataFrame,
    time_column: str | None,
    merge_overlapping: bool,
    value: str | None = None,
    sample_interval: timedelta | None = None,
) -> pd.DataFrame:
    """
    Gets the samples of observations: the `target_hash`, time (see
    `_get_observation_times`), longitude, and latitude of each (and a value,
    if named). With `merge_overlapping`, the overlapping observations of a
    target (by their `start` and `end`, if given, as by several platforms at
    once) are one sample, as `tatc.analysis.aggregate_observations` merges
    them: at the midpoint of the first and last of their times, with the
    mean of their values. With `sample_interval`, an observation (or merged
    observations) whose period (`start` to `end`) is longer is sampled every
    interval from its time, earlier and later, within its period.

    Args:
        observations (pandas.DataFrame): The observations.
        time_column (str | None): The column of the times, or `None` to find
            it.
        merge_overlapping (bool): `True` to merge overlapping observations.
        value (str | None): The column of a value of each observation.
        sample_interval (datetime.timedelta | None): The interval at which to
            sample long observations, or `None` for one sample each.

    Returns:
        pandas.DataFrame: the samples
    """
    if sample_interval is not None and sample_interval <= timedelta(0):
        raise ValueError(f"sample_interval must be positive, not {sample_interval}")
    longitude, latitude = _get_points(observations.geometry)
    samples = pd.DataFrame(
        {
            "target_hash": observations["target_hash"].to_numpy(),
            "time": _get_observation_times(observations, time_column).to_numpy(),
            "longitude": longitude,
            "latitude": latitude,
            **({} if value is None else {"value": observations[value].to_numpy()}),
        }
    )
    periods = (
        len(samples) > 0
        and "start" in observations.columns
        and "end" in observations.columns
    )
    if not periods:
        return samples
    samples["start"] = _get_naive_utc(observations["start"])
    samples["end"] = _get_naive_utc(observations["end"])
    if merge_overlapping:
        target, _ = pd.factorize(samples["target_hash"])
        order, group = _group_overlapping(
            target, samples["start"].to_numpy(), samples["end"].to_numpy()
        )
        grouped = samples.iloc[order].reset_index(drop=True).groupby(group, sort=False)
        first, last = grouped["time"].min(), grouped["time"].max()
        samples = pd.DataFrame(
            {
                "target_hash": grouped["target_hash"].first(),
                "time": first + (last - first) / 2,
                "longitude": grouped["longitude"].first(),
                "latitude": grouped["latitude"].first(),
                **({} if value is None else {"value": grouped["value"].mean()}),
                "start": grouped["start"].min(),
                "end": grouped["end"].max(),
            }
        ).reset_index(drop=True)
    if sample_interval is not None:
        samples = _sample_periods(samples, sample_interval)
    return samples.drop(columns=["start", "end"])


def _sample_periods(samples: pd.DataFrame, interval: timedelta) -> pd.DataFrame:
    """
    Samples the periods of samples (their naive UTC `start` and `end`) every
    interval from their times, earlier and later, within their periods (see
    `_get_samples`): each sample is repeated at each of those times.

    Args:
        samples (pandas.DataFrame): The samples, with a `time`, `start`, and
            `end` each.
        interval (datetime.timedelta): The interval.

    Returns:
        pandas.DataFrame: the samples
    """
    step = pd.Timedelta(interval).value
    time = pd.DatetimeIndex(pd.to_datetime(samples["time"], utc=True))
    epoch = time.tz_localize(None).as_unit("ns").asi8
    start = pd.DatetimeIndex(samples["start"]).as_unit("ns").asi8
    end = pd.DatetimeIndex(samples["end"]).as_unit("ns").asi8
    # the multiples of the interval from the time within the period (at
    # least the time itself), none for a sample without a period
    known = ~(pd.isna(samples["start"]).to_numpy() | pd.isna(samples["end"]).to_numpy())
    first = np.where(known, np.minimum(-((epoch - start) // step), 0), 0)
    last = np.where(known, np.maximum((end - epoch) // step, 0), 0)
    counts = (last - first + 1).astype(int)
    index = np.repeat(np.arange(len(samples)), counts)
    offset = np.arange(np.sum(counts)) - np.repeat(np.cumsum(counts) - counts, counts)
    multiple = first[index] + offset
    expanded = samples.iloc[index].reset_index(drop=True)
    expanded["time"] = time[index] + pd.to_timedelta(multiple * step, unit="ns")
    return expanded


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
    sample `j` with those before it), stepped for every group at once: step
    `k` takes the `k`-th sample of each group with more than `k` samples, so
    that the steps take each sample once, in memory proportional to the
    number of samples.
    """
    if length <= 0 or np.any(scales <= 0):
        raise ValueError("the length and decorrelation time must be positive")
    order = np.lexsort((seconds, codes))
    t, c = seconds[order], codes[order]
    counts = np.bincount(c, minlength=groups)
    first = np.concatenate([[0], np.cumsum(counts)[:-1]])
    position = np.arange(len(t)) - first[c]
    # the samples (indices in t) in order of their position in their group
    by_position = np.argsort(position, kind="stable")
    bounds = np.concatenate([[0], np.cumsum(np.bincount(position))])
    recent = np.zeros(groups)  # S of the latest sample of each group
    total = np.zeros(groups)  # sum of S over the samples of each group
    for k in range(1, len(bounds) - 1):
        step = by_position[bounds[k] : bounds[k + 1]]
        group = c[step]
        gap = t[step] - t[step - 1]
        recent[group] = np.exp(-gap / scales[group]) * (recent[group] + 1)
        total[group] += recent[group]
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
    time_column: str | None = None,
    merge_overlapping: bool = False,
    sample_interval: timedelta | None = None,
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
            (`time_column`) of each.
        start (datetime.datetime): The start of the period.
        end (datetime.datetime): The end of the period.
        decorrelation_time (datetime.timedelta | Mapping | pandas.Series): The
            decorrelation time of the field, or of each target (by
            `target_hash`), as for a field that decorrelates faster in some
            places (such as over land) than others.
        standard_deviation (float | Mapping | pandas.Series): The standard
            deviation of the field, or of each target (by `target_hash`).
        targets (geopandas.GeoDataFrame | geopandas.GeoSeries | None): The
            targets, including those not observed, with the other columns of
            a data frame of targets (such as a `point_id` or a weight,
            except those named as the results), or `None` for those
            observed.
        time_column (str | None): The column of the times of the samples,
            or `None` (default) for the first of `epoch`, `time`, or
            `closest_time`, or else the midpoint of `start` and `end` (as of
            the observations of regions).
        merge_overlapping (bool): `False` (default) to count each
            observation as a sample, as a study that averages every
            observation does (samples that repeat each other are accounted
            for by their correlation); `True` to count the overlapping
            observations of a target (as by several platforms at once) as one
            sample, at the midpoint of their times, as `aggregate_observations`
            merges them (so that `samples` counts the same looks as
            `reduce_observations`).
        sample_interval (datetime.timedelta | None): The interval at which
            to sample observations whose periods (`start` to `end`) are
            longer, such as those of a fixed platform that observes a
            target continuously (see `collect_footprint_observations`):
            every interval from the time of the observation, earlier and
            later, within its period; or `None` (default) for one sample
            per observation (or merged observations), however long.

    Returns:
        geopandas.GeoDataFrame: a row per target
    """
    length = (end - start).total_seconds()
    result = _get_targets(observations, targets)
    samples = _get_samples(
        observations, time_column, merge_overlapping, sample_interval=sample_interval
    )
    seconds = _get_seconds(samples["time"], start)
    codes, keys = pd.factorize(samples["target_hash"].to_numpy())
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

    def per_target(values: np.ndarray, missing: float) -> np.ndarray:
        # (without samples, there are no groups to pick from)
        if len(keys) == 0:
            return np.full(len(index), missing, dtype=float)
        return np.where(found, values[picked], missing)

    sample_counts = per_target(counts, 0).astype(int)
    relative_error = np.sqrt(per_target(variance, np.nan))
    with np.errstate(divide="ignore"):
        effective = per_target(1 / pairs, 0.0)
    if isinstance(standard_deviation, (Mapping, pd.Series)):
        sigma = result["target_hash"].map(standard_deviation).to_numpy(dtype=float)
    else:
        sigma = float(standard_deviation)
    return result.assign(
        samples=sample_counts,
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
    time_column: str | None = None,
    merge_overlapping: bool = False,
    sample_interval: timedelta | None = None,
) -> gpd.GeoDataFrame:
    """
    Computes statistics of the revisits of each target over a period: the
    number of `samples`, the mean and quantiles of the times between
    consecutive samples (`revisit`, and `revisit p50`, `revisit p90`, ...),
    and the longest time without a sample (`max_gap`), including the times
    from the start of the period to the first sample and from the last
    sample to the end, which a mean revisit hides. Targets without samples
    have a gap of the whole period.

    Each observation is a sample, by default, and its revisit is the time
    from the previous sample's time (such as its `epoch`) to its own, so that
    clusters of samples (as of satellites minutes apart) show in the
    quantiles. `reduce_observations` instead merges overlapping observations
    (see `aggregate_observations`) and measures each revisit from the end of
    one merged period to the start of the next; `merge_overlapping` counts
    samples as it does.

    Args:
        observations (geopandas.GeoDataFrame): The observations (see
            `collect_observations`), with the `target_hash` and the time
            (`time_column`) of each.
        start (datetime.datetime): The start of the period.
        end (datetime.datetime): The end of the period.
        quantiles (Sequence[float]): The quantiles of the revisits.
        targets (geopandas.GeoDataFrame | geopandas.GeoSeries | None): The
            targets, including those not observed, with the other columns of
            a data frame of targets (such as a `point_id` or a weight,
            except those named as the results), or `None` for those
            observed.
        time_column (str | None): The column of the times of the samples,
            or `None` (default) for the first of `epoch`, `time`, or
            `closest_time`, or else the midpoint of `start` and `end` (as of
            the observations of regions).
        merge_overlapping (bool): `False` (default) to count each
            observation as a sample, as a study that averages every
            observation does (samples that repeat each other are accounted
            for by their correlation); `True` to count the overlapping
            observations of a target (as by several platforms at once) as one
            sample, at the midpoint of their times, as `aggregate_observations`
            merges them (so that `samples` counts the same looks as
            `reduce_observations`).
        sample_interval (datetime.timedelta | None): The interval at which
            to sample observations whose periods (`start` to `end`) are
            longer, such as those of a fixed platform that observes a
            target continuously (see `collect_footprint_observations`):
            every interval from the time of the observation, earlier and
            later, within its period; or `None` (default) for one sample
            per observation (or merged observations), however long.

    Returns:
        geopandas.GeoDataFrame: a row per target
    """
    length = (end - start).total_seconds()
    result = _get_targets(observations, targets)
    samples = _get_samples(
        observations, time_column, merge_overlapping, sample_interval=sample_interval
    )
    seconds = _get_seconds(samples["time"], start)
    codes, keys = pd.factorize(samples["target_hash"].to_numpy())
    groups = len(keys)
    # the samples of each target in order of time, and the gaps between them
    order = np.lexsort((seconds, codes))
    t, c = seconds[order], codes[order]
    counts = np.bincount(c, minlength=groups)
    ends = np.cumsum(counts)
    same = c[1:] == c[:-1]
    gaps, gap_codes = (t[1:] - t[:-1])[same], c[1:][same]
    gap_counts = np.bincount(gap_codes, minlength=groups)
    with np.errstate(invalid="ignore", divide="ignore"):
        revisit = np.bincount(gap_codes, weights=gaps, minlength=groups) / gap_counts
    # the longest time without a sample, from the start to the end of the period
    longest = np.full(groups, -np.inf)
    if len(t) > 0:
        longest = np.maximum(t[ends - counts], length - t[ends - 1])
        np.maximum.at(longest, gap_codes, gaps)
    # quantiles of the gaps of each target (by linear interpolation, as numpy)
    gap_order = np.lexsort((gaps, gap_codes))
    ordered = gaps[gap_order]
    gap_starts = np.cumsum(gap_counts) - gap_counts
    has_gaps = gap_counts > 0
    revisits = {}
    for q in quantiles:
        position = q * np.maximum(gap_counts - 1, 0)
        below = np.floor(position).astype(int)
        above = np.minimum(below + 1, np.maximum(gap_counts - 1, 0))
        low = np.where(has_gaps, gap_starts + below, 0)
        high = np.where(has_gaps, gap_starts + above, 0)
        if len(ordered) > 0:
            value = ordered[low] + (ordered[high] - ordered[low]) * (position - below)
        else:
            value = np.zeros(groups)
        revisits[f"revisit p{round(100 * q)}"] = np.where(has_gaps, value, np.nan)
    # the statistics of every target (without samples: none, and the whole period)
    index = pd.Index(keys).get_indexer(result["target_hash"].to_numpy())
    found = index >= 0
    picked = np.where(found, index, 0)

    def per_target(values: np.ndarray, missing: float) -> np.ndarray:
        if groups == 0:
            return np.full(len(index), missing)
        return np.where(found, values[picked], missing)

    columns = {
        "samples": per_target(counts, 0).astype(int),
        "revisit": pd.to_timedelta(per_target(revisit, np.nan), unit="s"),
        **{
            name: pd.to_timedelta(per_target(values, np.nan), unit="s")
            for name, values in revisits.items()
        },
        "max_gap": pd.to_timedelta(per_target(longest, length), unit="s"),
    }
    return result.assign(**columns)


def compute_local_time_coverage(
    observations: gpd.GeoDataFrame,
    bins: int = 8,
    targets: gpd.GeoDataFrame | gpd.GeoSeries | None = None,
    time_column: str | None = None,
    merge_overlapping: bool = False,
    sample_interval: timedelta | None = None,
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
            (`time_column`) of each.
        bins (int): The number of intervals of the day.
        targets (geopandas.GeoDataFrame | geopandas.GeoSeries | None): The
            targets, including those not observed, with the other columns of
            a data frame of targets (such as a `point_id` or a weight,
            except those named as the results), or `None` for those
            observed.
        time_column (str | None): The column of the times of the samples,
            or `None` (default) for the first of `epoch`, `time`, or
            `closest_time`, or else the midpoint of `start` and `end` (as of
            the observations of regions).
        merge_overlapping (bool): `False` (default) to count each
            observation as a sample, as a study that averages every
            observation does (samples that repeat each other are accounted
            for by their correlation); `True` to count the overlapping
            observations of a target (as by several platforms at once) as one
            sample, at the midpoint of their times, as `aggregate_observations`
            merges them (so that `samples` counts the same looks as
            `reduce_observations`).
        sample_interval (datetime.timedelta | None): The interval at which
            to sample observations whose periods (`start` to `end`) are
            longer, such as those of a fixed platform that observes a
            target continuously (see `collect_footprint_observations`):
            every interval from the time of the observation, earlier and
            later, within its period; or `None` (default) for one sample
            per observation (or merged observations), however long.

    Returns:
        geopandas.GeoDataFrame: a row per target
    """
    if bins < 1:
        raise ValueError("bins must be at least 1")
    result = _get_targets(observations, targets)
    samples = _get_samples(
        observations, time_column, merge_overlapping, sample_interval=sample_interval
    )
    times = pd.to_datetime(samples["time"], utc=True)
    hours = (
        times.dt.hour + times.dt.minute / 60 + times.dt.second / 3600
    ).to_numpy() + samples["longitude"].to_numpy() / 15
    index = np.floor((hours % 24) / (24 / bins)).astype(int) % bins
    counts = (
        pd.crosstab(samples["target_hash"].to_numpy(), index)
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


def _coordinate_extent(values: np.ndarray) -> tuple[float, float]:
    """
    Gets the extent of a coordinate: half a spacing beyond its outer values
    (the cells of which they are the centers), or unbounded for one value.
    """
    values = np.unique(np.asarray(values, dtype=float))
    if len(values) < 2:
        return -np.inf, np.inf
    return (
        values[0] - (values[1] - values[0]) / 2,
        values[-1] + (values[-1] - values[-2]) / 2,
    )


def _outside_field(
    field: Any,
    lat_dim: str,
    lon_dim: str,
    latitude: np.ndarray,
    longitude: np.ndarray,
) -> np.ndarray:
    """
    Finds the places outside the extent of a field (see `_coordinate_extent`),
    as of a regional nature run, in latitude and, unless the field's
    longitudes go around the globe, in longitude (in the field's convention,
    see `_field_longitudes`).
    """
    low, high = _coordinate_extent(field[lat_dim].to_numpy())
    outside = (latitude < low) | (latitude > high)
    low, high = _coordinate_extent(field[lon_dim].to_numpy())
    if high - low < 360 - 1e-9:
        outside |= (longitude < low) | (longitude > high)
    return outside


def _sample_values(
    field: Any,
    times: pd.Series,
    longitude: np.ndarray,
    latitude: np.ndarray,
    method: str,
    dimensions: Mapping[str, str] | None,
) -> np.ndarray:
    """
    Samples a field at times and places (see `sample_field`): its nearest
    value (`NaN` outside the period or extent of the field) or its linear
    interpolation in time and space.

    Args:
        field (xarray.DataArray): The field.
        times (pandas.Series): The times of the samples.
        longitude (numpy.ndarray): The longitudes (degrees) of the samples.
        latitude (numpy.ndarray): The latitudes (degrees) of the samples.
        method (str): `"nearest"` or `"linear"`.
        dimensions (Mapping[str, str] | None): The names of the field's
            dimensions of `time`, `lat`, and `lon`, if not the defaults.

    Returns:
        numpy.ndarray: the values
    """
    # pylint: disable=import-outside-toplevel
    import xarray as xr

    if method not in ("nearest", "linear"):
        raise ValueError(f"method must be 'nearest' or 'linear', not {method!r}")
    dimensions = dict(dimensions or {})
    time_dim = _get_dimension(field, ["time"], dimensions.get("time"))
    lat_dim = _get_dimension(field, ["lat", "latitude"], dimensions.get("lat"))
    lon_dim = _get_dimension(field, ["lon", "longitude"], dimensions.get("lon"))
    where = {
        time_dim: xr.DataArray(
            _field_times(field, time_dim, pd.Series(times)), dims="sample"
        ),
        lat_dim: xr.DataArray(latitude, dims="sample"),
        lon_dim: xr.DataArray(
            _field_longitudes(field, lon_dim, longitude), dims="sample"
        ),
    }
    if method == "linear":
        return field.interp(where).to_numpy().astype(float)
    values = field.sel(where, method="nearest").to_numpy().astype(float)
    sample_times = where[time_dim].to_numpy()
    coordinate = field[time_dim].to_numpy()
    values[(sample_times < coordinate.min()) | (sample_times > coordinate.max())] = (
        np.nan
    )
    values[
        _outside_field(field, lat_dim, lon_dim, latitude, where[lon_dim].to_numpy())
    ] = np.nan
    return values


def sample_field(
    observations: gpd.GeoDataFrame,
    field: Any,
    name: str | None = None,
    method: str = "nearest",
    time_column: str | None = None,
    dimensions: Mapping[str, str] | None = None,
) -> gpd.GeoDataFrame:
    """
    Samples a field (such as a variable of a nature run) at the time and place
    of each observation: at its point (or the representative point of its
    region), by the nearest value or by linear interpolation in time and
    space. Observations outside the period of the field, or outside its
    extent (as of a regional nature run: by more than half a grid spacing
    beyond its outer latitudes or, unless they go around the globe,
    longitudes), have no value (`NaN`). Requires xarray.

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
        time_column (str | None): The column of the times of the samples,
            or `None` (default) for the first of `epoch`, `time`, or
            `closest_time`, or else the midpoint of `start` and `end`.
        dimensions (Mapping[str, str] | None): The names of the field's
            dimensions of `time`, `lat`, and `lon`, if not the defaults.

    Returns:
        geopandas.GeoDataFrame: the observations, with the sampled values
    """
    if method not in ("nearest", "linear"):
        raise ValueError(f"method must be 'nearest' or 'linear', not {method!r}")
    name = name or field.name or "value"
    result = observations.copy()
    if len(observations) == 0:
        result[name] = pd.Series(dtype=float)
        return result
    longitude, latitude = _get_points(observations.geometry)
    result[name] = _sample_values(
        field,
        _get_observation_times(observations, time_column),
        longitude,
        latitude,
        method,
        dimensions,
    )
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
    merge_overlapping: bool = False,
    time_column: str | None = None,
    sample_interval: timedelta | None = None,
) -> gpd.GeoDataFrame:
    """
    Compares the mean of each target's samples of a field (see
    `sample_field`) with the true mean of the field at the target over a
    period (the mean of its values at the target's place, nearest or
    interpolated, at every time of the field in the period): the number of
    `samples`, the `sampled_mean`, the `true_mean`, and the `error` (sampled
    less true), which measures the sampling error against a simulation of
    the field. Targets without samples have no sampled mean (`NaN`), and
    targets outside the extent of the field (see `sample_field`) no true
    mean. Requires xarray.

    Args:
        sampled (geopandas.GeoDataFrame): The observations with sampled
            values.
        field (xarray.DataArray): The field (see `sample_field`).
        start (datetime.datetime): The start of the period.
        end (datetime.datetime): The end of the period.
        name (str | None): The column of the sampled values, or `None` for the
            field's name (or `"value"`).
        targets (geopandas.GeoDataFrame | geopandas.GeoSeries | None): The
            targets, including those not observed, with the other columns of
            a data frame of targets (such as a `point_id` or a weight,
            except those named as the results), or `None` for those
            observed.
        method (str): `"nearest"` or `"linear"`, for the field at the targets.
        dimensions (Mapping[str, str] | None): The names of the field's
            dimensions of `time`, `lat`, and `lon`, if not the defaults.
        merge_overlapping (bool): `False` (default) to count each sample;
            `True` to count the overlapping samples of a target (as by
            several platforms at once) as one sample, with the mean of their
            values (see `estimate_sampling_error`).
        time_column (str | None): The column of the times of the samples,
            or `None` (default) for the first of `epoch`, `time`, or
            `closest_time`, or else the midpoint of `start` and `end`.
        sample_interval (datetime.timedelta | None): The interval at which
            to sample observations longer than it (as the continuous
            observations of fixed platforms), as for
            `estimate_sampling_error`, with the field sampled again (by
            `method`) at each of their times, or `None` (default) for one
            sample (the sampled value) of each observation.

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
    true_mean = at_targets.mean(dim=time_dim).to_numpy().astype(float)
    if method == "nearest":
        true_mean[
            _outside_field(field, lat_dim, lon_dim, latitude, where[lon_dim].to_numpy())
        ] = np.nan
    result["true_mean"] = true_mean
    samples = _get_samples(
        sampled, time_column, merge_overlapping, name, sample_interval
    )
    if sample_interval is not None and len(samples) > 0:
        # the field at each time of a long observation
        samples["value"] = _sample_values(
            field,
            samples["time"],
            samples["longitude"].to_numpy(),
            samples["latitude"].to_numpy(),
            method,
            dimensions,
        )
    groups = samples.groupby("target_hash")["value"]
    result["samples"] = result["target_hash"].map(groups.count()).fillna(0).astype(int)
    result["sampled_mean"] = result["target_hash"].map(groups.mean()).astype(float)
    result["error"] = result["sampled_mean"] - result["true_mean"]
    return result

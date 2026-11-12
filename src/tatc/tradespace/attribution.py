"""
Attribution of the variance of the outcomes of designs to their uncertain
parameters: which uncertainties spread the outcomes of a design between its
replicates, found by drawing each again while holding the others.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from itertools import combinations

import numpy as np
import pandas as pd

from .space import Design, DesignSpace, _check_uncertain, _designs, _generic
from .study import Study

_SHARES = ("first", "total")
"""The shares of variance of a parameter (see `attribute_variance`)."""


@dataclass
class Attribution:
    """
    The attribution of the variance of the metrics of designs between their
    replicates to their uncertain parameters (see `attribute_variance`).
    """

    results: pd.DataFrame
    """pandas.DataFrame: The results of every evaluation (see `Study.run`):
    the plain replicates, and the mixed replicates, whose `redrawn` column
    names the parameter (or group of parameters) drawn as in other
    replicates (`None` for plain replicates). Those of differences (see
    `differences`) are the differences of the metrics of two designs
    (`design` less `other`) in each replicate."""
    shares: pd.DataFrame
    """pandas.DataFrame: A row per design (and group of `by`), metric
    (`metric`), and parameter (`parameter`): the variance of the metric
    between the design's replicates (`variance`), the share of the variance
    that knowing the parameter would remove (`first`), the share that the
    parameter causes, alone or with others (`total`), and the standard
    errors of the shares (`first se` and `total se`, by a bootstrap). Those
    of differences (see `differences`) are of a pair of designs (`design`
    and `other`)."""
    by: tuple[str, ...] = ()
    """tuple[str, ...]: The columns that separate the records of an
    evaluation (see `attribute_variance`)."""
    samples: dict[str, np.ndarray] = field(default_factory=dict, repr=False)
    """dict[str, numpy.ndarray]: The variances that the parameter of each row
    of `shares` explains (`first` and `total`), and the variance of the
    metric (`variance`), in each resample of the bootstrap (a row per row of
    `shares`, and a column per resample, of the same replicates for every
    design), from which `summary` estimates the standard errors of shares
    over designs."""
    resamples: np.ndarray | None = field(default=None, repr=False)
    """numpy.ndarray | None: The rows of A, B, and the mixed replicates in
    each resample of the bootstrap (a row per resample), which `differences`
    resamples too."""

    def summary(self, share: str = "total") -> pd.DataFrame:
        """
        Gets the shares of the variance of each metric that each parameter
        explains over every design (and for each group of `by`): the
        variances that it explains, summed over the designs in which it is
        drawn, as a share of their variances (a mean of the designs' shares,
        weighted by their variances), with their standard errors (by the
        bootstrap of the shares, whose resamples are of the same replicates
        for every design).

        Args:
            share (str): `"total"` for the shares that the parameters cause,
                alone or with others, or `"first"` for the shares that knowing
                them would remove.

        Returns:
            pandas.DataFrame: a row per parameter (and group of `by`), and a
                column per metric, then a column of the standard errors of
                each (`<metric> se`)
        """
        if share not in _SHARES:
            raise ValueError(f"share must be 'first' or 'total', not {share!r}")
        frame = self.shares
        metrics = list(pd.unique(frame["metric"]))
        columns = [*metrics, *(f"{metric} se" for metric in metrics)]
        if frame.empty:
            return pd.DataFrame(columns=columns)
        index = [*self.by, "parameter"]
        usable = (frame[share].notna() & frame["variance"].notna()).to_numpy()
        explained = np.where(usable, frame[share] * frame["variance"], 0.0)
        variance = np.where(usable, frame["variance"], 0.0)
        sampled = np.where(usable[:, None], self.samples[share], np.nan)
        sampled_variance = np.where(usable[:, None], self.samples["variance"], np.nan)
        records = []
        with warnings.catch_warnings():
            # parameters drawn in no design, or without variance
            warnings.simplefilter("ignore", RuntimeWarning)
            groups = frame.groupby([*index, "metric"], sort=False).indices
            for key, rows in groups.items():
                # the share in each resample, of the designs it is estimated for
                known = ~np.isnan(sampled[rows]) & ~np.isnan(sampled_variance[rows])
                resampled = np.where(known, sampled[rows], 0).sum(axis=0) / np.where(
                    known, sampled_variance[rows], 0
                ).sum(axis=0)
                records.append(
                    {
                        **dict(zip([*index, "metric"], key)),
                        "share": explained[rows].sum() / variance[rows].sum(),
                        "se": np.nanstd(resampled, ddof=1),
                    }
                )
        table = pd.DataFrame(records)
        result = pd.concat(
            [
                table.pivot(index=index, columns="metric", values="share"),
                table.pivot(index=index, columns="metric", values="se").add_suffix(
                    " se"
                ),
            ],
            axis=1,
        ).reindex(columns=columns)
        result.columns.name = None
        if not self.by:
            result = result.reindex(index=list(pd.unique(frame["parameter"])))
        return result

    def differences(
        self, pairs: Sequence[tuple[Design | str, Design | str]] | None = None
    ) -> Attribution:
        """
        Attributes the variance of the differences of the metrics between
        pairs of designs (the first less the other, in the same replicates)
        to the uncertain parameters, from the evaluations of the designs,
        without evaluating them again: the plain and mixed replicates of both
        designs are the same replicates, so that their differences are those
        of the difference.

        With common random numbers (see `Design.realize`), the parameters that
        act alike on both designs cancel in their difference, and those that
        act differently remain: the parameters with large shares of the
        variance of a difference are those whose noise keeps the order of
        the designs unsettled (such as their order on a front, see
        `pareto_probability`), which more replicates, or less of that noise,
        would settle. A parameter may spread the metrics of both designs but
        not their difference, such as the growth of a cost that both share.
        Without common random numbers, nothing cancels: the variance of a
        difference is the sum of the designs' variances.

        A design that does not draw a parameter (such as one of a branch it
        does not choose) is not evaluated for it, as it does not depend on
        it: its values with the parameter drawn again are those of A. The
        shares of a parameter that neither design draws are `NaN`, and the
        replicates in which either design is infeasible are left out.

        Args:
            pairs (Sequence[tuple[Design | str, Design | str]] | None): The
                pairs of designs (or their keys), or `None` for every pair of
                the designs, in their order.

        Returns:
            Attribution: the differences (`results`) and their shares, whose
                rows give both designs (`design` and `other`)

        Raises:
            ValueError: if a design of a pair is not attributed, or a pair is
                of a design with itself.
        """
        if "other" in self.shares.columns or self.resamples is None:
            raise ValueError("differences need an attribution of designs")
        designs = list(pd.unique(self.results["design"]))
        pairs = list(combinations(designs, 2)) if pairs is None else pairs
        keyed = [(getattr(a, "key", a), getattr(b, "key", b)) for a, b in pairs]
        unknown = sorted({key for pair in keyed for key in pair} - set(designs))
        if unknown:
            raise ValueError(f"designs that are not attributed: {unknown}")
        same = [a for a, b in keyed if a == b]
        if same:
            raise ValueError(f"pairs of a design with itself: {same}")
        metrics = list(pd.unique(self.shares["metric"]))
        names = list(pd.unique(self.shares["parameter"]))
        n = self.resamples.shape[1]
        tables = {
            key: _tables(self.results, key, metrics, names, self.by)
            for key in dict.fromkeys(key for pair in keyed for key in pair)
        }
        frames = []
        for first, other in keyed:
            for name in [None, *names]:
                a, b = tables[first].get(name), tables[other].get(name)
                if a is None and b is None:
                    continue  # neither design draws the parameter
                # a design that does not draw a parameter is as in A
                a = _first_replicates(tables[first][None], n) if a is None else a
                b = _first_replicates(tables[other][None], n) if b is None else b
                frames.append(
                    _difference(a, b, metrics).assign(
                        design=first, other=other, redrawn=name
                    )
                )
        columns = [*self.by, "replicate", *metrics, "feasible"]
        differences = (
            pd.concat(frames, ignore_index=True)
            if frames
            else pd.DataFrame(columns=[*columns, "design", "other", "redrawn"])
        )
        shares, samples = _shares(
            differences,
            ("design", "other"),
            keyed,
            names,
            metrics,
            self.by,
            self.resamples,
        )
        return Attribution(
            results=differences,
            shares=shares,
            by=self.by,
            samples=samples,
            resamples=self.resamples,
        )


def _groups_of(
    space: DesignSpace, parameters: Sequence[str] | Mapping[str, Sequence[str]] | None
) -> dict[str, list[str]]:
    """
    Gets the groups of uncertain parameters to attribute variance to, by
    name: each parameter of the design space, each given parameter, or the
    given groups of parameters.
    """
    if not space.uncertain_paths():
        raise ValueError("the design space has no uncertain parameters")
    if parameters is None:
        groups = {path: [path] for path in space.uncertain_paths()}
    elif isinstance(parameters, str):
        groups = {parameters: [parameters]}
    elif isinstance(parameters, Mapping):
        groups = {
            name: [paths] if isinstance(paths, str) else list(paths)
            for name, paths in parameters.items()
        }
    else:
        groups = {path: [path] for path in parameters}
    if not groups:
        raise ValueError("no parameters to attribute variance to")
    paths = [path for group in groups.values() for path in group]
    _check_uncertain(space, paths)
    repeated = sorted({path for path in paths if paths.count(path) > 1})
    if repeated:
        raise ValueError(f"parameters in more than one group: {repeated}")
    empty = [name for name, group in groups.items() if not group]
    if empty:
        raise ValueError(f"groups without parameters: {empty}")
    return groups


def _draws(design: Design, paths: Sequence[str]) -> bool:
    """
    Whether a design draws any of the uncertain parameters (by their paths,
    see `DesignSpace.uncertain_paths`), which those of branches it does not
    choose (or of the items of an empty repetition) it does not.
    """
    drawn = {_generic(path) for path in design.flat()}
    return any(path in drawn for path in paths)


def _estimate(
    plain: np.ndarray,
    mixed: np.ndarray,
    n: int,
    resamples: np.ndarray,
) -> dict[str, np.ndarray]:
    """
    Estimates the variance of a metric and the first and total shares of each
    group of parameters, from the metric in the plain replicates (`2 n`: the
    rows of A, then of B) and the mixed replicates of each group (a row per
    group, of `n`: the rows of A with the group drawn as in B), which are
    `NaN` where infeasible (or not evaluated), with the variances that the
    groups explain and the variance in each resample of a bootstrap
    (`resamples`, indices of rows, resamples x n), and the standard errors
    of the shares.

    The estimators are Jansen's: the total effect of a group is half the
    mean squared difference between A and its mixed replicates (which differ
    in the group alone), and the first-order effect is the variance less
    half the mean squared difference between B and the mixed replicates
    (which share the group alone).
    """
    a, b = plain[:n], plain[n:]
    # the rows in which A, B, and the mixed replicates are all feasible
    joint = ~np.isnan(a) & ~np.isnan(b) & ~np.isnan(mixed)
    halves = {
        "first": np.where(joint, (b - mixed) ** 2, np.nan) / 2,
        "total": np.where(joint, (a - mixed) ** 2, np.nan) / 2,
    }
    values = plain[~np.isnan(plain)]
    variance = float(np.var(values, ddof=1)) if len(values) > 1 else np.nan
    # shares need a variance, and two rows feasible in every replicate
    valid = (joint.sum(axis=1) >= 2) & (variance > 0)
    estimates: dict = {"variance": variance}
    with warnings.catch_warnings():
        # designs that are infeasible (or groups not drawn) give empty slices
        warnings.simplefilter("ignore", RuntimeWarning)
        both = np.concatenate([a[resamples], b[resamples]], axis=1)
        sampled_variance = np.nanvar(both, axis=1, ddof=1)
        sampled_variance = np.where(sampled_variance > 0, sampled_variance, np.nan)
        estimates["variance samples"] = sampled_variance
        for share, terms in halves.items():
            explained = np.nanmean(terms, axis=1)
            sampled = np.nanmean(terms[:, resamples], axis=2)
            if share == "first":
                explained, sampled = variance - explained, sampled_variance - sampled
            estimates[share] = np.where(valid, explained / variance, np.nan)
            estimates[f"{share} samples"] = np.where(valid[:, None], sampled, np.nan)
            estimates[f"{share} se"] = np.where(
                valid, np.nanstd(sampled / sampled_variance, axis=1, ddof=1), np.nan
            )
    return estimates


def _tables(
    results: pd.DataFrame,
    key: str,
    metrics: Sequence[str],
    names: Sequence[str],
    by: Sequence[str],
) -> dict[str | None, pd.DataFrame]:
    """
    Gets the metrics of a design (`NaN` where infeasible) and whether it is
    feasible, by (the groups of `by` and) replicate: in its plain replicates
    (`None`), and in its mixed replicates of each group of parameters it was
    evaluated for.
    """
    rows = results[results["design"] == key]

    def table(part: pd.DataFrame) -> pd.DataFrame:
        indexed = part.set_index([*by, "replicate"])
        feasible = indexed["feasible"].astype(bool).to_numpy()
        values = indexed.reindex(columns=list(metrics)).to_numpy(dtype=float)
        values[~feasible] = np.nan
        frame = pd.DataFrame(values, index=indexed.index, columns=list(metrics))
        frame["feasible"] = feasible
        return frame

    tables: dict[str | None, pd.DataFrame] = {None: table(rows[rows["redrawn"].isna()])}
    for name in names:
        part = rows[rows["redrawn"] == name]
        if len(part):
            tables[name] = table(part)
    return tables


def _first_replicates(table: pd.DataFrame, n: int) -> pd.DataFrame:
    """
    Gets the rows of a table (see `_tables`) of the first `n` replicates (A).
    """
    return table[table.index.get_level_values("replicate") < n]


def _difference(
    a: pd.DataFrame, b: pd.DataFrame, metrics: Sequence[str]
) -> pd.DataFrame:
    """
    Gets the differences of the metrics of two tables (see `_tables`) in the
    same rows, which are feasible if both are.
    """
    index = a.index.union(b.index)
    difference = a.reindex(index)[list(metrics)] - b.reindex(index)[list(metrics)]
    difference["feasible"] = a["feasible"].reindex(index, fill_value=False) & b[
        "feasible"
    ].reindex(index, fill_value=False)
    return difference.reset_index()


def _by_replicate(
    group: pd.DataFrame, metric: str, rows: pd.Series, size: int
) -> np.ndarray:
    """
    Gets a metric of the rows of a group of results in the order of their
    replicates, up to `size` (`NaN` where a replicate is missing).
    """
    values = group.loc[rows, metric].astype(float)
    return values.set_axis(group.loc[rows, "replicate"]).reindex(range(size)).to_numpy()


def _shares(
    results: pd.DataFrame,
    identity: Sequence[str],
    keys: Sequence[tuple],
    names: Sequence[str],
    metrics: Sequence[str],
    by: Sequence[str],
    resamples: np.ndarray,
) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    """
    Estimates the shares of the variance of each metric of each design (or
    pair of designs, by the values of the `identity` columns of the results
    in `keys`, and of each group of `by`) of each group of parameters (see
    `attribute_variance`), with a bootstrap of the resamples of the rows of
    A, B, and the mixed replicates (`resamples x n`), whose samples it also
    returns (see `Attribution.samples`).
    """
    n = resamples.shape[1]
    feasible = results[results["feasible"].astype(bool)]
    grouped = {
        key if isinstance(key, tuple) else (key,): group
        for key, group in feasible.groupby(list(identity), sort=False)
    }
    rows = []
    samples: dict[str, list] = {"first": [], "total": [], "variance": []}
    for key in keys:
        if key not in grouped:
            continue  # infeasible in every replicate
        records = grouped[key]
        separated = records.groupby(list(by), sort=False) if by else [((), records)]
        for values, group in separated:
            labels = dict(zip(by, values if isinstance(values, tuple) else (values,)))
            redrawn = group["redrawn"]
            for metric in (metric for metric in metrics if metric in group.columns):
                estimates = _estimate(
                    _by_replicate(group, metric, redrawn.isna(), 2 * n),
                    np.array(
                        [
                            _by_replicate(group, metric, redrawn == name, n)
                            for name in names
                        ]
                    ),
                    n,
                    resamples,
                )
                rows.extend(
                    {
                        **dict(zip(identity, key)),
                        **labels,
                        "metric": metric,
                        "parameter": name,
                        "variance": estimates["variance"],
                        **{
                            share: estimates[share][i]
                            for share in ("first", "total", "first se", "total se")
                        },
                    }
                    for i, name in enumerate(names)
                )
                for share in _SHARES:
                    samples[share].extend(estimates[f"{share} samples"])
                samples["variance"].extend([estimates["variance samples"]] * len(names))
    columns = [*identity, *by, "metric", "parameter", "variance"]
    shares = pd.DataFrame(
        rows, columns=[*columns, "first", "total", "first se", "total se"]
    )
    return shares, {
        name: np.array(arrays, dtype=float).reshape(len(shares), len(resamples))
        for name, arrays in samples.items()
    }


def attribute_variance(
    study: Study,
    metrics: Sequence[str],
    designs: Sequence[Design] | pd.DataFrame,
    *,
    parameters: Sequence[str] | Mapping[str, Sequence[str]] | None = None,
    replicates: int | None = None,
    by: Sequence[str] = (),
    n_resamples: int = 500,
    seed: int = 0,
) -> Attribution:
    """
    Attributes the variance of metrics between the replicates of each design
    to its uncertain parameters (see `Uncertain`): which uncertainties spread
    its outcomes, and by how much. For example, the spread of a share of
    points sampled between replicates may be noise of the sample of points
    (which more points would reduce), and that of a cost, a risk of the
    growth of the cost.

    The effect of a parameter is measured by drawing it again while holding
    the others (a pick-freeze design, as of Sobol' indices), as its values
    need not relate to the outcomes in any order (as those of a seed do
    not). The designs are evaluated in `2 n` plain replicates (`n` of
    `replicates`): A, the first `n`, and B, the next `n` (those evaluated
    before are taken from the study's store). For each parameter (or group
    of parameters), they are also evaluated in `n` mixed replicates: those
    of A, with the parameter drawn as in B (see `Study.run`'s `shift`), which
    takes `n` evaluations of each design for each parameter.

    Two shares of the variance of a metric are estimated for each parameter,
    by Jansen's estimators (with common random numbers, the designs share
    the parameters of each replicate):

    * `total`: the share that the parameter causes, alone or with others
      (half the mean squared difference between A and the mixed replicates,
      which differ in the parameter alone, as a share of the variance). A
      parameter whose total share is near 0 could be fixed at any value
      without changing the spread of the metric.
    * `first`: the share that knowing the parameter would remove (the
      variance less half the mean squared difference between B and the
      mixed replicates, which share the parameter alone, as a share of the
      variance). The first shares of the parameters sum to at most 1; the
      rest is caused by the parameters together (their interactions).

    The shares are of the variance of a replicate, and so also of the
    variance of a mean of replicates. They are not clipped (an estimate of a
    share near 0 may be slightly negative), and their standard errors are
    estimated by a bootstrap of the rows of A, B, and the mixed replicates
    together. A share is `NaN` if the metric does not vary between the
    design's replicates, if the design does not draw the parameter (such as
    one of a branch it does not choose, which it is not evaluated for), or
    if fewer than two rows are feasible in A, B, and the mixed replicates.
    Infeasible evaluations are left out (and designs infeasible in every
    replicate have no shares).

    The evaluations must be determined by the uncertain parameters (as a
    study's store assumes): randomness from another source (such as a
    generator that is not seeded by a parameter, see `Uncertain.seed`) adds
    to the total share of every parameter.

    Args:
        study (Study): The study that evaluates designs.
        metrics (Sequence[str]): The metrics (columns of the results).
        designs (Sequence[Design] | pandas.DataFrame): The designs, or a data
            frame of their values (such as rows of a summary that keeps them,
            see `DesignSpace.designs_of`).
        parameters (Sequence[str] | Mapping[str, Sequence[str]] | None): The
            uncertain parameters to attribute variance to (by their paths,
            see `DesignSpace.uncertain_paths`), groups of them (by name), or
            `None` for each of the design space's. Parameters not given vary
            between replicates, but are not attributed variance.
        replicates (int | None): The replicates of A and of B (`n`), or
            `None` for the study's.
        by (Sequence[str]): The columns that separate the records of an
            evaluation, if it has several (such as a sensor).
        n_resamples (int): The number of resamples of the bootstrap.
        seed (int): The seed of the bootstrap.

    Returns:
        Attribution: the results, and the shares of each design, metric, and
            parameter
    """
    n = study.replicates if replicates is None else int(replicates)
    if n < 2:
        raise ValueError("replicates must be at least 2")
    if n_resamples < 2:
        raise ValueError("n_resamples must be at least 2")
    groups = _groups_of(study.space, parameters)
    metrics, by = list(metrics), tuple(by)
    unique = _designs(study.space, designs)
    plain = study.run(unique, replicates=range(2 * n))
    missing = [column for column in (*metrics, *by) if column not in plain.columns]
    if missing and plain["feasible"].any():
        raise ValueError(f"columns that are not in the results: {missing}")
    keys = ["design", *by, "replicate"]
    if not missing and plain.duplicated(keys).any():
        raise ValueError(
            "the results have several records of an evaluation (a design in a "
            + "replicate): give the columns that separate them as `by`"
        )
    frames = [plain.assign(redrawn=None)]
    for name, paths in groups.items():
        drawing = [design for design in unique if _draws(design, paths)]
        if drawing:
            mixed = study.run(
                drawing, replicates=range(n), shift={path: n for path in paths}
            )
            frames.append(mixed.assign(redrawn=name))
    results = pd.concat(frames, ignore_index=True)

    resamples = np.random.default_rng(seed).integers(n, size=(n_resamples, n))
    shares, samples = _shares(
        results,
        ("design",),
        [(design.key,) for design in unique],
        list(groups),
        metrics,
        by,
        resamples,
    )
    return Attribution(
        results=results, shares=shares, by=by, samples=samples, resamples=resamples
    )

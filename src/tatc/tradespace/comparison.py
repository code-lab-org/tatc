"""
Comparisons of designs by the replicates of a metric: the probability that
the metric of one design exceeds that of another in the same replicate, for
every pair of designs (such as the probability that one design costs more
than another).

With common random numbers (see `Design.realize`), designs share the draws
of their uncertain parameters in each replicate, so that their metrics are
correlated by what they share (such as the growth of the cost of a
component that both carry, see `CostModel.draw`): the order of two designs
may then be settled although the metric of each spreads over the replicates
by more than they differ, as the noise that they share cancels in their
difference. Comparisons of the spreads of designs (such as their quantiles)
cannot see this, and judge the order of such designs unsettled.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import numpy as np
import pandas as pd

from .criteria import _get_matrices, _groups
from .space import Design, _design_items


def _keys_of(designs: Iterable[Design | str] | pd.DataFrame) -> list[str]:
    """
    Gets the keys of designs given as designs, keys, or a data frame of them
    (with keys `design`, a column or the name of the index), each once.
    """
    keys = [d.key if isinstance(d, Design) else d for d in _design_items(designs)]
    return list(dict.fromkeys(keys))


def exceedance(
    results: pd.DataFrame,
    metric: str,
    designs: Iterable[Design | str] | pd.DataFrame | None = None,
    *,
    by: Sequence[str] | str = (),
) -> pd.DataFrame:
    """
    Estimates the probability that a metric of each design exceeds that of
    each other design (and of itself) in the same replicate, from the
    replicates in which both are feasible: the share of these replicates in
    which it is greater, counting ties as half, so that the probabilities of
    a pair and of its reverse add to 1 (and that of a design with itself is
    0.5). A probability near 0 or 1 is an order of the pair that the
    replicates settle, and one near 0.5 an order that they do not::

        results = study.run(designs)
        compared = exceedance(results, "cost")
        plot_exceedance(compared)

    With common random numbers (see `Design.realize`), the replicates of
    designs share the draws of their uncertain parameters, and the
    probability is that of the order of the designs given what they share
    (such as the growth of the cost of a component that both carry, see
    `CostModel.draw`). Without them, the draws are independent, and the
    probability is that of the order of independent outcomes of the designs.

    Args:
        results (pandas.DataFrame): The results of a study (see
            `Study.run`), with a record per replicate of each design (and
            group of `by`). Infeasible results, and those whose metric is
            missing (`NaN`), are ignored.
        metric (str): The metric.
        designs (Iterable[Design | str] | pandas.DataFrame | None): The
            designs (or their keys, or a data frame of them, such as a
            summary), in order, or `None` for every design of the results
            (of each group), by the mean of the metric (least first, and
            those without a feasible result last).
        by (Sequence[str] | str): The columns of groups of the results (such
            as `scenario`) in which designs are compared separately.

    Returns:
        pandas.DataFrame: a row per ordered pair of designs (and group), in
            the order of the designs: the `design` and the `other`, the
            columns of the group, the mean difference of the metric (the
            design's less the other's, `<metric> difference`), the
            probability that the design's metric exceeds the other's
            (`exceedance`), and the number of replicates in which both are
            feasible (`replicates`)

    Raises:
        ValueError: if the results lack the metric or a design, or have
            several records of a design in a replicate (of a group).
    """
    by = (by,) if isinstance(by, str) else tuple(by)
    if metric not in results.columns:
        raise ValueError(f"no metric {metric!r} in the results")
    keys = None if designs is None else _keys_of(designs)
    if keys is not None:
        present = set(results["design"])
        missing = [key for key in keys if key not in present]
        if missing:
            raise ValueError(f"no results of designs {missing[:5]}")
    frames = []
    for group_key, positions in _groups(results, by):
        group = results.iloc[positions]
        index = pd.Index(group["design"].unique() if keys is None else keys)
        if index.empty:
            continue
        if "feasible" in group:
            group = group[group["feasible"].astype(bool)]
        _, _, matrices = _get_matrices(group, [metric], index)
        values = matrices[metric]
        valid = ~np.isnan(values)
        if keys is None:
            # by the mean of the metric, designs without feasible results last
            with np.errstate(invalid="ignore", divide="ignore"):
                means = np.where(valid, values, 0.0).sum(axis=1) / valid.sum(axis=1)
            order = np.argsort(means, kind="stable")
            index, values, valid = index[order], values[order], valid[order]
        parts: dict[str, list[np.ndarray]] = {
            "design": [],
            "other": [],
            "difference": [],
            "exceedance": [],
            "replicates": [],
        }
        for i, design in enumerate(index):
            both = valid[i] & valid
            counts = both.sum(axis=1)
            greater = (both & (values[i] > values)).sum(axis=1)
            ties = (both & (values[i] == values)).sum(axis=1)
            with np.errstate(invalid="ignore", divide="ignore"):
                difference = np.where(both, values[i] - values, 0.0)
                parts["difference"].append(difference.sum(axis=1) / counts)
                parts["exceedance"].append((greater + 0.5 * ties) / counts)
            parts["design"].append(np.full(len(index), design, dtype=object))
            parts["other"].append(np.asarray(index, dtype=object))
            parts["replicates"].append(counts)
        frame = pd.DataFrame(
            {
                "design": np.concatenate(parts["design"]),
                "other": np.concatenate(parts["other"]),
                **dict(zip(by, group_key)),
                f"{metric} difference": np.concatenate(parts["difference"]),
                "exceedance": np.concatenate(parts["exceedance"]),
                "replicates": np.concatenate(parts["replicates"]),
            }
        )
        frames.append(frame)
    if not frames:
        return pd.DataFrame(
            columns=[
                "design",
                "other",
                *by,
                f"{metric} difference",
                "exceedance",
                "replicates",
            ]
        )
    return pd.concat(frames, ignore_index=True)

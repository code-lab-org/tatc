"""
Utility: the value of the metrics of designs to their stakeholders, as a
single-attribute utility of each metric (from 0 at its worst acceptable
value to 1 at its best) combined into a multi-attribute utility (by the
weights of the metrics, additively or in Keeney and Raiffa's multiplicative
form), for trade spaces of utility against cost (as in multi-attribute
tradespace exploration).

A utility is a metric computed from others in each replicate (or each row
of a summary), so that statistics of it judge designs as statistics of any
metric do (see `Statistic`): its mean over the replicates is the expected
utility, whose curvature expresses an attitude to the risks that the
replicates draw.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd


def _normalized(values: Any, worst: float, best: float) -> np.ndarray:
    """Gets values as shares of the way from the worst to the best, clipped
    to [0, 1] (`NaN` kept)."""
    share = (np.asarray(values, dtype=float) - worst) / (best - worst)
    return np.clip(share, 0.0, 1.0)


@dataclass(frozen=True)
class LinearUtility:
    """
    A single-attribute utility linear in the value of a metric: 0 at the worst
    value and 1 at the best (and beyond them), for a metric to maximize
    (`worst < best`) or to minimize (`worst > best`).

    Args:
        worst (float): The worst acceptable value (utility 0).
        best (float): The best value (utility 1).
    """

    worst: float
    best: float

    def __post_init__(self):
        if self.worst == self.best:
            raise ValueError("the worst and best values must differ")

    def __call__(self, values: Any) -> np.ndarray:
        return _normalized(values, self.worst, self.best)


@dataclass(frozen=True)
class ExponentialUtility:
    """
    A single-attribute utility exponential in the value of a metric (of
    constant absolute risk aversion): 0 at the worst value and 1 at the best,
    and `(1 - exp(-risk z)) / (1 - exp(-risk))` between, where `z` is the share
    of the way from the worst value to the best. A positive `risk` is averse
    (concave: the first gains matter most, and a sure value is preferred to
    a gamble of the same mean), a negative one is seeking (convex), and 0 is
    linear (see `LinearUtility`).

    Args:
        worst (float): The worst acceptable value (utility 0).
        best (float): The best value (utility 1).
        risk (float): The coefficient of risk aversion (of the shares `z`).
    """

    worst: float
    best: float
    risk: float = 0.0

    def __post_init__(self):
        if self.worst == self.best:
            raise ValueError("the worst and best values must differ")

    def __call__(self, values: Any) -> np.ndarray:
        share = _normalized(values, self.worst, self.best)
        if self.risk == 0:
            return share
        return np.expm1(-self.risk * share) / math.expm1(-self.risk)


@dataclass(frozen=True)
class PiecewiseUtility:
    """
    A single-attribute utility interpolated linearly between points (value,
    utility), and constant beyond the first and last, such as a utility
    elicited at a few values.

    Args:
        points (Sequence[tuple[float, float]]): The points: values (distinct)
            and their utilities (between 0 and 1).
    """

    points: Sequence[tuple[float, float]]

    def __post_init__(self):
        points = tuple(sorted((float(x), float(u)) for x, u in self.points))
        if len(points) < 2:
            raise ValueError("a piecewise utility needs at least two points")
        values = [x for x, _ in points]
        if len(set(values)) != len(values):
            raise ValueError("the values of the points must be distinct")
        if not all(0 <= u <= 1 for _, u in points):
            raise ValueError("utilities must be between 0 and 1")
        object.__setattr__(self, "points", points)

    def __call__(self, values: Any) -> np.ndarray:
        array = np.asarray(values, dtype=float)
        xs, us = zip(*self.points)
        result = np.interp(array, xs, us)
        return np.where(np.isnan(array), np.nan, result)


_FORMS = ("additive", "multiplicative")


def _scaling_constant(weights: Sequence[float]) -> float:
    """
    Gets the scaling constant `K` of a multiplicative utility, the root other
    than 0 of `1 + K = prod(1 + K k)` (0 if the weights sum to 1): between -1
    and 0 if they sum to more than 1, and positive if less.
    """
    total = sum(weights)
    if len(weights) < 2 or math.isclose(total, 1.0):
        return 0.0

    def excess(k: float) -> float:
        return math.prod(1 + k * w for w in weights) - (1 + k)

    if total > 1:
        low, high = -1.0 + 1e-12, -1e-12
    else:
        low, high = 1e-12, 1.0
        while excess(high) < 0:
            high *= 2
    for _ in range(200):
        middle = (low + high) / 2
        if (excess(middle) > 0) == (excess(low) > 0):
            low = middle
        else:
            high = middle
    return (low + high) / 2


@dataclass(frozen=True, eq=False)
class Utility:
    """
    A multi-attribute utility of the metrics of a design: the single-attribute
    utility of each metric (see `LinearUtility`, `ExponentialUtility`, and
    `PiecewiseUtility`, or any function of the values that gives utilities
    from 0 to 1), combined by weights:

    * `"additive"`: the weighted sum of the utilities (by weights normalized
      to sum to 1), for metrics that are preferentially independent and whose
      gains add.
    * `"multiplicative"`: Keeney and Raiffa's multiplicative form, `1 + K U =
      prod(1 + K k u)`, with the weights `k` (each between 0 and 1) the utility
      of each metric at its best with the others at their worst, and `K` from
      them: weights that sum to less than 1 make the metrics complements (a
      design must do well in all), and more than 1 substitutes (doing well in
      one suffices).

    A utility is a metric of statistics (such as `Mean(utility)`, its expected
    value over the replicates, see `Statistic`), computed in each replicate
    from its metrics. The expected utility suits replicates that draw risks
    of a mission (such as cost growth, or failures): a concave utility then
    penalizes the spread of outcomes, as a risk-averse stakeholder would.
    Where the spread between replicates is noise of an estimate (such as of
    a sample of points, see `attribute_variance`), the expected utility
    penalizes the noise rather than the design: judge the utility of the
    means instead (`utility(summary)`, of a summary of the means of the
    metrics, see `aggregate`).

    Args:
        attributes (Mapping[str, Callable]): The single-attribute utility of
            each metric.
        weights (Mapping[str, float]): The weight of each metric.
        form (str): `"additive"` or `"multiplicative"`.
    """

    attributes: Mapping[str, Callable[[Any], Any]]
    weights: Mapping[str, float]
    form: str = "additive"

    def __post_init__(self):
        attributes, weights = dict(self.attributes), dict(self.weights)
        if not attributes:
            raise ValueError("a utility needs attributes")
        if set(weights) != set(attributes):
            raise ValueError(
                "a utility needs a weight for each attribute, and no others"
            )
        if self.form not in _FORMS:
            raise ValueError(f"form must be one of {list(_FORMS)}, not {self.form!r}")
        if any(not weight > 0 for weight in weights.values()):
            raise ValueError("weights must be positive")
        if self.form == "multiplicative" and any(w > 1 for w in weights.values()):
            raise ValueError(
                "the weights of a multiplicative utility must be at most 1"
            )
        object.__setattr__(self, "attributes", attributes)
        object.__setattr__(
            self, "weights", {name: weights[name] for name in attributes}
        )

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, Utility)
            and self.attributes == other.attributes
            and self.weights == other.weights
            and self.form == other.form
        )

    def __hash__(self) -> int:
        return hash(
            (tuple(self.attributes.items()), tuple(self.weights.items()), self.form)
        )

    @property
    def columns(self) -> tuple[str, ...]:
        """tuple[str, ...]: The metrics of the utility."""
        return tuple(self.attributes)

    @property
    def scaling(self) -> float:
        """float: The scaling constant `K` of a multiplicative utility (0 for
        an additive one, or weights that sum to 1)."""
        if self.form == "additive":
            return 0.0
        return _scaling_constant(list(self.weights.values()))

    def evaluate(self, values: Mapping[str, Any]) -> np.ndarray:
        """
        Computes the utility of values of the metrics (arrays of the same
        shape, such as of designs by replicates).

        Args:
            values (Mapping[str, Any]): The values of each metric.

        Returns:
            numpy.ndarray: the utility of each value
        """
        utilities = {
            name: np.asarray(function(values[name]), dtype=float)
            for name, function in self.attributes.items()
        }
        weights = self.weights
        scaling = self.scaling
        if scaling == 0.0:
            total = sum(weights.values())
            return sum(weights[name] / total * u for name, u in utilities.items())
        product = np.ones_like(next(iter(utilities.values())))
        for name, u in utilities.items():
            product = product * (1 + scaling * weights[name] * u)
        return (product - 1) / scaling

    def __call__(self, frame: pd.DataFrame | Mapping[str, Any]) -> Any:
        """
        Computes the utility of each row of a data frame (such as the records
        of a study, or a summary of the means of designs), or of values of
        the metrics.

        Args:
            frame (pandas.DataFrame | Mapping[str, Any]): The metrics.

        Returns:
            pandas.Series | numpy.ndarray: the utility of each row (with the
                frame's index), or of the values
        """
        if isinstance(frame, pd.DataFrame):
            missing = [name for name in self.columns if name not in frame.columns]
            if missing:
                raise ValueError(
                    f"metrics of the utility that are not columns: {missing}"
                )
            values = {name: frame[name].to_numpy(dtype=float) for name in self.columns}
            return pd.Series(self.evaluate(values), index=frame.index, name="utility")
        return self.evaluate(frame)

    def __repr__(self) -> str:
        return f"Utility({list(self.attributes)}, {self.form})"

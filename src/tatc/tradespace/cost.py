"""
Life-cycle cost of constellations for trade studies: the development,
production, launch, and operations of the satellites of a deployment (see
`DeploymentSchedule`), with learning in production, replacements over the
mission, and spares, dated for their present value and annual costs.

TAT-C does not include cost-estimating relationships, whose common forms are
restricted or need calibration to an organization's experience: the costs of
each type of satellite (such as the satellites of each payload) and of
launches and operations are given, as values or functions of a design (such
as a mass-based relationship, as a `MissionPlan` of each design sizes them),
in a currency of the user's choice and constant-year amounts. Cost risk
enters as uncertain growth factors (see `CostModel.estimate`), drawn for
each replicate like other uncertain parameters, so that statistics of the
cost over replicates (such as a 0.7 `Quantile`) are objectives or
requirements: of every cost, of each element, or of each element of each
type of satellite.

The costs of a type of satellite may also be the sum of those of its
components (such as its bus and instruments, see `Component`), each with an
uncertain growth by its maturity (such as a `Triangular` distribution, wider
for new components than for existing ones), correlated between the
components (`CostModel.correlation`), as in formal risk assessments of cost
estimates. Their growth is drawn from scores of the components and a common
score (`CostModel.draw`), uncertain parameters of the design space
(`risk_space`) that common random numbers give alike to every design of a
replicate: designs that share a component share its growth, so that the
costs of designs are correlated by what they share (as in M. W. Maier, K.
L. Yeakel, and E. B. Wendoloski, "Cost Variance Analysis in Constellation
Architecture Studies," IEEE Systems Journal 17(2), 2023), and comparisons
between them (such as the probability that one costs more than another, see
`exceedance`) see the differences of their costs.

A deployment dates the costs: development of each type of satellite (and
of each later block of its production, see `SatelliteCost.block_size`),
production of its satellites, launches, replacements, and operations from
launch to the end of the mission, so that each design's costs follow from
its satellites and planes. The dated costs give the present value of the
cost at a discount rate, and the cost of each year (such as its peak, for
annual budgets)::

    def evaluate(built):
        d, satellites = built
        deployment = SCHEDULE.plan(d.layout, kind=d.payload)
        cost = MODEL.estimate(deployment, growth={"development": d.growth})
        return cost.metrics()

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import dataclasses
import math
import numbers
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from .deployment import Deployment
from .mission import MissionHistory
from .nodes import Uncertain
from .space import _RESERVED_NAMES, DesignSpace

ELEMENTS = ("development", "production", "launch", "operations")
"""The elements of the cost of a constellation."""

_PATH_CHARACTERS = (".", "[", "]")
"""Characters of the paths of values of designs, which names of scores lack."""

COMMON = "common"
"""The name of the common score of the growth of components (see
`CostModel.draw`), which no component may have."""


@dataclass(frozen=True)
class Triangular:
    """
    A triangular distribution of the growth of a cost, from `low` to `high`
    and most likely `mode` (by default 1, the estimate), such as of a
    component of a given maturity: wider and skewed further up for new
    components than for existing ones. For example, formal risk assessments
    of satellites have used growths from 0.827 to 2.274 for new components,
    0.904 to 1.708 for modified ones, and 0.981 to 1.140 for existing ones
    (M. W. Maier, K. L. Yeakel, and E. B. Wendoloski, "Cost Variance Analysis
    in Constellation Architecture Studies," IEEE Systems Journal 17(2),
    2023, Table I).

    Args:
        low (float): The lowest growth.
        high (float): The highest growth.
        mode (float): The most likely growth.
    """

    low: float
    high: float
    mode: float = 1.0

    def __post_init__(self):
        if not 0 <= self.low <= self.mode <= self.high or self.low == self.high:
            raise ValueError("growth must satisfy 0 <= low <= mode <= high, low < high")

    def mean(self) -> float:
        """
        Gets the mean growth.

        Returns:
            float: the mean
        """
        return (self.low + self.mode + self.high) / 3

    def std(self) -> float:
        """
        Gets the standard deviation of the growth.

        Returns:
            float: the standard deviation
        """
        a, b, c = self.low, self.mode, self.high
        return math.sqrt((a * a + b * b + c * c - a * b - a * c - b * c) / 18)

    def ppf(self, q: float) -> float:
        """
        Gets the quantile of the growth at a probability (the inverse of its
        cumulative distribution).

        Args:
            q (float): The probability, in [0, 1].

        Returns:
            float: the growth
        """
        low, mode, high = self.low, self.mode, self.high
        split = (mode - low) / (high - low)
        if q <= split:
            return low + math.sqrt(q * (high - low) * (mode - low))
        return high - math.sqrt((1 - q) * (high - low) * (high - mode))


def _is_number(value: Any) -> bool:
    return isinstance(value, numbers.Real) and not isinstance(value, bool)


@dataclass(frozen=True)
class Component:
    """
    A component of a type of satellite (such as its bus, an instrument, or
    its software), whose costs of development and of its first unit grow by
    the same factor: a number, or a distribution (an object with methods
    `ppf` and `mean`, such as `Triangular`). A component's growth is drawn
    from the distribution by its score (see `CostModel.draw`), and is
    otherwise the mean of the distribution, so that estimates without drawn
    growth are expected costs.

    Components are identified by name: those of the same name draw the same
    score, in every type of satellite and design (such as an instrument
    that several types carry), and the types of satellites of a cost model
    that carry a component develop it once (see `CostModel.estimate`).

    Args:
        name (str): The name.
        development (float): The nonrecurring cost.
        first_unit (float): The recurring cost of the first unit.
        growth (Any): The growth factor of its costs, or its distribution.
    """

    name: str
    development: float = 0.0
    first_unit: float = 0.0
    growth: Any = 1.0

    def __post_init__(self):
        if not isinstance(self.name, str) or not self.name:
            raise ValueError("a component needs a name")
        if self.development < 0 or self.first_unit < 0:
            raise ValueError("costs must not be negative")
        if _is_number(self.growth):
            if self.growth < 0:
                raise ValueError("growth must not be negative")
        elif not (hasattr(self.growth, "ppf") and hasattr(self.growth, "mean")):
            raise TypeError(
                "growth must be a number or a distribution with methods ppf and mean"
            )

    @property
    def uncertain(self) -> bool:
        """bool: Whether its growth is drawn from a distribution."""
        return not _is_number(self.growth)

    @property
    def factor(self) -> float:
        """float: Its growth factor (or the mean of its distribution)."""
        if self.uncertain:
            return float(self.growth.mean())
        return float(self.growth)


@dataclass(frozen=True)
class SatelliteCost:
    """
    The cost of a type of satellite (bus and payload): a nonrecurring
    `development` cost, and the recurring cost of its units, from that of
    the first (`first_unit`, the theoretical first unit) by a learning curve:
    unit `k` costs `first_unit * k ** log2(learning_rate)`, so that each
    doubling of the units built costs `learning_rate` times as much per unit
    (typically 0.85 to 0.95 for satellites, 1 for no learning).

    The costs may also be, or include, those of its components (such as its
    bus and instruments, see `Component`), each with its growth: the costs
    of the type are its own and those of its components, grown by their
    factors (`total_development` and `total_first_unit`, of the type
    alone). The growths of components may be correlated (when drawn, see
    `CostModel.correlation`).

    Satellites may be produced in blocks of `block_size` units, each block
    after the first a new design of the type (as its components become
    obsolete over a long mission of replacements): its development is
    repeated (a share of it, `redevelopment`, 1 for a full redesign) before
    its first unit is produced, and the learning curve restarts with its
    first unit, so that a type replaced over decades costs the development
    and first units of each block (as in M. W. Maier, K. L. Yeakel, and E.
    B. Wendoloski, "Cost Variance Analysis in Constellation Architecture
    Studies," IEEE Systems Journal 17(2), 2023). Smaller blocks cost more
    developments and first units for the same satellites.

    Args:
        development (float): The nonrecurring cost (besides that of its
            components).
        first_unit (float): The recurring cost of the first unit (besides
            that of its components).
        learning_rate (float): The learning rate, in (0, 1].
        mass (float | None): The mass (kg), for launches priced by mass.
        components (Sequence[Component]): The components, if any, of
            distinct names.
        block_size (int | None): The units of each block of production, or
            `None` for one block of every unit.
        redevelopment (float): The share of the development (with that of
            its components, `total_development`) repeated for each block
            after the first.
    """

    development: float = 0.0
    first_unit: float = 0.0
    learning_rate: float = 1.0
    mass: float | None = None
    components: Sequence[Component] = ()
    block_size: int | None = None
    redevelopment: float = 1.0

    def __post_init__(self):
        if not 0 < self.learning_rate <= 1:
            raise ValueError("learning_rate must be in (0, 1]")
        if self.development < 0 or self.first_unit < 0:
            raise ValueError("costs must not be negative")
        if self.block_size is not None and (
            isinstance(self.block_size, bool)
            or not isinstance(self.block_size, numbers.Integral)
            or self.block_size < 1
        ):
            raise ValueError("block_size must be a positive integer (or None)")
        if self.redevelopment < 0:
            raise ValueError("redevelopment must not be negative")
        components = tuple(self.components)
        if not all(isinstance(c, Component) for c in components):
            raise TypeError("components must be Components")
        names = [c.name for c in components]
        if len(set(names)) < len(names):
            raise ValueError(
                "the components of a type of satellite need distinct names"
            )
        object.__setattr__(self, "components", components)

    @property
    def total_development(self) -> float:
        """float: The nonrecurring cost, with those of its components grown
        by their factors (of the type alone: types that share a component
        develop it once, see `CostModel.estimate`)."""
        return self.development + sum(c.development * c.factor for c in self.components)

    @property
    def total_first_unit(self) -> float:
        """float: The recurring cost of the first unit, with those of its
        components grown by their factors."""
        return self.first_unit + sum(c.first_unit * c.factor for c in self.components)

    def production(self, units: int) -> float:
        """
        Gets the recurring cost of building units (the sum of the costs of
        units 1 to `units` on the learning curve, from `total_first_unit`,
        which restarts with each block).

        Args:
            units (int): The number of units.

        Returns:
            float: the cost
        """
        if units <= 0:
            return 0.0
        exponent = math.log2(self.learning_rate)
        size = units if self.block_size is None else int(self.block_size)
        full, rest = divmod(units, size)

        def curve(count: int) -> float:
            return float(np.sum(np.arange(1, count + 1) ** exponent))

        return float(self.total_first_unit * (full * curve(size) + curve(rest)))

    def blocks(self, units: int) -> int:
        """
        Gets the number of blocks of production that building units starts.

        Args:
            units (int): The number of units.

        Returns:
            int: the blocks (1 for any units without a block size)
        """
        if units <= 0:
            return 0
        if self.block_size is None:
            return 1
        return -(-units // int(self.block_size))


@dataclass(frozen=True)
class LaunchCost:
    """
    The cost of launches: of each launch of a deployment (dedicated
    launches, see `DeploymentSchedule` for their capacity), and prices per
    satellite and per kg (as of rideshares). Give the fields of the launch
    service of a study and leave the others zero.

    Args:
        per_launch (float): The cost of a launch.
        per_satellite (float): The cost of launching each satellite.
        per_kg (float): The cost of launching each kg of satellite.
    """

    per_launch: float = 0.0
    per_satellite: float = 0.0
    per_kg: float = 0.0


@dataclass(frozen=True)
class OperationsCost:
    """
    The cost of operating a constellation: a fixed cost per year (such as
    the mission operations center) and a cost per satellite and year.

    Args:
        per_year (float): The fixed cost of each year.
        per_satellite_year (float): The cost of each satellite in each year.
    """

    per_year: float = 0.0
    per_satellite_year: float = 0.0


def _spread(amount: float, begin: float, end: float) -> list[tuple[float, float]]:
    """
    Spreads an amount evenly from `begin` to `end` (years), dated at the
    middle of each of its years (or at `begin`, if it has no duration).
    """
    if end <= begin:
        return [(begin, amount)]
    flows = []
    time = begin
    while time < end - 1e-12:
        step = min(1.0, end - time)
        flows.append((time + step / 2, amount * step / (end - begin)))
        time += step
    return flows


@dataclass
class CostEstimate:
    """
    An estimate of the life-cycle cost of a constellation (see
    `CostModel.estimate`).
    """

    development: float = 0.0
    """float: The nonrecurring cost of the types of satellites."""
    production: float = 0.0
    """float: The recurring cost of the units built."""
    launch: float = 0.0
    """float: The cost of launches."""
    operations: float = 0.0
    """float: The cost of operations."""
    units: dict[str, int] = field(default_factory=dict)
    """dict: The units built of each type of satellite (including
    replacements and spares)."""
    launches: int = 0
    """int: The launches priced per launch."""
    flows: list[tuple[float, str, float]] = field(default_factory=list)
    """list: The dated costs: the time (years from the start of
    development), element, and amount of each."""
    present: dict[str, float] = field(default_factory=dict)
    """dict: The present value of each element."""
    blocks: dict[str, int] = field(default_factory=dict)
    """dict: The blocks of production of each type of satellite (see
    `SatelliteCost.block_size`)."""

    @property
    def total(self) -> float:
        """float: The total cost."""
        return self.development + self.production + self.launch + self.operations

    @property
    def present_value(self) -> float:
        """float: The present value of the total cost."""
        return float(sum(self.present.values()))

    def annual(self) -> pd.DataFrame:
        """
        Gets the cost of each year (from the start of development), by
        element and in total.

        Returns:
            pandas.DataFrame: a row per year, a column per element, and the
                `total`
        """
        frame = pd.DataFrame(self.flows, columns=["time", "element", "amount"])
        if len(frame) == 0:
            return pd.DataFrame(columns=[*ELEMENTS, "total"])
        frame["year"] = np.floor(frame["time"]).astype(int)
        table = frame.pivot_table(
            index="year", columns="element", values="amount", aggfunc="sum"
        )
        table = table.reindex(columns=list(ELEMENTS)).fillna(0.0)
        table = table.reindex(range(table.index.min(), table.index.max() + 1))
        table = table.fillna(0.0)
        table["total"] = table.sum(axis=1)
        return table

    def metrics(self, prefix: str = "cost") -> dict[str, float]:
        """
        Gets the cost as metrics of an evaluation (see `Study`): the total
        (`<prefix>`), each element (`<prefix> development`, `<prefix>
        production`, `<prefix> launch`, and `<prefix> operations`), the
        present value of the total (`<prefix> present`), and the cost of its
        most expensive year (`<prefix> peak annual`).

        Args:
            prefix (str): The name of the total, and prefix of the others.

        Returns:
            dict[str, float]: the metrics
        """
        annual = self.annual()
        return {
            prefix: self.total,
            **{f"{prefix} {element}": getattr(self, element) for element in ELEMENTS},
            f"{prefix} present": self.present_value,
            f"{prefix} peak annual": (
                float(annual["total"].max()) if len(annual) > 0 else 0.0
            ),
        }


Growth = float | Mapping[str, float | Mapping[Any, float]]
"""Growth factors of costs: one for every element, or by element (see
`ELEMENTS`, others 1), and of development, production, and launches, by
type of satellite (a mapping of each kind to its factor)."""


@dataclass(frozen=True)
class CostModel:
    """
    A model of the life-cycle cost of constellations: the costs of their
    types of satellites, launches, and operations, applied to their
    deployments (see `DeploymentSchedule`), and the discount rate of their
    present value.

    A cost at time `t` (years from the start of development, which begins
    `start` years after the present) has a present value of `cost / (1 +
    discount_rate) ** (start + t)`.

    The growths of the components of satellites (see `Component`) may be
    correlated (when drawn, see `draw`): costs estimated together tend to be
    underestimated together, as are those of the software and management of
    a satellite whose payload costs more than estimated. A correlation of
    0.5 between the components of a satellite is commonly assumed; here, it
    is of every pair of components of the model (such as of the instruments
    of different types of satellites, estimated alike).

    Args:
        satellites (SatelliteCost | Mapping[str, SatelliteCost]): The cost of
            the satellites, or of each type of satellite (by name). The
            components of the same name of different types must be the same
            (equal costs and growth).
        launch (LaunchCost): The cost of launches.
        operations (OperationsCost): The cost of operations.
        discount_rate (float): The annual discount rate (such as 0.03).
        start (float): The years from the present to the start of
            development.
        correlation (float): The correlation of the scores of the components
            (see `draw`), in [0, 1].
    """

    satellites: SatelliteCost | Mapping[str, SatelliteCost]
    launch: LaunchCost = LaunchCost()
    operations: OperationsCost = OperationsCost()
    discount_rate: float = 0.0
    start: float = 0.0
    correlation: float = 0.0

    def __post_init__(self):
        if self.discount_rate <= -1:
            raise ValueError("discount_rate must be greater than -1")
        if not 0 <= self.correlation <= 1:
            raise ValueError("correlation must be in [0, 1]")
        components: dict[str, Component] = {}
        for _, _, cost in self._types():
            for component in cost.components:
                if component.name == COMMON:
                    raise ValueError(
                        f"{COMMON!r} names the common score of cost risk, not a component"
                    )
                if components.setdefault(component.name, component) != component:
                    raise ValueError(
                        f"the components named {component.name!r} differ: types of "
                        + "satellites that share a component share its costs and growth"
                    )

    def discount(self, time: float) -> float:
        """
        Gets the factor of the present value of a cost at a time.

        Args:
            time (float): The time (years from the start of development).

        Returns:
            float: the factor
        """
        return float((1 + self.discount_rate) ** -(self.start + time))

    def _get_type(self, kind: Any) -> SatelliteCost:
        if isinstance(self.satellites, SatelliteCost):
            return self.satellites
        if kind is None:
            raise ValueError("a model of several types of satellites needs a kind")
        if kind not in self.satellites:
            raise ValueError(f"no cost of satellites of kind {kind!r}")
        return self.satellites[kind]

    def _types(self) -> list[tuple[str, Any, SatelliteCost]]:
        """
        Gets the types of satellites: the name of each (as of the units of
        estimates), its kind, and its cost.
        """
        if isinstance(self.satellites, SatelliteCost):
            return [(_type_name(None), None, self.satellites)]
        return [
            (_type_name(kind), kind, cost) for kind, cost in self.satellites.items()
        ]

    def score_names(self) -> list[str]:
        """
        Gets the names of the scores from which the growth of the components
        is drawn (see `draw`): those of the components with distributions of
        growth, after the common score (`"common"`) if they are correlated.

        Returns:
            list[str]: the names, in the order of the types and components
        """
        names = list(
            dict.fromkeys(
                component.name
                for _, _, cost in self._types()
                for component in cost.components
                if component.uncertain
            )
        )
        if names and self.correlation > 0:
            names.insert(0, COMMON)
        return names

    def draw(self, scores: Mapping[str, float]) -> CostModel:
        """
        Draws the growth of the costs of the components of its types of
        satellites (see `Component`) from their scores: standard normal
        values of the components and a common score (by name, see
        `score_names`), such as uncertain parameters of a design space drawn
        for each replicate (see `risk_space`). The growth of a component is
        the quantile of its distribution at the normal probability of its
        score, mixed with the common score by the correlation (ρ, see
        `correlation`): `sqrt(ρ) * score["common"] + sqrt(1 - ρ) *
        score[component]`. The growths of the components are then
        correlated (by about ρ), and a component that several types of
        satellites carry has the same growth in each.

        With common random numbers (see `Design.realize`), a component (and
        the common score) has the same score in the same replicate of every
        design, so that designs that share it share its growth: the costs of
        designs are correlated by what they share, and comparisons between
        them (such as `exceedance`) see the differences of their costs
        rather than of independent draws.

        Args:
            scores (Mapping[str, float]): The scores, by name (others are
                ignored).

        Returns:
            CostModel: the cost model, whose components have the drawn growths

        Raises:
            ValueError: if a score is missing.
        """
        names = self.score_names()
        missing = [name for name in names if name not in scores]
        if missing:
            raise ValueError(f"no scores of {missing}")
        common = float(scores[COMMON]) if COMMON in names else 0.0
        probabilities = {}
        for name in names:
            if name != COMMON:
                score = math.sqrt(self.correlation) * common + math.sqrt(
                    1 - self.correlation
                ) * float(scores[name])
                probabilities[name] = 0.5 * (1 + math.erf(score / math.sqrt(2)))
        drawn = {}
        for _, kind, cost in self._types():
            components = tuple(
                (
                    dataclasses.replace(
                        component,
                        growth=float(
                            component.growth.ppf(probabilities[component.name])
                        ),
                    )
                    if component.uncertain
                    else component
                )
                for component in cost.components
            )
            drawn[kind] = dataclasses.replace(cost, components=components)
        if isinstance(self.satellites, SatelliteCost):
            return dataclasses.replace(self, satellites=drawn[None])
        return dataclasses.replace(self, satellites=drawn)

    def estimate(
        self, deployment: Deployment | MissionHistory, growth: Growth = 1.0
    ) -> CostEstimate:
        """
        Estimates the life-cycle cost of a deployment (as planned, without
        failures) or of a simulated mission (see `Deployment.simulate`), with
        the satellites and launches its failures and restorations take,
        dating its costs: the development of each type of satellite, spread
        over its development; the production of its units (including spares,
        replacements, and restorations) along its learning curve, spread over
        their production; each launch (including those that fail), when it
        flies; a fixed cost of operations over the mission; and the
        operations of each satellite while it is in orbit (operating, or a
        spare waiting). The costs of types of satellites include those of
        their components, grown by their factors (drawn, see `draw`, or
        otherwise the means of their distributions, for expected costs). A
        component that several types of the deployment carry is developed
        once, with the type whose development ends first (and by its growth
        factor), and produced with each type (on its learning curve). Each
        block of production after the first (see `SatelliteCost.block_size`)
        develops its type again (a share of its development), dated over the
        type's development just before the production of its first unit.

        Growth factors multiply the estimate, such as uncertain parameters of
        cost growth drawn for each replicate (for example,
        `Uncertain(lambda rng: rng.triangular(0.9, 1.0, 1.5))`), whose
        statistics over the replicates give the cost at a confidence level
        (such as the 0.7 `Quantile`). Factors of each type of satellite (such
        as `{"development": {"narrow": "growth_narrow", "wide":
        "growth_wide"}}` in a `MissionPlan`) are independent between types
        but shared by the designs of each type of a replicate (with common
        random numbers, see `Design.realize`), so that the costs of designs
        are correlated by the types of satellites they share. Do not use
        uncertain growth in constraints of a design space, which see their
        nominal values (`NaN`).

        Args:
            deployment (Deployment | MissionHistory): The deployment (see
                `DeploymentSchedule.plan`), or a simulated mission.
            growth (float | Mapping): A growth factor of the cost, or of each
                element (`"development"`, `"production"`, `"launch"`, or
                `"operations"`; others 1), and of development, production, and
                launches, of each kind of satellite (a mapping of each kind of
                the deployment to its factor).

        Returns:
            CostEstimate: the estimate

        Raises:
            ValueError: if growth names an unknown element, is of operations
                by type of satellite, or lacks a kind of the deployment.
        """
        factors = _get_growth(growth)
        history = (
            deployment.simulate() if isinstance(deployment, Deployment) else deployment
        )
        estimate = CostEstimate()
        flows: list[tuple[float, str, float]] = []

        def add(
            element: str, amounts: list[tuple[float, float]], kind: Any = None
        ) -> None:
            factor = _get_factor(factors, element, kind)
            flows.extend((time, element, amount * factor) for time, amount in amounts)

        if self.launch.per_kg > 0 and any(
            self._get_type(kind).mass is None for kind in history.development
        ):
            raise ValueError("launches priced by mass require masses")
        developed_components: set[str] = set()
        for kind, developed in sorted(
            history.development.items(), key=lambda item: item[1]
        ):
            cost = self._get_type(kind)
            amount = cost.development + sum(
                c.development * c.factor
                for c in cost.components
                if c.name not in developed_components
            )
            developed_components.update(c.name for c in cost.components)
            add("development", _spread(amount, 0.0, developed), kind)
        built: dict[Any, int] = {}
        for batch in history.production:
            cost = self._get_type(batch.kind)
            before = built.get(batch.kind, 0)
            amount = cost.production(before + batch.satellites) - cost.production(
                before
            )
            built[batch.kind] = before + batch.satellites
            add("production", _spread(amount, batch.begin, batch.end), batch.kind)
            if cost.block_size is not None and cost.redevelopment > 0:
                # the development of each later block, before its first unit
                duration = history.development.get(batch.kind, 0.0)
                size = int(cost.block_size)
                for unit in range(before, before + batch.satellites):
                    if unit > 0 and unit % size == 0:
                        begins = (
                            batch.begin
                            + (unit - before)
                            * (batch.end - batch.begin)
                            / batch.satellites
                        )
                        add(
                            "development",
                            _spread(
                                cost.redevelopment * cost.total_development,
                                begins - duration,
                                begins,
                            ),
                            batch.kind,
                        )
        estimate.units = {_type_name(kind): units for kind, units in built.items()}
        estimate.blocks = {
            _type_name(kind): self._get_type(kind).blocks(units)
            for kind, units in built.items()
        }
        begin, end = history.begin, history.end
        for launch in history.launches:
            cost = self._get_type(launch.kind)
            estimate.launches += int(self.launch.per_launch > 0)
            price = self.launch.per_launch + launch.satellites * (
                self.launch.per_satellite + self.launch.per_kg * (cost.mass or 0.0)
            )
            add("launch", [(launch.time, price)], launch.kind)
        for start, stop, count in history.in_orbit():
            add(
                "operations",
                _spread(
                    self.operations.per_satellite_year * count * (stop - start),
                    start,
                    stop,
                ),
            )
        add("operations", _spread(self.operations.per_year * (end - begin), begin, end))
        flows = [flow for flow in flows if flow[2] != 0]
        estimate.flows = sorted(flows, key=lambda flow: (flow[0], flow[1]))
        for element in ELEMENTS:
            setattr(
                estimate,
                element,
                float(sum(amount for _, kind, amount in flows if kind == element)),
            )
        estimate.present = {
            element: float(
                sum(
                    amount * self.discount(time)
                    for time, kind, amount in flows
                    if kind == element
                )
            )
            for element in ELEMENTS
        }
        return estimate


def _type_name(kind: Any) -> str:
    """
    Gets the name of a type of satellite (as of the units of estimates): its
    kind, or `"satellites"` without one.
    """
    return "satellites" if kind is None else str(kind)


def _get_growth(growth: Growth) -> dict[str, float | dict[Any, float]]:
    """
    Gets the growth factor of each element of cost (or a mapping of each kind
    of satellite to its factor).
    """
    if not isinstance(growth, Mapping):
        return {element: float(growth) for element in ELEMENTS}
    unknown = set(growth) - set(ELEMENTS)
    if unknown:
        raise ValueError(f"growth of unknown elements: {sorted(unknown)}")
    factors: dict[str, float | dict[Any, float]] = {}
    for element in ELEMENTS:
        value = growth.get(element, 1.0)
        if isinstance(value, Mapping):
            if element == "operations":
                raise ValueError("the growth of operations is not by type of satellite")
            factors[element] = {kind: float(factor) for kind, factor in value.items()}
        else:
            factors[element] = float(value)
    return factors


def _get_factor(
    factors: Mapping[str, float | Mapping[Any, float]], element: str, kind: Any
) -> float:
    """
    Gets the growth factor of an element of the cost of a kind of satellite
    (`None` for costs of no type).
    """
    factor = factors[element]
    if not isinstance(factor, Mapping):
        return factor
    if kind is None:
        raise ValueError(f"growth of {element} by type of satellite needs kinds")
    if kind not in factor:
        raise ValueError(f"no growth of {element} of satellites of kind {kind!r}")
    return factor[kind]


def risk_space(*models: CostModel) -> DesignSpace:
    """
    Gets a group of a design space of the scores of the cost risk of cost
    models (from which the growth of their components is drawn, see
    `CostModel.draw`): an uncertain standard normal value of each name of
    their scores (see `CostModel.score_names`), drawn for each replicate::

        space = space.extend(cost_risk=risk_space(COST_MODEL))
        PLAN = MissionPlan(SCHEDULE, COST_MODEL, kind="payload", scores="cost_risk")

    With common random numbers (see `Design.realize`), each name has the same
    score in the same replicate of every design, so that the designs that
    share a component share its growth. The variance of costs (or of the
    differences of the costs of designs) may be attributed to each component
    by the paths of its score (such as `"cost_risk.imager"`, see
    `attribute_variance`), and to the risk that components share by that of
    the common score (`"cost_risk.common"`); a scenario may give a score a
    value (such as 1.28, for a growth at the 0.9 quantile of a component
    without correlation).

    Args:
        *models (CostModel): The cost models (such as the models of the
            designs of a study, whose scores are shared by name).

    Returns:
        DesignSpace: the group of scores, by name

    Raises:
        ValueError: if a name has a character of paths (`.`, `[`, or `]`), or
            is reserved (`index`, `count`, or `choice`).
    """
    names = list(
        dict.fromkeys(name for model in models for name in model.score_names())
    )
    invalid = [n for n in names if any(c in n for c in _PATH_CHARACTERS)]
    if invalid:
        raise ValueError(f"names of scores must not have '.', '[', or ']': {invalid}")
    reserved = sorted(_RESERVED_NAMES.intersection(names))
    if reserved:
        raise ValueError(
            f"names of scores must not be those of values of groups ({sorted(_RESERVED_NAMES)}): "
            + f"{reserved}"
        )
    return DesignSpace(**{name: Uncertain.normal(0.0, 1.0) for name in names})

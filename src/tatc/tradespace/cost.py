"""
Life-cycle cost of constellations for trade studies: the development,
production, launch, and operations of the satellites of a deployment (see
`DeploymentSchedule`), with learning in production, replacements over the
mission, and spares, dated for their present value and annual costs.

TAT-C does not include cost-estimating relationships, whose common forms are
restricted or need calibration to an organization's experience: the costs of
each type of satellite (such as the satellites of each payload) and of
launches and operations are given, as values or functions of a design (such
as a mass-based relationship), in a currency of the user's choice and
constant-year amounts. Cost risk enters as uncertain growth factors (see
`CostModel.estimate`), drawn for each replicate like other uncertain
parameters, so that statistics of the cost over replicates (such as a 0.7
`Quantile`) are objectives or requirements.

A deployment dates the costs: development of each type of satellite,
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

import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from .deployment import Deployment
from .mission import MissionHistory

ELEMENTS = ("development", "production", "launch", "operations")
"""The elements of the cost of a constellation."""


@dataclass(frozen=True)
class SatelliteCost:
    """
    The cost of a type of satellite (bus and payload): a nonrecurring
    `development` cost, and the recurring cost of its units, from that of
    the first (`first_unit`, the theoretical first unit) by a learning curve:
    unit `k` costs `first_unit * k ** log2(learning_rate)`, so that each
    doubling of the units built costs `learning_rate` times as much per unit
    (typically 0.85 to 0.95 for satellites, 1 for no learning).

    Args:
        development (float): The nonrecurring cost.
        first_unit (float): The recurring cost of the first unit.
        learning_rate (float): The learning rate, in (0, 1].
        mass (float | None): The mass (kg), for launches priced by mass.
    """

    development: float = 0.0
    first_unit: float = 0.0
    learning_rate: float = 1.0
    mass: float | None = None

    def __post_init__(self):
        if not 0 < self.learning_rate <= 1:
            raise ValueError("learning_rate must be in (0, 1]")
        if self.development < 0 or self.first_unit < 0:
            raise ValueError("costs must not be negative")

    def production(self, units: int) -> float:
        """
        Gets the recurring cost of building units (the sum of the costs of
        units 1 to `units` on the learning curve).

        Args:
            units (int): The number of units.

        Returns:
            float: the cost
        """
        if units <= 0:
            return 0.0
        exponent = math.log2(self.learning_rate)
        return float(self.first_unit * np.sum(np.arange(1, units + 1) ** exponent))


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


Growth = float | Mapping[str, float]
"""Growth factors of costs: one for every element, or by element (see
`ELEMENTS`, others 1)."""


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

    Args:
        satellites (SatelliteCost | Mapping[str, SatelliteCost]): The cost of
            the satellites, or of each type of satellite (by name).
        launch (LaunchCost): The cost of launches.
        operations (OperationsCost): The cost of operations.
        discount_rate (float): The annual discount rate (such as 0.03).
        start (float): The years from the present to the start of
            development.
    """

    satellites: SatelliteCost | Mapping[str, SatelliteCost]
    launch: LaunchCost = LaunchCost()
    operations: OperationsCost = OperationsCost()
    discount_rate: float = 0.0
    start: float = 0.0

    def __post_init__(self):
        if self.discount_rate <= -1:
            raise ValueError("discount_rate must be greater than -1")

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
        spare waiting).

        Growth factors multiply the estimate, such as uncertain parameters of
        cost growth drawn for each replicate (for example,
        `Uncertain(lambda rng: rng.triangular(0.9, 1.0, 1.5))`), whose
        statistics over the replicates give the cost at a confidence level
        (such as the 0.7 `Quantile`). Do not use uncertain growth in
        constraints of a design space, which see their nominal values
        (`NaN`).

        Args:
            deployment (Deployment | MissionHistory): The deployment (see
                `DeploymentSchedule.plan`), or a simulated mission.
            growth (float | Mapping[str, float]): A growth factor of the
                cost, or of each element (`"development"`, `"production"`,
                `"launch"`, or `"operations"`; others 1).

        Returns:
            CostEstimate: the estimate
        """
        factors = _get_growth(growth)
        history = (
            deployment.simulate() if isinstance(deployment, Deployment) else deployment
        )
        estimate = CostEstimate()
        flows: list[tuple[float, str, float]] = []

        def add(element: str, amounts: list[tuple[float, float]]) -> None:
            flows.extend(
                (time, element, amount * factors[element]) for time, amount in amounts
            )

        if self.launch.per_kg > 0 and any(
            self._get_type(kind).mass is None for kind in history.development
        ):
            raise ValueError("launches priced by mass require masses")
        for kind, developed in history.development.items():
            add(
                "development", _spread(self._get_type(kind).development, 0.0, developed)
            )
        built: dict[Any, int] = {}
        for batch in history.production:
            cost = self._get_type(batch.kind)
            before = built.get(batch.kind, 0)
            amount = cost.production(before + batch.satellites) - cost.production(
                before
            )
            built[batch.kind] = before + batch.satellites
            add("production", _spread(amount, batch.begin, batch.end))
        estimate.units = {
            "satellites" if kind is None else str(kind): units
            for kind, units in built.items()
        }
        begin, end = history.begin, history.end
        for launch in history.launches:
            cost = self._get_type(launch.kind)
            estimate.launches += int(self.launch.per_launch > 0)
            price = self.launch.per_launch + launch.satellites * (
                self.launch.per_satellite + self.launch.per_kg * (cost.mass or 0.0)
            )
            add("launch", [(launch.time, price)])
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


def _get_growth(growth: Growth) -> dict[str, float]:
    """
    Gets the growth factor of each element of cost.
    """
    if isinstance(growth, Mapping):
        unknown = set(growth) - set(ELEMENTS)
        if unknown:
            raise ValueError(f"growth of unknown elements: {sorted(unknown)}")
        return {element: float(growth.get(element, 1.0)) for element in ELEMENTS}
    return {element: float(growth) for element in ELEMENTS}

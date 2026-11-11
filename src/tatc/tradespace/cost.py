"""
Life-cycle cost of constellations for trade studies: the development,
production, launch, and operations of the satellites of their layouts (see
`walker_space`, `sso_planes_space`, and `planes_space`), with learning in
production, replenishment over the mission, and spares.

TAT-C does not include cost-estimating relationships, whose common forms are
restricted or need calibration to an organization's experience: the costs of
each type of satellite (such as the satellites of each payload) and of
launches and operations are given, as values or functions of a design (such
as a mass-based relationship), in a currency of the user's choice and
constant-year amounts. Cost risk enters as uncertain growth factors (see
`CostModel.estimate`), drawn for each replicate like other uncertain
parameters, so that statistics of the cost over replicates (such as a 0.7
`Quantile`) are objectives or requirements.

An optional schedule (`CostSchedule`) dates the costs from rates, so that
each design's schedule follows from its satellites and planes: development
of each type of satellite, the production of its satellites on a line at a
production rate, launches at a cadence as satellites are ready, replacements
at the end of their design life, and operations from launch to the end of
the mission. The dated costs give the present value of the cost at a
discount rate, and the cost of each year (such as its peak, for annual
budgets). A schedule that is itself a decision (such as the launch cadence)
is a design variable, with the model built for each design::

    def evaluate(built):
        d, satellites = built
        schedule = dataclasses.replace(SCHEDULE, launch_interval=d.launch_interval)
        model = dataclasses.replace(MODEL, schedule=schedule)
        return model.estimate(d.layout).metrics()

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from .constellations import PLANES, SSO_PLANES, WALKER, _kind, _layout, _layouts
from .values import Values

ELEMENTS = ("development", "production", "launch", "operations")
"""The elements of the cost of a constellation."""


def plane_sizes(layout: Any) -> list[int]:
    """
    Gets the number of satellites of each orbit plane of a layout (the
    values of `walker_space`, `sso_planes_space`, or `planes_space`, or of a
    choice of them), or of several (a list or a repetition of layouts).

    Args:
        layout (Any): The values of the layout (or layouts).

    Returns:
        list[int]: the number of satellites of each plane
    """
    several = _layouts(layout)
    if several is not None:
        return [size for item in several for size in plane_sizes(item)]
    layout = _layout(layout)
    kind = _kind(layout)
    if kind == WALKER:
        return [int(layout.per_plane)] * int(layout.planes)
    if kind in (SSO_PLANES, PLANES):
        return [int(plane.satellites) for plane in layout.planes]
    raise ValueError(f"unknown layout {kind!r}")


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
        development_years (float | None): The duration of its development,
            with a schedule (see `CostSchedule`), if not the schedule's.
    """

    development: float = 0.0
    first_unit: float = 0.0
    learning_rate: float = 1.0
    mass: float | None = None
    development_years: float | None = None

    def __post_init__(self):
        if not 0 < self.learning_rate <= 1:
            raise ValueError("learning_rate must be in (0, 1]")
        if self.development < 0 or self.first_unit < 0:
            raise ValueError("costs must not be negative")
        if self.development_years is not None and self.development_years < 0:
            raise ValueError("development_years must not be negative")

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
    The cost of launching satellites: dedicated launches of each plane (as
    many as it takes to carry its satellites, given a `capacity` in number
    or mass), and prices per satellite and per kg (as of rideshares). Give
    the fields of the launch service of a study and leave the others zero.

    Args:
        per_launch (float): The cost of a dedicated launch.
        capacity (int | None): The most satellites of a dedicated launch, if
            limited.
        mass_capacity (float | None): The most mass (kg) of a dedicated
            launch, if limited.
        per_satellite (float): The cost of launching each satellite.
        per_kg (float): The cost of launching each kg of satellite.
    """

    per_launch: float = 0.0
    capacity: int | None = None
    mass_capacity: float | None = None
    per_satellite: float = 0.0
    per_kg: float = 0.0

    def launches(self, satellites: int, mass: float | None) -> int:
        """
        Gets the dedicated launches of a plane of satellites (0 without a
        cost per launch).

        Args:
            satellites (int): The number of satellites of the plane.
            mass (float | None): The mass (kg) of each satellite.

        Returns:
            int: the number of launches
        """
        if self.per_launch == 0 or satellites <= 0:
            return 0
        count = 1
        if self.capacity is not None:
            count = max(count, math.ceil(satellites / self.capacity))
        if self.mass_capacity is not None:
            if mass is None:
                raise ValueError("a launch capacity in mass requires masses")
            per_launch = math.floor(self.mass_capacity / mass)
            if per_launch < 1:
                raise ValueError("a satellite exceeds the mass capacity of a launch")
            count = max(count, math.ceil(satellites / per_launch))
        return count


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


@dataclass(frozen=True)
class CostSchedule:
    """
    A schedule of the costs of a constellation, from rates, and the discount
    rate of their present value.

    Times are in years from the start of development, which begins `start`
    years after the present (the date of present values). The schedule of a
    design follows from its satellites and planes:

    * Each type of satellite is developed over `development_years` (or its
      own, see `SatelliteCost`), its cost spread evenly.
    * Its satellites are then produced on one line, in the order of the
      planes, at `production_rate` satellites a year (or at once, if
      `None`), each plane's paid for evenly while they are produced.
    * Planes are launched in the order they are ready, at most one launch
      every `launch_interval` years: a plane of several dedicated launches
      (see `LaunchCost`) takes several, and its satellites are spread over
      them. A launch is paid for when it flies.
    * The mission lasts `years` (see `CostModel`) from the first launch.
      Each satellite is replaced at the end of its design life (if before
      the end of the mission), its replacement produced at the production
      rate just before it is needed and launched then.
    * Operations cost a fixed amount in each year of the mission, and each
      satellite from its launch to the end of the mission.

    Costs spread over a period are dated at the middle of each of its
    years. A cost at time `t` has a present value of `cost / (1 +
    discount_rate) ** (start + t)`.

    Args:
        discount_rate (float): The annual discount rate (such as 0.03).
        development_years (float): The duration of development, before
            production.
        production_rate (float | None): The satellites of each type produced
            a year, or `None` for production that takes no time.
        launch_interval (float): The least time (years) between launches.
        start (float): The years from the present to the start of
            development.
    """

    discount_rate: float = 0.0
    development_years: float = 3.0
    production_rate: float | None = None
    launch_interval: float = 0.0
    start: float = 0.0

    def __post_init__(self):
        if self.discount_rate <= -1:
            raise ValueError("discount_rate must be greater than -1")
        for name in ("development_years", "launch_interval"):
            if getattr(self, name) < 0:
                raise ValueError(f"{name} must not be negative")
        if self.production_rate is not None and self.production_rate <= 0:
            raise ValueError("production_rate must be positive")

    def producing(self, satellites: int) -> float:
        """
        Gets the time to produce satellites on one line.

        Args:
            satellites (int): The number of satellites.

        Returns:
            float: the time (years)
        """
        return (
            0.0 if self.production_rate is None else satellites / self.production_rate
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
    """int: The dedicated launches."""
    flows: list[tuple[float, str, float]] = field(default_factory=list)
    """list: The dated costs, with a schedule (see `CostSchedule`): the time
    (years from the start of development), element, and amount of each."""
    present: dict[str, float] | None = None
    """dict | None: The present value of each element, with a schedule."""
    deployments: list[tuple[float, int]] = field(default_factory=list)
    """list: The launches of the first satellites of the planes, with a
    schedule: the time (years from the start of development) and number of
    satellites of each."""

    @property
    def deployed(self) -> float | None:
        """float | None: The time (years from the start of development) by
        which every plane is deployed, with a schedule."""
        return max((time for time, _ in self.deployments), default=None)

    @property
    def total(self) -> float:
        """float: The total cost."""
        return self.development + self.production + self.launch + self.operations

    @property
    def present_value(self) -> float | None:
        """float | None: The present value of the total cost, with a
        schedule."""
        return None if self.present is None else sum(self.present.values())

    def annual(self) -> pd.DataFrame:
        """
        Gets the cost of each year (from the start of development), by
        element and in total, with a schedule.

        Returns:
            pandas.DataFrame: a row per year, a column per element, and the
                `total`
        """
        if not self.flows:
            raise ValueError("annual costs require a schedule (see CostSchedule)")
        frame = pd.DataFrame(self.flows, columns=["time", "element", "amount"])
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
        (`<prefix>`) and each element (`<prefix> development`, `<prefix>
        production`, `<prefix> launch`, and `<prefix> operations`), and, with
        a schedule, the present value of the total (`<prefix> present`) and
        the cost of its most expensive year (`<prefix> peak annual`).

        Args:
            prefix (str): The name of the total, and prefix of the elements.

        Returns:
            dict[str, float]: the metrics
        """
        metrics = {
            prefix: self.total,
            **{f"{prefix} {element}": getattr(self, element) for element in ELEMENTS},
        }
        if self.present is not None:
            metrics[f"{prefix} present"] = self.present_value
            metrics[f"{prefix} peak annual"] = float(self.annual()["total"].max())
        return metrics


Growth = float | Mapping[str, float]
"""Growth factors of costs: one for every element, or by element (see
`ELEMENTS`, others 1)."""


@dataclass(frozen=True)
class CostModel:
    """
    A model of the life-cycle cost of constellations: the costs of their
    types of satellites, launches, and operations over a mission of `years`.

    Satellites are replaced at the end of their `design_life`: a mission of
    6 years with satellites of a design life of 3 years builds and launches
    two generations of each satellite, whose units continue the learning
    curve of the first.

    Args:
        satellites (SatelliteCost | Mapping[str, SatelliteCost]): The cost of
            the satellites, or of each type of satellite (by name).
        launch (LaunchCost): The cost of launches.
        operations (OperationsCost): The cost of operations.
        years (float): The duration of the mission.
        design_life (float | None): The design life of the satellites, or
            `None` if they last the mission.
        schedule (CostSchedule | None): The schedule of the costs (for their
            present value and annual costs), if any.
    """

    satellites: SatelliteCost | Mapping[str, SatelliteCost]
    launch: LaunchCost = LaunchCost()
    operations: OperationsCost = OperationsCost()
    years: float = 1.0
    design_life: float | None = None
    schedule: CostSchedule | None = None

    def __post_init__(self):
        if self.years <= 0:
            raise ValueError("years must be positive")
        if self.design_life is not None and self.design_life <= 0:
            raise ValueError("design_life must be positive")

    @property
    def generations(self) -> int:
        """int: The generations of satellites built over the mission."""
        if self.design_life is None:
            return 1
        return max(1, math.ceil(self.years / self.design_life - 1e-9))

    def _get_type(self, kind: str | None) -> SatelliteCost:
        if isinstance(self.satellites, SatelliteCost):
            return self.satellites
        if kind is None:
            raise ValueError("a model of several types of satellites needs a kind")
        if kind not in self.satellites:
            raise ValueError(f"no cost of satellites of kind {kind!r}")
        return self.satellites[kind]

    def estimate(
        self,
        layout: Any,
        kind: str | None = None,
        spares: int = 0,
        growth: Growth = 1.0,
    ) -> CostEstimate:
        """
        Estimates the life-cycle cost of a constellation: the layout (or
        layouts) of its satellites (see `plane_sizes`), of one type (`kind`,
        if the model has several), or a mapping of types to layouts (such as
        the shells of a mixed constellation, each of its own payload).

        Each plane has `spares` more satellites (spares in orbit, built,
        launched, and operated with it). Each type of satellite is developed
        once; its units (of every generation) follow its learning curve;
        each plane is launched (by dedicated launches, or per satellite or
        kg) in each generation; and every satellite is operated for the
        mission.

        With a schedule (see `CostSchedule`), the estimate dates its costs,
        and gives their present value, annual costs, and deployments:
        satellites launched after the first are operated, and may be
        replaced, for less of the mission, which lowers their costs.

        Growth factors multiply the estimate, such as uncertain parameters of
        cost growth drawn for each replicate (for example,
        `Uncertain(lambda rng: rng.triangular(0.9, 1.0, 1.5))`), whose
        statistics over the replicates give the cost at a confidence level
        (such as the 0.7 `Quantile`). Do not use uncertain growth in
        constraints of a design space, which see their nominal values
        (`NaN`).

        Args:
            layout (Any): The layout (or layouts) of the satellites, or a
                mapping of types of satellites to layouts.
            kind (str | None): The type of the satellites of a layout, if the
                model has several types.
            spares (int): The spares in orbit in each plane.
            growth (float | Mapping[str, float]): A growth factor of the
                cost, or of each element (`"development"`, `"production"`,
                `"launch"`, or `"operations"`; others 1).

        Returns:
            CostEstimate: the estimate
        """
        if spares < 0:
            raise ValueError("spares must not be negative")
        # a mapping of types to layouts (the values of a layout are a mapping too)
        if isinstance(layout, Mapping) and not isinstance(layout, Values):
            fleets = {name: plane_sizes(item) for name, item in layout.items()}
        else:
            fleets = {kind: plane_sizes(layout)}
        if self.schedule is not None:
            return self._estimate_scheduled(fleets, spares, _get_growth(growth))
        estimate = CostEstimate()
        operated = 0
        for name, sizes in fleets.items():
            cost = self._get_type(name)
            sizes = [size + spares for size in sizes if size > 0]
            if not sizes:
                continue
            units = self.generations * sum(sizes)
            estimate.units["satellites" if name is None else str(name)] = units
            estimate.development += cost.development
            estimate.production += cost.production(units)
            launches = self.generations * sum(
                self.launch.launches(size, cost.mass) for size in sizes
            )
            estimate.launches += launches
            if self.launch.per_kg > 0 and cost.mass is None:
                raise ValueError("launches priced by mass require masses")
            estimate.launch += launches * self.launch.per_launch + units * (
                self.launch.per_satellite + self.launch.per_kg * (cost.mass or 0.0)
            )
            operated += sum(sizes)
        estimate.operations = self.years * (
            self.operations.per_year + self.operations.per_satellite_year * operated
        )
        factors = _get_growth(growth)
        for element in ELEMENTS:
            setattr(estimate, element, getattr(estimate, element) * factors[element])
        return estimate

    def _estimate_scheduled(
        self, fleets: dict, spares: int, factors: dict[str, float]
    ) -> CostEstimate:
        """
        Estimates the dated costs of the planes of each type of satellite
        (see `estimate` and `CostSchedule`).
        """
        # pylint: disable=too-many-locals
        schedule = self.schedule
        estimate = CostEstimate()
        flows: list[tuple[float, str, float]] = []

        def add(element: str, amounts: list[tuple[float, float]]) -> None:
            flows.extend(
                (time, element, amount * factors[element]) for time, amount in amounts
            )

        if self.launch.per_kg > 0 and any(
            self._get_type(name).mass is None for name in fleets
        ):
            raise ValueError("launches priced by mass require masses")
        # development, then production of each type's planes on its line
        planes = []  # (ready, name, cost, satellites)
        batches: dict[Any, list[tuple[float, float, int]]] = {}
        for name, sizes in fleets.items():
            cost = self._get_type(name)
            sizes = [size + spares for size in sizes if size > 0]
            if not sizes:
                continue
            developed = (
                schedule.development_years
                if cost.development_years is None
                else cost.development_years
            )
            add("development", _spread(cost.development, 0.0, developed))
            line = developed
            for size in sizes:
                ready = line + schedule.producing(size)
                batches.setdefault(name, []).append((line, ready, size))
                planes.append((ready, name, cost, size))
                line = ready
        if not planes:
            estimate.flows, estimate.present = [], dict.fromkeys(ELEMENTS, 0.0)
            return estimate
        # launches in the order planes are ready, at the launch cadence
        launches = []  # (time, name, cost, satellites, dedicated)
        cadence = -math.inf
        for ready, name, cost, size in sorted(planes, key=lambda plane: plane[0]):
            count = self.launch.launches(size, cost.mass) or 1
            first = max(ready, cadence)
            for j in range(count):
                carried = size // count + (1 if j < size % count else 0)
                launches.append(
                    (
                        first + j * schedule.launch_interval,
                        name,
                        cost,
                        carried,
                        self.launch.per_launch > 0,
                    )
                )
            cadence = first + count * schedule.launch_interval
            estimate.deployments.append((first, size))
        originals = list(launches)
        begin = min(time for time, *_ in launches)
        end = begin + self.years
        # replacements at the end of each design life, before the end
        for time, name, cost, carried, dedicated in originals:
            if self.design_life is None:
                break
            need = time + self.design_life
            while need < end - 1e-9:
                produced = need - schedule.producing(carried)
                batches.setdefault(name, []).append((produced, need, carried))
                launches.append((need, name, cost, carried, dedicated))
                need += self.design_life
        # production in order along each type's learning curve
        for name, items in batches.items():
            cost = self._get_type(name)
            built = 0
            for produced, ready, size in sorted(items, key=lambda item: item[0]):
                amount = cost.production(built + size) - cost.production(built)
                built += size
                add("production", _spread(amount, produced, ready))
            estimate.units["satellites" if name is None else str(name)] = built
        for time, name, cost, carried, dedicated in launches:
            estimate.launches += int(dedicated)
            price = (self.launch.per_launch if dedicated else 0.0) + carried * (
                self.launch.per_satellite + self.launch.per_kg * (cost.mass or 0.0)
            )
            add("launch", [(time, price)])
        # operations of the mission, and of each satellite from its launch
        add("operations", _spread(self.operations.per_year * self.years, begin, end))
        for time, _, _, carried, _ in originals:
            add(
                "operations",
                _spread(
                    self.operations.per_satellite_year * carried * (end - time),
                    time,
                    end,
                ),
            )
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
                    amount * schedule.discount(time)
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

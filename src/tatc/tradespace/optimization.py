"""
Searches of design spaces for Pareto-efficient designs with optimizers:
Optuna (which searches the hierarchical design space directly, see
`DesignSpace.suggest`) and pymoo (which searches the unit hypercube that
encodes it, see `DesignSpace.decode`). Both are optional dependencies.

The objectives of a design are statistics over its replicates of the
metrics of its evaluations (see `Study.run`, which must give one record per
evaluation): their means, or other statistics (see `aggregate`), such as a
quantile or the probability of meeting a threshold. Requirements on
statistics (see `Requirement`) are output constraints, whose violations the
optimizers minimize along with those of the constraints of the design space.
Each design is evaluated once, however often an optimizer proposes it, and
designs proposed together are evaluated together, in parallel if the study
is (see `Study`).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd

from .criteria import (
    Requirement,
    Statistic,
    _check_names,
    _evaluate,
    _get_matrices,
    _get_statistics,
)
from .pareto import _check_directions
from .space import Design, Infeasible
from .study import Study


class Evaluations:
    """
    The evaluations of designs for a search: the objectives of each design
    (by its key), the violations of its requirements, and the records of the
    study's evaluations.

    Args:
        study (Study): The study that evaluates designs.
        objectives (Mapping[str, str]): The objectives: statistics (or
            metrics, whose means they are) and their directions (`"min"` or
            `"max"`).
        statistics (Mapping[str, Statistic | str] | None): Statistics, by
            name, of objectives and requirements that are not means.
        requirements (Sequence[Requirement]): Requirements on statistics.
    """

    def __init__(
        self,
        study: Study,
        objectives: Mapping[str, str],
        statistics: Mapping[str, Statistic | str] | None = None,
        requirements: Sequence[Requirement] = (),
    ):
        _check_directions(objectives)
        self.study, self.objectives = study, dict(objectives)
        self.requirements = list(requirements)
        names = list(objectives) + [r.name for r in self.requirements]
        self.statistics = _get_statistics(names, statistics)
        _check_names(self.statistics, self.requirements)
        self.values: dict[str, tuple[float, ...] | None] = {}
        """dict: The objectives of each evaluated design (by key), or `None`
        if it is infeasible."""
        self.violations: dict[str, tuple[float, ...] | None] = {}
        """dict: The violation of each requirement by each evaluated design
        (by key), or `None` if it is infeasible."""
        self._records: list[pd.DataFrame] = []

    @property
    def results(self) -> pd.DataFrame:
        """pandas.DataFrame: The results of every evaluation (see
        `Study.run`)."""
        if not self._records:
            return pd.DataFrame()
        return pd.concat(self._records, ignore_index=True)

    def evaluate(self, designs: list[Design]) -> list[tuple[float, ...] | None]:
        """
        Gets the objectives of designs, evaluating those not yet evaluated
        together (see also `violations`).

        Args:
            designs (list[Design]): The designs.

        Returns:
            list[tuple[float, ...] | None]: the objectives of each design, in
                the order of the objectives, or `None` if it is infeasible
        """
        new = list({d.key: d for d in designs if d.key not in self.values}.values())
        if new:
            records = self.study.run(new)
            self._records.append(records)
            if records.groupby(["design", "replicate"]).size().max() > 1:
                raise ValueError(
                    "objectives need one record per evaluation: make the "
                    + "study's evaluate function return a mapping of metrics"
                )
            infeasible = set(records.loc[~records["feasible"], "design"])
            for key in infeasible:
                self.values[key] = self.violations[key] = None
            feasible = records[~records["design"].isin(infeasible)]
            if len(feasible) > 0:
                metrics = sorted({s.column for s in self.statistics.values()})
                keys, _, matrices = _get_matrices(feasible, metrics)
                values, _ = _evaluate(self.statistics, [], matrices)
                for i, key in enumerate(keys):
                    self.values[key] = tuple(
                        float(values[name][i]) for name in self.objectives
                    )
                    self.violations[key] = tuple(
                        float(requirement.violation(values[requirement.name][i]))
                        for requirement in self.requirements
                    )
        return [self.values[d.key] for d in designs]

    def worst(self) -> list[float]:
        """
        Gets the worst value of each objective of the feasible designs
        evaluated (0 if none is), as for infeasible designs.

        Returns:
            list[float]: the worst values, in the order of the objectives
        """
        feasible = [
            values
            for values in self.values.values()
            if values is not None and not np.isnan(values).any()
        ]
        if not feasible:
            return [0.0] * len(self.objectives)
        table = np.array(feasible)
        return [
            float(table[:, i].max() if direction == "min" else table[:, i].min())
            for i, direction in enumerate(self.objectives.values())
        ]


def optimize_optuna(
    study: Study,
    objectives: Mapping[str, str],
    n_trials: int,
    sampler: Any = None,
    optuna_study: Any = None,
    batch_size: int | None = None,
    constraints: bool = True,
    statistics: Mapping[str, Statistic | str] | None = None,
    requirements: Sequence[Requirement] = (),
) -> tuple[Any, pd.DataFrame]:
    """
    Searches a study's design space with Optuna: each trial suggests a design
    (see `DesignSpace.suggest`) and is told the means of its objectives over
    its replicates (or other statistics, see `statistics`).

    With `constraints`, each trial records the violation of the constraints
    of its design (`trial.set_constraint`, under the name `"feasibility"`):
    0 for a feasible design, the amount by which it violates a constraint of
    the design space (see `DesignSpace`), or 1 if the study finds it
    infeasible. Optuna's constraint-aware samplers (such as TPE and NSGA-II)
    then prefer feasible designs, and infeasible designs that violate their
    constraints less. An infeasible trial is told the worst objectives of the
    feasible designs evaluated before it, which these samplers do not use to
    rank it. Without `constraints` (or with a version of Optuna before 5,
    which lacks `set_constraint`), an infeasible design prunes its trial.

    Each requirement is a constraint of the trials (named by its statement,
    such as `"sufficient p10 >= 80"`), whose value is its violation (see
    `Requirement.violation`), and which requires `constraints`.

    Trials are asked in batches, whose designs are evaluated together (in
    parallel, if the study is) before they are told. A batch defaults to one
    trial per worker of the study.

    Args:
        study (Study): The study that evaluates designs.
        objectives (Mapping[str, str]): The objectives: columns of the
            results and their directions (`"min"` or `"max"`).
        n_trials (int): The number of trials.
        sampler (optuna.samplers.BaseSampler | None): The sampler, if not
            Optuna's default (TPE with one objective, NSGA-II with several);
            ignored with an `optuna_study`.
        optuna_study (optuna.Study | None): A study to continue, whose
            directions are those of the objectives.
        batch_size (int | None): The number of trials evaluated together.
        constraints (bool): `True` to record the violations of constraints,
            `False` to prune infeasible trials.
        statistics (Mapping[str, Statistic | str] | None): Statistics, by
            name, of objectives and requirements that are not means.
        requirements (Sequence[Requirement]): Requirements on statistics.

    Returns:
        tuple[optuna.Study, pandas.DataFrame]: the Optuna study, and the
            results of every evaluation (see `Study.run`)
    """
    # pylint: disable=import-outside-toplevel
    try:
        import optuna
    except ImportError as error:
        raise ImportError(
            "optimize_optuna requires optuna (pip install optuna)"
        ) from error
    evaluations = Evaluations(study, objectives, statistics, requirements)
    directions = [
        "minimize" if direction == "min" else "maximize"
        for direction in evaluations.objectives.values()
    ]
    if optuna_study is None:
        optuna_study = optuna.create_study(directions=directions, sampler=sampler)
    batch_size = batch_size or study.workers
    constrained = constraints and hasattr(optuna.trial.Trial, "set_constraint")
    if requirements and not constrained:
        raise ValueError("requirements need constraints (and Optuna 5 or later)")
    names = [f"{r.name} {r.operator} {r.value}" for r in evaluations.requirements]

    def set_constraints(trial: Any, feasibility: float, violations) -> None:
        trial.set_constraint("feasibility", feasibility)
        for name, violation in zip(names, violations):
            trial.set_constraint(name, violation)

    def infeasible(trial: Any, violation: float) -> None:
        if constrained:
            set_constraints(trial, max(violation, 1e-9), [0.0] * len(names))
            optuna_study.tell(trial, evaluations.worst())
        else:
            optuna_study.tell(trial, state=optuna.trial.TrialState.PRUNED)

    remaining = n_trials
    while remaining > 0:
        trials, designs = [], []
        for _ in range(min(batch_size, remaining)):
            trial = optuna_study.ask()
            try:
                design = study.space.suggest(trial)
            except Infeasible as error:
                trial.set_user_attr("error", str(error))
                infeasible(trial, error.violation)
                continue
            trial.set_user_attr("design", design.key)
            trials.append(trial)
            designs.append(design)
        remaining -= min(batch_size, remaining)
        for trial, design, values in zip(
            trials, designs, evaluations.evaluate(designs)
        ):
            if values is None:
                infeasible(trial, 1.0)
            else:
                if constrained:
                    violations = evaluations.violations[design.key]
                    set_constraints(trial, 0.0, np.nan_to_num(violations, nan=1.0))
                optuna_study.tell(
                    trial, list(np.where(np.isnan(values), evaluations.worst(), values))
                )
    return optuna_study, evaluations.results


def make_pymoo_problem(
    study: Study,
    objectives: Mapping[str, str],
    statistics: Mapping[str, Statistic | str] | None = None,
    requirements: Sequence[Requirement] = (),
) -> Any:
    """
    Makes a pymoo problem of a study's design space: its variables are the
    unit hypercube that encodes the design space (see `DesignSpace.decode`
    and `DesignSpace.unit_paths`), its objectives are minimized (maximized
    objectives are negated), and its constraint (`G <= 0`) is violated by
    designs that violate a constraint of the design space (by the amount of
    the violation, see `DesignSpace`) or that the study finds infeasible (by
    1). Each requirement is a further constraint, whose value is its
    violation (see `Requirement.violation`). The designs of a population are
    evaluated together.

    The problem keeps its evaluations (`evaluations`, see `Evaluations`),
    whose `results` are those of every evaluation (see `Study.run`).

    Args:
        study (Study): The study that evaluates designs.
        objectives (Mapping[str, str]): The objectives: statistics (or
            metrics) and their directions (`"min"` or `"max"`).
        statistics (Mapping[str, Statistic | str] | None): Statistics, by
            name, of objectives and requirements that are not means.
        requirements (Sequence[Requirement]): Requirements on statistics.

    Returns:
        pymoo.core.problem.Problem: the problem
    """
    # pylint: disable=import-outside-toplevel
    try:
        from pymoo.core.problem import Problem
    except ImportError as error:
        raise ImportError(
            "make_pymoo_problem requires pymoo (pip install pymoo)"
        ) from error

    evaluations = Evaluations(study, objectives, statistics, requirements)
    n_requirements = len(evaluations.requirements)
    signs = np.array(
        [1.0 if direction == "min" else -1.0 for direction in objectives.values()]
    )
    n_var = len(study.space.unit_paths())

    class DesignSpaceProblem(Problem):
        """
        A pymoo problem of a design space (see `make_pymoo_problem`).
        """

        def __init__(self):
            super().__init__(
                n_var=n_var,
                n_obj=len(signs),
                n_ieq_constr=1 + n_requirements,
                xl=0.0,
                xu=1.0,
            )
            self.evaluations = evaluations

        def _evaluate(self, x, out, *args, **kwargs):
            designs: list[Design | None] = []
            violation = np.ones(len(np.atleast_2d(x)))
            for i, row in enumerate(np.atleast_2d(x)):
                try:
                    designs.append(study.space.decode(row))
                except Infeasible as error:
                    designs.append(None)
                    violation[i] = max(error.violation, 1e-9)
            feasible = [design for design in designs if design is not None]
            values = iter(evaluations.evaluate(feasible))
            costs = np.zeros((len(designs), len(signs)))
            unmet = np.zeros((len(designs), n_requirements))
            for i, design in enumerate(designs):
                value = None if design is None else next(values)
                if value is not None:
                    value = np.where(np.isnan(value), evaluations.worst(), value)
                    costs[i] = signs * np.asarray(value)
                    violation[i] = 0.0
                    unmet[i] = np.nan_to_num(
                        evaluations.violations[design.key], nan=1.0
                    )
            out["F"] = costs
            out["G"] = np.column_stack([violation, unmet])

    return DesignSpaceProblem()


def optimize_pymoo(
    study: Study,
    objectives: Mapping[str, str],
    algorithm: Any = None,
    termination: Any = ("n_gen", 10),
    seed: int | None = None,
    statistics: Mapping[str, Statistic | str] | None = None,
    requirements: Sequence[Requirement] = (),
    **kwargs: Any,
) -> tuple[Any, pd.DataFrame]:
    """
    Searches a study's design space with a pymoo algorithm (see
    `make_pymoo_problem`).

    Args:
        study (Study): The study that evaluates designs.
        objectives (Mapping[str, str]): The objectives: columns of the
            results and their directions (`"min"` or `"max"`).
        algorithm (pymoo.core.algorithm.Algorithm | None): The algorithm, or
            `None` for NSGA-II with populations of 20.
        termination (Any): The termination criterion of `pymoo.optimize.minimize`.
        seed (int | None): The seed of the algorithm.
        statistics (Mapping[str, Statistic | str] | None): Statistics, by
            name, of objectives and requirements that are not means.
        requirements (Sequence[Requirement]): Requirements on statistics.
        **kwargs (Any): Other arguments of `pymoo.optimize.minimize`.

    Returns:
        tuple[pymoo.core.result.Result, pandas.DataFrame]: the result of the
            search, and the results of every evaluation (see `Study.run`)
    """
    # pylint: disable=import-outside-toplevel
    from pymoo.optimize import minimize

    if algorithm is None:
        from pymoo.algorithms.moo.nsga2 import NSGA2

        algorithm = NSGA2(pop_size=20)
    problem = make_pymoo_problem(study, objectives, statistics, requirements)
    result = minimize(problem, algorithm, termination, seed=seed, **kwargs)
    return result, problem.evaluations.results

"""
Searches of design spaces for Pareto-efficient designs with optimizers:
Optuna (which searches the hierarchical design space directly, suggesting
each decision in turn) and pymoo (which searches the unit hypercube that
encodes it, a dimension per decision). Both are optional dependencies.

Designs are judged by criteria (see `Criteria`, without `by`): their
objectives are statistics over their replicates of the metrics of their
evaluations (see `Study.run`, which must give one record per evaluation):
their means, or other statistics, such as a quantile or the probability of
meeting a threshold. Requirements on statistics (see `Requirement`) are
output constraints, whose violations the optimizers minimize along with
those of the constraints of the design space. The tolerances of the
criteria do not apply to the optimizers' own ranking of designs, but they
do to the front of the designs they evaluate (see `Search`).

A search of a study with several scenarios (see `Study`) evaluates each
design in every scenario, and judges it by criteria `across` the scenarios
(see `Across`): by its worst or its mean over them, such as for designs
robust to scenarios that cannot be told apart in advance.
Each design is evaluated once, however often an optimizer proposes it, and
designs proposed together are evaluated together, in parallel if the study
is (see `Study`).

Designs are judged by their feasible replicates, as summaries judge them
(see `summarize`): a design infeasible in some replicates is judged by the
others (a requirement on `feasible`, the share of its replicates that are
feasible, such as `Requirement("feasible", ">=", 1)`, rejects it), and a
design without a feasible replicate is infeasible, which the optimizers'
constraint of feasibility rejects.

A search chooses its front from one evaluation of each design (with the
study's replicates), which favors the designs whose replicates were lucky:
refine the designs it evaluated (see `refine`) and confirm its front (see
`confirm`, of the designs of the search, see `Search`), as for any study.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import copy
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from .criteria import Criteria, Requirement, _as_criteria, _signs
from .space import Design, Infeasible
from .study import Study
from .summary import _summarize
from .tables import FEASIBLE, check_records, efficient_keys


@dataclass
class Search:
    """
    The outcome of a search with an optimizer (see `optimize_optuna` and
    `optimize_pymoo`): as the outcomes of other steps of a trade study (see
    `Refinement`), its `results`, the `summary` of its designs, the `designs`
    it selects for the next step (the efficient designs it found), and its
    `criteria`, which refine and confirm its front as for any study::

        search = optimize_optuna(study, criteria, n_trials=200)
        refinement = refine(study, criteria, search.results)
        confirmation = confirm(study, criteria, refinement.results)
    """

    results: pd.DataFrame
    """pandas.DataFrame: The results of every evaluation (see `Study.run`)."""
    summary: pd.DataFrame
    """pandas.DataFrame: A row per design that the search found feasible
    (with a feasible replicate, in every scenario with criteria `across`
    scenarios, as designs are judged by their feasible replicates, see
    `summarize`): the summary of the criteria over its replicates (with the
    spreads of means, as `summarize` gives them), with whether it is
    efficient among these designs (`pareto`, with the tolerances of the
    criteria)."""
    designs: list[Design]
    """list[Design]: The efficient designs that the search found (`pareto` of
    `summary`), for the next step."""
    criteria: Criteria
    """Criteria: The criteria of the search."""
    evaluated: list[Design] = field(default_factory=list)
    """list[Design]: The designs evaluated, in the order of their evaluation."""
    state: Any = None
    """Any: The optimizer's own outcome: the Optuna study (`optuna.Study`) or
    the result of pymoo (`pymoo.core.result.Result`)."""


class Evaluations:
    """
    The evaluations of designs for a search: the objectives of each design
    (by its key), the violations of its requirements, and the records of the
    study's evaluations. Designs are judged by their feasible replicates (see
    `summarize`), and a design without a feasible replicate is infeasible.

    Args:
        study (Study): The study that evaluates designs.
        criteria (Criteria | Mapping[str, str]): The criteria, without `by`
            (see `Criteria`), or the objectives (metrics whose means they
            are) and their directions (`"min"` or `"max"`).
    """

    def __init__(self, study: Study, criteria: Criteria | Mapping[str, str]):
        self.study, self.criteria = study, _as_criteria(criteria)
        if self.criteria.by:
            raise ValueError(
                "searches need one record per evaluation: give criteria without `by`"
            )
        self.values: dict[str, tuple[float, ...] | None] = {}
        """dict: The objectives of each evaluated design (by key), or `None`
        if it is infeasible."""
        self.violations: dict[str, tuple[float, ...] | None] = {}
        """dict: The violation of each requirement by each evaluated design
        (by key), or `None` if it is infeasible."""
        self.infeasibility: dict[str, float] = {}
        """dict: The violation of the feasibility of each evaluated design (by
        key), for constraint-aware optimizers: 0 for a design that has
        feasible replicates (by which it is judged), and otherwise the share
        of its evaluations that are infeasible (1 for a design infeasible in
        every replicate, and less, with criteria `across` scenarios, for one
        without a feasible replicate in some of them)."""
        self.designs: dict[str, Design] = {}
        """dict: The evaluated designs (by key), in the order of their
        evaluation."""
        self._records: list[pd.DataFrame] = []
        self._checked = False

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
        together (see also `violations`). A design is infeasible if none of
        its replicates is feasible (see `infeasibility`). The first design
        evaluated is evaluated alone, to check that the study's evaluations
        give one record each before evaluating others.

        Args:
            designs (list[Design]): The designs.

        Returns:
            list[tuple[float, ...] | None]: the objectives of each design, in
                the order of the objectives, or `None` if it is infeasible
        """
        new = list({d.key: d for d in designs if d.key not in self.values}.values())
        if new and not self._checked:
            self._run(new[:1])
            new = new[1:]
        if new:
            self._run(new)
        return [self.values[d.key] for d in designs]

    def _run(self, designs: list[Design]) -> None:
        """
        Evaluates designs, recording their objectives and violations (of
        their summary, see `summarize`).
        """
        self.designs.update((design.key, design) for design in designs)
        records = self.study.run(designs)
        self._records.append(records)
        across = self.criteria.across
        keys = ["design", "replicate", *([] if across is None else [across.column])]
        if across is None and "scenario" in records.columns:
            hint = (
                "a study with several scenarios gives a record in each: give "
                + "criteria `across` the scenarios (see Across), or search a study "
                + "with one scenario"
            )
        else:
            hint = (
                "searches need one record per evaluation: make the study's evaluate "
                + "function return a mapping of metrics"
            )
        check_records(records, keys, hint, self.criteria.metrics)
        self._checked = self._checked or bool(records[FEASIBLE].any())
        summary = _summarize(records, self.criteria, design_values=False)
        summary = summary.set_index("design")
        shares = records.groupby("design")[FEASIBLE].mean()
        for design in designs:
            key = design.key
            if key not in summary.index:
                self.values[key] = self.violations[key] = None
                self.infeasibility[key] = 1.0 - float(shares[key])
                continue
            row = summary.loc[key]
            self.values[key] = tuple(
                float(row[name]) for name in self.criteria.objectives
            )
            self.violations[key] = tuple(
                float(requirement.violation(row[requirement.name]))
                for requirement in self.criteria.requirements
            )
            self.infeasibility[key] = 0.0

    def search(self, state: Any) -> Search:
        """
        Gets the outcome of a search of these evaluations.

        Args:
            state (Any): The optimizer's own outcome.

        Returns:
            Search: the search
        """
        results = self.results
        summary = _summarize(results, self.criteria, spread=True)
        return Search(
            results=results,
            criteria=self.criteria,
            summary=summary,
            designs=[self.designs[key] for key in efficient_keys(summary)],
            evaluated=list(self.designs.values()),
            state=state,
        )

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
            return [0.0] * len(self.criteria.objectives)
        table = np.array(feasible)
        return [
            float(table[:, i].max() if direction == "min" else table[:, i].min())
            for i, direction in enumerate(self.criteria.objectives.values())
        ]


def optimize_optuna(
    study: Study,
    criteria: Criteria | Mapping[str, str],
    *,
    n_trials: int = 100,
    seed: int | None = 0,
    sampler: Any = None,
    optuna_study: Any = None,
    batch_size: int | None = None,
    constraints: bool = True,
    max_designs: int | None = None,
) -> Search:
    """
    Searches a study's design space with Optuna: each trial suggests a design
    (each of its decisions in turn) and is told the objectives of the criteria:
    statistics (such as means) of its metrics over its replicates.

    With `constraints`, each trial records the violation of the constraints
    of its design (`trial.set_constraint`, under the name `"feasibility"`):
    0 for a design with feasible replicates (by which it is judged, see
    `summarize`), the amount by which it violates a constraint of the design
    space (see `DesignSpace`), or, for a design without a feasible replicate,
    the share of its evaluations that the study finds infeasible (see
    `Evaluations.infeasibility`). Optuna's constraint-aware samplers (such
    as TPE and NSGA-II) then prefer feasible designs, and infeasible designs
    that violate their constraints less. An infeasible trial is told the
    worst objectives of the feasible designs evaluated before it, which
    these samplers do not use to rank it. Without `constraints` (or with a
    version of Optuna before 5, which lacks `set_constraint`), an infeasible
    design prunes its trial.

    Each requirement of the criteria is a constraint of the trials (named by
    its statement, `str(requirement)`, such as `"sufficient p10 >= 80"`, or,
    if requirements state themselves alike, by their values in full and
    their order), whose value is its violation (see `Requirement.violation`),
    and which requires `constraints`. A requirement on `feasible` (the share
    of a design's replicates that are feasible) constrains the designs that
    are infeasible in some of their replicates.

    Trials are asked in batches, whose designs are evaluated together (in
    parallel, if the study is) before they are told. A batch defaults to one
    trial per worker of the study. The default sampler is TPE, which, in
    batches, tells its pending trials the worst values evaluated so far
    (`constant_liar`), so that a batch spreads its trials rather than
    suggesting nearly the same design several times.

    Samplers suggest designs evaluated before, especially in small design
    spaces: such trials are told the objectives of the design without
    evaluating it again. The trials (`n_trials`) are then a poor measure of
    the cost of a search; `max_designs` limits the designs evaluated
    instead (the search stops once it has evaluated as many, or after
    `n_trials`).

    Args:
        study (Study): The study that evaluates designs.
        criteria (Criteria | Mapping[str, str]): The criteria, without `by`
            (see `Criteria`), or the objectives (metrics whose means they
            are) and their directions (`"min"` or `"max"`).
        n_trials (int): The number of trials.
        seed (int | None): The seed of the default sampler (`None` for one
            drawn anew); a sampler given is seeded by its own arguments.
        sampler (optuna.samplers.BaseSampler | None): The sampler, or `None`
            for TPE (with a constant liar in batches); ignored with an
            `optuna_study`.
        optuna_study (optuna.Study | None): A study to continue, whose
            directions are those of the objectives.
        batch_size (int | None): The number of trials evaluated together.
        constraints (bool): `True` to record the violations of constraints,
            `False` to prune infeasible trials.
        max_designs (int | None): The most designs to evaluate, if limited.

    Returns:
        Search: the results and designs of the search, and the Optuna study

    Raises:
        ValueError: if the criteria have requirements without `constraints`
            (or with a version of Optuna before 5), or if `max_designs` is
            less than 1.
    """
    # pylint: disable=import-outside-toplevel
    try:
        import optuna
    except ImportError as error:
        raise ImportError(
            "optimize_optuna requires optuna (pip install optuna)"
        ) from error
    if max_designs is not None and max_designs < 1:
        raise ValueError("max_designs must be at least 1")
    evaluations = Evaluations(study, criteria)
    requirements = evaluations.criteria.requirements
    directions = [
        "minimize" if direction == "min" else "maximize"
        for direction in evaluations.criteria.objectives.values()
    ]
    batch_size = batch_size or study.workers
    if optuna_study is None:
        if sampler is None:
            sampler = optuna.samplers.TPESampler(
                constant_liar=batch_size > 1, seed=seed
            )
        optuna_study = optuna.create_study(directions=directions, sampler=sampler)
    constrained = constraints and hasattr(optuna.trial.Trial, "set_constraint")
    if requirements and not constrained:
        raise ValueError("requirements need constraints (and Optuna 5 or later)")
    names = _constraint_names(requirements)

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
        size = min(batch_size, remaining)
        if max_designs is not None:
            size = min(size, max_designs - len(evaluations.values))
            if size <= 0:
                break
        trials, designs = [], []
        for _ in range(size):
            trial = optuna_study.ask()
            try:
                design = study.space._suggest(trial)  # pylint: disable=protected-access
            except Infeasible as error:
                trial.set_user_attr("error", str(error))
                infeasible(trial, error.violation)
                continue
            trial.set_user_attr("design", design.key)
            trials.append(trial)
            designs.append(design)
        remaining -= size
        for trial, design, values in zip(
            trials, designs, evaluations.evaluate(designs)
        ):
            if values is None:
                infeasible(trial, evaluations.infeasibility[design.key])
            else:
                if constrained:
                    violations = evaluations.violations[design.key]
                    set_constraints(trial, 0.0, np.nan_to_num(violations, nan=1.0))
                optuna_study.tell(
                    trial, list(np.where(np.isnan(values), evaluations.worst(), values))
                )
    return evaluations.search(optuna_study)


def _constraint_names(requirements: Sequence[Requirement]) -> list[str]:
    """
    Gets distinct names of the constraints of requirements (other than
    `"feasibility"`): their statements (see `Requirement`), with their values
    in full where statements are alike, and their order where these are too
    (such as of requirements that differ only in their scales).
    """
    names = [str(requirement) for requirement in requirements]
    names = [
        (f"{r.name} {r.operator} {r.value!r}" if names.count(name) > 1 else name)
        for r, name in zip(requirements, names)
    ]
    return [
        name if names[:i].count(name) == 0 else f"{name} ({names[:i].count(name) + 1})"
        for i, name in enumerate(names)
    ]


def make_pymoo_problem(study: Study, criteria: Criteria | Mapping[str, str]) -> Any:
    """
    Makes a pymoo problem of a study's design space: its variables are the
    unit hypercube that encodes the design space (a dimension per decision,
    whose levels divide it equally), its objectives are minimized (maximized
    objectives are negated), and its constraint (`G <= 0`) is violated by
    designs that violate a constraint of the design space (by the amount of
    the violation, see `DesignSpace`) or that have no feasible replicate (by
    the share of their infeasible evaluations, see
    `Evaluations.infeasibility`). Each requirement of the
    criteria is a further constraint, whose value is its violation (see
    `Requirement.violation`).
    The designs of a population are evaluated together.

    The problem keeps its evaluations (`evaluations`, see `Evaluations`),
    whose `results` are those of every evaluation (see `Study.run`), and an
    elimination of duplicates by design (`duplicates`, a pymoo
    `DuplicateElimination` for the `eliminate_duplicates` of an algorithm):
    many points of the unit hypercube decode to the same design (those in the
    same level of each variable, or that differ only in the variables of the
    branches the design does not choose), which pymoo's elimination of equal
    points does not see. A proposal is a duplicate if it decodes to the same
    design (by key) as another proposal, a member of the population, or a
    design evaluated before; proposals that violate a constraint of the
    design space are not duplicates.

    Args:
        study (Study): The study that evaluates designs.
        criteria (Criteria | Mapping[str, str]): The criteria, without `by`
            (see `Criteria`), or the objectives (metrics whose means they
            are) and their directions (`"min"` or `"max"`).

    Returns:
        pymoo.core.problem.Problem: the problem
    """
    # pylint: disable=import-outside-toplevel
    try:
        from pymoo.core.duplicate import DuplicateElimination
        from pymoo.core.problem import Problem
    except ImportError as error:
        raise ImportError(
            "make_pymoo_problem requires pymoo (pip install pymoo)"
        ) from error

    evaluations = Evaluations(study, criteria)
    n_requirements = len(evaluations.criteria.requirements)
    signs = _signs(evaluations.criteria)
    space = study.space
    decode = space._decode  # pylint: disable=protected-access
    n_var = len(space._unit_paths())  # pylint: disable=protected-access

    def key_of(individual: Any) -> str | None:
        """The key of the design of a proposal, or `None` if it violates a
        constraint of the design space."""
        try:
            return decode(individual.X).key
        except Infeasible:
            return None

    class DesignDuplicates(DuplicateElimination):
        """
        An elimination of the proposals of designs proposed or evaluated
        before (see `make_pymoo_problem`).
        """

        def _do(self, pop: Any, other: Any, is_duplicate: np.ndarray) -> np.ndarray:
            seen = set(evaluations.designs)
            if other is not None:
                seen.update(key for key in map(key_of, other) if key is not None)
            for i, individual in enumerate(pop):
                key = key_of(individual)
                if key is None:
                    continue
                if key in seen:
                    is_duplicate[i] = True
                else:
                    seen.add(key)
            return is_duplicate

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
            self.duplicates = DesignDuplicates()

        def _evaluate(self, x, out, *args, **kwargs):
            designs: list[Design | None] = []
            violation = np.ones(len(np.atleast_2d(x)))
            for i, row in enumerate(np.atleast_2d(x)):
                try:
                    designs.append(decode(row))
                except Infeasible as error:
                    designs.append(None)
                    violation[i] = max(error.violation, 1e-9)
            feasible = [design for design in designs if design is not None]
            values = iter(evaluations.evaluate(feasible))
            costs = np.zeros((len(designs), len(signs)))
            unmet = np.zeros((len(designs), n_requirements))
            for i, design in enumerate(designs):
                value = None if design is None else next(values)
                if design is not None and value is None:
                    violation[i] = max(evaluations.infeasibility[design.key], 1e-9)
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


def _eliminate_duplicates(algorithm: Any, duplicates: Any) -> None:
    """
    Makes a pymoo algorithm eliminate duplicates by design: pymoo's genetic
    algorithms eliminate duplicates in their initialization and mating.
    """
    parts = [
        algorithm,
        getattr(algorithm, "initialization", None),
        getattr(algorithm, "mating", None),
    ]
    if any(part is None or not hasattr(part, "eliminate_duplicates") for part in parts):
        raise ValueError(
            "distinct searches need an algorithm that eliminates duplicates in its "
            + "initialization and mating (as pymoo's genetic algorithms do): "
            + "pass distinct=False for others"
        )
    for part in parts:
        part.eliminate_duplicates = duplicates


def _designs_termination(
    evaluations: Evaluations, max_designs: int, termination: Any
) -> Any:
    """
    Gets a pymoo termination of a search once it has evaluated `max_designs`
    designs, or by another termination (if any), whichever is first.
    """
    # pylint: disable=import-outside-toplevel
    from pymoo.core.termination import Termination
    from pymoo.util.misc import termination_from_tuple

    class DesignsTermination(Termination):
        """
        A termination by the number of designs evaluated (see
        `optimize_pymoo`), or by another termination.
        """

        def __init__(self, other: Any):
            super().__init__()
            self.other = other

        def _update(self, algorithm: Any) -> float:
            progress = len(evaluations.values) / max_designs
            if self.other is not None:
                progress = max(progress, self.other.update(algorithm))
            return progress

    return DesignsTermination(
        None if termination is None else termination_from_tuple(termination)
    )


def optimize_pymoo(
    study: Study,
    criteria: Criteria | Mapping[str, str],
    *,
    algorithm: Any = None,
    termination: Any = ("n_gen", 10),
    seed: int | None = 0,
    distinct: bool = True,
    max_designs: int | None = None,
    **kwargs: Any,
) -> Search:
    """
    Searches a study's design space with a pymoo algorithm (see
    `make_pymoo_problem`). The algorithm eliminates proposals of designs
    proposed or evaluated before (its initialization and mating, see
    `make_pymoo_problem`'s `duplicates`), unless `distinct` is `False`, so
    that each generation proposes as many new designs as it can: as many
    generations then evaluate more designs, and find more of the front. The
    algorithm given is not changed: a distinct search changes a copy.

    The generations of a search (of the `termination`) are a poor measure of
    its cost if they propose designs evaluated before (without `distinct`),
    whose objectives are known without evaluating them again: `max_designs`
    limits the designs evaluated instead, as of `optimize_optuna` (the search
    stops after the generation in which it has evaluated as many, or by the
    `termination`, whichever is first). As a generation evaluates its new
    designs together, the last may evaluate more designs than are left.

    Args:
        study (Study): The study that evaluates designs.
        criteria (Criteria | Mapping[str, str]): The criteria, without `by`
            (see `Criteria`), or the objectives (metrics whose means they
            are) and their directions (`"min"` or `"max"`).
        algorithm (pymoo.core.algorithm.Algorithm | None): The algorithm, or
            `None` for NSGA-II with populations of 20.
        termination (Any): The termination criterion of `pymoo.optimize.minimize`.
        seed (int | None): The seed of the algorithm (`None` for one drawn
            anew, as pymoo does).
        distinct (bool): `True` to eliminate the proposals of designs proposed
            or evaluated before, `False` for the algorithm's own elimination
            of duplicates (of equal points).
        max_designs (int | None): The most designs to evaluate (but those of
            the last generation), if limited.
        **kwargs (Any): Other arguments of `pymoo.optimize.minimize`.

    Returns:
        Search: the results and designs of the search, and the result of pymoo

    Raises:
        ValueError: if the search is distinct but the algorithm does not
            eliminate duplicates as pymoo's genetic algorithms do, or if
            `max_designs` is less than 1.
    """
    # pylint: disable=import-outside-toplevel
    try:
        from pymoo.optimize import minimize
    except ImportError as error:
        raise ImportError(
            "optimize_pymoo requires pymoo (pip install pymoo)"
        ) from error

    if max_designs is not None and max_designs < 1:
        raise ValueError("max_designs must be at least 1")
    if algorithm is None:
        from pymoo.algorithms.moo.nsga2 import NSGA2

        algorithm = NSGA2(pop_size=20)
    problem = make_pymoo_problem(study, criteria)
    if distinct:
        # a copy, so that the caller's algorithm keeps its own elimination
        algorithm = copy.deepcopy(algorithm)
        _eliminate_duplicates(algorithm, problem.duplicates)
    if max_designs is not None:
        termination = _designs_termination(
            problem.evaluations, max_designs, termination
        )
    result = minimize(problem, algorithm, termination, seed=seed, **kwargs)
    evaluations = problem.evaluations
    return evaluations.search(result)

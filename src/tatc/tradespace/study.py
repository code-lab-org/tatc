"""
Trade studies: the evaluation of designs from a design space, with
replicates that draw their uncertain parameters, in sequence or in parallel.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import hashlib
import inspect
import json
import math
import os
import time
import warnings
from collections.abc import Callable, Iterable, Mapping, Sequence
from concurrent.futures import Executor, as_completed
from typing import Any, NamedTuple

import pandas as pd
from pydantic import ValidationError

from .fingerprint import (
    _fingerprint_of,
    _without_code,
    compare_definitions,
    get_definition,
)
from .nodes import Infeasible
from .space import (
    Design,
    DesignSpace,
    _check_scenario,
    _check_uncertain,
    _describe_scenario,
    _designs,
)
from .store import ResultStore

_STATUS_COLUMNS = ["design", "replicate", "feasible", "error", "seconds"]
"""Columns of the results of a study that are neither design values nor metrics."""

_UNSTORED = "_unstored"
"""A key of the records of evaluations that are not stored: those that failed
(with `errors="record"`) and designs that failed validation, which may be
mistakes in a study's functions that its fingerprint cannot see."""


_WORKER_TIMEOUT = 120
"""The seconds that an idle worker process waits for evaluations before it
exits: long enough to outlast the pauses of a study between its runs (such
as the bootstrap of each round of `refine`), which would otherwise start new
workers (with a warning of loky's, if evaluations arrive as a worker exits)."""


def _get_executor(n_jobs: int) -> Executor:
    """
    Gets a reusable executor of `n_jobs` worker processes (from joblib's
    loky, which serializes functions defined in notebooks and lambdas).
    """
    try:
        # pylint: disable=import-outside-toplevel
        from joblib.externals.loky import get_reusable_executor
    except ImportError as error:
        raise ImportError(
            "parallel trade studies require joblib (pip install joblib), "
            + "or pass an executor"
        ) from error
    # one thread per worker for the numerical libraries, whose thread pools
    # (one thread per processor in each worker) otherwise exhaust memory
    threads = {
        name: "1"
        for name in (
            "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS",
            "NUMBA_NUM_THREADS",
        )
    }
    return get_reusable_executor(
        max_workers=n_jobs, env=threads, timeout=_WORKER_TIMEOUT
    )


class Study:
    """
    A trade study: designs from a design space, each built into an
    architecture (or any object to evaluate) and evaluated in replicates::

        study = Study(space, build, evaluate, replicates=5, n_jobs=-1, store=store)
        results = study.run()  # every design, a row per replicate
        more = study.replace(replicates=20).run(designs)

    Each replicate realizes a design (see `Design.realize`), builds it with
    `build`, and evaluates the result with `evaluate`, which returns its
    outcomes: a mapping (or `pandas.Series`) of metrics, or a
    `pandas.DataFrame` of several records of metrics (such as one per
    sensor), whose named index is kept as columns. An `evaluate` function
    with a parameter named `design` is also given the realized design (with
    the uncertain parameters of the replicate), as for an evaluation that
    depends on values of the design other than those it builds::

        def build(d):
            return build_constellation(d.layout, payload, epoch=START)

        def evaluate(satellites, design):
            deployment = SCHEDULE.plan(design.layout)
            ...

    A design is recorded as infeasible, with the error, rather than raised,
    if `build` or `evaluate` rejects it by raising `Infeasible`, or if
    `build` fails the validation of TAT-C's schemas (a pydantic
    `ValidationError`, such as of an orbit too high to be sun-synchronous).
    Other errors (including other `ValueError`s, which are as likely
    mistakes in the functions as properties of designs) are raised, which
    stops the study, unless `errors` is `"record"`: then the evaluation is
    recorded as infeasible, with the type of the error and its message.
    Evaluations that fail and designs that fail validation are not kept in
    the study's store, so that a later run evaluates them again (as after
    fixing a mistake in code that the study's fingerprint cannot see, see
    `fingerprint`). A warning is given if no design of a run (of several)
    is feasible.

    The metrics must not be named as the columns of the results that
    identify an evaluation (`design`, `replicate`, `feasible`, `error`, and
    `seconds`), and a metric named as a value of the design (see
    `Design.flat`) must equal it (restating it, as a design variable that is
    also an objective): metrics that would replace these columns raise a
    `ValueError`, whatever `errors` is.

    A study may evaluate designs in `scenarios`: conditions under which
    designs are compared rather than averaged over (as uncertain parameters
    are), such as the requirements of a mission or the field it must sample.
    A scenario gives fixed values or uncertain parameters of the design space
    other values (by path, see `Design.realize`), which the functions read
    from the design as any other value (so that they need not change), and
    each design is evaluated in each scenario in each replicate (with the
    same uncertain parameters in every scenario), recorded with the name of
    its scenario (`scenario`, so that the design space may not have a value
    of that name). Judge designs in each scenario with criteria `by` the
    scenario (see `Criteria`), and across them with `robustness`. Scenarios
    are not part of the study's fingerprint: the evaluations of each are
    stored by its values (as a variant of the replicates, see
    `ResultStore`), so that a scenario added later is evaluated alone, and
    one without values shares the evaluations of the study without
    scenarios.

    Evaluations run in parallel with `n_jobs` worker processes (or an
    `executor` with the interface of `concurrent.futures.Executor`). They are
    sent to the workers in chunks of consecutive evaluations of the same
    replicates, so that a cache of each worker (see `SatelliteCache`) is
    shared by the designs of a replicate, which share their uncertain
    parameters (with common random numbers). The functions, and the objects
    they refer to, are copied to each worker for each chunk, so keep large
    inputs out of them (for example, load them from files in the workers) if
    copying them takes long. Each evaluation is added to the study's store
    (if any) as it completes, by the worker that evaluates it, so that a run
    that is interrupted loses only the evaluations in progress; chunks that
    have not begun are then cancelled.

    Args:
        space (DesignSpace): The design space.
        build (Callable): A function of a realized design that builds the
            object to evaluate (for example, a list of satellites).
        evaluate (Callable): A function of a built object that evaluates it
            (and of the realized design, if it has a parameter named
            `design`).
        replicates (int): The number of replicates of each design.
        seed (int): The seed of the uncertain parameters.
        common_random_numbers (bool): `True` to draw the same uncertain
            parameters in the same replicate of every design (see
            `Design.realize`).
        n_jobs (int): The number of worker processes (1 to evaluate in this
            process, -1 for one per processor), each limited to one thread of
            the numerical libraries (OpenMP, OpenBLAS, MKL, and numba).
        executor (concurrent.futures.Executor | None): An executor to use
            instead of worker processes.
        chunks_per_job (int): The number of chunks of evaluations per worker.
        errors (str): `"raise"` to raise the errors of evaluations (other than
            those of designs that `build` rejects), or `"record"` to record
            them as infeasible evaluations.
        store (ResultStore | None): A store of the records of evaluations, from
            which evaluations stored before (by a study of the same
            fingerprint, see `fingerprint`) are taken rather than evaluated
            again, and to which new ones are added as they complete (see
            `ResultStore`), by the workers that evaluate them (so its file
            must be reachable from them). Studies of other definitions may
            share a store (and its namespace): their fingerprints keep their
            evaluations apart. A warning is given if the store has
            evaluations of the same study with other code (of its functions,
            of TAT-C, or of Python), which are not reused.
        scenarios (Mapping[str, Mapping[str, Any]] | None): The scenarios
            in which to evaluate designs, if any: the values of fixed values
            and uncertain parameters of each (by path), by name.
    """

    def __init__(
        self,
        space: DesignSpace,
        build: Callable[[Design], Any],
        evaluate: Callable[[Any], Any],
        replicates: int = 1,
        seed: int = 0,
        common_random_numbers: bool = True,
        n_jobs: int = 1,
        executor: Executor | None = None,
        chunks_per_job: int = 2,
        store: ResultStore | None = None,
        errors: str = "raise",
        scenarios: Mapping[str, Mapping[str, Any]] | None = None,
    ):
        if errors not in ("raise", "record"):
            raise ValueError(f"errors must be 'raise' or 'record', not {errors!r}")
        if replicates < 1:
            raise ValueError("replicates must be at least 1")
        if n_jobs == 0 or n_jobs < -1:
            raise ValueError("n_jobs must be positive or -1")
        if chunks_per_job < 1:
            raise ValueError("chunks_per_job must be at least 1")
        self.space, self.build, self.evaluate = space, build, evaluate
        self.replicates, self.seed = replicates, seed
        self.common_random_numbers = common_random_numbers
        self.n_jobs, self.executor = n_jobs, executor
        self.chunks_per_job = chunks_per_job
        self.store = store
        self.errors = errors
        self.scenarios = _check_scenarios(space, scenarios)

    def __getstate__(self) -> dict:
        # an executor is not sent to workers (which add evaluations to the store)
        state = self.__dict__.copy()
        state["executor"] = None
        return state

    def _evaluate(
        self,
        design: Design,
        replicate: int,
        shift: Mapping[str, int] | None = None,
        case: _Case | None = None,
    ) -> list[dict]:
        """
        Evaluates a replicate of a design (with parameters drawn as in other
        replicates, if shifted, and in a scenario, if any), returning its
        records.
        """
        realized = design.realize(
            replicate,
            self.seed,
            self.common_random_numbers,
            shift,
            None if case is None else case.values,
        )
        base: dict = {"design": design.key, "replicate": replicate}
        if case is not None and case.name is not None:
            base["scenario"] = case.name
        base.update(realized.flat())
        began = time.perf_counter()
        try:
            built = self.build(realized)
        except Infeasible as error:
            return [{**base, "feasible": False, "error": str(error)}]
        except ValidationError as error:
            return [{**base, "feasible": False, "error": str(error), _UNSTORED: True}]
        except Exception as error:  # pylint: disable=broad-exception-caught
            if self.errors == "raise":
                raise
            return [self._failure(base, error)]
        try:
            if _takes_design(self.evaluate):
                outcome = self.evaluate(built, design=realized)
            else:
                outcome = self.evaluate(built)
        except Infeasible as error:
            return [{**base, "feasible": False, "error": str(error)}]
        except Exception as error:  # pylint: disable=broad-exception-caught
            if self.errors == "raise":
                raise
            return [self._failure(base, error)]
        base.update(feasible=True, seconds=time.perf_counter() - began)
        if isinstance(outcome, pd.DataFrame):
            named = any(name is not None for name in outcome.index.names)
            frame = outcome.reset_index() if named else outcome
            records = frame.to_dict("records")
        else:
            records = [dict(outcome)]
        for record in records:
            _check_metrics(record, base)
        return [{**base, **record} for record in records]

    @staticmethod
    def _failure(base: dict, error: Exception) -> dict:
        """
        Gets the record of an evaluation that failed.
        """
        message = f"{type(error).__name__}: {error}"
        return {**base, "feasible": False, "error": message, _UNSTORED: True}

    def _evaluate_and_store(
        self,
        design: Design,
        replicate: int,
        context: str,
        shift: Mapping[str, int] | None = None,
        case: _Case | None = None,
    ) -> list[dict]:
        """
        Evaluates a replicate of a design (shifted, if given, and in a
        scenario, if any), adding its records to the store (if any, and if
        they are stored, see `_UNSTORED`) in a context.
        """
        records = self._evaluate(design, replicate, shift, case)
        if self.store is not None and not _unstored(records):
            variant = _variant_of(shift) if case is None else case.variant
            if case is not None and case.name is not None:
                # the evaluations of a scenario are stored by its values, not
                # its name (which the results give, see `_in_scenario`)
                kept = [
                    {k: v for k, v in r.items() if k != "scenario"} for r in records
                ]
            else:
                kept = records
            self.store.put([(design.key, replicate, kept)], context, variant)
        return records

    def _run_chunk(
        self,
        tasks: list[tuple[int, Design, int, int]],
        context: str,
        shift: Mapping[str, int] | None,
        cases: Sequence[_Case],
    ) -> list[list[dict]]:
        """
        Evaluates a chunk of tasks (index, design, replicate, and case) in
        turn (shifted, if given), adding each to the store (if any) in a
        context.
        """
        return [
            self._evaluate_and_store(design, replicate, context, shift, cases[case])
            for _, design, replicate, case in tasks
        ]

    def _cases(self, shift: Mapping[str, int]) -> list[_Case]:
        """
        Gets the cases of a run: its scenarios (or none), and the variants of
        their evaluations in a store.
        """
        if not self.scenarios:
            return [_Case(None, {}, _variant_of(shift))]
        return [
            _Case(name, values, _variant_of(shift, _describe_scenario(values)))
            for name, values in self.scenarios.items()
        ]

    def replace(self, **changes: Any) -> Study:
        """
        Gets a study like this one with some of its arguments changed (by
        name, as of its class), such as a variant with more replicates, in
        scenarios, or of another design space::

            confirming = study.replace(replicates=20)
            policies = study.replace(space=space.extend(policy=POLICY))

        Args:
            **changes (Any): The arguments to change.

        Returns:
            Study: the study, of the same class

        Raises:
            TypeError: if a change names no argument of the study's class.
        """
        arguments = self._arguments()
        unknown = sorted(set(changes) - set(arguments))
        if unknown:
            raise TypeError(f"{type(self).__name__} has no arguments {unknown}")
        return type(self)(**{**arguments, **changes})

    def _arguments(self) -> dict[str, Any]:
        """
        Gets the arguments of the study, by name, from which `replace` makes
        another (a subclass with other arguments gives its own).
        """
        if type(self).__init__ is not Study.__init__:
            raise TypeError(
                f"{type(self).__name__} does not give its arguments to replace"
            )
        return {
            "space": self.space,
            "build": self.build,
            "evaluate": self.evaluate,
            **self._settings(),
        }

    def _settings(self) -> dict[str, Any]:
        """Gets the settings of the study (the arguments of every study but
        its design space and functions), by name."""
        return {
            "replicates": self.replicates,
            "seed": self.seed,
            "common_random_numbers": self.common_random_numbers,
            "n_jobs": self.n_jobs,
            "executor": self.executor,
            "chunks_per_job": self.chunks_per_job,
            "store": self.store,
            "errors": self.errors,
            "scenarios": self.scenarios,
        }

    @property
    def definition(self) -> dict:
        """dict: The definition of the study (see `get_definition`), which
        its fingerprint hashes: its seed, its use of common random numbers,
        its design space (including fixed values and the distributions of
        uncertain parameters), its build and evaluation functions (their
        code, the values they close over, and the plain values and functions
        of their module they refer to), and the source of the modules of
        TAT-C that these, and studies, depend on."""
        return get_definition(
            study=Study,
            seed=self.seed,
            common_random_numbers=self.common_random_numbers,
            space=self.space,
            build=self.build,
            evaluate=self.evaluate,
        )

    @property
    def fingerprint(self) -> str:
        """str: A key of the definition of the study (see `definition`), which
        changes if it does. It keys the records of a store (see
        `ResultStore`); data that it cannot describe (such as data frames)
        must be identified by the store's namespace."""
        return _fingerprint_of(self.definition)

    @property
    def workers(self) -> int:
        """int: The number of workers that evaluate designs (1 in this process)."""
        if self.executor is not None:
            return max(getattr(self.executor, "_max_workers", 1) or 1, 1)
        if self.n_jobs == -1:
            return os.cpu_count() or 1
        return self.n_jobs

    def run(
        self,
        designs: Iterable[Design] | pd.DataFrame | None = None,
        progress: Callable[[int, int], None] | None = None,
        replicates: Iterable[int] | None = None,
        shift: Mapping[str, int] | None = None,
    ) -> pd.DataFrame:
        """
        Evaluates designs in every replicate (or in given replicates, such as
        additional ones, see `refine`). Repeated designs (with the same key)
        are evaluated once.

        With a `shift`, each replicate draws the named uncertain parameters
        as in another replicate (see `Design.realize`): for example, with a
        shift of `{"points_seed": 10}`, replicate 3 draws `points_seed` as
        replicate 13 does, and its other parameters as replicate 3 does.
        Such mixed replicates are kept apart from plain ones in the store
        (by their variant, see `ResultStore`), and are recorded with the
        replicate they shift (and the values they draw, see `Design.flat`).

        With scenarios (see `Study`), each design is evaluated in each
        scenario in each replicate.

        Args:
            designs (Iterable[Design] | pandas.DataFrame | None): The designs,
                a data frame of their values (such as a summary that keeps
                them, see `DesignSpace.designs_of`), or `None` for every
                design of the design space (see `DesignSpace.enumerate`).
            progress (Callable | None): A function called as evaluations
                complete with the number completed and the number of
                evaluations (not counting those taken from the store).
            replicates (Iterable[int] | None): The replicates (their
                indices, which seed their uncertain parameters), or `None` for
                the study's (`range(replicates)`).
            shift (Mapping[str, int] | None): The offsets of the replicates as
                in which the named uncertain parameters are drawn (by their
                paths, see `DesignSpace.uncertain_paths`), if any.

        Returns:
            pandas.DataFrame: The outcomes, a record per evaluation (or per
                record of an evaluation), in the order of the designs (and
                scenarios) and replicates, with the design's key (`design`),
                the `replicate` (and the `scenario`, with scenarios), the
                values of the design by path (see `Design.flat`), whether it
                is `feasible` (and the `error` if not), the time to build and
                evaluate it (`seconds`), and the metrics.
        """
        shift = {name: int(offset) for name, offset in (shift or {}).items() if offset}
        _check_uncertain(self.space, shift)
        cases = self._cases(shift)
        unique = (
            self.space.enumerate() if designs is None else _designs(self.space, designs)
        )
        # tasks in the order of replicates (and scenarios), whose designs share
        # uncertain values
        tasks = [
            (i, design, replicate, case)
            for replicate in (
                range(self.replicates)
                if replicates is None
                else sorted(set(replicates))
            )
            for case in range(len(cases))
            for i, design in enumerate(unique)
        ]
        results: list[tuple[int, int, int, list[dict]]] = []
        context = ""
        if self.store is not None:
            # evaluations stored before (by an interrupted or earlier run)
            definition = self.definition
            context = _fingerprint_of(definition)
            stored = {}
            for k, case in enumerate(cases):
                found = self.store.get(
                    ((d.key, r) for _, d, r, c in tasks if c == k),
                    context,
                    case.variant,
                )
                stored.update(
                    {(key, r, k): records for (key, r), records in found.items()}
                )
            if tasks and not stored:
                _warn_other_definitions(self.store, context, definition)
            self.store.define(context, definition)
            results.extend(
                (i, case, replicate, stored[design.key, replicate, case])
                for i, design, replicate, case in tasks
                if (design.key, replicate, case) in stored
            )
            tasks = [t for t in tasks if (t[1].key, t[2], t[3]) not in stored]
        total = len(tasks)
        if self.executor is None and self.n_jobs == 1:
            for done, (i, design, replicate, case) in enumerate(tasks, 1):
                records = self._evaluate_and_store(
                    design, replicate, context, shift, cases[case]
                )
                results.append((i, case, replicate, records))
                if progress is not None:
                    progress(done, total)
        elif tasks:
            executor = self.executor or _get_executor(self.workers)
            n_chunks = min(total, self.workers * self.chunks_per_job)
            size = math.ceil(total / max(n_chunks, 1))
            chunks = [tasks[k : k + size] for k in range(0, total, size)]
            futures = {
                executor.submit(self._run_chunk, chunk, context, shift, cases): chunk
                for chunk in chunks
            }
            done = 0
            try:
                for future in as_completed(futures):
                    chunk = futures[future]
                    evaluated = future.result()
                    for (i, _, replicate, case), records in zip(chunk, evaluated):
                        results.append((i, case, replicate, records))
                    done += len(chunk)
                    if progress is not None:
                        progress(done, total)
            except BaseException:
                # interrupted (or an error): chunks that have not begun are
                # cancelled, and those in progress add their evaluations to
                # the store as they complete
                for future in futures:
                    future.cancel()
                raise
        results.sort(key=lambda result: result[:3])
        frame = pd.DataFrame.from_records(
            [
                _in_scenario(record, cases[case].name)
                for _, case, _, records in results
                for record in records
            ]
        )
        _check_feasible(frame)
        return frame


def _in_scenario(record: dict, name: str | None) -> dict:
    """
    Gets a record of an evaluation for the results of a run: without the
    marker of unstored records (see `_UNSTORED`), and with the name of its
    scenario (if any) after its replicate (as records of a scenario are
    stored without it).
    """
    items = [(k, v) for k, v in record.items() if k != _UNSTORED]
    if name is None:
        return dict(items)
    items = [(k, v) for k, v in items if k != "scenario"]
    return dict([*items[:2], ("scenario", name), *items[2:]])


class _Case(NamedTuple):
    """
    A scenario of a run (or none, without scenarios): its name, its values,
    and the variant of its evaluations in a store (see `_variant_of`).
    """

    name: str | None
    values: Mapping[str, Any]
    variant: str


def _check_scenarios(
    space: DesignSpace, scenarios: Mapping[str, Mapping[str, Any]] | None
) -> dict[str, dict[str, Any]] | None:
    """
    Checks the scenarios of a study (see `Study`), returning them as plain
    mappings (or `None` without scenarios).
    """
    if not scenarios:
        return None
    if not isinstance(scenarios, Mapping):
        raise TypeError("scenarios must be a mapping of names to values")
    if "scenario" in space.children:
        raise ValueError(
            "a design space with a value named 'scenario' cannot have scenarios, "
            + "whose names the results give as `scenario`"
        )
    checked = {}
    for name, values in scenarios.items():
        if not isinstance(name, str):
            raise TypeError(f"the name of a scenario must be a string, not {name!r}")
        if not isinstance(values, Mapping):
            raise TypeError(f"the values of scenario {name!r} must be a mapping")
        _check_scenario(space, values)
        checked[name] = dict(values)
    return checked


def _variant_of(shift: Mapping[str, int] | None, scenario: str = "") -> str:
    """
    Gets the variant of replicates in a store (see `ResultStore`): `""` for
    plain replicates, or (as JSON) a hash of the description of the values of
    a scenario (see `_describe_scenario`) and the shift of each parameter.
    """
    shift = {name: int(offset) for name, offset in (shift or {}).items() if offset}
    parts: dict[str, Any] = {}
    if scenario:
        parts["scenario"] = hashlib.sha1(scenario.encode()).hexdigest()[:16]
    if shift:
        parts["shift"] = shift
    return json.dumps(parts, sort_keys=True) if parts else ""


def _takes_design(function: Callable) -> bool:
    """
    Whether a function has a parameter named `design`, which a study gives
    the realized design (see `Study`).
    """
    try:
        parameter = inspect.signature(function).parameters.get("design")
    except (TypeError, ValueError):  # such as of some built-in functions
        return False
    return parameter is not None and parameter.kind in (
        inspect.Parameter.POSITIONAL_OR_KEYWORD,
        inspect.Parameter.KEYWORD_ONLY,
    )


def _warn_other_definitions(store: ResultStore, context: str, definition: dict) -> None:
    """
    Warns if a store has evaluations of its namespace by the same study with
    other code (its functions', TAT-C's, or Python's: the same values but for
    the code of functions, see `_without_code`), which a study does not
    reuse, with the differences (see `compare_definitions`). Studies of other
    values (such as other settings) share a namespace without a warning, as
    their fingerprints keep their evaluations apart.
    """
    contexts = store.contexts()
    others = [other for other in contexts if other != context]
    if not others:
        return
    values = _without_code(definition.get("values", {}))
    for other in reversed(others):
        before = store.definition(other)
        if before is None or _without_code(before.get("values", {})) != values:
            continue
        differences = compare_definitions(before, definition)
        warnings.warn(
            f"the store has {contexts[other]} evaluations of namespace "
            + f"{store.namespace!r} by this study with other code, which are not "
            + "reused (see Study.fingerprint, and ResultStore.prune to remove "
            + "them); it differs in: "
            + "; ".join(differences),
            stacklevel=3,
        )
        return


def _unstored(records: list[dict]) -> bool:
    """
    Whether the records are of an evaluation that is not stored (see
    `_UNSTORED`).
    """
    return any(record.get(_UNSTORED, False) for record in records)


def _same(value: Any, other: Any) -> bool:
    """
    Whether two values are equal (`NaN` equal to `NaN`).
    """
    try:
        # pylint: disable-next=comparison-with-itself
        return bool(value == other) or bool(value != value and other != other)
    except (TypeError, ValueError):  # such as arrays
        return False


def _check_metrics(metrics: dict, base: dict) -> None:
    """
    Checks that the metrics of a record do not replace the columns of the
    results of a study: they must not be named as its status columns, and
    a metric named as a value of the design (a key of `base`) must equal it.
    """
    clashes = [
        str(name)
        for name, value in metrics.items()
        if name in _STATUS_COLUMNS or (name in base and not _same(value, base[name]))
    ]
    if clashes:
        raise ValueError(
            f"metrics named as columns of the results of a study: {sorted(clashes)}; "
            + f"rename them (the results have the columns {_STATUS_COLUMNS}, and "
            + "the values of each design, see Design.flat, which a metric of the "
            + "same name must equal)"
        )


def _check_feasible(results: pd.DataFrame) -> None:
    """
    Warns if none of several designs of the results of a study is feasible,
    which is as likely a mistake in the study's functions as a property of
    its designs.
    """
    if len(results) == 0 or results["feasible"].any():
        return
    designs = results["design"].nunique()
    if designs < 2:
        return
    errors = results["error"].value_counts()
    warnings.warn(
        f"none of the {designs} designs is feasible; the most common error "
        + f"({errors.iloc[0]} of {len(results)} evaluations): {errors.index[0]}",
        stacklevel=3,
    )

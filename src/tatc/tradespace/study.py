"""
Trade studies: the evaluation of designs from a design space, with
replicates that draw their uncertain parameters, in sequence or in parallel.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import math
import os
import time
import warnings
from collections.abc import Callable, Iterable, Sequence
from concurrent.futures import Executor, as_completed
from typing import Any

import numpy as np
import pandas as pd

from .fingerprint import get_fingerprint
from .space import Design, DesignSpace
from .store import ResultStore

_STATUS_COLUMNS = ["design", "replicate", "feasible", "error", "seconds"]
"""Columns of the results of a study that are neither design values nor metrics."""

_FAILED = "_failed"
"""A key of the records of evaluations that failed (with `errors="record"`), which
are not stored."""


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
    return get_reusable_executor(max_workers=n_jobs, env=threads)


class Study:
    """
    A trade study: designs from a design space, each built into an
    architecture (or any object to evaluate) and evaluated in replicates.

    Each replicate realizes a design (see `Design.realize`), builds it with
    `build`, and evaluates the result with `evaluate`, which returns its
    outcomes: a mapping (or `pandas.Series`) of metrics, or a
    `pandas.DataFrame` of several records of metrics (such as one per
    sensor), whose named index is kept as columns.

    A design that `build` rejects with a `ValueError` (which includes the
    validation errors of TAT-C's schemas) is recorded as infeasible, with
    the error, rather than raised. Other errors of `build` or `evaluate`
    are raised, which stops the study, unless `errors` is `"record"`: then
    the evaluation is recorded as infeasible, with the type of the error and
    its message, and is not kept in the study's store, so that a later run
    evaluates it again.

    Evaluations run in parallel with `n_jobs` worker processes (or an
    `executor` with the interface of `concurrent.futures.Executor`). They are
    sent to the workers in chunks of consecutive evaluations of the same
    replicates, so that a cache of each worker (see `SatelliteCache`) is
    shared by the designs of a replicate, which share their uncertain
    parameters (with common random numbers). The functions, and the objects
    they refer to, are copied to each worker for each chunk, so keep large
    inputs out of them (for example, load them from files in the workers) if
    copying them takes long.

    Args:
        space (DesignSpace): The design space.
        build (Callable): A function of a realized design that builds the
            object to evaluate (for example, a list of satellites).
        evaluate (Callable): A function of a built object that evaluates it.
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
            `ResultStore`).
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

    def __getstate__(self) -> dict:
        # an executor and a store are not sent to workers
        state = self.__dict__.copy()
        state["executor"] = state["store"] = None
        return state

    def _evaluate(self, design: Design, replicate: int) -> list[dict]:
        """
        Evaluates a replicate of a design, returning its records.
        """
        realized = design.realize(replicate, self.seed, self.common_random_numbers)
        base = {"design": design.key, "replicate": replicate, **realized.flat()}
        began = time.perf_counter()
        try:
            built = self.build(realized)
        except ValueError as error:
            return [{**base, "feasible": False, "error": str(error)}]
        except Exception as error:  # pylint: disable=broad-exception-caught
            if self.errors == "raise":
                raise
            return [self._failure(base, error)]
        try:
            outcome = self.evaluate(built)
        except Exception as error:  # pylint: disable=broad-exception-caught
            if self.errors == "raise":
                raise
            return [self._failure(base, error)]
        base.update(feasible=True, seconds=time.perf_counter() - began)
        if isinstance(outcome, pd.DataFrame):
            named = any(name is not None for name in outcome.index.names)
            frame = outcome.reset_index() if named else outcome
            return [{**base, **record} for record in frame.to_dict("records")]
        return [{**base, **dict(outcome)}]

    @staticmethod
    def _failure(base: dict, error: Exception) -> dict:
        """
        Gets the record of an evaluation that failed.
        """
        message = f"{type(error).__name__}: {error}"
        return {**base, "feasible": False, "error": message, _FAILED: True}

    def _run_chunk(self, tasks: list[tuple[int, Design, int]]) -> list[list[dict]]:
        """
        Evaluates a chunk of tasks (index, design, replicate) in turn.
        """
        return [self._evaluate(design, replicate) for _, design, replicate in tasks]

    @property
    def fingerprint(self) -> str:
        """str: A key of the definition of the study, which changes with its
        seed, its use of common random numbers, its design space (including
        fixed values and the distributions of uncertain parameters), its
        build and evaluation functions (their code, the values they close
        over, and the plain values and functions of their module they refer
        to), and the version and source of TAT-C (see `get_source_hash`). It
        keys the records of a store (see
        `ResultStore`); data that it cannot describe (such as data frames)
        must be identified by the store's namespace."""
        return get_fingerprint(
            self.seed,
            self.common_random_numbers,
            self.space,
            self.build,
            self.evaluate,
        )

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
        designs: Iterable[Design] | None = None,
        progress: Callable[[int, int], None] | None = None,
        replicates: Iterable[int] | None = None,
    ) -> pd.DataFrame:
        """
        Evaluates designs in every replicate (or in given replicates, such as
        additional ones, see `refine`). Repeated designs (with the same key)
        are evaluated once.

        Args:
            designs (Iterable[Design] | None): The designs, or `None` for
                every design of the design space (see
                `DesignSpace.enumerate`).
            progress (Callable | None): A function called as evaluations
                complete with the number completed and the number of
                evaluations (not counting those taken from the store).
            replicates (Iterable[int] | None): The replicates (their
                indices, which seed their uncertain parameters), or `None` for
                the study's (`range(replicates)`).

        Returns:
            pandas.DataFrame: The outcomes, a record per evaluation (or per
                record of an evaluation), in the order of the designs and
                replicates, with the design's key (`design`), the
                `replicate`, the values of the design by path (see
                `Design.flat`), whether it is `feasible` (and the `error` if
                not), the time to build and evaluate it (`seconds`), and the
                metrics.
        """
        designs = self.space.enumerate() if designs is None else designs
        unique = list({design.key: design for design in designs}.values())
        # tasks in the order of replicates, whose designs share uncertain values
        tasks = [
            (i, design, replicate)
            for replicate in (
                range(self.replicates)
                if replicates is None
                else sorted(set(replicates))
            )
            for i, design in enumerate(unique)
        ]
        results: list[tuple[int, int, list[dict]]] = []
        context = ""
        if self.store is not None:
            # evaluations stored before (by an interrupted or earlier run)
            context = self.fingerprint
            stored = self.store.get(
                ((d.key, replicate) for _, d, replicate in tasks), context
            )
            others = {
                other: count
                for other, count in self.store.contexts().items()
                if other != context
            }
            if tasks and not stored and others:
                warnings.warn(
                    f"the store has {sum(others.values())} evaluations of namespace "
                    + f"{self.store.namespace!r} by studies of other definitions "
                    + "(seed, design space, functions, or version or source of TAT-C), "
                    + "which are not reused (see Study.fingerprint)",
                    stacklevel=2,
                )
            results.extend(
                (i, replicate, stored[design.key, replicate])
                for i, design, replicate in tasks
                if (design.key, replicate) in stored
            )
            tasks = [t for t in tasks if (t[1].key, t[2]) not in stored]
        total = len(tasks)
        if self.executor is None and self.n_jobs == 1:
            for done, (i, design, replicate) in enumerate(tasks, 1):
                records = self._evaluate(design, replicate)
                results.append((i, replicate, records))
                if self.store is not None and not _failed(records):
                    self.store.put([(design.key, replicate, records)], context)
                if progress is not None:
                    progress(done, total)
        elif tasks:
            executor = self.executor or _get_executor(self.workers)
            n_chunks = min(total, self.workers * self.chunks_per_job)
            size = math.ceil(total / max(n_chunks, 1))
            chunks = [tasks[k : k + size] for k in range(0, total, size)]
            futures = {
                executor.submit(self._run_chunk, chunk): chunk for chunk in chunks
            }
            done = 0
            for future in as_completed(futures):
                chunk = futures[future]
                evaluated = future.result()
                for (i, _, replicate), records in zip(chunk, evaluated):
                    results.append((i, replicate, records))
                if self.store is not None:
                    self.store.put(
                        (
                            (design.key, replicate, records)
                            for (_, design, replicate), records in zip(chunk, evaluated)
                            if not _failed(records)
                        ),
                        context,
                    )
                done += len(chunk)
                if progress is not None:
                    progress(done, total)
        results.sort(key=lambda result: (result[0], result[1]))
        return pd.DataFrame.from_records(
            [
                {k: v for k, v in record.items() if k != _FAILED}
                for _, _, records in results
                for record in records
            ]
        )


def _failed(records: list[dict]) -> bool:
    """
    Whether the records are of an evaluation that failed.
    """
    return any(record.get(_FAILED, False) for record in records)


def summarize(
    results: pd.DataFrame,
    metrics: Sequence[str],
    by: Sequence[str] = (),
    design_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """
    Summarizes the metrics of the feasible results of a study over their
    replicates: for each design (and each group of the columns `by`, such as
    a sensor), the mean of each metric, its standard deviation (`<metric>
    std`) and standard error (`<metric> sem`) between replicates, and the
    number of replicates (`replicates`), with the values of the design.

    Args:
        results (pandas.DataFrame): The results (see `Study.run`).
        metrics (Sequence[str]): The columns of the metrics.
        by (Sequence[str]): Other columns that group records of an
            evaluation.
        design_columns (Sequence[str] | None): The columns of the values of
            designs, or `None` for those that are the same in every
            replicate of a design (except the metrics).

    Returns:
        pandas.DataFrame: the summary, a row per design (and group)
    """
    feasible = results[results["feasible"]] if "feasible" in results else results
    keys = ["design", *by]
    if design_columns is None:
        excluded = set(_STATUS_COLUMNS) | set(metrics) | set(by)
        candidates = [column for column in feasible.columns if column not in excluded]
        constant = feasible.groupby("design")[candidates].nunique(dropna=False)
        design_columns = [
            column for column in candidates if (constant[column] <= 1).all()
        ]
    # metrics (such as a number of satellites) may also be values of designs
    design_columns = [
        column
        for column in design_columns
        if column not in metrics and column not in by
    ]
    grouped = feasible.groupby(keys, sort=False)
    summary = grouped[list(design_columns)].first() if design_columns else None
    statistics = grouped[list(metrics)].agg(["mean", "std", "count"])
    columns = {}
    for metric in metrics:
        count = statistics[(metric, "count")]
        columns[metric] = statistics[(metric, "mean")]
        columns[f"{metric} std"] = statistics[(metric, "std")]
        columns[f"{metric} sem"] = statistics[(metric, "std")] / np.sqrt(count)
    columns["replicates"] = grouped["replicate"].nunique()
    frame = pd.DataFrame(columns)
    if summary is not None:
        frame = summary.join(frame)
    return frame.reset_index()

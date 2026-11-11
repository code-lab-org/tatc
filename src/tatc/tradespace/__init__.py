"""
Trade space analysis: design spaces of (hierarchical) design variables and
uncertain parameters, building blocks of constellations,
trade studies that evaluate their designs (with caches
of analyses by satellite, in parallel), statistics of replicates and
requirements on them, Pareto efficiency (with the uncertainty of replicates),
adaptive replication, confirmation of selected designs with new replicates,
persistent stores of results (keyed by fingerprints of studies), sensitivity analysis,
screening with analyses of increasing fidelity,
and searches with optimizers (Optuna and pymoo, optional dependencies).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from .cache import SatelliteCache, get_satellite_key
from .confirmation import Confirmation, confirm
from .constellations import (
    build_constellation,
    count_planes,
    count_satellites,
    hours_to_time,
    planes_space,
    remove_failures,
    sso_planes_space,
    walker_space,
)
from .criteria import (
    Maximum,
    Mean,
    Minimum,
    Probability,
    Quantile,
    Requirement,
    StandardDeviation,
    Statistic,
    aggregate,
)
from .optimization import (
    Evaluations,
    make_pymoo_problem,
    optimize_optuna,
    optimize_pymoo,
)
from .fingerprint import get_fingerprint
from .pareto import hypervolume, pareto_front, pareto_probability, pareto_ranks
from .refinement import Refinement, refine
from .screening import Screening, screen
from .sensitivity import main_effects, sensitivity
from .nodes import (
    Categorical,
    Choice,
    Derived,
    Fixed,
    Infeasible,
    Integer,
    Real,
    Repeat,
    Uncertain,
    Variable,
)
from .space import Design, DesignSpace
from .store import ResultStore
from .study import Study, summarize
from .values import ChoiceValues, RepeatValues, Scope, Values

__all__ = [
    "Categorical",
    "Choice",
    "ChoiceValues",
    "Confirmation",
    "Derived",
    "Design",
    "DesignSpace",
    "Evaluations",
    "Fixed",
    "Infeasible",
    "Integer",
    "Maximum",
    "Mean",
    "Minimum",
    "Probability",
    "Quantile",
    "Real",
    "Refinement",
    "Repeat",
    "RepeatValues",
    "Requirement",
    "ResultStore",
    "SatelliteCache",
    "Scope",
    "Screening",
    "StandardDeviation",
    "Statistic",
    "Study",
    "Uncertain",
    "Values",
    "Variable",
    "aggregate",
    "build_constellation",
    "confirm",
    "count_planes",
    "count_satellites",
    "get_fingerprint",
    "get_satellite_key",
    "hours_to_time",
    "hypervolume",
    "main_effects",
    "make_pymoo_problem",
    "optimize_optuna",
    "optimize_pymoo",
    "pareto_front",
    "pareto_probability",
    "pareto_ranks",
    "planes_space",
    "refine",
    "remove_failures",
    "screen",
    "sensitivity",
    "sso_planes_space",
    "summarize",
    "walker_space",
]

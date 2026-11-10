"""
Trade space analysis: design spaces of (hierarchical) design variables and
uncertain parameters, trade studies that evaluate their designs (with caches
of analyses by satellite, in parallel), Pareto efficiency (with the
uncertainty of replicates), screening with analyses of increasing fidelity,
and searches with optimizers (Optuna and pymoo, optional dependencies).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from .cache import SatelliteCache, get_satellite_key
from .optimization import (
    Evaluations,
    make_pymoo_problem,
    optimize_optuna,
    optimize_pymoo,
)
from .pareto import hypervolume, pareto_front, pareto_probability, pareto_ranks
from .screening import Screening, screen
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
from .study import Study, summarize
from .values import ChoiceValues, RepeatValues, Scope, Values

__all__ = [
    "Categorical",
    "Choice",
    "ChoiceValues",
    "Derived",
    "Design",
    "DesignSpace",
    "Evaluations",
    "Fixed",
    "Infeasible",
    "Integer",
    "Real",
    "Repeat",
    "RepeatValues",
    "SatelliteCache",
    "Screening",
    "Scope",
    "Study",
    "Uncertain",
    "Values",
    "Variable",
    "get_satellite_key",
    "hypervolume",
    "make_pymoo_problem",
    "optimize_optuna",
    "optimize_pymoo",
    "pareto_front",
    "pareto_probability",
    "pareto_ranks",
    "screen",
    "summarize",
]

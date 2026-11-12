"""
Trade space analysis: design spaces of (hierarchical) design variables and
uncertain parameters, building blocks of constellations (their deployment and life-cycle cost,
and capability over the phases of deployment),
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
from .cost import (
    CostEstimate,
    CostModel,
    LaunchCost,
    OperationsCost,
    SatelliteCost,
    plane_sizes,
)
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
from .deployment import (
    Deployment,
    DeploymentSchedule,
    Launch,
    Production,
)
from .histories import (
    allocate_histories,
    evaluate_missions,
    expected_missions,
    summarize_missions,
)
from .mission import (
    MissionHistory,
    Phase,
    PhaseCache,
    evaluate_phases,
    summarize_phases,
)
from .failures import (
    Bathtub,
    Exponential,
    FailureModel,
    Lifetime,
    Piecewise,
    RestorationPolicy,
    Weibull,
    policy_space,
)
from .criteria import (
    Criteria,
    Maximum,
    Mean,
    Minimum,
    Probability,
    Quantile,
    Requirement,
    StandardDeviation,
    StandardError,
    Statistic,
    aggregate,
    summarize,
)
from .optimization import (
    Evaluations,
    Search,
    make_pymoo_problem,
    optimize_optuna,
    optimize_pymoo,
)
from .fingerprint import compare_definitions, get_definition, get_fingerprint
from .pareto import (
    hypervolume,
    hypervolume_progress,
    pareto_front,
    pareto_probability,
    pareto_ranks,
)
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
from .study import Study
from .values import (
    ChoiceValues,
    RepeatValues,
    Scope,
    Values,
    configure,
)

__all__ = [
    "Bathtub",
    "Categorical",
    "Choice",
    "ChoiceValues",
    "Confirmation",
    "CostEstimate",
    "CostModel",
    "Criteria",
    "Deployment",
    "DeploymentSchedule",
    "Derived",
    "Design",
    "DesignSpace",
    "Evaluations",
    "Exponential",
    "FailureModel",
    "Fixed",
    "Infeasible",
    "Integer",
    "Launch",
    "LaunchCost",
    "Lifetime",
    "Maximum",
    "Mean",
    "Minimum",
    "MissionHistory",
    "OperationsCost",
    "Phase",
    "PhaseCache",
    "Piecewise",
    "Probability",
    "Production",
    "Quantile",
    "Real",
    "Refinement",
    "Repeat",
    "RepeatValues",
    "Requirement",
    "RestorationPolicy",
    "ResultStore",
    "SatelliteCache",
    "SatelliteCost",
    "Scope",
    "Screening",
    "Search",
    "StandardDeviation",
    "StandardError",
    "Statistic",
    "Study",
    "Uncertain",
    "Values",
    "Variable",
    "Weibull",
    "aggregate",
    "allocate_histories",
    "build_constellation",
    "compare_definitions",
    "configure",
    "confirm",
    "count_planes",
    "count_satellites",
    "evaluate_missions",
    "evaluate_phases",
    "expected_missions",
    "get_definition",
    "get_fingerprint",
    "get_satellite_key",
    "hours_to_time",
    "hypervolume",
    "hypervolume_progress",
    "main_effects",
    "make_pymoo_problem",
    "optimize_optuna",
    "optimize_pymoo",
    "pareto_front",
    "pareto_probability",
    "pareto_ranks",
    "plane_sizes",
    "planes_space",
    "policy_space",
    "refine",
    "remove_failures",
    "screen",
    "sensitivity",
    "sso_planes_space",
    "summarize",
    "summarize_missions",
    "summarize_phases",
    "walker_space",
]

"""
Trade space analysis: design spaces of (hierarchical) design variables and
uncertain parameters, building blocks of constellations (their deployment and life-cycle cost,
and capability over the phases of deployment),
trade studies that evaluate their designs (in scenarios, with caches
of analyses by satellite, in parallel), studies of the observations of
constellations (their targets, deployment, costs, failures, and phases), criteria that judge designs by
statistics of their replicates (of metrics or their utility) and requirements on them, Pareto efficiency
(with the uncertainty of replicates),
adaptive replication, confirmation of selected designs with new replicates,
persistent stores of results (keyed by fingerprints of studies), sensitivity analysis
(to design variables, and attribution of variance to uncertain parameters),
screening with analyses of increasing fidelity,
searches with optimizers (Optuna and pymoo, optional dependencies),
and plots of results (matplotlib, an optional dependency).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from .attribution import Attribution, attribute_variance
from .cache import SatelliteCache, get_satellite_key
from .confirmation import Confirmation, confirm
from .cost import (
    CostEstimate,
    CostModel,
    LaunchCost,
    OperationsCost,
    SatelliteCost,
)
from .constellations import (
    build_constellation,
    count_planes,
    count_satellites,
    Layout,
    Plane,
    hours_to_time,
    layout_of,
    plane_sizes,
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
    Stage,
)
from .histories import (
    allocate_histories,
    evaluate_missions,
    expected_missions,
    summarize_missions,
)
from .mission import (
    Decision,
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
from .observations import (
    Failures,
    MissionPlan,
    ObservationStudy,
    Observations,
    Targets,
    value_of,
)
from .pareto import (
    hypervolume,
    hypervolume_progress,
    pareto_front,
    pareto_probability,
    pareto_ranks,
)
from .plotting import plot_designs, plot_means, plot_progress, style_axes
from .refinement import Refinement, refine
from .screening import Screening, screen
from .scenarios import robustness
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
from .utility import ExponentialUtility, LinearUtility, PiecewiseUtility, Utility
from .values import ChoiceValues, RepeatValues, Scope, Values, configure

__all__ = [
    "Attribution",
    "Bathtub",
    "Categorical",
    "Choice",
    "ChoiceValues",
    "Confirmation",
    "CostEstimate",
    "CostModel",
    "Criteria",
    "Decision",
    "Deployment",
    "DeploymentSchedule",
    "Derived",
    "Design",
    "DesignSpace",
    "Evaluations",
    "Exponential",
    "ExponentialUtility",
    "FailureModel",
    "Failures",
    "Fixed",
    "Infeasible",
    "Integer",
    "Launch",
    "LaunchCost",
    "Layout",
    "Lifetime",
    "LinearUtility",
    "Maximum",
    "Mean",
    "Minimum",
    "MissionHistory",
    "MissionPlan",
    "ObservationStudy",
    "Observations",
    "OperationsCost",
    "Phase",
    "PhaseCache",
    "Piecewise",
    "PiecewiseUtility",
    "Plane",
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
    "Stage",
    "StandardDeviation",
    "StandardError",
    "Statistic",
    "Study",
    "Targets",
    "Uncertain",
    "Utility",
    "Values",
    "Variable",
    "Weibull",
    "aggregate",
    "allocate_histories",
    "attribute_variance",
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
    "layout_of",
    "main_effects",
    "make_pymoo_problem",
    "optimize_optuna",
    "optimize_pymoo",
    "pareto_front",
    "pareto_probability",
    "pareto_ranks",
    "plane_sizes",
    "planes_space",
    "plot_designs",
    "plot_means",
    "plot_progress",
    "policy_space",
    "refine",
    "remove_failures",
    "robustness",
    "screen",
    "sensitivity",
    "sso_planes_space",
    "style_axes",
    "summarize",
    "summarize_missions",
    "summarize_phases",
    "value_of",
    "walker_space",
]

"""
Trade space analysis of constellations.

A trade study declares a design space (`DesignSpace`), evaluates its designs
in replicates that draw its uncertain parameters (`Study`, or
`ObservationStudy` for the observations, deployment, costs, and failures of
constellations) into results (a row per replicate), judges designs by
criteria (`Criteria`) on summaries of the results (`aggregate`, a row per
design) and their Pareto efficiency (`pareto_front`), and settles the
uncertainty of its replicates (`refine`, `confirm`), or finds efficient
designs of spaces too large to enumerate (`screen`, `optimize_optuna`,
`optimize_pymoo`).

This package exports the functions and classes of trade studies. The parts
from which they are made (such as the layouts of constellations, the
simulation of missions, and fingerprints) are in its modules.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from .attribution import attribute_variance
from .cache import SatelliteCache
from .comparison import exceedance
from .confirmation import confirm
from .constellations import (
    build_constellation,
    count_planes,
    count_satellites,
    planes_space,
    sso_planes_space,
    walker_space,
)
from .cost import (
    Component,
    CostModel,
    LaunchCost,
    OperationsCost,
    SatelliteCost,
    Triangular,
    risk_space,
)
from .criteria import (
    Across,
    Criteria,
    Maximum,
    Mean,
    Minimum,
    Probability,
    Quantile,
    Requirement,
    StandardDeviation,
    StandardError,
    aggregate,
    summarize,
)
from .deployment import DeploymentSchedule, Stage
from .failures import (
    Bathtub,
    Exponential,
    FailureModel,
    Piecewise,
    RestorationPolicy,
    Weibull,
    policy_space,
)
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
)
from .observations import Failures, MissionPlan, ObservationStudy, Observations, Targets
from .optimization import optimize_optuna, optimize_pymoo
from .pareto import (
    hypervolume,
    hypervolume_progress,
    pareto_front,
    pareto_probability,
    pareto_ranks,
)
from .plotting import (
    plot_designs,
    plot_exceedance,
    plot_means,
    plot_progress,
    style_axes,
)
from .refinement import refine
from .scenarios import robustness
from .screening import screen
from .sensitivity import main_effects, sensitivity
from .space import Design, DesignSpace
from .store import ResultStore
from .study import Study
from .utility import ExponentialUtility, LinearUtility, PiecewiseUtility, Utility

__all__ = [
    # design spaces
    "Categorical",
    "Choice",
    "Derived",
    "Design",
    "DesignSpace",
    "Fixed",
    "Infeasible",
    "Integer",
    "Real",
    "Repeat",
    "Uncertain",
    "build_constellation",
    "count_planes",
    "count_satellites",
    "planes_space",
    "policy_space",
    "sso_planes_space",
    "walker_space",
    # studies
    "Bathtub",
    "Component",
    "CostModel",
    "DeploymentSchedule",
    "Exponential",
    "FailureModel",
    "Failures",
    "LaunchCost",
    "MissionPlan",
    "ObservationStudy",
    "Observations",
    "OperationsCost",
    "Piecewise",
    "RestorationPolicy",
    "ResultStore",
    "SatelliteCache",
    "SatelliteCost",
    "Stage",
    "Study",
    "Targets",
    "Triangular",
    "Weibull",
    "risk_space",
    # criteria and summaries
    "Across",
    "Criteria",
    "ExponentialUtility",
    "LinearUtility",
    "Maximum",
    "Mean",
    "Minimum",
    "PiecewiseUtility",
    "Probability",
    "Quantile",
    "Requirement",
    "StandardDeviation",
    "StandardError",
    "Utility",
    "aggregate",
    "exceedance",
    "hypervolume",
    "hypervolume_progress",
    "main_effects",
    "pareto_front",
    "pareto_probability",
    "pareto_ranks",
    "robustness",
    "sensitivity",
    "summarize",
    # refinement and searches
    "attribute_variance",
    "confirm",
    "optimize_optuna",
    "optimize_pymoo",
    "refine",
    "screen",
    # plots
    "plot_designs",
    "plot_exceedance",
    "plot_means",
    "plot_progress",
    "style_axes",
]

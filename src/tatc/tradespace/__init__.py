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

This package exports the functions and classes of trade studies, and the
outcomes of their steps (such as `Refinement`), for type hints. The parts
from which they are made (such as the layouts of constellations, the
simulation of missions, and fingerprints) are in its modules.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from .attribution import Attribution, attribute_variance
from .cache import SatelliteCache
from .comparison import exceedance
from .confirmation import Confirmation, confirm
from .constellations import (
    build_constellation,
    count_planes,
    count_satellites,
    layout_of,
    mog_space,
    planes_space,
    soc_space,
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
from .criteria import Across, Criteria, Requirement
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
from .observations import ObservationStudy, Observations, Targets
from .optimization import Search, optimize_optuna, optimize_pymoo
from .pareto import hypervolume, hypervolume_progress, pareto_front, pareto_ranks
from .plans import Failures, MissionPlan
from .plotting import (
    COLORS,
    GRAY,
    INK,
    plot_designs,
    plot_exceedance,
    plot_means,
    plot_progress,
    plot_shares,
    style_axes,
    use_style,
)
from .references import configure, value_of
from .refinement import Refinement, refine
from .scenarios import robustness
from .screening import Screening, screen
from .sensitivity import main_effects, sensitivity
from .space import Design, DesignSpace
from .statistics import (
    Maximum,
    Mean,
    Minimum,
    Probability,
    Quantile,
    StandardDeviation,
    StandardError,
    Statistic,
)
from .store import ResultStore
from .study import Study
from .summary import aggregate, pareto_probability, summarize
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
    "configure",
    "count_planes",
    "count_satellites",
    "layout_of",
    "mog_space",
    "planes_space",
    "policy_space",
    "soc_space",
    "sso_planes_space",
    "value_of",
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
    "Statistic",
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
    "Attribution",
    "Confirmation",
    "Refinement",
    "Screening",
    "Search",
    "attribute_variance",
    "confirm",
    "optimize_optuna",
    "optimize_pymoo",
    "refine",
    "screen",
    # plots
    "COLORS",
    "GRAY",
    "INK",
    "plot_designs",
    "plot_exceedance",
    "plot_means",
    "plot_progress",
    "plot_shares",
    "style_axes",
    "use_style",
]

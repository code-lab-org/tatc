"""
The trade space of the trade space tutorial, which its notebooks import: the
parts that a notebook defines and later notebooks use, under the name of the
notebook that defines and explains them.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import tempfile
from datetime import datetime, timedelta, timezone
from pathlib import Path

from tatc.analysis import estimate_sampling_error
from tatc.generation import generate_points_rotated_lattice
from tatc.schemas import Instrument
from tatc.tradespace import (
    Categorical,
    Choice,
    CostModel,
    DeploymentSchedule,
    DesignSpace,
    Integer,
    LaunchCost,
    MissionPlan,
    ObservationStudy,
    Observations,
    OperationsCost,
    Real,
    ResultStore,
    SatelliteCost,
    Targets,
    Uncertain,
    build_constellation,
    count_satellites,
    sso_planes_space,
    walker_space,
)

# EvaluateTradeSpace: the design space, and the study of sufficient sampling

START = datetime(2026, 7, 1, tzinfo=timezone.utc)
DURATION = timedelta(days=1)  # the period of the analysis
N_POINTS = 100  # points sampled in each replicate
N_REPLICATES = 5  # replicates (samples of points) of each design
SAMPLES_PER_DAY = 4  # observations a day of a sufficiently sampled point
MAX_SATELLITES = 10

space = DesignSpace(
    altitude=Real(500e3, 800e3, step=100e3),
    field_of_regard=Categorical([30, 60]),
    layout=Choice(
        walker=walker_space(
            inclination=Real(45, 75, step=15),
            planes=Integer(1, 4),
            per_plane=Integer(1, 4),
        ),
        sso=sso_planes_space(
            crossing=Real(6, 18, step=3),
            satellites=Integer(1, 2),
            count=Integer(1, 2),
        ),
    ),
    points_seed=Uncertain.seed(),
    constraints=[lambda d: count_satellites(d.layout) - MAX_SATELLITES],
)
designs = space.enumerate()


def build(d):
    """The satellites of a design: an imager in each orbit of its layout."""
    imager = Instrument(name="imager", field_of_regard=d.field_of_regard)
    return build_constellation(d.layout, imager, epoch=START, altitude=d.altitude)


def count_samples(observations, duration=DURATION):
    """The share (%) of points observed at least SAMPLES_PER_DAY times a day."""
    if len(observations) == 0:  # no satellites
        return {"sufficient": 0.0}
    samples = observations.groupby("target_hash").size()
    threshold = SAMPLES_PER_DAY * duration / timedelta(days=1)
    return {"sufficient": 100 * (samples >= threshold).sum() / N_POINTS}


CACHE_DIRECTORY = Path(tempfile.gettempdir()) / "tatc-trade-space"
# the points of each replicate: a lattice in a random orientation, from its seed
TARGETS = Targets(generate_points_rotated_lattice, N_POINTS, seed="points_seed")
OBSERVATIONS = Observations(TARGETS, START, DURATION, directory=CACHE_DIRECTORY)

study = ObservationStudy(
    space,
    build,
    OBSERVATIONS,
    count_samples,
    layout="layout",
    replicates=N_REPLICATES,
    n_jobs=-1,
    store=ResultStore(CACHE_DIRECTORY / "results.sqlite"),
)
OBJECTIVES = {"satellites": "min", "planes": "min", "sufficient": "max"}


def describe(d):
    """A short description of a design (a row of its values)."""
    if d.layout == "walker":
        layout = (
            f"Walker {d['layout.walker.inclination']:.0f} deg: "
            f"{d['layout.walker.planes']:.0f} x {d['layout.walker.per_plane']:.0f} "
            f"(spacing {d['layout.walker.relative_spacing']:.0f})"
        )
    else:
        planes = [
            f"{d[f'layout.sso.planes[{k}].satellites']:.0f} at "
            f"{d[f'layout.sso.planes[{k}].crossing']:02.0f}:00"
            for k in range(int(d["layout.sso.planes.count"]))
        ]
        layout = "SSO: " + ", ".join(planes)
    return f"{d.altitude / 1e3:.0f} km, {d.field_of_regard} deg FOR, {layout}"


def planes_name(planes):
    """The name of a number of planes."""
    return f"{planes:.0f} plane{'s' if planes > 1 else ''}"


# ExplainTradeSpace: the study of accurate sampling as well

DECORRELATION_TIME = timedelta(hours=3)  # of the field whose daily mean is estimated
MAX_RELATIVE_ERROR = 0.5  # of an accurate sample, in standard deviations of the field


def sample(observations, decorrelation_time=DECORRELATION_TIME):
    """The shares (%) of points sampled sufficiently and accurately."""
    if len(observations) == 0:  # no satellites
        return {"sufficient": 0.0, "accurate": 0.0}
    samples = observations.groupby("target_hash").size()
    end = START + DURATION
    error = estimate_sampling_error(observations, START, end, decorrelation_time)
    return {
        "sufficient": 100 * (samples >= SAMPLES_PER_DAY).sum() / N_POINTS,
        "accurate": 100 * (error.relative_error <= MAX_RELATIVE_ERROR).sum() / N_POINTS,
    }


accuracy_study = study.replace(metrics=sample)
ACCURACY = {"satellites": "min", "planes": "min", "accurate": "max"}

# CompareScenarios: fields of other decorrelation times

FIELDS = {
    "1 hour": {"decorrelation_time": timedelta(hours=1)},
    "3 hours": {},
    "12 hours": {"decorrelation_time": timedelta(hours=12)},
}


def sample_field(observations, design):
    """The shares (%) of points sampled sufficiently and accurately (of its field)."""
    return sample(observations, design.decorrelation_time)


# TradeCostAndUtility: the deployment and cost of each design

SCHEDULE = DeploymentSchedule(
    years=5,
    development_years=3,
    production_rate=4,
    launch_interval=0.5,
    launch_capacity=3,
)
COST_MODEL = CostModel(
    satellites={
        30: SatelliteCost(development=40, first_unit=12, learning_rate=0.9, mass=60),
        60: SatelliteCost(development=60, first_unit=18, learning_rate=0.9, mass=90),
    },
    launch=LaunchCost(per_launch=8),
    operations=OperationsCost(per_year=3, per_satellite_year=0.5),
    discount_rate=0.03,
)
PLAN = MissionPlan(
    SCHEDULE, COST_MODEL, kind="field_of_regard", growth={"development": "cost_growth"}
)
cost_space = space.extend(
    cost_growth=Uncertain(lambda rng: float(rng.triangular(0.9, 1.0, 1.6)))
)
cost_study = accuracy_study.replace(space=cost_space, plan=PLAN)

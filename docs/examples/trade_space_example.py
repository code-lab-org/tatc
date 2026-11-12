"""
The design space and study of the example `SearchTradeSpace`, which the
examples that follow it (`RankUncertainDesigns`, `CompareScenarios`,
`TradeCostAndUtility`, and `PlanMissions`) import. See `SearchTradeSpace`
for what each part means.

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
    Fixed,
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

START = datetime(2026, 7, 1, tzinfo=timezone.utc)
DURATION = timedelta(days=1)  # the period of the analysis
N_POINTS = 100  # points sampled in each replicate
N_REPLICATES = 5  # replicates (samples of points) of each design
SAMPLES_PER_DAY = (
    4  # a point is sufficiently sampled if it is observed at least this often
)
MAX_SATELLITES = 10
DECORRELATION_TIME = timedelta(hours=3)  # of the field whose daily mean is estimated
MAX_RELATIVE_ERROR = 0.5  # a point is accurately sampled within this share of the field's standard deviation

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
            crossing=Real(6, 18, step=3), satellites=Integer(1, 2), count=Integer(1, 2)
        ),
    ),
    points_seed=Uncertain.seed(),
    failure_seed=Uncertain.seed(),
    cost_growth=Uncertain(
        lambda rng: float(rng.triangular(0.9, 1.0, 1.6))
    ),  # of development, for cost risk
    decorrelation_time=Fixed(
        DECORRELATION_TIME
    ),  # of the field, which scenarios may change
    epoch=Fixed(
        START
    ),  # the start of the analysis, which may be drawn for each replicate
    constraints=[lambda d: count_satellites(d.layout) - MAX_SATELLITES],
)
designs = space.enumerate()

# the deployment of each design: 3 years of development, a production line of 4 satellites a year, launches of up to 3
# satellites at most every half year, and 5 years of operations from the first launch
SCHEDULE = DeploymentSchedule(
    years=5,
    development_years=3,
    production_rate=4,
    launch_interval=0.5,
    launch_capacity=3,
)
# illustrative costs (millions of dollars), not estimates of any mission: satellites with an imager of each field of
# regard, small launchers, and operations, discounted at 3% a year
COST_MODEL = CostModel(
    satellites={
        30: SatelliteCost(development=40, first_unit=12, learning_rate=0.9, mass=60),
        60: SatelliteCost(development=60, first_unit=18, learning_rate=0.9, mass=90),
    },
    launch=LaunchCost(per_launch=8),
    operations=OperationsCost(per_year=3, per_satellite_year=0.5),
    discount_rate=0.03,
)

# the deployment and cost of each design: satellites of the kind of its field of regard, with its growth of
# development costs
PLAN = MissionPlan(
    SCHEDULE, COST_MODEL, kind="field_of_regard", growth={"development": "cost_growth"}
)
# the points of each replicate: a randomly oriented lattice, generated from the replicate's seed
TARGETS = Targets(generate_points_rotated_lattice, N_POINTS, seed="points_seed")
CACHE_DIRECTORY = Path(tempfile.gettempdir()) / "tatc-search-trade-space"


def sample(observations, design, duration=DURATION):
    """The shares (%) of points sampled sufficiently and accurately by observations over a duration from the design's
    epoch (of the field of the design's decorrelation time)."""
    if len(observations) == 0:  # no satellites
        return {"sufficient": 0.0, "accurate": 0.0}
    samples = observations.groupby("target_hash").size()
    error = estimate_sampling_error(
        observations, design.epoch, design.epoch + duration, design.decorrelation_time
    )
    threshold = SAMPLES_PER_DAY * duration / timedelta(days=1)
    return {
        "sufficient": 100 * (samples >= threshold).sum() / N_POINTS,
        "accurate": 100 * (error.relative_error <= MAX_RELATIVE_ERROR).sum() / N_POINTS,
    }


def build(d):
    """The satellites of a design: an imager of its field of regard in each orbit of its layout."""
    imager = Instrument(name="imager", field_of_regard=d.field_of_regard)
    return build_constellation(d.layout, imager, epoch=START, altitude=d.altitude)


# the study of the shares (%) of points sampled sufficiently and accurately over the period by the satellites of each
# design, and of its deployment and cost (by the mission plan), whose evaluations are kept in a store
study = ObservationStudy(
    space,
    build,
    Observations(TARGETS, "epoch", DURATION, directory=CACHE_DIRECTORY),
    sample,
    layout="layout",
    plan=PLAN,
    replicates=N_REPLICATES,
    n_jobs=-1,
    store=ResultStore(CACHE_DIRECTORY / "results.sqlite"),
)

OBJECTIVES = {"satellites": "min", "planes": "min", "sufficient": "max"}
ACCURACY = {"satellites": "min", "planes": "min", "accurate": "max"}
REFERENCE = {
    "satellites": MAX_SATELLITES + 1,
    "planes": 5,
    "sufficient": 0,
}  # worse than every design, for hypervolumes


def describe(d):
    """A short description of a design (a row of its values)."""
    if d.layout == "walker":
        layout = (
            f"Walker {d['layout.walker.inclination']:.0f} deg: {d['layout.walker.planes']:.0f} x "
            f"{d['layout.walker.per_plane']:.0f} (spacing {d['layout.walker.relative_spacing']:.0f})"
        )
    else:
        planes = [
            f"{d[f'layout.sso.planes[{k}].satellites']:.0f} at {d[f'layout.sso.planes[{k}].crossing']:02.0f}:00"
            for k in range(int(d["layout.sso.planes.count"]))
        ]
        layout = "SSO: " + ", ".join(planes)
    return f"{d.altitude / 1e3:.0f} km, {d.field_of_regard} deg FOR, {layout}"


def planes_name(planes):
    """The name of a number of planes."""
    return f"{planes:.0f} plane{'s' if planes > 1 else ''}"

"""
Defines analysis functions.
"""

from .coverage import (
    aggregate_observations,
    collect_multi_observations,
    collect_observations,
    grid_observations,
    reduce_observations,
)
from .dop import DopMethod, compute_dop
from .latency import (
    collect_downlinks,
    compute_latencies,
    grid_latencies,
    reduce_latencies,
)
from .ro_coverage import (
    collect_ro_observations,
)
from .track import (
    OrbitCoordinate,
    OrbitOutput,
    collect_ground_pixels,
    collect_ground_track,
    collect_orbit_track,
    compute_ground_track,
)

__all__ = [
    "DopMethod",
    "OrbitCoordinate",
    "OrbitOutput",
    "aggregate_observations",
    "collect_downlinks",
    "collect_ground_pixels",
    "collect_ground_track",
    "collect_multi_observations",
    "collect_observations",
    "collect_orbit_track",
    "collect_ro_observations",
    "compute_dop",
    "compute_ground_track",
    "compute_latencies",
    "grid_latencies",
    "grid_observations",
    "reduce_latencies",
    "reduce_observations",
]

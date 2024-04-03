"""
Defines analysis functions.
"""

from .coverage import (
    collect_observations,
    collect_multi_observations,
    aggregate_observations,
    reduce_observations,
    grid_observations,
)
from .track import collect_orbit_track, collect_ground_track, compute_ground_track, OrbitCoordinate, OrbitOutput
from .latency import (
    collect_downlinks,
    compute_latencies,
    reduce_latencies,
    grid_latencies,
    compute_latencies_former
)

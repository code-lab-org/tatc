"""
Defines analysis functions.
"""

from .dop import DopMethod, compute_dop
from .ground_track import (
    collect_ground_pixels,
    collect_ground_track,
    compute_ground_track,
)
from .latency import (
    collect_downlinks,
    compute_latencies,
    grid_latencies,
    reduce_latencies,
)
from .limb_coverage import (
    ScanDirection,
    collect_limb_observations,
)
from .orbit_track import (
    OrbitCoordinate,
    OrbitOutput,
    collect_orbit_track,
)
from .point_coverage import (
    aggregate_observations,
    collect_multi_observations,
    collect_observations,
    grid_observations,
    reduce_observations,
)
from .radar import (
    collect_radar_track,
    compute_radar_track,
)
from .region_coverage import (
    collect_multi_region_observations,
    collect_region_observations,
)
from .ro_coverage import (
    collect_ro_observations,
)

__all__ = [
    "DopMethod",
    "OrbitCoordinate",
    "OrbitOutput",
    "ScanDirection",
    "aggregate_observations",
    "collect_downlinks",
    "collect_ground_pixels",
    "collect_ground_track",
    "collect_limb_observations",
    "collect_multi_observations",
    "collect_multi_region_observations",
    "collect_observations",
    "collect_orbit_track",
    "collect_radar_track",
    "collect_region_observations",
    "collect_ro_observations",
    "compute_dop",
    "compute_ground_track",
    "compute_latencies",
    "compute_radar_track",
    "grid_latencies",
    "grid_observations",
    "reduce_latencies",
    "reduce_observations",
]

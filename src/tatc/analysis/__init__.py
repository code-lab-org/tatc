"""
Defines analysis functions.
"""

from .coverage_metrics import (
    aggregate_observations,
    grid_observations,
    reduce_observations,
)
from .dop_sampling import DopMethod, compute_dop
from .ground_track import (
    collect_ground_pixels,
    collect_ground_track,
    compute_ground_track,
)
from .latency_sampling import (
    collect_downlinks,
    compute_latencies,
    grid_latencies,
    reduce_latencies,
)
from .limb_sampling import (
    ScanDirection,
    collect_limb_observations,
)
from .orbit_track import (
    OrbitCoordinate,
    OrbitOutput,
    collect_orbit_track,
)
from .point_sampling import (
    collect_multi_observations,
    collect_observations,
    compute_access_periods,
)
from .radar_track import (
    collect_radar_track,
    compute_radar_track,
)
from .region_sampling import (
    collect_multi_region_observations,
    collect_region_observations,
    compute_region_access_periods,
)
from .ro_sampling import (
    collect_ro_observations,
)
from .space_sampling import (
    collect_space_observations,
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
    "collect_space_observations",
    "compute_access_periods",
    "compute_dop",
    "compute_ground_track",
    "compute_latencies",
    "compute_radar_track",
    "compute_region_access_periods",
    "grid_latencies",
    "grid_observations",
    "reduce_latencies",
    "reduce_observations",
]

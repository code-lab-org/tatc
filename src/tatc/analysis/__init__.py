"""
Defines analysis functions.
"""

from .airborne_track import (
    collect_airborne_ground_track,
    compute_airborne_ground_track,
)
from .architecture import (
    ArchitectureEvaluation,
    collect_platform_observations,
    evaluate_architecture,
)
from .coincidence import (
    find_coincident_ground_track,
    find_coincident_observations,
)
from .coverage_metrics import (
    aggregate_observations,
    count_observations,
    grid_observations,
    reduce_observations,
)
from .data_flow import (
    DataFlowResult,
    analyze_data_flow,
    compute_data_flow,
    compute_generation,
    find_data_flow_downlinks,
    is_data_flow_modeled,
)
from .dop_sampling import DopMethod, compute_dop
from .feature_sampling import collect_feature_observations
from .field_sampling import sample_field, sample_profiles
from .footprints import (
    collect_footprint_observations,
    collect_footprints,
    dissolve_footprints,
)
from .ground_track import (
    collect_ground_pixels,
    collect_ground_track,
    compute_ground_track,
)
from .latency_sampling import (
    collect_downlinks,
    collect_relay_contacts,
    compute_latencies,
    grid_latencies,
    reduce_latencies,
)
from .link_sampling import (
    compute_link_rates,
    get_generation_rate,
    sample_link_geometry,
)
from .limb_sampling import (
    ScanDirection,
    collect_limb_observations,
    collect_limb_soundings,
)
from .observation_quality import (
    compute_observation_quality,
    simulate_observation_errors,
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
    collect_occultations,
    collect_ro_observations,
)
from .sample_sizes import compute_sample_sizes
from .space_sampling import (
    collect_space_observations,
)
from .sufficiency import (
    compare_sampled_means,
    compute_local_time_coverage,
    compute_revisit_statistics,
    estimate_sampling_error,
    sampling_error_variance,
)
from .surface_track import (
    collect_platform_track,
    collect_surface_track,
    compute_platform_track,
    compute_surface_track,
)

__all__ = [
    "DataFlowResult",
    "analyze_data_flow",
    "compute_data_flow",
    "compute_generation",
    "find_data_flow_downlinks",
    "is_data_flow_modeled",
    "compute_link_rates",
    "get_generation_rate",
    "sample_link_geometry",
    "ArchitectureEvaluation",
    "DopMethod",
    "OrbitCoordinate",
    "OrbitOutput",
    "ScanDirection",
    "aggregate_observations",
    "collect_airborne_ground_track",
    "collect_downlinks",
    "collect_feature_observations",
    "collect_footprint_observations",
    "collect_footprints",
    "collect_ground_pixels",
    "collect_ground_track",
    "collect_limb_observations",
    "collect_limb_soundings",
    "collect_multi_observations",
    "collect_multi_region_observations",
    "collect_observations",
    "collect_orbit_track",
    "collect_platform_observations",
    "collect_platform_track",
    "collect_radar_track",
    "collect_region_observations",
    "collect_relay_contacts",
    "collect_occultations",
    "collect_ro_observations",
    "collect_space_observations",
    "collect_surface_track",
    "compare_sampled_means",
    "compute_access_periods",
    "compute_airborne_ground_track",
    "compute_dop",
    "compute_ground_track",
    "compute_latencies",
    "compute_local_time_coverage",
    "compute_observation_quality",
    "compute_platform_track",
    "compute_radar_track",
    "compute_region_access_periods",
    "compute_revisit_statistics",
    "compute_sample_sizes",
    "compute_surface_track",
    "count_observations",
    "dissolve_footprints",
    "estimate_sampling_error",
    "evaluate_architecture",
    "find_coincident_ground_track",
    "find_coincident_observations",
    "grid_latencies",
    "grid_observations",
    "reduce_latencies",
    "reduce_observations",
    "sample_field",
    "sample_profiles",
    "sampling_error_variance",
    "simulate_observation_errors",
]

"""
Utility functions for the TATC library.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from .formatting import zero_pad
from .geometry import normalize_geometry, project_polygon_to_elevation, split_polygon
from .observation import (
    compute_field_of_regard,
    compute_max_access_time,
    compute_max_transit_time,
    compute_min_along_track_distance,
    compute_min_elevation_angle,
    field_of_regard_to_swath_width,
    swath_width_to_field_of_regard,
    swath_width_to_field_of_view,
)
from .orbital import (
    compute_ground_inertial_velocity,
    compute_ground_surface_velocity,
    compute_j2_aop_rate,
    compute_j2_raan_rate,
    compute_orbit_inertial_velocity,
    mean_anomaly_to_true_anomaly,
    mean_motion_to_orbit_period,
    mean_motion_to_semimajor_axis,
    semimajor_axis_to_mean_motion,
    semimajor_axis_to_orbit_period,
    true_anomaly_to_mean_anomaly,
)
from .projection import (
    buffer_footprint,
    buffer_target,
    compute_footprint,
    compute_limb,
    compute_projected_ray_position,
)
from .surface import compute_number_samples
from .time import to_datetime64_ns

__all__ = [
    "buffer_footprint",
    "buffer_target",
    "compute_field_of_regard",
    "compute_footprint",
    "compute_ground_inertial_velocity",
    "compute_ground_surface_velocity",
    "compute_j2_aop_rate",
    "compute_j2_raan_rate",
    "compute_limb",
    "compute_max_access_time",
    "compute_max_transit_time",
    "compute_min_along_track_distance",
    "compute_min_elevation_angle",
    "compute_number_samples",
    "compute_orbit_inertial_velocity",
    "compute_projected_ray_position",
    "field_of_regard_to_swath_width",
    "mean_anomaly_to_true_anomaly",
    "mean_motion_to_orbit_period",
    "mean_motion_to_semimajor_axis",
    "normalize_geometry",
    "project_polygon_to_elevation",
    "semimajor_axis_to_mean_motion",
    "semimajor_axis_to_orbit_period",
    "split_polygon",
    "swath_width_to_field_of_regard",
    "swath_width_to_field_of_view",
    "to_datetime64_ns",
    "true_anomaly_to_mean_anomaly",
    "zero_pad",
]

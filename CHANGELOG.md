# TAT-C Change Log

## 3.4.5

Added:
 - Support for Python 3.13.
 - Analysis function `collect_orbit_track` arguments `sat_sunlit` and `solar_altaz` optionally report satellite sunlit and solar altitude/azimuth metrics.
 
Changed:
 - Fixes bug in analysis function `compute_ground_track` when specifying a `mask`.

## 3.4.4
 
Changed:
 - Fixes bug in analysis function `collect_orbit_track` when specifying a `mask`.
 - Allows `mask` argument to analysis functions `collect_orbit_track`, `collect_ground_track`, `compute_ground_track`, and `collect_ground_pixels` to have a GeoDataFrame, GeoSeries, or list-like type.

## 3.4.3

Added:
 - Analysis function `collect_ground_track` arguments `sat_altaz` and `solar_altaz` optionally report satellite and solar altitude/azimuth metrics.
 - Analysis function `collect_ground_pixels` allows pixel-level analysis similar to ground tracks.
 - Analysis function `compute_limb` computes the observable limb (maximum extent of viewable Earth).
 - Schema `PointedInstrument` variables `cross_track_pixels`, `along_track_pixels`, `cross_track_oversampling`, `along_track_oversampling` enable pixel-level analysis.
 - Schema `PointedInstrument` methods `compute_footprint`, `compute_footprint_center`, `compute_projected_pixel_position`, `get_pixel_cone_and_clock_angle`, and `compute_footprint_pixel_array` refactor utility functions and add pixel-level capabilities.
 
Changed:
 - Updates analysis function `collect_orbit_track` to check if a masked area contains each subsatellite point prior to additional computation.
 - Refactors analysis function `collect_orbit_track` and utility function `compute_projected_ray_position` to use more vectorized calculations.
 - Updates analysis function `collect_ground_track` to check if a masked area contains each footprint center point prior to additional computation.
 - Updates analysis function `collect_ground_track` to compute ground tracks at a designated elevation.
 - Renames utility function `_get_footprint_point` to `compute_projected_ray_position`.
 - Refactors `Instrument` method `is_valid_observation` argument from `targets` to `target`.
 - Refactors utility function `compute_footprint_center` return type from a Skyfield GeographicPosition to a Shapely Geometry.
 - Replaced memory-intensive cross-product in analysis function `collect_ro_observations` and utility function `compute_projected_ray_position`.
 - Fixes a bug in utility function `_split_polygon_antimeridian` where a ground track footprint could circle the wrong pole if the smallest longitude coordinate is in the opposite north/south hemisphere.

## 3.4.2

Changed:
 - Analysis functions `collect_observations` and `collect_multi_observations` constrain observations to be contained within the ground track when using an instrument that inherits from `PointedInstrument` (address off-nadir pointing).

## 3.4.1

Changed:
 - Use default runtime configuration settings from constructor if default config file cannot be loaded.

## 3.4.0

Added:
 - Dependency `spiceypy` to access SPICE routines.
 - Dependency `pyyaml` to read YAML files.
 - Argument `dissolve_orbits` to analysis function `compute_ground_track` that, if set to `False`, optionally preserves individual orbits.
 - Argument `crs="spice"` to analysis functions `collect_ground_track` and `compute_ground_track` to leverage the SPICE routines to quickly compute footprints.
 - Off-nadir pointing instruments for analysis functions `collect_ground_track` and `compute_ground_track` by instantiating `PointedInstrument` instruments.
 - Schema for off-nadir pointing instruments `PointedInstruments`.
 - Utility function `swath_width_to_field_of_view` to compute along/cross track swath width for off-nadir pointing instruments.
 - Utility function `compute_footprint` to compute footprints using SPICE routines.
 - Utility function `buffer_footprint` to compute footprints using shapely buffering (refactored).
 - Utility function `project_polygon_to_elevation` to assign a constant elevation to polygons (refactored).
 - Configuration file to store/load runtime configurations. Defaults located at `resources/defaults.yaml` which are loaded to `tatc.config.rc` with schema `tatc.config.RuntimeConfiguration`.
 - Member function `as_skyfield` to `TwoLineElements` to easily construct a Skyfield `EarthSatellite` object.
 - Member function `get_repeat_cycle` to `TwoLineElements` to lazy-load a calculated repeat cycle.
 - Member function `get_orbit_track` to `TwoLineElements` to easily compute a Skyfield `Geocentric` object optionally leveraging a repeat cycle.
 - Member function `get_observation_events` to `TwoLineElements` to replicate the Skyfield `find_events` method optionally leveraging a repeat cycle.

Changed:
 - Refactored analysis function `collect_observations` to propagate with new `TwoLineElements` methods capable of leveraging a computed repeat cycle.
 - Refactored analysis function `collect_ro_observations` to propagate with new `TwoLineElements` methods capable of leveraging a computed repeat cycle.
 - Refactored private analysis functions `_get_visible_interval_series`, `_get_satellite_altaz_series`, `_get_satellite_sunlit_series`, `_get_solar_altaz_series`, and `_get_solar_time_series` to use TAT-C objects as input arguments.
 - Refactored `Instrument` member function `is_valid_observation` to reference a pre-computed Skyfield `Geocentric` input argument.
 - Refactored utility function `_split_polygon_antimeridian` within `split_polygon` to check if a footprint contains a pole before splitting along the antimeridian.

Removed:
 - Analysis function `_get_access_series` (refactored).
 - Analysis function `_get_revisit_series` (refactored).


## 3.3.0

Added:
- Schema for streets-of-coverage constellation `SOCConstellation`.
- Schemas for eccentric high-altitude orbits: `MolniyaOrbit` and `TundraOrbit`.
- Optional dependency `cartopy` for examples.

Changed:
- Signature for functions `collect_orbit_track`, `collect_ground_track`, `compute_ground_track`, and `collect_observations` passes in `instrument_index` rather than the `instrument` object, with a default value of 0 (first instrument).
- Signature for function `collect_ro_observations` replaces `max_azimuth` with `max_yaw` constraint.
- Signature for private function `_collect_ro_series` passes in receiver vertical, normal, binormal (VNB) unit vectors.
- Examples use `cartopy` for coastlines and stock imagery.
# TAT-C Change Log

## 3.6.0

Added new observing capabilities for limb sounding, ground-based radar, and regions (polygons), with support for shapely points. Changed how general perturbations orbits are propagated, with options to model orbits maintained on a repeat ground track. Also changed the tangent point for GNSS radio occultation (RO) and limb sounding to the point of minimum WGS 84 geodetic altitude along the line of sight.

Added:
 - Analysis method `collect_limb_observations` (with enumeration `ScanDirection`) for limb sounding observations from periodic scans at a constant angular rate.
 - Analysis method `collect_space_observations` for the periods when satellites can observe other satellites (as by a sensor or a radio frequency link), subject to slant range limits and Earth occlusion (a minimum WGS 84 geodetic altitude of the line of sight). Each observation records the lines of sight at its start and end as its geometry.
 - Object schemas `RadarStation` (with `RadarBand` and `RadarStation.from_band`) and `TerrainMask`, and analysis methods `collect_radar_track` and `compute_radar_track`, for ground-based radar coverage. `AllSurfaceObjects` includes `RadarStation`.
 - Module `preprocess` to derive inputs from external datasets, including terrain masks from Copernicus DEM data (`compute_terrain_mask`, `compute_terrain_mask_for_station`, `get_copernicus_dem_tile_urls`, and `sample_dem_elevation`).
 - Analysis methods `collect_region_observations` and `collect_multi_region_observations` to collect the periods (start and end) when an instrument observes any part of a region (shapely `Polygon` or `MultiPolygon`), using the footprints of pointed and conical instruments and the field of regard of other instruments. Each observation records its swath (the part of the region swept by the footprint) as its geometry and identifies the region by the hash of its geometry (`target_hash`). Methods `aggregate_observations` and `reduce_observations` combine the swaths of a region into the part of it observed.
 - Analysis methods `compute_access_periods`, the periods when a satellite is in view of a point above a minimum elevation angle, and `compute_region_access_periods`, a conservative superset of the periods when an instrument's field of regard may observe any part of a region (for example, to cull analysis times).
 - Support for shapely `Point` objects (longitude, latitude, and optional elevation) wherever a TAT-C `Point` is accepted (`collect_observations`, `collect_multi_observations`, `compute_dop`, and `GeneralPerturbationsOrbit.get_observation_events`), in anticipation of replacing TAT-C points.
 - Object schema `ConicalInstrument` (included in `AllInstruments`) for conically scanning instruments, observed when a point crosses the scanned cone, and utility method `compute_cone_and_azimuth`.
 - Instrument pointing options: fields `PointedInstrument.velocity_frame` (enumeration `VelocityFrame`, for spacecraft without yaw steering), `Instrument.nadir_reference` (enumeration `NadirReference`, for a geocentric nadir), `PointedInstrument.view_geometry` (enumeration `ViewGeometry`, for cross-track scanners), and `PointedInstrument.roll_angle_profile` and `pitch_angle_profile` (with methods `get_roll_angle` and `get_pitch_angle`) to vary pointing around the orbit. Projection functions accept matching arguments, including `tilt_angle`.
 - Fields `Instrument.min_target_solar_elevation` and `Instrument.max_target_solar_elevation` to require a range of solar elevation angles at the target.
 - Fields `GeneralPerturbationsOrbit.remove_drag`, to propagate without drag, and `GeneralPerturbationsOrbit.repeat_cycle`, to propagate the elements directly (`None`, the default), maintained on a repeat ground track with a repeat cycle found from the elements (`"auto"`), or with a declared repeat cycle (a duration), also accepted by `from_tle`, `from_omm_csv`, and `from_omm_json`. Related methods `GeneralPerturbationsOrbit.get_repeat_element` and `GeneralPerturbationsElements.get_repeat_element`, `refine_repeat_cycle`, `without_drag`, `is_sun_synchronous`, and property `has_drag`.
 - Field `inclination` of `MolniyaOrbit` and `TundraOrbit` (by default, the critical inclination), and methods `from_apogee_longitude` and `get_apogee_longitude` to place an orbit by the longitude of its apogee.
 - Field `WalkerConstellation.seam_spacing` (and method `get_seam_spacing`) to set the seam of a star configuration.
 - Argument `during_contact` of `compute_latencies` to choose when an observation that ends during a downlink is downlinked (`"end"`, `"next"`, or `"immediate"`).
 - Argument `systems` of `compute_dop` to estimate a receiver clock bias for each navigation system.
 - Utility methods `compute_view_tangents`, `compute_view_angles`, `compute_along_track_field_of_view`, `compute_argument_of_latitude`, `geodesic_destination`, `pressure_to_altitude`, and `altitude_to_pressure`, radar geometry methods (`compute_radar_beam_height`, `compute_radar_ground_range`, `compute_radar_ground_range_bounds`, `compute_radar_slant_range`, `compute_terrain_elevation_angle`, `compute_radar_footprint`, and `compute_radar_footprint_profile`), and method `PointedInstrument.is_in_field_of_view`.
 - Utility methods for the WGS 84 ellipsoid (module `tatc.utils.ellipsoid`): `geodetic_to_rectangular`, `rectangular_to_geodetic`, `compute_ellipsoid_intersection` (intersections of rays with the ellipsoid), and `compute_tangent_point` (the point of minimum geodetic altitude on a line), utility method `compute_vnb_frame` (a satellite's velocity, normal, and binormal frame), and utility method `hash_geometry` (a compact hash identifying a geometry).
 - Runtime configurations `footprint_points_radar_azimuthal`, `repeat_cycle_delta_semimajor_axis_m`, and `nutation_interpolation_minutes` (the step of a cached table from which nutation angles are interpolated, or `null` to compute them for every time), and constants `EFFECTIVE_EARTH_RADIUS_FACTOR`, `EARTH_ROTATION_RATE`, and `TROPICAL_YEAR_S`.
 - Validation notebooks in `docs/validation` comparing TAT-C with data from operational missions.
 - Optional dependencies `preprocess` and `validation`, and development dependency `pytest-xdist` to run tests in parallel (e.g. `pytest -n 8`).

Changed:
 - `GeneralPerturbationsOrbit` propagates orbit tracks and observation events consistently, as selected by its `repeat_cycle`: if set, times before the first and after the last element's epoch are propagated with that element maintained on its repeat ground track; other times use the closest element. Orbits are now propagated directly with SGP4 by default; use `remove_drag=True` and `repeat_cycle="auto"` for the former default. Removed the runtime settings `repeat_cycle_for_orbit_track` and `repeat_cycle_for_observation_events`, the `try_repeat` arguments, and methods `get_observation_repeat_cycle`, `get_repeat_shifts`, and `get_repeat_orbit_track`. A warning is issued when elements with drag are propagated far from their epoch.
 - `get_repeat_cycle` verifies repeat cycles without drag, reports a whole number of nodal days (or mean solar days for a sun-synchronous orbit), and also confirms repeat cycles for elements within `repeat_cycle_delta_semimajor_axis_m` of an exact repeat. Cached repeat cycles are recomputed when the search options change. Elements with a perigee below the Earth's surface (for example, from a mean motion in revolutions per day rather than radians per minute) have no repeat cycle, with a warning, rather than a search lasting minutes.
 - Removed the runtime settings `repeat_cycle_lazy_load` and `gp_orbit_lazy_load` and the `lazy_load` arguments of `to_gp_orbit` and `get_repeat_cycle`: converted general perturbations orbits and repeat cycles are always cached, and recomputed when the fields or search options they depend on change.
 - Guards against common input errors: orbit and general perturbations element epochs must be timezone-aware, orbits with a perigee below the Earth's surface (for example, from a semimajor axis in kilometers) raise a validation error, and analyses warn of satellites with a perigee below 100 km (for example, from an altitude in kilometers). Analysis time windows must have timezone-aware `start` and `end` datetimes with `end` no earlier than `start` (a `ValueError`), and an `instrument_index` out of range raises an `IndexError` naming the satellite.
 - `get_observation_events` refines rise and set times to a millisecond, and finds passes that culminate outside the search period or span a switch between elements.
 - `get_observation_events` completes passes whose rise or set Skyfield misses (grazing passes, which barely exceed the minimum elevation angle), finding the missing event from the culmination and dropping culminations below the minimum. `compute_access_periods` and `collect_observations` no longer treated such a pass as visible from the start (or to the end) of the analysis period, which made a conical instrument observe points through the Earth.
 - `collect_observations` refines access periods to the instrument's field of regard, uses the time the view sweeps over the point as the epoch of a `PointedInstrument` observation, and evaluates orbits with a repeat cycle consistently.
 - `collect_observations` and `collect_region_observations` accept one or more points or regions and one or more satellites, and `instrument_index=None` for every instrument of each satellite, computing all of their observations together (concatenated and sorted by start). `collect_multi_observations` and `collect_multi_region_observations` are kept as wrappers for backwards compatibility.
 - `collect_orbit_track`, `collect_ground_track`, `collect_ground_pixels`, `collect_limb_observations`, and `collect_downlinks` accept one or more satellites (and the track methods `instrument_index=None` for every instrument), and `collect_ro_observations` one or more receivers, with the results of several concatenated and sorted by time (or start). Downlinks and RO observations of several satellites are computed together.
 - `PointedInstrument` roll and pitch angles rotate the view rigidly (roll, then pitch about the rolled cross-track axis), and `compute_view_tangents` accepts them. Documented the pointing conventions (positive roll looks left, positive pitch looks forward, and cross-track pixel indices run from right to left).
 - `compute_latencies` downlinks an observation that ends during a downlink at the end of that downlink by default (`during_contact="end"`; formerly `"next"`).
 - `collect_ro_observations` uses the point of minimum WGS 84 geodetic altitude on the receiver-transmitter line as the tangent point.
 - `SOCConstellation` uses the polar streets-of-coverage pattern for near-polar orbits (new field `polar` and methods `is_polar`, `get_footprint_angle`, `get_polar_design`, `get_satellites_per_plane`, and `get_number_planes`; `generate_walker` raises a `ValueError` for polar designs) and a continuously covering hexagonal lattice for inclined orbits, which requires more satellites at a given packing distance.
 - `collect_ground_track` and `collect_ground_pixels` cull times with a mask using the periods when the instrument's field of regard may observe it, which no longer misses footprints near the poles or across the anti-meridian. Masks are split along the anti-meridian and poles before use, and results are sorted by time.
 - Observations of points (`collect_observations`, `collect_multi_observations`, and `compute_latencies`) identify each point by the hash of its geometry (`target_hash`, see `hash_geometry`) rather than a `point_id`, and `aggregate_observations`, `reduce_observations`, and `reduce_latencies` group observations by `target_hash`. The `id` of a TAT-C `Point` is no longer used.
 - Generated points and cells are identified by the hashes of their geometries (`point_id` and `cell_id`) rather than integer indices, so that a generated point's `point_id` equals the `target_hash` of its observations.
 - Generating points and cells is about 3 to 15 times faster (vectorized geometry construction, hashing, and clipping). Points and cells clipped to a mask keep the grid or lattice order rather than the arbitrary order of `geopandas.clip`.
 - Analysis methods that accept satellites raise a `TypeError` for a constellation (use its `generate_members` method).
 - Deprecated utility method `buffer_target`.
 - Fixed `SunSynchronousOrbit` to place the ascending node in local mean solar time and to keep the local time of its equator crossings constant as propagated.
 - Fixed `TrainConstellation` for long intervals between members, and `MolniyaOrbit` and `TundraOrbit` periods to repeat the ground track as propagated.
 - Fixed rectangular footprints of wide, elongated views.
 - Fixed `split_polygon` for longitudes beyond 180 degrees, for polygons spanning -180 to 180 degrees of longitude, and for holes across the anti-meridian or around a pole (holes are now split like exteriors), and to drop degenerate (zero-area) parts left by repairing invalid polygons.
 - Fixed cached orbit computations for copies with changed fields, and orbits and satellites to be picklable after propagation.
 - Improved the performance of footprints (rays projected together from one view frame, vectorized ray and limb intersections, polygon assembly, splitting along the anti-meridian, and `split_polygon`), masked ground tracks, multi-element orbit propagation, `collect_ro_observations` (observation periods of all transmitters searched and sampled together, with sample times and profiles built at once), `collect_limb_observations` (all scans sampled together, with sample times built at once), radar footprints (terrain masks interpolated at all azimuths at once, and a closed-form local projection instead of PROJ transformers for each station), `collect_observations` (linear-time pairing of rise and set events, and validity checks of all epochs at once), `collect_region_observations` (validity checks and swaths of all periods at once), `collect_orbit_track` (masks tested and points built for all times at once), observations of several points, regions, or satellites (Earth orientation computed together at each step of their searches), rise and set times (each bracketed near Skyfield's estimate and bisected only until it is narrower than a millisecond), orbit propagation (IAU 2000A nutation angles interpolated from a cached table, to within about a microarcsecond, rather than computed for every time), instrument validity checks (without iterating over Skyfield times), and the searches that refine observation periods (times built from offsets without `datetime` objects). `to_datetime64_ns` and `GeneralPerturbationsOrbit.get_closest_element_index` accept Skyfield `Time` objects.
 - Reorganized modules: `tatc.analysis.coverage` is split into `point_sampling`, `region_sampling`, `sampling` (shared by both), and `coverage_metrics` (the aggregation, reduction, and gridding of observations), and `tatc.analysis.track` into `orbit_track` and `ground_track`; `tatc.analysis.ro_coverage`, `latency`, and `dop` are renamed `ro_sampling`, `latency_sampling`, and `dop_sampling`; tangent point geometry moved from `tatc.analysis` to `tatc.utils.ellipsoid`, and the radar footprint functions from `tatc.utils.projection` to `tatc.utils.radar`. Public functions remain available from `tatc.analysis` and `tatc.utils`. Orbit propagation algorithms moved to `tatc.utils.propagation`, with the conversion of Skyfield times in `tatc.utils.time`, interpolated nutation angles in `tatc.utils.earth_orientation`, and computations run together in `tatc.utils.computation`.
 - Updated example notebooks.

## 3.5.1

Minor refactoring to improve backwards compatibility.

Added:
 - Added a simplified `TwoLineElements` object schema as a thin wrapper around the new `GeneralPerturbationsOrbit`.
 - Added `altitude` as an alias for `mean_altitude` for circular orbits.

Changed:
 - Fixed an overly-sensitive unit test that failed on different platforms due to numerical differences.

## 3.5.0

Major refactoring that focused on completing unit tests to approach full code coverage. Drops support for Python < 3.10 to improve compatibility with modern libraries. Replaces the `TwoLineElements` schema with `GeneralPerturbationsOrbit` to support modern orbit specifications including OMM CSV and JSON. Drops the CRS-based buffering approach to determine ground track in favor of geometric projections using the SPICE library. During refactoring, a few bug fixes and breaking changes were also made.

Added:
 - Added `GeosynchronousOrbit` orbit schema.
 - Support for Python 3.14.

Changed:
 - Changed the default orbit epoch from `datetime.now()` to `2020-01-01T00:00:00Z`.
 - Refactored all orbit schemas to inherit uniform getters from `OrbitBase`.
 - Replaced `TwoLineElements` with a more general `GeneralPerturbationsOrbit` to accommodate post-TLE GP data formats.
 - Improved `TundraOrbit` and `MolniyaOrbit` schemas to consider J2 perturbations when calculating orbit period.
 - Fixed a bug where SPICE projected instrument footprints around the geocentric, rather than geodetic, pointing vector, leading to ~2 km positioning errors.
 - Fixed `TrainConstellation` member right ascension of ascending node spacing to be based on a sidereal day, rather than a solar day.
 - Fixed `SOCConstellation` member generation to use hexagonal spacing with a non-zero `relative_spacing` value.
 - Fixed `MOGConstellation` member generation to use the mean anomaly of the reference orbit.
 - Fixed a bug in `config.py` where default configurations (`defaults.yml`) were never provided in wheels.
 - Fixed `collect_orbit_track` to only assign the `EPSG:4326` CRS to output when requesting WGS84 coordinates, rather than for all coordinate systems.
 - Fixed `collect_orbit_track` to use a proper East/North/Up velocity when requesting WGS84 coordinates.
 - Improved `collect_observations` to use the apogee altitude, rather than the initial altitude, to determine the maximum access duration.
 - Fixed a bug where `collect_multi_observations` could crash on an empty satellite list.
 - Fixed a bug in `grid_observations` and `grid_latencies` where spatial aggregation never actually worked.
 - Fixed a bug in `compute_dop` where the latitude and longitude were reversed in output geometry.
 - Improved `compute_dop` to use the nearest GP element to each time rather than the first one specified for an orbit.
 - Fixed the definition of binormal unit vector in `ro_coverage.py` to accommodate eccentric orbits.
 - Improved `collect_ro_observations` to interpolate among samples closest to target elevation.
 - Improved performance of `collect_ro_observations` through vectorized profile sampling, a more direct interface to Skyfield for orbit track computation, and making tangent point velocity calculation optional.
 - Fixed a bug where `collect_ro_observations` could use incorrect inertial positions for repeat track orbits more than 1 cycle after epoch.
 - Added utility methods: `compute_apoapsis_radius` and `geodesic_distance`.
 - Requires `setuptools >= 77.0.0` and switches `project.license` to an SPDX expression string to resolve a build metadata deprecation warning (issue #124).

Removed
 - Support for Python < 3.10.
 - Removed the `TwoLineElements` orbit schema.
 - Removed the legacy `crs` and `method` parameters from `collect_ground_track` and `compute_ground_track`; all instrument projection now uses the SPICE library.

## 3.4.10

Added:
 - Argument `solar_beta` to `collect_orbit_track` to optionally compute solar beta angle.
 - `CITATION.cff` file to provide guidance on citing this project.

## 3.4.9

Added:
 - Utility function `buffer_target` to buffer a target region based on orbit and instrument geometry to help with culling.

Changed:
 - Analysis functions `collect_ground_track` and `collect_ground_pixels` check whether the sub-satellite point is inside an expanded mask region (using `buffer_target`) for culling, rather a slower check for nonzero intersections between the mask region and the viewable limb projection.


## 3.4.8

Changed:
 - Improves performance for `TwoLineElements` orbits having multiple (perhaps many) TLE pairs.
 - Coerces `collect_observations` observation intervals constrained by `start` and `end` times to use the `datetime.timezone.utc` timezone value to address inability for pandas to combine `TzInfo(UTC)` and `datetime.UTC` representations.
 - Fixes a geometry indexing bug in `collect_ground_track` with option `crs="spice"` when `times` is a list of length 1.

## 3.4.7

Added:
 - Allows `TwoLineElements` orbits to accept multiple TLE pairs to more accurately reproduce historical orbits.
 
Changed:
 - Improves analysis function `collect_ground_pixels` when specifying a `mask` to include pixels when the footprint center falls outside the masked domain.
 - Fixes bugs in utility functions `compute_projected_ray_position` and `compute_limb` if the requested times are an array of length 1.
 - Fixes a bug in orbit determination for orbits with a detected repeat cycle where large state errors occurred at simulation times prior to the epoch date.

## 3.4.6
 
Changed:
 - Improves analysis function `collect_ground_track` when specifying a `mask` to include ground track area when the footprint center falls outside the masked domain.
 - Fixes bugs in analysis functions `collect_orbit_track`, `collect_ground_track`, and `collect_ground_pixels` when specifying a `mask` with multiple geometries in a GeoSeries or GeoDataFrame.

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
# TAT-C Change Log

## 3.6.0

Added new observing capabilities for limb sounding, ground-based radar, and regions (polygons), with support for shapely points. Also changed the tangent point definition for GNSS radio occultation (RO) and limb sounding from the line's closest approach to Earth's center to its minimum WGS 84 geodetic altitude.

Added:
 - Analysis method `collect_limb_observations` to identify limb sounding observations based on periodic scans at a constant angular rate, with a `ScanDirection` enumeration to specify scan direction.
 - Object schema `RadarStation` to model ground-based radar, including beam width and an optional `RadarBand` with `RadarStation.from_band` to apply a nominal maximum range typical of the band. `AllSurfaceObjects` now includes `RadarStation`.
 - Object schema `TerrainMask` to model terrain features that impede radar observation.
 - Analysis methods `collect_radar_track` and `compute_radar_track` to determine observable radar geometries.
 - Module `preprocess` with methods that derive TAT-C inputs from external datasets. For example, `compute_terrain_mask` (and `compute_terrain_mask_for_station`) creates a `TerrainMask` from Digital Elevation Model (DEM) data, `get_copernicus_dem_tile_urls` locates Copernicus DEM GLO-30 tiles, and `sample_dem_elevation` samples DEM elevation at a point.
 - Utility methods `pressure_to_altitude` and `altitude_to_pressure` based on the US Standard Atmosphere (1976).
 - Utility methods for radar geometry: `compute_radar_beam_height`, `compute_radar_ground_range`, `compute_radar_ground_range_bounds`, `compute_radar_slant_range`, `compute_terrain_elevation_angle`, `compute_radar_footprint`, and `compute_radar_footprint_profile`.
 - Utility method `geodesic_destination`.
 - Runtime configuration `footprint_points_radar_azimuthal` and constant `EFFECTIVE_EARTH_RADIUS_FACTOR` (4/3 Earth radius refraction model).
 - Enumeration `VelocityFrame` and field `PointedInstrument.velocity_frame` to orient an instrument's view relative to the inertial (orbital) velocity, as for a spacecraft without yaw steering, instead of the Earth-fixed velocity (the default, as for a yaw-steered spacecraft). The two differ by up to about 4 degrees in low Earth orbit, which shifts the edges of a wide view by up to about 100 km along track. Projection functions accept a matching `velocity_frame` argument.
 - Utility method `compute_view_tangents` to compute a target's along-track and cross-track view angles from a satellite, method `PointedInstrument.is_in_field_of_view`, and constant `EARTH_ROTATION_RATE`.
 - Object schema `ConicalInstrument` to model conically scanning instruments (such as microwave imagers), specified by a cone angle, a scan sector (center azimuth and half width, for forward, aft, or full-rotation scans), a velocity frame, and an along-track field of view. `collect_observations` reports the time a point crosses the cone within the scan sector as the observation epoch (up to two per pass), and its footprint is the scanned arc swept along track by the distance the along-track field of view subtends at nadir. `AllInstruments` now includes `ConicalInstrument`, and `collect_orbit_track` reports its swath width from the scan sector.
 - Utility method `compute_cone_and_azimuth` to compute a target's cone angle and scan azimuth from a satellite.
 - Enumeration `NadirReference` and field `Instrument.nadir_reference` to orient an instrument's view relative to the geocentric nadir (toward the Earth's center), as for spacecraft such as Sentinel-2, instead of the geodetic nadir (the WGS 84 ellipsoid normal, the default). The two differ by up to about 0.19 degrees at middle latitudes, which shifts a projected view by up to about 2.7 km from 800 km altitude. Projection functions accept a matching `nadir_reference` argument.
 - Utility method `compute_along_track_field_of_view` to tailor an instrument's along-track field of view to a time step (the angle subtended at nadir by the distance the ground track advances in one time step), by default at the fastest ground track velocity over the orbit, so that consecutive footprints tile without gaps.
 - Method `GeneralPerturbationsOrbit.get_repeat_element` to get an element maintained on its repeat ground track, with which times before the orbit's first and after its last element's epoch are propagated. Methods `GeneralPerturbationsElements.without_drag`, `is_sun_synchronous`, and property `has_drag`, and argument `remove_drag` of `GeneralPerturbationsElements.to_skyfield`, whose result is now cached.
 - Fields `Instrument.min_target_solar_elevation` and `Instrument.max_target_solar_elevation` to require a range of solar elevation angles at the target for valid observations (for example, a minimum for daylight optical imaging or a maximum for night-time imaging), extending the target sunlit requirement (`req_target_sunlit`, equivalent to a minimum or maximum of 0 degrees).
 - Argument `during_contact` of `compute_latencies` to set when an observation that ends while a downlink is in progress is downlinked: `"end"` (the default) downlinks it at the end of the downlink in progress (as for stored data played back after the data recorded before the contact, such as NOAA-20's), `"next"` (the previous behavior) waits for the next downlink, and `"immediate"` downlinks it as it is observed (real-time downlink).
 - Field `inclination` of `MolniyaOrbit` and `TundraOrbit`, defaulting to the critical inclination (about 63.4 degrees, the previously fixed value), to model orbits such as QZSS's quasi-zenith orbits (about 40 degrees). Away from the critical inclination, the argument of perigee precesses as propagated, and the orbit period accounts for its effect on the apogee's longitude so that the apogee longitudes still repeat.
 - Argument `systems` of `compute_dop` to label each satellite's navigation system (for example, "GPS" or "Galileo") and estimate a receiver clock bias per system with visible satellites, as multi-GNSS receivers do, instead of a single clock bias (the default). The system of the first satellite is the reference: TDOP is its clock's dilution of precision and GDOP is the square root of PDOP squared plus TDOP squared, both undefined (NaN) at times when no satellite of the reference system is visible. Times with fewer visible satellites than 3 plus the number of visible systems return NaN. A clock per system raises the PDOP of combined constellations by about 1%.
 - Fields `GeneralPerturbationsOrbit.remove_drag`, to propagate an orbit without drag (ignoring the elements' B* and mean motion derivatives, while the Earth's oblateness still perturbs the orbit), as for an orbit maintained against drag, and `GeneralPerturbationsOrbit.repeat_cycle`, to select how the orbit is propagated: `None` (the default) to propagate the elements directly with SGP4; `"auto"` to propagate it as maintained on a repeat ground track with the repeat cycle found from the elements, if any; or a duration to declare the repeat cycle of such an orbit instead of searching for it. A repeat cycle (found or declared) requires `remove_drag`, as for an orbit maintained against drag. The elements are kept as given, and propagated maintained on the repeat ground track (with method `GeneralPerturbationsElements.get_repeat_element`): without drag, with the mean motion adjusted so that a whole number of orbits spans the repeat cycle, refined to the nearest whole number of nodal days or of mean solar days for a sun-synchronous orbit (with method `GeneralPerturbationsElements.refine_repeat_cycle`), and, for a sun-synchronous orbit, with the inclination adjusted (by a small fraction of a degree) so that the nodal day is a mean solar day, as for an orbit maintained at a constant local time of ascending node. The satellite then returns to its initial position after each cycle: with its own mean motion, an ICESat-2 element set (67 m of semimajor axis from the exact repeat) returns 874 km (about 113 s) away after 91 days, a discontinuity at the end of each repeated cycle, and within 1 km once adjusted (the remainder from the eccentricity terms as the argument of perigee precesses). Long repeat cycles cannot be found from a single element set: for ICESat-2 (91 days), drag and small errors of mean motion move the propagated satellite farther from its initial position than chance near-repeats, whereas the declared repeat cycle is within 25 s of the measured one. `from_tle`, `from_omm_csv`, and `from_omm_json` accept both options.
 - Field `PointedInstrument.roll_angle_profile` to vary an instrument's roll angle around the orbit, as (argument of latitude, roll angle) pairs interpolated linearly and periodically (for example, the roll steering of Sentinel-1, whose look angles vary by up to 2.5 degrees: a fitted profile reduces swath edge errors from up to 20 km to 7 km), method `PointedInstrument.get_roll_angle`, and utility method `compute_argument_of_latitude`.
 - Field `PointedInstrument.pitch_angle_profile` to vary an instrument's pitch angle around the orbit, as (argument of latitude, pitch angle) pairs interpolated linearly and periodically like `roll_angle_profile` (for example, the tilt of PACE OCI, which changes once per orbit from 20 degrees aft to 20 degrees forward: with each day's profile, observation times are within 1 s of OCI's, instead of 38 s without the tilt), and method `PointedInstrument.get_pitch_angle`.
 - Enumeration `ViewGeometry` and field `PointedInstrument.view_geometry` to define an instrument's fields of view and pixels in cross-track and along-track angles (`scan`, as for a cross-track scanner, whose along-track angular extent is constant across the scan) instead of in the plane perpendicular to the boresight (`frame`, the default, as for a framing camera or pushbroom array). In scan geometry, the pitch angle tilts the scan plane (a rotation about the cross-track axis applied before the roll angle), as for a scanner that tilts fore or aft such as PACE OCI, whose pixels a tilted scan reproduces within 0.8 km near nadir and 9 km at the swath edges (a constant along-track angle would instead sweep a cone that misplaces the edges by 210 km). For VIIRS, the along-track length of a scan footprint at the swath edges is 25.9 km in scan geometry (25 km observed) and 14.5 km in frame geometry. Projection functions accept a matching `view_geometry` argument, `compute_projected_ray_position` accepts a `tilt_angle` argument, and utility method `compute_view_angles` gives a target's roll and pitch angles (in a scan plane tilted by its `tilt_angle` argument).
 - Method `from_apogee_longitude` of `MolniyaOrbit` and `TundraOrbit` to place an orbit by the longitude of its apogee (the center of a Tundra or quasi-zenith orbit's figure-8 ground track), solving for the right ascension of ascending node so that the first apogee after the epoch, as propagated, is over the longitude, and method `get_apogee_longitude` to get that longitude.
 - Field `WalkerConstellation.seam_spacing` to narrow (or widen) the seam of a star configuration, between the counter-rotating first and last planes, with the other planes equally spaced over the remainder of 180 degrees (for example, Iridium's six planes 31.6 degrees apart with a 22 degree seam), and method `WalkerConstellation.get_seam_spacing`. By default, star planes remain equally spaced. Satellites moving in opposite directions across the seam cover a narrower street, so equally spaced planes leave coverage gaps at the seam (up to 1.9 minutes for an Iridium-like 66/6 constellation).
 - Analysis methods `collect_region_observations` and `collect_multi_region_observations` (module `region_coverage`) to collect observations of a region: a shapely `Polygon` or `MultiPolygon`, at the elevation of the mean of its z coordinates (if any). A region is observed while any part of it lies within the instrument's field of regard and above the satellite's horizon (only the field of regard is considered, for any instrument), and each observation's epoch is the period's midpoint, with its validity (illumination) and satellite and solar angles evaluated for the region's point nearest to nadir. Observation periods are found from conservative overpass periods (see below), refined so that a pass over separate parts of a region is divided into separate observations. A region about 100 m across is observed within 0.03 s of the point at its center.
 - Support for shapely `Point` objects (longitude, latitude, and optional elevation in meters) wherever a TAT-C `Point` is accepted: `collect_observations`, `collect_multi_observations`, `compute_dop`, and `GeneralPerturbationsOrbit.get_observation_events`, in anticipation of replacing TAT-C points with shapely points. Observations record a TAT-C point's `id` as their `point_id`, or 0 for a shapely point or a region.
 - Validation notebooks in `docs/validation` for the MLS and SABER limb sounders, COSMIC-2 and PlanetiQ GNSS RO, NEXRAD radar, the NOAA-20 ATMS and VIIRS imagers, the Sentinel-1C and NISAR synthetic aperture radars, the Sentinel-2 MSI pushbroom imager, Sentinel-2 constellation revisit, and Sentinel-2 region coverage, the AMSR2 and GMI conically scanning radiometers, orbit and constellation generation, GNSS dilution of precision, NOAA-20 downlink latency, Landsat 8 and 9 solar geometry and acquisition, the SWOT wide-swath altimeter, the ICESat-2 narrow-beam lidar, the Iridium streets-of-coverage constellation, the GOES-19 geostationary imager, the International Space Station's low orbit (ECOSTRESS and EMIT), Landsat 9's maintained orbit over 19 months, the PACE OCI tilted scanner, Terra MODIS and Terra's drifting orbit, OCO-2's nadir and glint pointing, and GPM's radar and microwave imager across a yaw flip.
 - Optional dependencies `preprocess` (for the `preprocess` module) and `validation` (for validation notebooks; includes `examples` and `preprocess`).

Changed:
 - Changed `collect_observations` for a `PointedInstrument` to check the field of view at the time the view sweeps over the point (when the point's along-track view angle equals the pitch angle), which is now the observation epoch, rather than at the midpoint of the field of regard access period. Views that are narrow along track (for example, tailored to a scan period) previously missed observations, especially near swath edges.
 - Changed `PointedInstrument` roll and pitch angles (and the `roll_angle` and `pitch_angle` arguments of projection functions) to rotate the view rigidly: first by the roll angle about the along-track axis, then by the pitch angle about the rolled cross-track axis. Previously, views were offset from nadir in tangent space, so a rolled view was narrower than its field of view on the side away from nadir (for example, a 40 degree view rolled by 10 degrees spanned -10.8 to 28.4 degrees instead of -10 to 30 degrees), and combined roll and pitch angles also shifted the view center. Views with only a roll or only a pitch angle keep the same center. `compute_view_tangents` accepts roll and pitch angles to express a target relative to a pointed view's center.
 - Corrected documentation of pointing conventions: positive roll looks to the left of the direction of motion, positive pitch looks forward, and cross-track pixel indices run from right to left (index 0 is on the right of the direction of motion).
 - Fixed rectangular footprints for wide, elongated views: corners are now located using the ratio of half-width tangents (previously the ratio of fields of view, which overshot the long sides near the corners) and polygon points are evenly spaced along each side (previously clustered near the middle of the long sides, cutting off parts of the footprint).
 - Fixed `SunSynchronousOrbit` to place the ascending node at the equator crossing time in local mean solar time (universal time plus longitude), using the Greenwich mean sidereal time of the epoch. It previously used the apparent right ascension of the Sun in the J2000 frame, which shifted the ascending node by the equation of time and precession since 2000 (for example, 12.7 minutes early in October 2026, and from about -15 to +18 minutes over a year).
 - Fixed `SunSynchronousOrbit` inclination to keep the local time of the equator crossings constant as propagated by SGP4: the classical J2 estimate is refined so that SGP4's secular precession of the ascending node (which includes the J4 term) equals the mean Sun's (360 degrees per tropical year). The inclination increases by about 0.012 degrees (for example, 98.192 instead of 98.180 degrees at 705 km); previously, the local time drifted by about 2 minutes per year. Added constant `TROPICAL_YEAR_S`.
 - Fixed cached orbit computations (`to_gp_orbit`, the `MolniyaOrbit` and `TundraOrbit` periods, and repeat cycles) to be recomputed for copies with changed fields: `model_copy(update=...)` copies the cache, so that, for example, a `SunSynchronousOrbit` copied with a new altitude was converted with its original altitude.
 - Fixed `TrainConstellation` for intervals longer than a few minutes: member satellites now trail along the orbit at the lead orbit's propagated (SGP4) rate, and, when repeating the ground track, their ascending nodes account for the precession of the orbit plane during the interval. Previously, a sun-synchronous train with an interval of days was placed out of the lead's plane (by about 1 degree per day) and out of phase (for example, 1,800 km from Sentinel-2C trailing Sentinel-2B by five days).
 - Fixed `MolniyaOrbit` and `TundraOrbit` periods to repeat the ground track as propagated, accounting for the precession of the ascending node (computed with SGP4's secular rates). Previously, the apogee longitudes of a Molniya orbit drifted westward by about 0.1 degrees per day.
 - Changed `SOCConstellation` to use the polar streets-of-coverage pattern (Rider, 1985; Adams and Rider, 1987) for near-polar orbits (inclination within 10 degrees of 90 degrees, or as set by the new field `polar`): planes spanning 180 degrees of right ascension, with co-rotating planes and the counter-rotating planes across the seam spaced for continuous single coverage, adjacent planes offset by half the in-plane spacing, and the number of satellites per plane that minimizes the total. The packing distance reduces the footprint used in the design as a coverage margin. Previously, near-polar designs used the Walker delta pattern of inclined orbits, with planes spanning 360 degrees and footprints that only touched along each plane: for Iridium's altitude and minimum elevation angle, 99 satellites that still left coverage gaps, instead of Iridium's 66 (6 planes of 11), which the polar pattern reproduces. Added methods `is_polar`, `get_footprint_angle`, `get_polar_design`, `get_satellites_per_plane`, and `get_number_planes`; `generate_walker` raises a `ValueError` for polar designs, whose planes are unequally spaced.
 - Changed `SOCConstellation` for inclined orbits (the Walker delta pattern) to place footprints on the hexagonal lattice that just covers continuously at a packing distance of 1: footprint centers are spaced by the square root of 3 footprint radii along each plane and 1.5 footprint radii between planes, instead of 2 and the square root of 3 footprint radii (Eqs. (23)-(24) in Anderson et al. (2022)), with which footprints only touched along each plane and left gaps (for example, 0.4% of the area below 46 degrees latitude uncovered at any time for 1,500 km swaths from 550 km altitude at 53 degrees inclination). Designs at a given packing distance have about 30% more satellites; the previous design corresponds to a packing distance about 15% larger (2 divided by the square root of 3).
 - Changed `GeneralPerturbationsOrbit.get_repeat_cycle` (and `GeneralPerturbationsElements.get_repeat_cycle`) to verify candidate repeat cycles without drag (as for an orbit maintained on its repeat ground track) and to report a whole number of nodal days (of the element maintained on its repeat ground track), or of mean solar days for a sun-synchronous orbit, rather than the time of the propagated satellite's closest return to its initial position. Detected repeat cycles change by seconds (for example, from 16 days minus 12.5 s to 16 days for a Landsat 9 TLE), which previously accumulated over the repeated cycles of long analyses (about 7.5 minutes of timing error over 19 months for Landsat 9). A repeat cycle is also confirmed if the elements' semimajor axis is within the new runtime configuration `repeat_cycle_delta_semimajor_axis_m` (100 m by default) of the exact repeat and the elements adjusted to it return to their initial state, as for elements taken anywhere within an orbit's maintenance band (for example, within about 100 m rather than 70 m of a 16-day repeat), for repeat cycles short enough (up to about 32 days in low Earth orbit) that this tolerance does not admit chance near-repeats. Previously, detection failed (and long analyses fell back to propagation with drag) for elements far enough from the exact repeat for the ground track to drift by more than 10 km per cycle. Cached repeat cycles are recomputed for searches with other options, including changes of the runtime configuration.
 - Changed `compute_latencies` to downlink an observation that ends while a downlink is in progress at the end of that downlink (`during_contact="end"`) instead of at the midpoint of the next downlink. Latencies of such observations are shorter; for NOAA-20, the previous behavior (`during_contact="next"`) overestimated the mean latency from observation to downlink by about a third (31 versus 23 minutes).
 - Changed `collect_observations` to refine access periods to the times when a point's angle from nadir equals half the instrument's field of regard. Access periods were previously derived from a minimum elevation angle for a spherical Earth and the orbit's apogee altitude above the mean Earth radius, which did not account for the Earth's oblateness (periods could be too short at high latitudes) and only approximated the field of regard. The start, end, and (for nadir instruments) epoch of observations change by up to a few seconds, and passes that only graze the field of regard may no longer be reported.
 - Changed `collect_observations`, when the orbit is propagated with a repeat cycle, to model an orbit maintained on its repeat ground track consistently: epochs, fields of view, and satellite angles are evaluated with the same repeated orbit track as the observation events. Previously, epochs were solved on the directly propagated orbit, which drifts from the repeated access periods (by about a minute after 90 days for a sun-synchronous orbit), so that observations by narrow views could be missed in long analyses.
 - Changed `GeneralPerturbationsOrbit` to propagate orbit tracks (`get_orbit_track`, `get_orbit_track_at_time`, `get_geographic_position`) and observation events (`get_observation_events`) consistently with repeat cycles: if the orbit is propagated with a repeat cycle (`GeneralPerturbationsOrbit.repeat_cycle`), declared or found, times before the first element's epoch and after the last element's epoch (for a single element, all times) are propagated with that element maintained on its repeat ground track (`get_repeat_element`), shifted by whole repeat cycles to within one repeat cycle of its epoch, with the satellite's Earth-fixed position and velocity expressed in the inertial frame at the time itself; times between the first and last epochs are propagated with the closest element. So that consecutive repeated cycles join, the maintained element is used whether the repeat cycle is declared or detected: with the element's own mean motion, the repeated orbit tracks of Landsat 9 element sets jumped by a median of 65 km (up to about 200 km, or about 9 s along track) at the end of each repeat cycle, and with the maintained element, by about 1 km. Previously, observation events were repeated with the first cycle of the analysis period, propagated (with drag) to its start however far from the epoch (for a Sentinel-2A TLE, observation times were off by a median of 9 minutes for an analysis starting 6 months after the epoch, and 34 minutes after 12 months); orbit tracks, and so ground tracks and footprints (`collect_orbit_track`, `collect_ground_track`, `collect_ground_pixels`, `compute_ground_track`), were always propagated directly, with drag, so that far from the epoch they did not match observations (by about 1,200 km after 100 days for a Landsat 8 element set) and a `mask`, culled with repeated positions, dropped footprints; and orbits with multiple elements (such as an element set history) were never repeated. Methods `get_observation_repeat_cycle`, `get_repeat_shifts`, and `get_repeat_orbit_track` were removed.
 - Removed the runtime settings `repeat_cycle_for_orbit_track` and `repeat_cycle_for_observation_events` and the `try_repeat` arguments of `GeneralPerturbationsOrbit.get_geographic_position`, `get_geographic_position_at_time`, and `get_observation_events` in favor of `GeneralPerturbationsOrbit.repeat_cycle`, so that each orbit selects how it is propagated, consistently in every analysis; the settings are ignored if present in a configuration file. Orbits are now propagated directly with SGP4 by default (as with `try_repeat=False`), including those converted from other orbit types (`to_gp_orbit`); to propagate an orbit as maintained on its repeat ground track (formerly the default), use `remove_drag=True` and `repeat_cycle="auto"`. A warning is issued when elements with drag are propagated directly farther from their epoch than the `repeat_cycle_search_duration_days` setting.
 - Fixed `GeneralPerturbationsOrbit.get_observation_events` (and so the access periods of `collect_observations`) to refine rise and set times by bisection to within a millisecond. Skyfield's `find_events` stops refining once the first of its unequal search brackets converges, so over long analysis periods some rise and set times were several seconds early or late (for example, 5 s for a 35 degree minimum elevation from a 440 km orbit over one day), which also caused `ConicalInstrument` cone crossings near the start of an access period to be missed.
 - Fixed `GeneralPerturbationsOrbit.get_observation_events` (and so `collect_observations` and `collect_downlinks`) to find the set of a pass that culminates before the start of a search period and the rise of one that culminates after its end, which Skyfield's `find_events` misses. Passes at the start or end of an analysis period were missed or truncated, as were passes spanning the switch between elements of an orbit with multiple elements (for a Landsat 9 element set history, 2 of 8,241 observations over 19 months). A pass also rises or sets at a switch between elements (or repeat cycles) if the satellite crosses the minimum elevation angle there.
 - Fixed orbits (and satellites) to be picklable, and deep copies and derived orbits to succeed, after propagation: cached values, such as Skyfield satellites, which cannot be pickled, are neither pickled nor copied, so that, for example, `joblib` can send satellites to other processes.
 - Changed the tangent point in `collect_ro_observations` to the point of minimum WGS 84 geodetic altitude on the receiver-transmitter line, matching operational RO geolocation. Tangent point locations shift by up to about 20 km horizontally at middle latitudes relative to 3.5.x; altitudes change by only meters.
 - Improved `GeneralPerturbationsOrbit` multi-TLE propagation by sharing expensive Skyfield quantities across time slices.
 - Improved `collect_ro_observations` performance by computing transmitter azimuth from existing positions rather than re-propagating.
 - `to_datetime64_ns` and `GeneralPerturbationsOrbit.get_closest_element_index` accept Skyfield `Time` objects using a vectorized conversion.
 - Changed `collect_ground_track` and `collect_ground_pixels` to cull times with a mask by first finding the periods when the instrument's field of regard may observe any part of the mask, and propagating the orbit only at the requested times within them, instead of propagating every requested time and testing each sub-satellite point against a mask buffered in an equidistant cylindrical projection (`buffer_target`). The periods are found from the angular distance between the satellite's direction and the mask (on a sphere, with conservative margins for the Earth's oblateness and the orbit's eccentricity), whose rate of change is bounded by the satellite's angular rate, so that no observation is missed however brief, including over the poles and across the anti-meridian. Previously, footprints over polar regions were dropped (about 31% of those over Antarctica south of 70 degrees, for NOAA-20 VIIRS) and masks across the anti-meridian could be buffered incorrectly. Masks are now split along the anti-meridian and poles (`split_polygon`) before use, so that masks given with longitudes beyond 180 degrees are matched with footprints, and results are sorted by time (clipping to the mask previously reordered rows). For example, a 7-day ground track sampled every 10 s over Arizona takes about 2.4 s.
 - Deprecated utility method `buffer_target`, which does not conservatively bound observations near the poles or across the anti-meridian; it is no longer used and will be removed in a future release.
 - Fixed `split_polygon` for polygons with longitudes continuing beyond 180 (or -180) degrees, such as 170 to 190 degrees, which were returned unchanged, and for polygons spanning exactly -180 to 180 degrees of longitude, such as a polar cap given as a box, which became empty. Polar caps are now detected by the winding of their longitudes around a pole, and edges of exactly 180 degrees of longitude are no longer treated as crossing the anti-meridian.
 - Changed analysis methods that accept satellites (`collect_observations`, `collect_multi_observations`, `collect_downlinks`, `collect_limb_observations`, `collect_ro_observations`, `collect_orbit_track`, `collect_ground_track`, `compute_ground_track`, `collect_ground_pixels`, and `compute_dop`) to raise a `TypeError` for a constellation, which was previously analyzed silently as its lead satellite (whose name, orbit, and instruments it shares); use the constellation's `generate_members` method.
 - Changed `aggregate_observations`, `reduce_observations`, and `reduce_latencies` to group observations by target, its `point_id` and geometry, rather than by `point_id` alone, so that targets with distinct geometries that share an identifier (such as shapely points with the default identifier of 0) are no longer merged: previously, into a single row with a `MultiPoint` geometry and statistics mixed across the targets (for example, a revisit of 4.4 hours for two points revisited every 7.1 and 13.2 hours), which `grid_observations` and `grid_latencies` then dropped, since no cell contained it. Results are unchanged for targets with unique identifiers.
 - Reorganized analysis modules: `tatc.analysis.coverage` is split into `point_coverage` and `region_coverage`, and `tatc.analysis.track` into `orbit_track` and `ground_track`; all public functions remain available from `tatc.analysis`. Orbit propagation algorithms (observation event search, repeat ground track propagation, and repeat cycle search) moved from the general perturbations schema modules to the new internal module `tatc.utils.propagation`.
 - Updated example notebooks: corrected descriptions of observation periods, latency, terrain-masked radar coverage, and dilution of precision, among others; the constellation example now passes the constellation's member satellites; and the limb sounding example plots each scan's position.

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
==================
Validation Summary
==================

The validation notebooks compare TAT-C's analyses with data from operational
missions: recorded orbits, geolocated samples, acquisition catalogs, and
archives of observations, mostly from 2025 and 2026. This page summarizes
their headline results by capability, the modeling choices they found to
matter, and the limits of what TAT-C models. Each result links to the notebook
with its full method, data, and figures.

Unless noted, orbits are propagated with SGP4 from public two-line element
(TLE) sets near the analysis period, and percentiles are of the absolute
differences from the reference data.

Orbits and Propagation
======================

.. list-table::
   :header-rows: 1
   :widths: 20 25 55

   * - Capability
     - Reference
     - Result
   * - TLE propagation (`GeneralPerturbationsOrbit`)
     - Recorded or precise orbits of NOAA-20, Sentinel-1C, Sentinel-2, NISAR, Terra, GCOM-W, Aura, PACE, GRACE-FO
     - 0.5 to 1.7 km along track (95th percentile) within days of the TLE epoch, under 0.5 km across track and radially
       (:doc:`ValidateImagerATMS`, :doc:`ValidatePushbroomMSI`, :doc:`ValidateSpaceGRACEFO`); errors grow with drag:
       29 km after 20 days for Sentinel-2A (:doc:`ValidateRevisitSentinel2`), 300 km after 26 days for GPM
       (:doc:`ValidateYawFlipGPM`), and about 250 km per day before the ISS's last reboost (:doc:`ValidateOrbitISS`)
   * - Maintained orbits (`remove_drag`, `repeat_cycle`)
     - Landsat 9 element set history, 19 months
     - Repeating the first 16-day cycle keeps the ground track within 3 km and the timing within 50 s, and predicts
       99.4% of 8,241 observations; direct propagation diverges within months (:doc:`ValidateMaintainedOrbitLandsat9`)
   * - Repeat cycles (`get_repeat_cycle`)
     - SWOT, OCO-2, ICESat-2
     - SWOT's 20.86-day cycle within 3.1 s (:doc:`ValidateSwathSWOT`); OCO-2's 16 days (:doc:`ValidatePointingOCO2`);
       ICESat-2's 91 days not detected, but refined to within 25 s when declared (:doc:`ValidateLidarICESat2`)
   * - Drifting orbits
     - Terra, since its inclination adjustments ended
     - A single TLE reproduces the local time drift (39.5 vs. 37.3 minutes per year) within 0.3 minutes over three
       months, but departs by 4.4 minutes after a year (:doc:`ValidateDriftMODIS`)
   * - Orbit generators (`SunSynchronousOrbit`, `WalkerConstellation`, `TrainConstellation`, `MolniyaOrbit`,
       `TundraOrbit`, `GeosynchronousOrbit`)
     - 13 sun-synchronous missions; Galileo, BeiDou-3, GLONASS, OneWeb; Sentinel-2, Landsat, NOAA-20/21; Meridian,
       EKS, QZSS, GOES
     - Ascending node within 0.02 minutes of local time; Walker slots within 1 to 3 deg; trains within 29 to 45 km;
       Molniya apogee drift under 0.07 deg per day; Tundra apogee latitudes within 0.03 deg (:doc:`ValidateOrbitGeneration`)
   * - Streets of coverage (`SOCConstellation`)
     - Iridium
     - Reproduces Iridium's 6 planes of 11 satellites and its continuous coverage (:doc:`ValidateConstellationIridium`)

Instruments and Pointing
========================

.. list-table::
   :header-rows: 1
   :widths: 20 25 55

   * - Capability
     - Reference
     - Result
   * - Cross-track scanners (`PointedInstrument`, `ViewGeometry.SCAN`)
     - NOAA-20 ATMS, CrIS, and VIIRS; Terra MODIS
     - More than 99% of passes predicted; epochs within 0.26 s (ATMS) and 2 s (VIIRS) at the 95th percentile, 0.22 s
       (MODIS) median; swath width within 2 to 3 km with measured fields of regard; CrIS footprint centers within
       0.7 km (95th percentile) (:doc:`ValidateImagerATMS`, :doc:`ValidateImagerCrIS`, :doc:`ValidateImagerVIIRS`,
       :doc:`ValidateDriftMODIS`)
   * - Pushbroom imagers (`PointedInstrument`)
     - Sentinel-2 MSI
     - Footprint IoU 0.965 (median, 122 datastrips); epochs within 2 to 3 s with one view, 0.26 s (95th percentile)
       with a view per detector module and a geocentric nadir (:doc:`ValidatePushbroomMSI`)
   * - Conical scanners (`ConicalInstrument`)
     - GCOM-W AMSR2, GPM GMI, SMAP radiometer
     - 99.95 to 100% of passes, epochs within 0.27 s (AMSR2) and 0.69 s (GMI) at the 95th percentile; swath IoU 0.998
       to 0.999; SMAP's forward and aft looks within 0.2 s (median) (:doc:`ValidateConicalAMSR2`,
       :doc:`ValidateConicalGMI`, :doc:`ValidateGridSMAP`)
   * - Side-looking SAR (`PointedInstrument` with roll, `roll_angle_profile`)
     - NISAR, Sentinel-1C
     - Epochs within 0.3 and 0.44 s (95th percentile) of zero-Doppler times; footprint IoU 0.93 to 0.96 (median); roll
       steering profiles bring swath edges within 6 to 7 km (:doc:`ValidateSARNISAR`, :doc:`ValidateSARSentinel1`)
   * - Wide-swath interferometer (two rolled instruments)
     - SWOT KaRIn
     - Every pixel 11 to 59 km from the ground track observed by the instrument on its side, none in the nadir gap;
       epochs within 0.04 s after removing the orbit offset (:doc:`ValidateSwathSWOT`)
   * - Narrow-beam lidar
     - ICESat-2 ATLAS
     - Beams located within about 2 km with a TLE, so a view about 4 km wide predicts 95% of their observations
       (:doc:`ValidateLidarICESat2`)
   * - Geostationary imagers
     - GOES-19 ABI
     - A 17.4 deg field of regard observes the full disk; elevation angles within 0.0004 deg; product coverage limits
       within about 0.3 deg of local zenith angle (:doc:`ValidateGeostationaryGOES`)
   * - Tilted scanners (`pitch_angle_profile`), sun glint (`compute_glint_angle`)
     - PACE OCI
     - Pixels within 0.8 km near nadir and 9 km at the edges (median); epochs within 0.15 s (median) with a daily tilt
       profile; glint angles within 0.26 deg (:doc:`ValidateTiltPACE`)
   * - Attitudes (`Attitude`, `AttitudeSchedule`, `TargetAttitude`, `GlintAttitude`)
     - OCO-2 modes, GPM yaw flips, ISS attitude (EMIT)
     - OCO-2 views within 1.5 km (glint), 2.1 km (nadir), and 5.5 km (target) median under one schedule
       (:doc:`ValidatePointingOCO2`); GMI times relative to the radar within 0.6 s (median) across a yaw flip
       (:doc:`ValidateYawFlipGPM`); EMIT epochs within 1.1 s (median) with the recorded ISS attitude
       (:doc:`ValidateOrbitISS`)
   * - Sample sizes (`compute_sample_sizes`)
     - ATMS, CrIS, VIIRS, MODIS, AMSR2, GMI
     - Sample spacing within 1% across track and 0.2 to 0.8% along track; 55 of 58 published footprint dimensions
       within 3% (:doc:`ValidateSampleSizes`)

Coverage Analyses
=================

.. list-table::
   :header-rows: 1
   :widths: 20 25 55

   * - Capability
     - Reference
     - Result
   * - Revisit (`collect_observations`, `aggregate_observations`, `reduce_observations`)
     - Sentinel-2A, 2B, and 2C acquisitions over 20 days
     - Every observed pass of 2B and 2C predicted; median 4 samples and 5.0-day mean revisit, observed and predicted;
       samples equal at 88% of points (:doc:`ValidateRevisitSentinel2`)
   * - Regions (`collect_region_observations`)
     - Sentinel-2 over six countries; Sentinel-1C over polar regions
     - All 114 and all 210 acquired passes predicted, including regions around both poles and across the
       anti-meridian; sample counts and mean revisit equal to observed for five of six countries; swaths reach 87.58 deg
       north (products: 87.57 deg) (:doc:`ValidateRegionSentinel2`, :doc:`ValidateRegionPolarSentinel1`)
   * - Features (`collect_feature_observations`)
     - TC PRIMED tropical cyclone overpasses, 2025 Atlantic season, 10 satellites
     - All 1,433 archived overpasses predicted; viewed fraction within 0.2 percentage points (median)
       (:doc:`ValidateFeaturesTCPRIMED`)
   * - Gridded statistics (`grid_observations`)
     - SMAP radiometer, one day
     - Same number of passes at 90% of points; gridded samples within 2.4% and revisit equal (median)
       (:doc:`ValidateGridSMAP`)
   * - Illumination (`min_target_solar_elevation`, `omit_solar=False`)
     - Landsat 8 and 9 scenes
     - Solar elevation and azimuth within 0.025 and 0.06 deg, epochs within 0.7 s; a 5 deg minimum solar elevation
       predicts 94 to 95% of acquisitions, 91 to 93% of predictions acquired (:doc:`ValidateSolarLandsat`)
   * - Coincidences (`find_coincident_observations`, `find_coincident_ground_track`)
     - GMI and AMSR2 within 15 minutes
     - With `instantaneous=True`, 97.7% of true coincidences found and 99.9% of those found true, within 2.4 s (95th
       percentile); 98.2% of coincident points within coincident footprints (:doc:`ValidateCoincidenceGMIAMSR2`)

Limb, Occultation, and Satellite-to-Satellite Analyses
======================================================

.. list-table::
   :header-rows: 1
   :widths: 20 25 55

   * - Capability
     - Reference
     - Result
   * - Limb sounding (`collect_limb_observations`)
     - Aura MLS, TIMED SABER
     - MLS tangent points within 1.8 km (median) at matched times and altitudes; SABER scan sample times within
       0.05 s (median); daily latitude ranges within 0.2 deg and equator local times within 0.02 h
       (:doc:`ValidateLimbMLS`, :doc:`ValidateLimbSABER`)
   * - Radio occultation (`collect_ro_observations`)
     - COSMIC-2 and PlanetiQ, GPS, GLONASS, and Galileo
     - All 5,699 and 3,838 recorded occultations predicted with the correct direction; times within 0.25 to 0.30 s
       (95th percentile); tangent points 0.33 to 0.38 km (median) from truth (:doc:`ValidateROCOSMIC2`,
       :doc:`ValidateROPlanetiQ`)
   * - Earth occlusion (`collect_space_observations`)
     - COSMIC-2 occultations; GRACE-FO precise orbits
     - Window edges within 0.25 to 0.30 s (95th percentile) of recorded occultations and within 2 ms of
       `collect_ro_observations`; all GRACE-FO occlusion edges predicted (:doc:`ValidateSpaceCOSMIC2`,
       :doc:`ValidateSpaceGRACEFO`)
   * - Slant range (`collect_space_observations`)
     - GRACE-FO precise orbits; SOCRATES close approaches
     - Range between formation satellites within 21 m (standard deviation); closest approaches within 1 m and 0.3 ms
       (median) of SOCRATES; every window at least `min_duration` long found, down to 0.1 s
       (:doc:`ValidateSpaceGRACEFO`, :doc:`ValidateSpaceSOCRATES`)

Ground Systems and Navigation
=============================

.. list-table::
   :header-rows: 1
   :widths: 20 25 55

   * - Capability
     - Reference
     - Result
   * - Ground radar (`RadarStation`, `compute_terrain_mask`)
     - Four NEXRAD radars in Arizona, Level II volume scans
     - DEM terrain masks agree with observed blockage at 90 to 92% of azimuths; observable area IoU 0.77 to 0.92 per
       station and 0.84 to 0.93 for the network (:doc:`ValidateRadarNEXRAD`)
   * - Downlink latency (`collect_downlinks`, `compute_latencies`)
     - NOAA-20 product granules
     - Mean latency 21 vs. 23 minutes, within 9 minutes per granule (95th percentile), the difference being processing
       time (:doc:`ValidateLatencyNOAA20`)
   * - Dilution of precision (`compute_dop`)
     - GPS, Galileo, BeiDou, and GLONASS broadcast ephemerides at six IGS stations
     - Within 0.01% (median) of an independent computation for the same satellites, including a clock per system
       (:doc:`ValidateDopGNSS`)

Modeling Guidance
=================

The validations found that a few modeling choices matter more than any
numerical error in TAT-C:

* **Velocity frame.** Scanners on spacecraft without yaw steering (NOAA-20,
  Terra, GCOM-W, GPM) scan perpendicular to the *inertial* velocity, which
  `Attitude(velocity_frame=VelocityFrame.INERTIAL)` models; TAT-C's default
  Earth-fixed frame misplaces ATMS's swath edges by up to 90 km along track.
  SAR images (zero-Doppler) and yaw-steered spacecraft (Sentinel-2, Metop)
  follow the default Earth-fixed frame.
* **Measured instrument parameters.** Nominal fields of regard gave swaths 11 km
  too narrow (ATMS) and 47 km too wide (VIIRS); parameters fitted to geolocated
  samples close the gap.
* **Element set age.** Use element sets near the analysis period, especially
  for low orbits with strong drag or frequent maneuvers (ISS, GPM), a
  multi-element `GeneralPerturbationsOrbit` for past periods, and a maintained
  orbit (`remove_drag=True` with a `repeat_cycle`) for analyses of months or
  more. Align the element set switches of formations with
  `align_element_switches`.
* **Predicted versus acquired.** TAT-C predicts what a system *can* observe.
  Acquisition plans, tracking limits, and archive gaps make the recorded
  observations a subset: Sentinel-2A's reduced plan, Landsat's low-Sun and
  night requests, 37 to 81% of predicted COSMIC-2 occultations recorded, and
  missing SSMIS data in TC PRIMED. A minimum target solar elevation of about
  5 deg reproduces the acquisitions of sun-synchronous land imagers.
* **Coincidences of samples.** Observation periods of scanning instruments last
  minutes; compare their epochs with `instantaneous=True`.
* **Navigation constellations.** Model the operational satellites (those with
  healthy broadcast ephemerides); public TLE groups include others, which make
  DOP 4 to 16% optimistic.

Not Modeled
===========

* Refraction: occultations and limb views use straight lines of sight.
* Scan duration: every pixel of a view is observed at one instant, which moves
  the pixels at the edges of a sweep by 1.5 to 4.5 km along track (MODIS, ATMS).
* Non-uniform scan rates: limb scans are at a constant angular rate, so MLS's
  two-speed scan differs by several seconds (tens of kilometers).
* Pointing limits in `collect_space_observations`, such as an antenna's field
  of view, and slews between attitudes, which are instantaneous.
* Acquisition plans, downlink processing time, and data gaps.

Changes to TAT-C from Validation
================================

Several validations found errors or missing capabilities, which were corrected
or added before the results above, including:

* the velocity frame and nadir reference options, scan view geometry (the bow
  tie of cross-track scanners), and `ConicalInstrument`;
* roll and pitch angle profiles, attitudes and attitude schedules, and target
  and glint pointing;
* maintained orbits, declared repeat cycles, and element set switch alignment;
* sun-synchronous ascending nodes (from the apparent Sun to sidereal time),
  train spacing, Molniya and Tundra periods and inclinations, and the polar
  streets-of-coverage pattern;
* refined access periods (rise and set times to a millisecond, grazing passes,
  and fields of regard on the ellipsoid);
* `split_polygon` for regions around the poles, `grid_observations` for points
  with a single sample, and the direction of grazing occultations;
* minimum and maximum target solar elevations, instantaneous coincidences,
  closest approaches, the `during_contact` option of `compute_latencies`, and a
  clock per navigation system in `compute_dop`.

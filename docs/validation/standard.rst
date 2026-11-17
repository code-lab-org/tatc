===================
Validation Standard
===================

The entries of :mod:`tatc.catalog` are held to this standard. It states what an
entry must show before it is called validated, so that entries can be compared
and the claim means the same for each. An entry's tier is computed from its
results and these thresholds, and cannot be declared.

This is version 6 of the standard. The thresholds of version 1 were drafted from
the earlier results below, so they can show little about those instruments. They
were frozen before validating instruments new to TAT-C, which are their first real
test. Changing a threshold requires a new version and a stated reason here.

Versions
========

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - Version
     - Change and reason
   * - 1
     - The first thresholds, for cross-track, pushbroom, conical, SAR, and
       geostationary instruments.
   * - 2
     - Adds the lidar class, for ICESat-2 ATLAS, whose beams fit none of the other
       classes. No threshold of version 1 changed, so results recorded under
       version 1 still apply. The lidar gates were drafted knowing an earlier result
       (:doc:`instruments/ValidateLidarICESat2`: the central beam pair within about 2 km of the
       ground track over Antarctica, and 95% of beam segments within views about
       4 km wide), and frozen before any data were analyzed for the catalog.
   * - 3
     - Adds the altimeter class, for nadir radar altimeters (Sentinel-6 Michael
       Freilich's Poseidon-4), whose measurements lie along the ground track rather
       than in a swath. No threshold of earlier versions changed, so their results
       still apply. The altimeter gates mirror the lidar gates and were frozen
       before any altimeter data were analyzed.
   * - 4
     - Adds the scatterometer class, for real-aperture scatterometers (Metop ASCAT),
       which observe swaths beside a nadir gap but are not SAR. Its gates are those
       of the SAR class, so that a class filter does not mix the two. No threshold
       of earlier versions changed, so their results still apply. The gates were
       adopted before any scatterometer data were analyzed.
   * - 5
     - Adds the regional geostationary class, for geostationary instruments that
       scan a fixed region rather than the Earth's disk (TEMPO), whose fields of
       regard a local zenith angle does not bound. Its gate is the intersection over
       union of the predicted field of regard with the area imaged in a scan, at
       least 0.95, as for swaths. No threshold of earlier versions changed, so their
       results still apply. The gate was frozen before any TEMPO data were analyzed.
   * - 6
     - Adds the profiler class, for single-beam atmospheric profilers near nadir
       (EarthCARE's ATLID lidar and CPR radar), which observe a narrow track, as
       altimeters do, but are neither altimeters nor multi-beam lidars, whose beam
       spacing gate needs several beams. Its gates are those of the altimeter
       class, so that a class filter does not mix the two. No threshold of earlier
       versions changed, so their results still apply. The gates were adopted
       before any profiler data were analyzed.

Rules
=====

1. **Public methods, exact values.** A result counts only if it comes from public
   TAT-C methods using the entry's exact values. A record identifies the
   configuration it tested by a digest of its values, so results of a
   configuration that changed after its test cannot make it `validated`.
2. **Holdout.** Values fitted to reference data are tested on data at least 30
   days from their fitting data. The gap is computed, not declared: it is the
   days between the nearest days of the period of each fitted value's data,
   which its provenance states, and of the test's data, which its record states.
   Values from published specifications or computed from other values need no
   holdout. Results on another platform of the same instrument are recorded as
   transfers, which do not count: using a configuration on an untested platform
   warns.
3. **Geometry, not acquisition.** Passes are scored before solar limits or
   acquisition plans. TAT-C predicts what a system can observe, not what it chose
   to observe.
4. **Provenance.** Each parameter states whether it is `nominal`, `measured`,
   `derived`, or a `placeholder`, and its source. A published specification is
   never reported as measured.
5. **Version.** Each configuration's record states the version of the standard
   of its test. Its results apply under a compatible version that has the gates
   of its class.

Tiers
=====

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Tier
     - Meaning
   * - `validated`
     - All applicable gates pass under the rules above. An entry is validated only
       if every configuration is, and it has configurations at the grades its
       class requires.
   * - `partial`
     - Some gates are unmet or untested, or the evidence is incomplete. The entry
       lists each gap.
   * - `nominal`
     - Published values, without a comparison with reference data.

Gates
=====

A configuration's grade states its use. Coverage and timing differ only in the
epoch gate. Cross-track and conical instruments also need a configuration at the
sample size grade. Orbit error against recorded positions is reported but not
gated, because it grades the element set rather than the instrument. What each
metric measures, and how, is stated under `Metrics`_.

Lidar, altimeter, and profiler gates are separate, because they observe narrow
tracks rather than a swath. Profilers take the altimeter's gates, with the track
offset measured from the predicted track of the view's center, which is the ground
track for a view at nadir; a profiler pointed off nadir (ATLID, a few degrees behind)
observes a track displaced from the ground track by the Earth's rotation:

.. list-table::
   :header-rows: 1
   :widths: 50 25 25

   * - Metric
     - Lidar
     - Altimeter
   * - Segment recall: share of beam segments (or measurements) within the views
       at their times, where the beams follow their reference tracks
     - ≥ 99%
     - ≥ 99%
   * - Beam spacing error, median: distance across track between beams, against
       the spacing of the views
     - ≤ 0.2 km
     -
   * - Track offset, median: distance across track of measurements from the
       predicted ground track
     -
     - ≤ 0.5 km
   * - Epoch, 95th percentile (coverage)
     - ≤ 3 s
     - ≤ 3 s
   * - Epoch, 95th percentile (timing)
     - ≤ 1 s
     - ≤ 1 s

.. list-table::
   :header-rows: 1
   :widths: 34 13 13 13 13 14

   * - Metric
     - Cross-track
     - Pushbroom
     - Conical
     - SAR, scatterometer
     - Geostationary
   * - Pass recall
     - ≥ 99%
     - ≥ 99%
     - ≥ 99%
     - ≥ 99%
     -
   * - Pass precision
     - ≥ 95%
     -
     - ≥ 95%
     -
     -
   * - Epoch, 95th percentile (coverage)
     - ≤ 3 s
     - ≤ 3 s
     - ≤ 3 s
     - ≤ 3 s
     -
   * - Epoch, 95th percentile (timing)
     - ≤ 1 s
     - ≤ 1 s
     - ≤ 1 s
     - ≤ 1 s
     -
   * - Swath width, median
     - ≤ 5 km
     - ≤ 5 km
     - ≤ 5 km
     -
     -
   * - Swath edges, 95th percentile
     -
     -
     -
     - ≤ 10 km
     -
   * - Footprint or swath IoU, median
     - ≥ 0.95
     - ≥ 0.95
     - ≥ 0.95
     - ≥ 0.90
     -
   * - Footprint IoU, 5th percentile
     -
     -
     -
     - ≥ 0.90
     -
   * - Pixel position, median (sample size)
     - ≤ 6 km
     -
     -
     -
     -
   * - Dimensions within 3% (sample size)
     - ≥ 90%
     -
     - ≥ 90%
     -
     -
   * - Product limit of local zenith angle
     -
     -
     -
     -
     - ≤ 0.5 deg
   * - Field of regard IoU with the area imaged in a scan (regional geostationary)
     -
     -
     -
     -
     - ≥ 0.95

Pushbroom imagers, SAR, and lidar have no precision gate because they acquire or
point according to plans, so a predicted observation that was not made is not an
error of geometry. Nor do altimeters, whose products omit measurements by surface
and quality, profilers, whose products omit profiles by quality and mode, or
scatterometers, whose wind products omit land and sea ice.

A swath's IoU is measured at several time steps, up to 120 s, and gated at the
lowest. Configurations whose instruments are swept (`is_swept`), as every
pointed configuration with a view tailored to the time step is, take as their
swath the strip their cross-track line sweeps, as observations of points do, so
the IoU at each time step measures that line, not a view stretched along track
to cover the step. Swept swaths changed no threshold. Records were rescored with
them on the same held-out data, and configurations whose tier changed were tested
again on a later held-out day (:doc:`summary`).

Metrics
=======

Each metric compares a prediction made with public TAT-C methods and the entry's
values against reference data: the instrument's geolocated products, whose positions
and times are taken as true. The catalog's validation scripts compute the metrics,
with shared computations in :mod:`tatc.catalog.metrics`. Unit tests recompute them
from a reduced fixture of each record. A test may score several parts or conditions
separately (instruments, platforms, time steps, or the two edges of a swath); each
is recorded with its scope, and a gate takes the worst value over scopes.

Coverage
--------

**Pass recall** and **pass precision** test whether TAT-C predicts the right
overpasses. Points are spread evenly over the Earth, a few hundred kilometers apart,
on a Fibonacci lattice. A point is observed on a pass when it lies inside the cells
formed by the product's geolocated samples, with one pass for each visit (samples
more than 20 minutes apart, or in another acquisition, belong to separate visits).
TAT-C predicts the observations of the same points over the same period with
:func:`~tatc.analysis.collect_observations`, leaving out the periods without data.
Each observed pass is matched with the nearest prediction for its point within 10
minutes.

* *Recall* is the share of observed passes that were predicted. A missed pass
  means TAT-C would undercount coverage, so every class has this gate.
* *Precision* is the share of predicted passes that were observed. A predicted
  pass that was not observed means TAT-C would overcount coverage. The gate applies
  only to instruments that observe whatever is in view, not to those whose
  acquisitions follow plans or whose products omit observations (see `Gates`_).

**Swath width, median** tests the width of the predicted swath across track. At a
sample of scans, the width the instrument predicts (`get_swath_width`) at the
satellite's altitude at the time of the scan is compared with the width of the
data:

* for cross-track scanners, the length of the scan line through the centers of
  its samples, from one end to the other;
* for conical scanners, the chord between the first and last samples of the
  scanned arc;
* for pushbroom imagers, the width of the valid data across track (of a tile or
  through the center of a scene).

The gate takes the median absolute difference, in km.

**Swath edges, 95th percentile** tests the edges of side-looking swaths (SAR and
scatterometers), whose roll angle may vary around the orbit. At each line across
track of a product's geolocation grid, the rays at the near and far edges of the
view (the roll angle plus or minus half the field of view, at the time of the line)
are projected to the ground. The metric is the horizontal distance from each edge
of the line to its ray. The gate takes the 95th percentile, in km, at whichever of
the near and far edges is worse.

**Footprint or swath IoU** tests the area covered, as the intersection over union
of the predicted and observed areas, weighted by area: 1 when they coincide, and 0
when they are disjoint. The predicted area is the swath of
:func:`~tatc.analysis.compute_ground_track` over the same interval.

* For global swaths (cross-track and conical scanners), the observed area is the
  set of cells of a 0.1 deg grid inside the product's geolocation over one orbit.
  Only cells within 300 km of the data, and away from the start and end of the
  orbit, are compared, so the ends of the interval do not count as errors.
* For pushbroom imagers, SAR, and scatterometers, whose products are tiles, scenes,
  or frames of an acquisition, the observed area is the union of a chunk of
  consecutive products, cut at lines across the track in the middle of the first
  and last product, so that the ragged ends of products do not count as errors. The
  median is taken over the chunks; for SAR and scatterometers, so is the 5th
  percentile, so that a few poor chunks cannot hide behind a good median.

Each is scored at several time steps, as stated above, and gated at the lowest.

**Segment recall** tests narrow tracks (lidar, altimeters, and profilers), the
counterpart of pass recall for views too narrow for a lattice of points. For a
sample of measurements (beam segments, altimeter measurements, or profiles), TAT-C
is asked with :func:`~tatc.analysis.collect_observations` whether the view of the
measurement's beam observes its location within 2 minutes of its time. The metric
is the share observed, at each time step, and the gate takes the lowest. For
ICESat-2, only measurements over Antarctica count, where the beams follow their
reference tracks rather than being pointed off them.

**Beam spacing error, median** tests the distance across track between the beams of
a multi-beam lidar. Each segment's offset across track from the predicted ground
track is measured, and the median offset of each beam is taken in bins of 10 s.
In each bin, the distance between adjacent beams is compared with the distance
between the centers of their views on the ground. The gate takes the median
absolute difference, in km. As a difference of offsets, it leaves out errors of
the element set that move all beams together.

**Track offset, median** tests where a narrow track lies. Each measurement's
distance across track is measured from the predicted track: the ground track for an
altimeter, and the track of the view's center (from `compute_footprint_center`
within 30 s of the measurement) for a profiler. The gate takes the median absolute
distance, in km.

**Product limit of local zenith angle** tests the field of regard of geostationary
imagers that observe the Earth's disk. Their products state a bound of local zenith
angle beyond which retrievals are not of quantitative quality, and the field of
regard is set to that bound. A sample of retrieved pixels within 5 deg of the bound
is checked for whether TAT-C predicts each as observed. The metric is the larger of
the distances from the bound of the greatest zenith angle predicted as observed and
of the smallest predicted as not observed, in degrees: how far the predicted edge
lies from the product's.

**Field of regard IoU with the area imaged in a scan** tests geostationary
instruments that scan a fixed region (TEMPO). The area imaged is the outline of the
outer corners of the outermost pixels of a full scan. On a grid of 0.05 deg around
it, each point is predicted as observed if it is in the field of view from the
satellite's position at the middle of the scan and above its horizon. The metric is
the intersection over union of the two areas, weighted by area.

Timing
------

**Epoch, 95th percentile** tests when TAT-C predicts an observation. It is the
difference between the predicted time of observation of a point and the time the
reference shows the point was imaged, over matched passes or a sample of points.
The reference time is that of the instrument's own geolocation: the time of the
sample nearest the point (interpolated between the scans before and after it, for
scanners), the zero-Doppler time of the grid line for SAR, the imaging time of the
line for pushbroom imagers, and the time of the measurement for lidar, altimeters,
and profilers. The gate takes the 95th percentile of the absolute differences, in
seconds: up to 3 s at the coverage grade, and up to 1 s at the timing grade, for
uses that need the time of each observation.

Sample size
-----------

**Dimensions within 3%** tests the size of samples (footprints) on the ground.
Published footprint dimensions, along and across track at stated scan angles, are
compared with the sample sizes that `compute_sample_sizes` predicts for a satellite
in a circular orbit at the published altitude and inclination, observing at those
scan angles. Each footprint gives two dimensions, and the metric is the share whose
relative error is within 3%.

**Pixel position, median** tests where samples fall within the swath. At a sample of
scans, the position that the instrument projects for each pixel at the middle of
the scan (`compute_projected_pixel_position`) is compared with the geolocated center
of the same sample in the product. The gate takes the median distance, in km.

Back-Test
=========

Results from :doc:`summary` against the thresholds. These come from each notebook
as it was run when this standard was drafted, generally with fitted values and
without holdout, so passing here does not make an instrument `validated`. Some of
these tests have since moved from the notebooks to the catalog's validation
scripts, whose results are on the :doc:`catalog` page. Failures are listed, not explained away:
no threshold was loosened to pass one of them.

.. list-table::
   :header-rows: 1
   :widths: 24 34 22 20

   * - Result
     - Value
     - Gate
     - Outcome
   * - ATMS, VIIRS, MODIS: passes
     - More than 99% predicted
     - Recall ≥ 99%
     - Pass
   * - AMSR2, GMI, SMAP: passes
     - 99.95 to 100%
     - Recall ≥ 99%
     - Pass
   * - ATMS: epochs
     - 0.26 s
     - ≤ 1 s
     - Pass (timing)
   * - VIIRS: epochs
     - 2 s
     - ≤ 3 s, ≤ 1 s
     - Pass (coverage); **fail** (timing)
   * - MODIS: epochs
     - 0.22 s median, no 95th percentile
     - ≤ 1 s
     - Not comparable
   * - AMSR2, GMI: epochs
     - 0.27 s, 0.69 s
     - ≤ 1 s
     - Pass (timing)
   * - MSI: epochs, one view
     - 2 to 3 s
     - ≤ 3 s, ≤ 1 s
     - Pass at the limit (coverage); **fail** (timing)
   * - MSI: epochs, detector views
     - 0.26 s
     - ≤ 1 s
     - Pass (timing)
   * - NISAR, Sentinel-1C: epochs
     - 0.3 s, 0.44 s
     - ≤ 1 s
     - Pass (timing)
   * - Scanners: swath width, measured fields of regard
     - 2 to 3 km
     - ≤ 5 km
     - Pass
   * - ATMS, VIIRS: swath width, nominal fields of regard
     - 11 km, 47 km
     - ≤ 5 km
     - **Fail**, so these need measured values
   * - MSI: footprint IoU
     - 0.965
     - ≥ 0.95
     - Pass
   * - NISAR, Sentinel-1C: footprint IoU
     - 0.93 to 0.96
     - ≥ 0.90
     - Pass, with a lower gate than other classes
   * - Sentinel-1C, NISAR: swath edges, roll profiles
     - 6 to 7 km
     - ≤ 10 km
     - Pass
   * - Conical scanners: swath IoU
     - 0.998 to 0.999
     - ≥ 0.95
     - Pass
   * - CrIS: footprint centers, public methods
     - 19 km (central), 43 km (all), 95th percentile
     - Median ≤ 6 km
     - **Likely fail**; the 0.7 km in the summary needs a custom loop, so it cannot
       be claimed
   * - OCI: pixels
     - 0.8 km near nadir, 9 km at the edges (median)
     - Median ≤ 6 km
     - Pass near nadir; **fail** at the edges
   * - ATMS, CrIS, VIIRS, MODIS, AMSR2, GMI: sample dimensions
     - 55 of 58 within 3% (95%)
     - ≥ 90%
     - Pass
   * - ABI: product limits
     - About 0.3 deg
     - ≤ 0.5 deg
     - Pass
   * - Landsat: passes, 5 deg solar limit
     - 94 to 95%
     - Recall ≥ 99%
     - Not comparable: rule 3 scores passes before solar limits

The metric definitions of some earlier results differ from the gates (for example,
in the percentile reported). Their entries will state values computed to the
definitions here.

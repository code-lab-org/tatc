===================
Validation Standard
===================

The entries of :mod:`tatc.catalog` are held to this standard. It states what an
entry must show before it is called validated, so that entries can be compared
and the claim means the same for each. An entry's tier is computed from its
results and these thresholds, and cannot be declared.

This is version 1 of the standard. Its thresholds were drafted from the earlier
results below, so they can show little about those instruments. They were frozen
before validating instruments new to TAT-C, which are their first real test.
Changing a threshold requires a new version and a stated reason here.

Rules
=====

1. **Public methods, exact values.** A result counts only if it comes from public
   TAT-C methods using the entry's exact values. Results from notebook-only
   refinements, or from rounded or other values, are recorded as `reported` and
   cannot make a configuration `validated`.
2. **Holdout.** Values fitted to reference data are tested on data at least 30
   days from their fitting data. Values from published specifications or computed
   from other values need no holdout. Results on another platform of the same
   instrument do not count: using a configuration on an untested platform warns.
3. **Geometry, not acquisition.** Passes are scored before solar limits or
   acquisition plans. TAT-C predicts what a system can observe, not what it chose
   to observe.
4. **Provenance.** Each parameter states whether it is `nominal`, `measured`,
   `derived`, or a `placeholder`, and its source. A published specification is
   never reported as measured.

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
gated, because it grades the element set rather than the instrument.

.. list-table::
   :header-rows: 1
   :widths: 34 13 13 13 13 14

   * - Metric
     - Cross-track
     - Pushbroom
     - Conical
     - SAR
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

Pushbroom imagers and SAR have no precision gate because they acquire according to
plans, so a predicted pass that was not acquired is not an error of geometry.

Back-Test
=========

Results from :doc:`summary` against the thresholds. These come from each notebook
as it was run, generally with fitted values and without holdout, so passing here
does not make an instrument `validated`. Failures are listed, not explained away:
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

Ground and Orbit Track
======================

Track functions collect the footprints of instruments on observing platforms
at times (one row per platform, instrument, and time), and compute their
aggregated (dissolved) coverage. Each footprint names its platform in the
column of its kind: ``satellite``, ``aircraft``, ``surface_platform``, or
``radar_station``.

.. list-table:: Analysis functions by kind of platform
  :header-rows: 1
  :stub-columns: 1

  * - Analysis
    - Satellite
    - Aircraft
    - Surface platform
    - Radar station
    - Any platforms
  * - Footprints
    - ``collect_ground_track``
    - ``collect_airborne_ground_track``
    - ``collect_surface_track``
    - ``collect_radar_track``
    - ``collect_footprints``
  * - Coverage
    - ``compute_ground_track``
    - ``compute_airborne_ground_track``
    - ``compute_surface_track``
    - ``compute_radar_track``
    - ``dissolve_footprints``
  * - Pixels
    - ``collect_ground_pixels``
    -
    -
    -
    - ``collect_footprints(pixels=True)``
  * - Observations of points
    - ``collect_observations``
    -
    -
    -
    - ``collect_footprint_observations``
  * - Observations of features
    -
    -
    -
    -
    - ``collect_feature_observations``
  * - Coincident footprints
    -
    -
    -
    -
    - ``find_coincident_ground_track``
  * - Coincident observations
    -
    -
    -
    -
    - ``find_coincident_observations``

``collect_observations`` searches a satellite's orbit for the periods when it
observes points; ``collect_footprint_observations`` samples the observations
of points from footprints collected at times, for any platforms. Fixed
platforms (surface platforms and radar stations) give static footprints
without times, or the same footprint at each of the times given.

Any Platforms
-------------

.. autofunction:: tatc.analysis.collect_footprints

.. autofunction:: tatc.analysis.dissolve_footprints

.. autofunction:: tatc.analysis.collect_footprint_observations

Satellites
----------

.. autoclass:: tatc.analysis.OrbitCoordinate
  :members:
  :undoc-members:

.. autoclass:: tatc.analysis.OrbitOutput
  :members:
  :undoc-members:

.. autofunction:: tatc.analysis.collect_orbit_track

.. autofunction:: tatc.analysis.collect_ground_track

.. autofunction:: tatc.analysis.compute_ground_track

.. autofunction:: tatc.analysis.collect_ground_pixels

Aircraft
--------

.. autofunction:: tatc.analysis.collect_airborne_ground_track

.. autofunction:: tatc.analysis.compute_airborne_ground_track

Surface Platforms
-----------------

.. autofunction:: tatc.analysis.collect_surface_track

.. autofunction:: tatc.analysis.compute_surface_track

Radar Stations
--------------

.. autofunction:: tatc.analysis.collect_radar_track

.. autofunction:: tatc.analysis.compute_radar_track

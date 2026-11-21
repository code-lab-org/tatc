=======================
Preprocessing Functions
=======================

Preprocessing functions derive TAT-C model inputs from external datasets and
tools. The terrain and raster functions require the optional `preprocess`
dependencies (``pip install tatc[preprocess]``); the airborne functions read
the outputs of HyPlan, installed with the optional `airborne` dependencies
(``pip install tatc[airborne]``).

Terrain Masks
=============

Terrain masks (:class:`tatc.schemas.TerrainMask`) for radar stations and
surface platforms are derived by ray-tracing a digital elevation model (DEM),
such as the public Copernicus DEM GLO-30.

Copernicus DEM Tiles
--------------------

.. autofunction:: tatc.preprocess.get_copernicus_dem_tile_urls

Sample DEM Elevation
--------------------

.. autofunction:: tatc.preprocess.sample_dem_elevation

Compute Terrain Mask
--------------------

.. autofunction:: tatc.preprocess.compute_terrain_mask

Compute Terrain Mask for Station
--------------------------------

.. autofunction:: tatc.preprocess.compute_terrain_mask_for_station

Compute Terrain Mask for Platform
---------------------------------

.. autofunction:: tatc.preprocess.compute_terrain_mask_for_platform

Terrain Slopes and Incidence Angles
===================================

The terrain slope and aspect at observed targets, from a DEM, give the local
incidence angle of observations on sloped terrain (see
:func:`tatc.utils.compute_local_incidence_angle`), for example, to constrain
side-looking observations of mountainous regions.

Compute Terrain Slope and Aspect
--------------------------------

.. autofunction:: tatc.preprocess.compute_terrain_slope_aspect

Add Local Incidence Angles
--------------------------

.. autofunction:: tatc.preprocess.add_local_incidence_angles

Rasters
=======

Read Raster Weights
-------------------

.. autofunction:: tatc.preprocess.read_raster_weights

Aircraft
========

The states of an aircraft's trajectory (read by
:meth:`tatc.schemas.Trajectory.from_dataframe`) and the instruments it carries
are derived from flight plans, flight recordings, and sensor models, such as
those of `HyPlan <https://github.com/ryanpavlick/hyplan>`_. Altitudes above
mean sea level are converted to altitudes above the WGS 84 ellipsoid with the
EGM2008 geoid (:func:`tatc.utils.compute_geoid_height`) or a constant geoid
height.

Sample Flight Plan
------------------

.. autofunction:: tatc.preprocess.sample_flight_plan

Estimate Attitude
-----------------

.. autofunction:: tatc.preprocess.estimate_attitude

States from IWG1
----------------

.. autofunction:: tatc.preprocess.states_from_iwg1

Instrument from Sensor
----------------------

.. autofunction:: tatc.preprocess.instrument_from_sensor

Trajectory from Waypoints
-------------------------

.. autofunction:: tatc.preprocess.trajectory_from_waypoints

Aircraft from HyPlan
--------------------

.. autofunction:: tatc.preprocess.aircraft_from_hyplan

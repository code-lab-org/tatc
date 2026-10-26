=======================
Preprocessing Functions
=======================

Preprocessing functions derive TAT-C model inputs from external geospatial
datasets. They require the optional `preprocess` dependencies
(``pip install tatc[preprocess]``).

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

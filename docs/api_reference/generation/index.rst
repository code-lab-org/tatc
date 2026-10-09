====================
Generation Functions
====================

Generation functions create distributions of points or cells for spatial analysis.

Points
======

Fibonacci Lattice
-----------------

.. autofunction:: tatc.generation.generate_points_fibonacci_lattice

Equally Spaced
--------------

.. autofunction:: tatc.generation.generate_points_uniform_spacing

Equal Angular Distance
----------------------

.. autofunction:: tatc.generation.generate_points_uniform_angular_distance

Random
------

.. autofunction:: tatc.generation.generate_points_random

Weights for random points can be read from a GeoTIFF raster (requires the
optional `preprocess` dependencies, ``pip install tatc[preprocess]``):

.. autofunction:: tatc.preprocess.read_raster_weights

Cells
=====

Equally Spaced
--------------

.. autofunction:: tatc.generation.generate_cells_uniform_spacing

Equal Angular Spacing
---------------------

.. autofunction:: tatc.generation.generate_cells_uniform_angular_spacing

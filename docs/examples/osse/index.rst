=======================================
Observing System Simulation Experiments
=======================================

These notebooks run observing system simulation experiments (OSSEs): they
sample features of simulated weather (the 7-km GEOS-5 Nature Run, and a 250-m
GEM test frame of the EarthCARE mission) with the swaths of the atmospheric
sounders of current polar-orbiting satellites (ATMS and CrIS) and of notional
future ones. The hurricane and cold pool OSSEs use the observing systems of the
monsoon OSSE, so read that one first:

#. :doc:`RunOSSEMonsoon`: the cold cloud tops of North American Monsoon
   convection over the American Southwest, the horizontal resolution of the
   sounders, and a trade of the satellites and orbit planes of a
   constellation.
#. :doc:`RunOSSEHurricane`: the cold cloud tops and heavy precipitation of a
   simulated North Atlantic hurricane, with the same sensors, and a trade of
   the satellites and orbit planes of a constellation over the storm's two days.
#. :doc:`RunOSSEColdPools`: the convective cold pools of a snapshot of the
   tropical Pacific (the EarthCARE Hawaii test frame), the images of them that
   each sensor would record, and a trade of the satellites and orbit planes of
   a constellation for this scene.

The notebook :doc:`VisualizeSensingGeometry` draws the sensing geometry of
the sounders of these OSSEs: their scans, the size and spacing of their
footprints from nadir to the edge of the swath, and a day of their swaths.

Each OSSE downloads the part of its nature run that it needs once, to a
``data`` directory beside it. The notebooks require::

  pip install tatc[osse]

.. toctree::
  :hidden:

  RunOSSEMonsoon.ipynb
  RunOSSEHurricane.ipynb
  RunOSSEColdPools.ipynb
  VisualizeSensingGeometry.ipynb

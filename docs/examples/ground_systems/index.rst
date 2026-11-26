=============================
Ground Systems and Navigation
=============================

These notebooks analyze observing systems on the surface, navigation
constellations, and the ground systems that receive satellites' data:

#. :doc:`CollectRadarTrack`: the coverage of ground-based radars (NEXRAD,
   around Phoenix, Arizona), with masks of terrain blockage, and a comparison
   of radar bands.
#. :doc:`CollectSurfaceTrack`: the coverage of instruments on towers and
   tethered balloons.
#. :doc:`ComputeDOP`: the dilution of precision of a GPS-like constellation
   over a day, compared with STK.
#. :doc:`ComputeDataFlow`: the latency of a constellation's observations
   through the flow of its data, with link budgets, shared antennas of KSAT's
   stations, and a TDRS relay.

The masks of terrain blockage that the radar and platform examples read are
beside them.

.. toctree::
  :hidden:

  CollectRadarTrack.ipynb
  CollectSurfaceTrack.ipynb
  ComputeDOP.ipynb
  ComputeDataFlow.ipynb

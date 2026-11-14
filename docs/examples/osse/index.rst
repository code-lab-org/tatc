=======================================
Observing System Simulation Experiments
=======================================

These notebooks run observing system simulation experiments (OSSEs): they
sample features of the simulated weather of the 7-km GEOS-5 Nature Run with
the swaths of the atmospheric sounders of current polar-orbiting satellites
(ATMS and CrIS) and of notional future ones. The hurricane OSSE uses the
observing systems of the monsoon OSSE, so read that one first:

#. :doc:`RunOSSEMonsoon`: the cold cloud tops of North American Monsoon
   convection over the American Southwest, the horizontal resolution of the
   sounders, and a trade of the satellites and orbit planes of a
   constellation.
#. :doc:`RunOSSEHurricane`: the cold cloud tops and heavy precipitation of a
   simulated North Atlantic hurricane, with the same sensors.

Each notebook downloads the part of the Nature Run that it needs once, to a
``data`` directory beside it. The notebooks require::

  pip install tatc[osse]

.. toctree::
  :hidden:

  RunOSSEMonsoon.ipynb
  RunOSSEHurricane.ipynb

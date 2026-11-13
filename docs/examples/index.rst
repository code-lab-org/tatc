========
Examples
========

The following provide examples of TAT-C use in a Jupyter notebook environment.

.. toctree::
  :maxdepth: 1

  CollectObservations.ipynb
  GenerateRandomPoints.ipynb
  CollectOrbitGroundTrack.ipynb
  CollectRadarTrack.ipynb
  CollectPlatformTrack.ipynb
  CollectAirborneGroundTrackP3B.ipynb
  ComputeCoverage.ipynb
  FindCoincidentObservations.ipynb
  ComputeLatency.ipynb
  ComputeDOP.ipynb
  CollectRO.ipynb
  CollectLimb.ipynb
  CollectSpace.ipynb
  RunOSSEMonsoon.ipynb
  RunOSSEHurricane.ipynb

Note that running the above examples in a local Python environment requires additional dependencies which can be installed via::

  pip install tatc[examples]

using the pip package manager. The airborne example (``CollectAirborneGroundTrackP3B.ipynb``) additionally requires the HyPlan flight planning library, installed via::

  pip install hyplan

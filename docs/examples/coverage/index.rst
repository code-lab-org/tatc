=================
Coverage Analyses
=================

These notebooks choose where to sample and analyze what observations
sample:

#. :doc:`GenerateRandomPoints`: generate random points, uniformly or weighted
   by a population grid, and find when an imaging satellite observes them.
#. :doc:`FindCoincidentObservations`: find when two or more satellites
   observe the same place at nearly the same time, as points and as ground
   tracks.
#. :doc:`CombineObservingPlatforms`: plan an aircraft underflight of a
   satellite's overpass, and find where and when the satellite, the
   aircraft, and a ground-based weather radar observe the same places and a
   storm cell, from the footprints of all three together.
#. :doc:`AssessSamplingSufficiency`: assess whether a constellation samples a
   varying field sufficiently, from the error of its estimates, its revisits,
   and the times of day it sees, and test the estimates with a synthetic
   nature run.

The AssessSamplingSufficiency example also requires xarray
(``pip install xarray``) for its nature run.

.. toctree::
  :hidden:

  GenerateRandomPoints.ipynb
  FindCoincidentObservations.ipynb
  CombineObservingPlatforms.ipynb
  AssessSamplingSufficiency.ipynb

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
  AssessSamplingSufficiency.ipynb

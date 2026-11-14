========
Examples
========

The following provide examples of TAT-C use in a Jupyter notebook environment.

.. toctree::
  :maxdepth: 1

  CollectObservations.ipynb
  UseInstrumentCatalog.ipynb
  GenerateRandomPoints.ipynb
  CollectOrbitGroundTrack.ipynb
  CollectRadarTrack.ipynb
  CollectPlatformTrack.ipynb
  ComputeCoverage.ipynb
  FindCoincidentObservations.ipynb
  ComputeLatency.ipynb
  ComputeDOP.ipynb
  CollectRO.ipynb
  CollectLimb.ipynb
  CollectSpace.ipynb
  AssessSamplingSufficiency.ipynb
  osse/index
  tradespace/index

Note that running the above examples in a local Python environment requires additional dependencies which can be installed via::

  pip install tatc[examples]

using the pip package manager. The AssessSamplingSufficiency example also requires xarray (``pip install xarray``) for its nature run. The trade space tutorial also requires::

  pip install tatc[tradespace]

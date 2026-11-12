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
  RunOSSEMonsoon.ipynb
  RunOSSEHurricane.ipynb
  SearchTradeSpace.ipynb
  RankUncertainDesigns.ipynb
  CompareScenarios.ipynb
  TradeCostAndUtility.ipynb
  PlanMissions.ipynb

Note that running the above examples in a local Python environment requires additional dependencies which can be installed via::

  pip install tatc[examples]

using the pip package manager. The trade space examples also require::

  pip install tatc[tradespace]

The trade space examples after the first (SearchTradeSpace) import its design space and study from
``trade_space_example.py``, beside the notebooks.

===============
Getting Started
===============

These notebooks introduce the analysis functions of TAT-C with one satellite,
NOAA-20, and its Visible Infrared Imaging Radiometer Suite (VIIRS). Read them
in order:

#. :doc:`CollectObservations`: define a satellite and its instrument, and
   collect its observations of a point.
#. :doc:`UseInstrumentCatalog`: take instruments of past and present missions
   from the catalog, with the attitude they were tested with, and read how
   well each was tested.
#. :doc:`CollectOrbitGroundTrack`: collect the orbit and ground track of a
   satellite, with a requirement that targets be sunlit, and compare a
   detailed model of VIIRS with the simple one.
#. :doc:`ComputeCoverage`: compute the coverage of a region of interest.
#. :doc:`ComputeLatency`: compute the latency of observations downlinked to
   ground stations.

.. toctree::
  :hidden:

  CollectObservations.ipynb
  UseInstrumentCatalog.ipynb
  CollectOrbitGroundTrack.ipynb
  ComputeCoverage.ipynb
  ComputeLatency.ipynb

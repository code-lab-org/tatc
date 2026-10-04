==========
Validation
==========

The following compare TAT-C analysis results with reference data from operational missions.

.. toctree::
  :maxdepth: 1

  ValidateLimbMLS.ipynb
  ValidateLimbSABER.ipynb
  ValidateROCOSMIC2.ipynb
  ValidateROPlanetiQ.ipynb
  ValidateRadarNEXRAD.ipynb

These notebooks download reference data on first run (the MLS data require a free `NASA Earthdata Login <https://urs.earthdata.nasa.gov/>`_ account) and require additional dependencies, which can be installed via::

  pip install tatc[validation]

using the pip package manager.

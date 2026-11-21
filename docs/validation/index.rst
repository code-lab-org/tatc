==========
Validation
==========

The following compare TAT-C analysis results with reference data from operational missions, or with published studies. The :doc:`summary` collects their headline results by capability, with modeling guidance and the limits of what TAT-C models. The :doc:`catalog` gives instruments of operational missions, ready to use, with their results against the :doc:`standard`, and the :doc:`systems` gives the satellite programs, radar networks, and surface networks of baseline systems, with their results against the system standard.

.. toctree::
  :maxdepth: 1

  summary
  catalog
  systems
  standard

The validation notebooks are grouped by capability, as in the :doc:`summary`. Those on instruments and pointing derive the configurations of several catalog entries and explain the modeling choices behind them; the entries' results on data held out from their fits are on the :doc:`catalog` page.

.. toctree::
  :maxdepth: 2
  :caption: Notebooks

  orbits/index
  instruments/index
  coverage/index
  limb_and_space/index
  ground_systems/index
  trade_studies/index

These notebooks download reference data on first run, except ``ValidateSampleSizes``, which reuses the data of the imager and conical scanner notebooks, and ``ValidateSpaceCOSMIC2``, which reuses the data of ``ValidateROCOSMIC2``. All of them cache their downloads in one ``data`` directory in ``docs/validation``, which they share (for example, element sets retrieved from CelesTrak). ``ValidateSpaceSOCRATES`` includes a sample of a SOCRATES run (data from CelesTrak). ``ValidateCostVarianceNSOSA`` needs no data: its reference values are those published in the paper it reproduces. The MLS, NISAR, AMSR2, GMI, SWOT, ICESat-2, ECOSTRESS, PACE, MODIS, OCO-2, GPM, GRACE-FO, SMAP, and IMERG data require a free `NASA Earthdata Login <https://urs.earthdata.nasa.gov/>`_ account, and the orbit histories of the Landsat 9, TC PRIMED, and GRACE-FO notebooks a free `Space-Track.org <https://www.space-track.org/>`_ account. The notebooks require additional dependencies, which can be installed via::

  pip install tatc[validation]

using the pip package manager.

The scripts that validate catalog entries and write their records are in the ``validation/catalog`` folder of the repository; see its `README <https://github.com/code-lab-org/tatc/blob/main/validation/README.md>`_.

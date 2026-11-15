==========
Validation
==========

The following compare TAT-C analysis results with reference data from operational missions, or with published studies. The :doc:`summary` collects their headline results by capability, with modeling guidance and the limits of what TAT-C models. The :doc:`catalog` gives instruments of operational missions, ready to use, with their results against the :doc:`standard`.

.. toctree::
  :maxdepth: 1

  summary
  catalog
  standard

The notebooks below are grouped as in the :doc:`summary`. Those on instruments and pointing derive the configurations of several catalog entries and explain the modeling choices behind them; the entries' results on data held out from their fits are on the :doc:`catalog` page.

.. toctree::
  :maxdepth: 1
  :caption: Orbits and Propagation

  ValidateOrbitGeneration.ipynb
  ValidateMaintainedOrbitLandsat9.ipynb
  ValidateOrbitISS.ipynb
  ValidateConstellationIridium.ipynb

.. toctree::
  :maxdepth: 1
  :caption: Instruments and Pointing

  ValidateImagerATMS.ipynb
  ValidateImagerCrIS.ipynb
  ValidateImagerVIIRS.ipynb
  ValidateDriftMODIS.ipynb
  ValidatePushbroomMSI.ipynb
  ValidateConicalAMSR2.ipynb
  ValidateConicalGMI.ipynb
  ValidateSARNISAR.ipynb
  ValidateSARSentinel1.ipynb
  ValidateSwathSWOT.ipynb
  ValidateLidarICESat2.ipynb
  ValidateTiltPACE.ipynb
  ValidatePointingOCO2.ipynb
  ValidateYawFlipGPM.ipynb
  ValidateSampleSizes.ipynb

.. toctree::
  :maxdepth: 1
  :caption: Coverage Analyses

  ValidateRevisitSentinel2.ipynb
  ValidateRegionSentinel2.ipynb
  ValidateRegionPolarSentinel1.ipynb
  ValidateFeaturesTCPRIMED.ipynb
  ValidateGridSMAP.ipynb
  ValidateSamplingErrorIMERG.ipynb
  ValidateSolarLandsat.ipynb
  ValidateCoincidenceGMIAMSR2.ipynb

.. toctree::
  :maxdepth: 1
  :caption: Limb, Occultation, and Satellite-to-Satellite

  ValidateLimbMLS.ipynb
  ValidateLimbSABER.ipynb
  ValidateROCOSMIC2.ipynb
  ValidateROPlanetiQ.ipynb
  ValidateSpaceCOSMIC2.ipynb
  ValidateSpaceGRACEFO.ipynb
  ValidateSpaceSOCRATES.ipynb

.. toctree::
  :maxdepth: 1
  :caption: Ground Systems and Navigation

  ValidateRadarNEXRAD.ipynb
  ValidateLatencyNOAA20.ipynb
  ValidateDopGNSS.ipynb

.. toctree::
  :maxdepth: 1
  :caption: Trade Studies

  ValidateCostVarianceNSOSA.ipynb

These notebooks download reference data on first run, except ``ValidateSampleSizes``, which reuses the data of the imager and conical scanner notebooks, and ``ValidateSpaceCOSMIC2``, which reuses the data of ``ValidateROCOSMIC2``. ``ValidateSpaceSOCRATES`` includes a sample of a SOCRATES run (data from CelesTrak). ``ValidateCostVarianceNSOSA`` needs no data: its reference values are those published in the paper it reproduces. The MLS, NISAR, AMSR2, GMI, SWOT, ICESat-2, ECOSTRESS, PACE, MODIS, OCO-2, GPM, GRACE-FO, SMAP, and IMERG data require a free `NASA Earthdata Login <https://urs.earthdata.nasa.gov/>`_ account, and the orbit histories of the Landsat 9, TC PRIMED, and GRACE-FO notebooks a free `Space-Track.org <https://www.space-track.org/>`_ account. The notebooks require additional dependencies, which can be installed via::

  pip install tatc[validation]

using the pip package manager.

The scripts that validate catalog entries and write their records are in the ``validation/catalog`` folder of the repository; see its `README <https://github.com/code-lab-org/tatc/blob/main/validation/README.md>`_.

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
  ValidateSpaceCOSMIC2.ipynb
  ValidateSpaceGRACEFO.ipynb
  ValidateRadarNEXRAD.ipynb
  ValidateImagerATMS.ipynb
  ValidateImagerCrIS.ipynb
  ValidateImagerVIIRS.ipynb
  ValidateSARSentinel1.ipynb
  ValidateSARNISAR.ipynb
  ValidatePushbroomMSI.ipynb
  ValidateRevisitSentinel2.ipynb
  ValidateRegionSentinel2.ipynb
  ValidateConicalAMSR2.ipynb
  ValidateConicalGMI.ipynb
  ValidateOrbitGeneration.ipynb
  ValidateDopGNSS.ipynb
  ValidateLatencyNOAA20.ipynb
  ValidateSolarLandsat.ipynb
  ValidateSwathSWOT.ipynb
  ValidateLidarICESat2.ipynb
  ValidateConstellationIridium.ipynb
  ValidateGeostationaryGOES.ipynb
  ValidateOrbitISS.ipynb
  ValidateMaintainedOrbitLandsat9.ipynb
  ValidateTiltPACE.ipynb
  ValidateDriftMODIS.ipynb
  ValidatePointingOCO2.ipynb
  ValidateYawFlipGPM.ipynb
  ValidateSampleSizes.ipynb
  ValidateFeaturesTCPRIMED.ipynb

These notebooks download reference data on first run, except ``ValidateSampleSizes``, which reuses the data of the imager and conical scanner notebooks, and ``ValidateSpaceCOSMIC2``, which reuses the data of ``ValidateROCOSMIC2``. The MLS, NISAR, AMSR2, GMI, SWOT, ICESat-2, ECOSTRESS, PACE, MODIS, OCO-2, GPM, and GRACE-FO data require a free `NASA Earthdata Login <https://urs.earthdata.nasa.gov/>`_ account, and the orbit histories of the Landsat 9, TC PRIMED, and GRACE-FO notebooks a free `Space-Track.org <https://www.space-track.org/>`_ account. The notebooks require additional dependencies, which can be installed via::

  pip install tatc[validation]

using the pip package manager.

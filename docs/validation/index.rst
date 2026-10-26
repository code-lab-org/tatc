==========
Validation
==========

The following compare TAT-C analysis results with reference data from operational missions. The :doc:`summary` collects their headline results by capability, with modeling guidance and the limits of what TAT-C models.

.. toctree::
  :maxdepth: 1

  summary
  ValidateLimbMLS.ipynb
  ValidateLimbSABER.ipynb
  ValidateROCOSMIC2.ipynb
  ValidateROPlanetiQ.ipynb
  ValidateSpaceCOSMIC2.ipynb
  ValidateSpaceGRACEFO.ipynb
  ValidateSpaceSOCRATES.ipynb
  ValidateRadarNEXRAD.ipynb
  ValidateImagerATMS.ipynb
  ValidateImagerCrIS.ipynb
  ValidateImagerVIIRS.ipynb
  ValidateSARSentinel1.ipynb
  ValidateSARNISAR.ipynb
  ValidatePushbroomMSI.ipynb
  ValidateRevisitSentinel2.ipynb
  ValidateRegionSentinel2.ipynb
  ValidateRegionPolarSentinel1.ipynb
  ValidateConicalAMSR2.ipynb
  ValidateConicalGMI.ipynb
  ValidateCoincidenceGMIAMSR2.ipynb
  ValidateGridSMAP.ipynb
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

These notebooks download reference data on first run, except ``ValidateSampleSizes``, which reuses the data of the imager and conical scanner notebooks, and ``ValidateSpaceCOSMIC2``, which reuses the data of ``ValidateROCOSMIC2``. ``ValidateSpaceSOCRATES`` includes a sample of a SOCRATES run (data from CelesTrak). The MLS, NISAR, AMSR2, GMI, SWOT, ICESat-2, ECOSTRESS, PACE, MODIS, OCO-2, GPM, GRACE-FO, and SMAP data require a free `NASA Earthdata Login <https://urs.earthdata.nasa.gov/>`_ account, and the orbit histories of the Landsat 9, TC PRIMED, and GRACE-FO notebooks a free `Space-Track.org <https://www.space-track.org/>`_ account. The notebooks require additional dependencies, which can be installed via::

  pip install tatc[validation]

using the pip package manager.

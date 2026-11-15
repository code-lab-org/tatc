=============================================
Limb, Occultation, and Satellite-to-Satellite
=============================================

These notebooks analyze observations that look through the atmosphere or
across space rather than down at the surface. Each builds on the one before
it:

#. :doc:`CollectRO`: GNSS radio occultations, at the tangent point between a
   receiver and a transmitter. The element sets of the GNSS constellations
   it reads are in the ``elements`` folder beside it.
#. :doc:`CollectLimb`: the tangent points of a limb sounder, which scans in
   elevation from its own spacecraft.
#. :doc:`CollectSpace`: the periods when one satellite can observe another,
   with the line of sight between them clear of the Earth, including
   satellites observing geostationary ones.

.. toctree::
  :hidden:

  CollectRO.ipynb
  CollectLimb.ipynb
  CollectSpace.ipynb

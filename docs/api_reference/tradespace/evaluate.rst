=======
Studies
=======

Trade Studies
=============

.. autoclass:: tatc.tradespace.Study
   :members: run, realize, replace, prune, workers, definition, fingerprint

Observation Studies
===================

.. automodule:: tatc.tradespace.observations
   :no-members:

.. autoclass:: tatc.tradespace.ObservationStudy
   :members: architecture_of, satellites_of, targets_of, observations_of, evaluation_of, mission_of, prune

.. autoclass:: tatc.tradespace.Observations
   :members: start_of, end_of, evaluate, prune, __call__

.. autoclass:: tatc.tradespace.Targets

The parts of observation studies read the values of designs by references
(see :doc:`define` for :py:data:`tatc.tradespace.references.Reference` and
:py:func:`tatc.tradespace.value_of`).

Missions of Designs
===================

.. automodule:: tatc.tradespace.plans
   :no-members:

.. autoclass:: tatc.tradespace.MissionPlan
   :members: plan_mission

.. autoclass:: tatc.tradespace.Failures

.. autoclass:: tatc.tradespace.mission.Mission
   :members: phases, estimate, histories, capability, requirements

Stores of Results
=================

.. autoclass:: tatc.tradespace.ResultStore
   :members: clear

Caches of Satellites
====================

.. automodule:: tatc.tradespace.cache
   :no-members:

.. autoclass:: tatc.tradespace.SatelliteCache
   :members: __call__, prune, clear_memory, folder, fingerprint, hits, misses

.. autofunction:: tatc.tradespace.cache.get_satellite_key

Fingerprints
============

.. automodule:: tatc.tradespace.fingerprint
   :no-members:

.. autofunction:: tatc.tradespace.fingerprint.get_fingerprint

.. autofunction:: tatc.tradespace.fingerprint.get_definition

.. autofunction:: tatc.tradespace.fingerprint.compare_definitions

.. autofunction:: tatc.tradespace.fingerprint.describe

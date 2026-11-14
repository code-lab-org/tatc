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
   :members: satellites_of, targets_of, observations_of, deployment_of, phases_of, estimate_of

.. autoclass:: tatc.tradespace.Observations
   :members: end, start_of, end_of, cache_of, context_of, prune, __call__

.. autoclass:: tatc.tradespace.Targets
   :members: of, seed_of

.. autoclass:: tatc.tradespace.MissionPlan
   :members: schedule_of, cost_model_of, masses_of, growth_of, stages_of

.. autoclass:: tatc.tradespace.Failures

.. autodata:: tatc.tradespace.observations.Reference
   :annotation:

.. autofunction:: tatc.tradespace.observations.value_of

Stores of Results
=================

.. autoclass:: tatc.tradespace.ResultStore
   :members: get, put, define, definition, contexts, remove, clear

Caches of Satellites
====================

.. automodule:: tatc.tradespace.cache
   :no-members:

.. autoclass:: tatc.tradespace.SatelliteCache
   :members: get, prune, clear_memory, folder, __call__

.. autofunction:: tatc.tradespace.cache.get_satellite_key

Fingerprints
============

.. automodule:: tatc.tradespace.fingerprint
   :no-members:

.. autofunction:: tatc.tradespace.fingerprint.get_fingerprint

.. autofunction:: tatc.tradespace.fingerprint.get_definition

.. autofunction:: tatc.tradespace.fingerprint.compare_definitions

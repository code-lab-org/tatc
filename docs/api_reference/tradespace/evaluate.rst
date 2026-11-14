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
   :members: satellites_of, targets_of, observations_of, mission_of, prune

.. autoclass:: tatc.tradespace.Observations
   :members: end, start_of, end_of, cache_of, context_of, prune, __call__

.. autoclass:: tatc.tradespace.Targets
   :members: of, seed_of

.. autodata:: tatc.tradespace.references.Reference
   :annotation:

.. autofunction:: tatc.tradespace.references.value_of

Missions of Designs
===================

.. automodule:: tatc.tradespace.plans
   :no-members:

.. autoclass:: tatc.tradespace.MissionPlan
   :members: schedule_of, cost_model_of, masses_of, growth_of, stages_of, policy_of, mission_of

.. autoclass:: tatc.tradespace.Failures
   :members: seed_of, simulate, evaluate, summarize

.. autoclass:: tatc.tradespace.plans.Mission
   :members: planned, phases, estimate, histories, metrics, summarize

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

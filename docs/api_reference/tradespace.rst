Trade Space
===========

.. automodule:: tatc.tradespace

Design Spaces
-------------

.. autoclass:: tatc.tradespace.DesignSpace
   :members: enumerate, size, sample, design, suggest, unit_paths, encode, decode

.. autoclass:: tatc.tradespace.Variable
   :members: domain

.. autoclass:: tatc.tradespace.Real

.. autoclass:: tatc.tradespace.Integer

.. autoclass:: tatc.tradespace.Categorical

.. autoclass:: tatc.tradespace.Fixed

.. autoclass:: tatc.tradespace.Derived

.. autoclass:: tatc.tradespace.Uncertain
   :members: uniform, normal, choice, seed

.. autoclass:: tatc.tradespace.Choice

.. autoclass:: tatc.tradespace.Repeat

.. autoexception:: tatc.tradespace.Infeasible

Designs
-------

.. autoclass:: tatc.tradespace.Design
   :members: flat, realize, values

.. autoclass:: tatc.tradespace.Values

.. autoclass:: tatc.tradespace.ChoiceValues
   :members: choice, value

.. autoclass:: tatc.tradespace.RepeatValues
   :members: count

.. autoclass:: tatc.tradespace.Scope

Constellations
--------------

.. automodule:: tatc.tradespace.constellations
   :no-members:

.. autofunction:: tatc.tradespace.walker_space

.. autofunction:: tatc.tradespace.sso_planes_space

.. autofunction:: tatc.tradespace.planes_space

.. autofunction:: tatc.tradespace.build_constellation

.. autofunction:: tatc.tradespace.count_satellites

.. autofunction:: tatc.tradespace.count_planes

.. autofunction:: tatc.tradespace.hours_to_time

.. autofunction:: tatc.tradespace.remove_failures

Trade Studies
-------------

.. autoclass:: tatc.tradespace.Study
   :members: run, workers, fingerprint

.. autofunction:: tatc.tradespace.summarize

Statistics and Requirements
---------------------------

.. automodule:: tatc.tradespace.criteria
   :no-members:

.. autofunction:: tatc.tradespace.aggregate

.. autoclass:: tatc.tradespace.Statistic
   :members: reduce

.. autoclass:: tatc.tradespace.Mean

.. autoclass:: tatc.tradespace.StandardDeviation

.. autoclass:: tatc.tradespace.Minimum

.. autoclass:: tatc.tradespace.Maximum

.. autoclass:: tatc.tradespace.Quantile

.. autoclass:: tatc.tradespace.Probability

.. autoclass:: tatc.tradespace.Requirement
   :members: violation

Stores of Results
-----------------

.. autoclass:: tatc.tradespace.ResultStore
   :members: get, put, contexts, clear

.. autofunction:: tatc.tradespace.get_fingerprint

Caches
------

.. automodule:: tatc.tradespace.cache
   :no-members:

.. autoclass:: tatc.tradespace.SatelliteCache
   :members: get, clear, __call__

.. autofunction:: tatc.tradespace.get_satellite_key

Pareto Efficiency
-----------------

.. autofunction:: tatc.tradespace.pareto_front

.. autofunction:: tatc.tradespace.pareto_ranks

.. autofunction:: tatc.tradespace.hypervolume

.. autofunction:: tatc.tradespace.pareto_probability

Adaptive Replication
--------------------

.. autofunction:: tatc.tradespace.refine

.. autoclass:: tatc.tradespace.Refinement
   :members: uncertain

Confirmation
------------

.. autofunction:: tatc.tradespace.confirm

.. autoclass:: tatc.tradespace.Confirmation
   :members: optimism

Sensitivity
-----------

.. autofunction:: tatc.tradespace.sensitivity

.. autofunction:: tatc.tradespace.main_effects

Screening
---------

.. autofunction:: tatc.tradespace.screen

.. autoclass:: tatc.tradespace.Screening
   :members: final

Optimizers
----------

.. automodule:: tatc.tradespace.optimization
   :no-members:

.. autofunction:: tatc.tradespace.optimize_optuna

.. autofunction:: tatc.tradespace.optimize_pymoo

.. autofunction:: tatc.tradespace.make_pymoo_problem

.. autoclass:: tatc.tradespace.Evaluations
   :members: evaluate, results, worst

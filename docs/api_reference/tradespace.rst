Trade Space
===========

.. automodule:: tatc.tradespace

Design Spaces
-------------

.. autoclass:: tatc.tradespace.DesignSpace
   :members: enumerate, size, sample, design, suggest, unit_paths, encode, decode

.. autoclass:: tatc.tradespace.Real

.. autoclass:: tatc.tradespace.Integer

.. autoclass:: tatc.tradespace.Categorical

.. autoclass:: tatc.tradespace.Fixed

.. autoclass:: tatc.tradespace.Derived

.. autoclass:: tatc.tradespace.Uncertain
   :members: uniform, normal, choice

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

Trade Studies
-------------

.. autoclass:: tatc.tradespace.Study
   :members: run, workers

.. autofunction:: tatc.tradespace.summarize

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
   :members: evaluate, results

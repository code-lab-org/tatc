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

.. autofunction:: tatc.tradespace.configure

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

Deployment
----------

.. automodule:: tatc.tradespace.deployment
   :no-members:

.. autoclass:: tatc.tradespace.DeploymentSchedule
   :members: plan, launches_of

.. autoclass:: tatc.tradespace.Deployment
   :members: simulate, simulate_many, expected_metrics, planned, phases, begin, end, deployed, units

.. autoclass:: tatc.tradespace.Launch

.. autoclass:: tatc.tradespace.Production

Failures and Missions
---------------------

.. automodule:: tatc.tradespace.failures
   :no-members:

.. autoclass:: tatc.tradespace.FailureModel
   :members: life_of

.. autoclass:: tatc.tradespace.Lifetime
   :members: survival, hazard, cumulative_hazard, sample, from_uniform, dimensions

.. autoclass:: tatc.tradespace.Exponential

.. autoclass:: tatc.tradespace.Weibull

.. autoclass:: tatc.tradespace.Bathtub

.. autoclass:: tatc.tradespace.Piecewise

.. autoclass:: tatc.tradespace.RestorationPolicy
   :members: ground_spares_of

.. autofunction:: tatc.tradespace.policy_space

.. automodule:: tatc.tradespace.mission
   :no-members:

.. autoclass:: tatc.tradespace.MissionHistory
   :members: phases, metrics, in_orbit

.. autoclass:: tatc.tradespace.Phase
   :members: duration

.. autofunction:: tatc.tradespace.evaluate_phases

.. autofunction:: tatc.tradespace.summarize_phases

.. autoclass:: tatc.tradespace.PhaseCache
   :members: evaluate, clear

.. autofunction:: tatc.tradespace.evaluate_missions

.. autofunction:: tatc.tradespace.summarize_missions

.. autofunction:: tatc.tradespace.allocate_histories

Costs
-----

.. automodule:: tatc.tradespace.cost
   :no-members:

.. autoclass:: tatc.tradespace.CostModel
   :members: estimate, discount

.. autoclass:: tatc.tradespace.SatelliteCost
   :members: production

.. autoclass:: tatc.tradespace.LaunchCost

.. autoclass:: tatc.tradespace.OperationsCost

.. autoclass:: tatc.tradespace.CostEstimate
   :members: total, present_value, annual, metrics

.. autofunction:: tatc.tradespace.plane_sizes

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

.. autofunction:: tatc.tradespace.get_source_hash

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

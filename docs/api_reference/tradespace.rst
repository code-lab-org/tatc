Trade Space
===========

.. automodule:: tatc.tradespace

Design Spaces
-------------

.. autoclass:: tatc.tradespace.DesignSpace
   :members: enumerate, size, sample, design, suggest, unit_paths, uncertain_paths, encode, decode

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

.. autofunction:: tatc.tradespace.plane_sizes

.. autoclass:: tatc.tradespace.Layout
   :members: planes, orbits, names, satellites, build

.. autoclass:: tatc.tradespace.Plane

.. autofunction:: tatc.tradespace.layout_of

.. autoclass:: tatc.tradespace.constellations.WalkerLayout

.. autoclass:: tatc.tradespace.constellations.SSOPlanesLayout

.. autoclass:: tatc.tradespace.constellations.PlanesLayout

.. autoclass:: tatc.tradespace.constellations.Layouts

.. autofunction:: tatc.tradespace.hours_to_time

.. autofunction:: tatc.tradespace.remove_failures

Deployment
----------

.. automodule:: tatc.tradespace.deployment
   :no-members:

.. autoclass:: tatc.tradespace.DeploymentSchedule
   :members: plan, launches_of, mass_capacity_of

.. autoclass:: tatc.tradespace.Deployment
   :members: simulate, simulate_many, expected_metrics, planned, phases, begin, end, deployed, units, stages, capability

.. autoclass:: tatc.tradespace.Launch

.. autoclass:: tatc.tradespace.Production

.. autoclass:: tatc.tradespace.Stage

.. autoclass:: tatc.tradespace.Decision
   :members: time, stage, operating, history, deployed, capability

Failures and Missions
---------------------

.. automodule:: tatc.tradespace.failures
   :no-members:

.. autoclass:: tatc.tradespace.FailureModel
   :members: life_of

.. autoclass:: tatc.tradespace.Lifetime
   :members: survival, hazard, cumulative_hazard, sample, from_uniform, dimensions, stressed

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

.. autoclass:: tatc.tradespace.Capability
   :members: at
   :special-members: __call__

.. autofunction:: tatc.tradespace.evaluate_missions

.. autofunction:: tatc.tradespace.summarize_missions

.. autofunction:: tatc.tradespace.allocate_histories

.. autofunction:: tatc.tradespace.expected_missions

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

Trade Studies
-------------

.. autoclass:: tatc.tradespace.Study
   :members: run, workers, definition, fingerprint

Observation Studies
-------------------

.. automodule:: tatc.tradespace.observations
   :no-members:

.. autoclass:: tatc.tradespace.ObservationStudy
   :members: satellites_of, targets_of, observations_of

.. autoclass:: tatc.tradespace.Targets
   :members: of, seed_of

.. autoclass:: tatc.tradespace.Observations
   :members: end, start_of, end_of, cache_of, context_of, __call__

.. autoclass:: tatc.tradespace.MissionPlan
   :members: schedule_of, cost_model_of, masses_of, growth_of, stages_of

.. autoclass:: tatc.tradespace.Failures

.. autofunction:: tatc.tradespace.value_of

Criteria, Statistics, and Requirements
--------------------------------------

.. automodule:: tatc.tradespace.criteria
   :no-members:

.. autoclass:: tatc.tradespace.Criteria
   :members: names, metrics, statistic

.. autofunction:: tatc.tradespace.aggregate

.. autofunction:: tatc.tradespace.summarize

.. autoclass:: tatc.tradespace.Statistic
   :members: reduce, columns, values_of

.. autoclass:: tatc.tradespace.Mean

.. autoclass:: tatc.tradespace.StandardDeviation

.. autoclass:: tatc.tradespace.StandardError

.. autoclass:: tatc.tradespace.Minimum

.. autoclass:: tatc.tradespace.Maximum

.. autoclass:: tatc.tradespace.Quantile

.. autoclass:: tatc.tradespace.Probability

.. autoclass:: tatc.tradespace.Requirement
   :members: violation

.. autoclass:: tatc.tradespace.Across
   :members: weights_of

Utility
-------

.. automodule:: tatc.tradespace.utility
   :no-members:

.. autoclass:: tatc.tradespace.Utility
   :members: columns, scaling, evaluate, __call__

.. autoclass:: tatc.tradespace.LinearUtility

.. autoclass:: tatc.tradespace.ExponentialUtility

.. autoclass:: tatc.tradespace.PiecewiseUtility

Stores of Results
-----------------

.. autoclass:: tatc.tradespace.ResultStore
   :members: get, put, define, definition, contexts, prune, clear

Fingerprints
------------

.. automodule:: tatc.tradespace.fingerprint
   :no-members:

.. autofunction:: tatc.tradespace.get_fingerprint

.. autofunction:: tatc.tradespace.get_definition

.. autofunction:: tatc.tradespace.compare_definitions

Caches
------

.. automodule:: tatc.tradespace.cache
   :no-members:

.. autoclass:: tatc.tradespace.SatelliteCache
   :members: get, prune, clear, __call__

.. autofunction:: tatc.tradespace.get_satellite_key

Pareto Efficiency
-----------------

.. autofunction:: tatc.tradespace.pareto_front

.. autofunction:: tatc.tradespace.pareto_ranks

.. autofunction:: tatc.tradespace.hypervolume

.. autofunction:: tatc.tradespace.hypervolume_progress

.. autofunction:: tatc.tradespace.pareto_probability

Adaptive Replication
--------------------

.. autofunction:: tatc.tradespace.refine

.. autoclass:: tatc.tradespace.Refinement
   :members: uncertain, indifferent

Confirmation
------------

.. autofunction:: tatc.tradespace.confirm

.. autoclass:: tatc.tradespace.Confirmation
   :members: optimism

Sensitivity
-----------

.. autofunction:: tatc.tradespace.sensitivity

.. autofunction:: tatc.tradespace.main_effects

Scenarios
---------

.. automodule:: tatc.tradespace.scenarios
   :no-members:

.. autofunction:: tatc.tradespace.robustness

Attribution of Variance
-----------------------

.. automodule:: tatc.tradespace.attribution
   :no-members:

.. autofunction:: tatc.tradespace.attribute_variance

.. autoclass:: tatc.tradespace.Attribution
   :members: summary, differences

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

.. autoclass:: tatc.tradespace.Search

.. autoclass:: tatc.tradespace.Evaluations
   :members: evaluate, results, worst, designs

Plots
-----

.. automodule:: tatc.tradespace.plotting
   :no-members:

.. autofunction:: tatc.tradespace.plot_designs

.. autofunction:: tatc.tradespace.plot_progress

.. autofunction:: tatc.tradespace.plot_means

.. autofunction:: tatc.tradespace.style_axes

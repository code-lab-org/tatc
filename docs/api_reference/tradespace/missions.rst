===========================
Constellations and Missions
===========================

The parts of observation studies (see `ObservationStudy`): the layouts of
constellations, their deployment, the failures of their satellites, the
phases and histories of their missions, and their costs. Studies use them
through a `MissionPlan` and `Failures`; the rest are for analyses of their
own.

Constellations
==============

.. automodule:: tatc.tradespace.constellations
   :no-members:

.. autoclass:: tatc.tradespace.constellations.Layout
   :members: planes, orbits, names, describe, satellites, build

.. autoclass:: tatc.tradespace.constellations.Plane

.. autofunction:: tatc.tradespace.layout_of

.. autoclass:: tatc.tradespace.constellations.WalkerLayout

.. autoclass:: tatc.tradespace.constellations.SSOPlanesLayout

.. autoclass:: tatc.tradespace.constellations.PlanesLayout

.. autoclass:: tatc.tradespace.constellations.SOCLayout

.. autoclass:: tatc.tradespace.constellations.MOGLayout

.. autoclass:: tatc.tradespace.constellations.Layouts

.. autofunction:: tatc.tradespace.constellations.plane_sizes

.. autofunction:: tatc.tradespace.constellations.hours_to_time

.. autofunction:: tatc.tradespace.constellations.remove_failures

Deployment
==========

.. automodule:: tatc.tradespace.deployment
   :no-members:

.. autoclass:: tatc.tradespace.DeploymentSchedule
   :members: plan, launches_of

.. autoclass:: tatc.tradespace.Stage

.. autoclass:: tatc.tradespace.deployment.Deployment
   :members: simulate, simulate_many, planned, phases, begin, end, deployed, units, stages, policy, capability

.. automodule:: tatc.tradespace.events
   :no-members:

.. autoclass:: tatc.tradespace.events.PlannedPlane

.. autoclass:: tatc.tradespace.events.Launch

.. autoclass:: tatc.tradespace.events.Production

.. autoclass:: tatc.tradespace.mission.Decision
   :members: time, since_launch, stage, satellites_operating, history, stages_deployed, capability

Failures
========

.. automodule:: tatc.tradespace.failures
   :no-members:

.. autoclass:: tatc.tradespace.FailureModel
   :members: life_of

.. autoclass:: tatc.tradespace.failures.Lifetime
   :members: survival, hazard, cumulative_hazard, sample, from_uniform, dimensions, stressed

.. autoclass:: tatc.tradespace.Exponential

.. autoclass:: tatc.tradespace.Weibull

.. autoclass:: tatc.tradespace.Bathtub

.. autoclass:: tatc.tradespace.Piecewise

.. autoclass:: tatc.tradespace.RestorationPolicy
   :members: orbit_spares_of, ground_spares_of

Missions
========

.. automodule:: tatc.tradespace.mission
   :no-members:

.. autoclass:: tatc.tradespace.mission.MissionHistory
   :members: phases, metrics, in_orbit, deployed, decided

.. autoclass:: tatc.tradespace.mission.Phase
   :members: duration

.. autoclass:: tatc.tradespace.mission.Failure

.. autoclass:: tatc.tradespace.mission.Restoration

.. autofunction:: tatc.tradespace.mission.evaluate_phases

.. autofunction:: tatc.tradespace.mission.summarize_phases

.. autoclass:: tatc.tradespace.mission.PhaseCache
   :members: clear

.. autoclass:: tatc.tradespace.mission.Capability
   :members: evaluate, names
   :special-members: __call__

.. automodule:: tatc.tradespace.histories
   :no-members:

.. autofunction:: tatc.tradespace.histories.evaluate_histories

.. autofunction:: tatc.tradespace.histories.summarize_histories

.. autofunction:: tatc.tradespace.histories.allocate_histories

.. automodule:: tatc.tradespace.expectations
   :no-members:

.. autofunction:: tatc.tradespace.expectations.exact_expectations

.. autoexception:: tatc.tradespace.expectations.ExactExpectationsError

.. autodata:: tatc.tradespace.expectations.MAX_SETS
   :annotation:

Costs
=====

.. automodule:: tatc.tradespace.cost
   :no-members:

.. autoclass:: tatc.tradespace.CostModel
   :members: estimate, satellite_cost_of, draw, score_names, discount

.. autoclass:: tatc.tradespace.SatelliteCost
   :members: production, blocks, total_development, total_first_unit

.. autoclass:: tatc.tradespace.Component
   :members: uncertain, factor

.. autoclass:: tatc.tradespace.Triangular
   :members: mean, std, ppf

.. autofunction:: tatc.tradespace.risk_space

.. autoclass:: tatc.tradespace.LaunchCost

.. autoclass:: tatc.tradespace.OperationsCost

.. autoclass:: tatc.tradespace.cost.CostEstimate
   :members: total, present_value, annual, metrics

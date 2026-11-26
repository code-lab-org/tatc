Simulated Observations
======================

.. automodule:: tatc.analysis.field_sampling
   :no-members:

Functions that simulate the observations of an observing system simulation
experiment (OSSE) at the level of sampling: the quality of observations, a
field of a nature run sampled at each observation (at its point, over its
footprint, inside its region, at a level, or by weights of levels), profiles
at heights, and observation errors of a size that the user gives. Forward
models (such as radiative transfer) and retrievals are outside their scope.
See :func:`tatc.preprocess.find_features` for the features of a nature run,
and :mod:`tatc.export` for files of data assimilation systems.

Observation Quality
-------------------
.. autofunction:: tatc.analysis.compute_observation_quality

Sample Fields
-------------
.. autofunction:: tatc.analysis.sample_field

.. autofunction:: tatc.analysis.sample_profiles

Observation Errors
------------------
.. autofunction:: tatc.analysis.simulate_observation_errors

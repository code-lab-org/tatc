=======================
Criteria and Summaries
=======================

Criteria, Statistics, and Requirements
======================================

.. automodule:: tatc.tradespace.criteria
   :no-members:

.. autoclass:: tatc.tradespace.Criteria
   :members: names, metrics, statistic

.. autofunction:: tatc.tradespace.aggregate

.. autofunction:: tatc.tradespace.summarize

.. autoclass:: tatc.tradespace.Mean

.. autoclass:: tatc.tradespace.StandardDeviation

.. autoclass:: tatc.tradespace.StandardError

.. autoclass:: tatc.tradespace.Minimum

.. autoclass:: tatc.tradespace.Maximum

.. autoclass:: tatc.tradespace.Quantile

.. autoclass:: tatc.tradespace.Probability

.. autoclass:: tatc.tradespace.criteria.Statistic
   :members: reduce, columns, values_of

.. autoclass:: tatc.tradespace.Requirement
   :members: violation

.. autoclass:: tatc.tradespace.Across
   :members: weights_of

Pareto Efficiency
=================

.. automodule:: tatc.tradespace.pareto
   :no-members:

.. autofunction:: tatc.tradespace.pareto_front

.. autofunction:: tatc.tradespace.pareto_ranks

.. autofunction:: tatc.tradespace.pareto_probability

.. autofunction:: tatc.tradespace.hypervolume

.. autofunction:: tatc.tradespace.hypervolume_progress

Comparisons
===========

.. automodule:: tatc.tradespace.comparison
   :no-members:

.. autofunction:: tatc.tradespace.exceedance

Scenarios
=========

.. automodule:: tatc.tradespace.scenarios
   :no-members:

.. autofunction:: tatc.tradespace.robustness

Sensitivity
===========

.. autofunction:: tatc.tradespace.sensitivity

.. autofunction:: tatc.tradespace.main_effects

Utility
=======

.. automodule:: tatc.tradespace.utility
   :no-members:

.. autoclass:: tatc.tradespace.Utility
   :members: columns, scaling, evaluate, __call__

.. autoclass:: tatc.tradespace.LinearUtility

.. autoclass:: tatc.tradespace.ExponentialUtility

.. autoclass:: tatc.tradespace.PiecewiseUtility

Plots
=====

.. automodule:: tatc.tradespace.plotting
   :no-members:

.. autofunction:: tatc.tradespace.plot_designs

.. autofunction:: tatc.tradespace.plot_means

.. autofunction:: tatc.tradespace.plot_exceedance

.. autofunction:: tatc.tradespace.plot_progress

.. autofunction:: tatc.tradespace.style_axes

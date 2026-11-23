=======================
Criteria and Summaries
=======================

Criteria, Statistics, and Requirements
======================================

.. automodule:: tatc.tradespace.criteria
   :no-members:

.. autoclass:: tatc.tradespace.Criteria
   :members: names, metrics

.. autoclass:: tatc.tradespace.Requirement
   :members: met, violation

.. autoclass:: tatc.tradespace.Across

.. autodata:: tatc.tradespace.criteria.Objectives

.. autodata:: tatc.tradespace.criteria.Direction

.. autodata:: tatc.tradespace.criteria.Rule

.. automodule:: tatc.tradespace.statistics
   :no-members:

.. autodata:: tatc.tradespace.statistics.Operator

.. autoclass:: tatc.tradespace.Mean

.. autoclass:: tatc.tradespace.StandardDeviation

.. autoclass:: tatc.tradespace.Minimum

.. autoclass:: tatc.tradespace.Maximum

.. autoclass:: tatc.tradespace.Quantile

.. autoclass:: tatc.tradespace.Probability

.. autoclass:: tatc.tradespace.Statistic
   :members: name, metric, reduce

Summaries
=========

.. automodule:: tatc.tradespace.summary
   :no-members:

.. autofunction:: tatc.tradespace.summarize

.. autofunction:: tatc.tradespace.pareto_probability

Pareto Efficiency
=================

.. automodule:: tatc.tradespace.pareto
   :no-members:

.. autofunction:: tatc.tradespace.pareto_front

.. autofunction:: tatc.tradespace.pareto_ranks

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

.. autofunction:: tatc.tradespace.plot_shares

.. autofunction:: tatc.tradespace.style_axes

.. autofunction:: tatc.tradespace.use_style

.. autodata:: tatc.tradespace.plotting.STYLE
   :annotation:

.. autodata:: tatc.tradespace.plotting.COLORS
   :annotation:

.. autodata:: tatc.tradespace.plotting.INK
   :annotation:

.. autodata:: tatc.tradespace.plotting.GRAY
   :annotation:

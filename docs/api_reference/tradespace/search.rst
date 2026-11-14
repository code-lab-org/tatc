========================
Refinement and Searches
========================

Each of these steps takes a study and criteria (or metrics), and returns its
`results`, the `summary` of its designs, and the `designs` themselves, which
the next step takes.

Adaptive Replication
====================

.. autofunction:: tatc.tradespace.refine

.. autoclass:: tatc.tradespace.refinement.Refinement
   :members: results, summary, history, designs, uncertain, indifferent

Confirmation
============

.. autofunction:: tatc.tradespace.confirm

.. autoclass:: tatc.tradespace.confirmation.Confirmation
   :members: results, summary, criteria, designs, optimism

Screening
=========

.. autofunction:: tatc.tradespace.screen

.. autoclass:: tatc.tradespace.screening.Screening
   :members: results, summary, summaries, designs

Attribution of Variance
=======================

.. automodule:: tatc.tradespace.attribution
   :no-members:

.. autofunction:: tatc.tradespace.attribute_variance

.. autoclass:: tatc.tradespace.attribution.Attribution
   :members: summary, differences

Optimizers
==========

.. automodule:: tatc.tradespace.optimization
   :no-members:

.. autofunction:: tatc.tradespace.optimize_optuna

.. autofunction:: tatc.tradespace.optimize_pymoo

.. autoclass:: tatc.tradespace.optimization.Search
   :members: results, summary, designs, state

.. autofunction:: tatc.tradespace.optimization.make_pymoo_problem

.. autoclass:: tatc.tradespace.optimization.Evaluations
   :members: evaluate, results, worst, designs, search

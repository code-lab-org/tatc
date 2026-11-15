========================
Refinement and Searches
========================

Each of these steps takes a study and criteria (or metrics), and returns an
outcome (exported for type hints, such as `Refinement`) of its `results`, the
`summary` of its designs, the `designs` themselves, which the next step
takes, and its `criteria`, in that order (an `Attribution` gives the
`shares` of the variance of each design, and their `summary`).

Adaptive Replication
====================

.. autofunction:: tatc.tradespace.refine

.. autoclass:: tatc.tradespace.Refinement
   :members: results, summary, designs, criteria, rounds, candidates, evaluations

Confirmation
============

.. autofunction:: tatc.tradespace.confirm

.. autoclass:: tatc.tradespace.Confirmation
   :members: results, summary, designs, criteria, optimism

Screening
=========

.. autofunction:: tatc.tradespace.screen

.. autoclass:: tatc.tradespace.Screening
   :members: results, summary, designs, criteria, history, summaries

Attribution of Variance
=======================

.. automodule:: tatc.tradespace.attribution
   :no-members:

.. autofunction:: tatc.tradespace.attribute_variance

.. autoclass:: tatc.tradespace.Attribution
   :members: results, shares, designs, by, summary, summary_of, differences

Optimizers
==========

.. automodule:: tatc.tradespace.optimization
   :no-members:

.. autofunction:: tatc.tradespace.optimize_optuna

.. autofunction:: tatc.tradespace.optimize_pymoo

.. autoclass:: tatc.tradespace.Search
   :members: results, summary, designs, criteria, evaluated, state

.. autofunction:: tatc.tradespace.optimization.make_pymoo_problem

.. autoclass:: tatc.tradespace.optimization.Evaluations
   :members: evaluate, results, worst, designs, search

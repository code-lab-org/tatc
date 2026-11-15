========================
Refinement and Searches
========================

Each of these steps takes a study and criteria (or metrics), and returns an
outcome (exported for type hints, such as `Refinement`) of its `results`, the
`summary` of its designs, the `designs` themselves, which the next step
takes, and its `criteria`.

Adaptive Replication
====================

.. autofunction:: tatc.tradespace.refine

.. autoclass:: tatc.tradespace.Refinement
   :members: results, summary, rounds, designs, candidates, criteria, evaluations, uncertain, indifferent

Confirmation
============

.. autofunction:: tatc.tradespace.confirm

.. autoclass:: tatc.tradespace.Confirmation
   :members: results, summary, criteria, designs, optimism

Screening
=========

.. autofunction:: tatc.tradespace.screen

.. autoclass:: tatc.tradespace.Screening
   :members: results, summary, summaries, history, designs, criteria

Attribution of Variance
=======================

.. automodule:: tatc.tradespace.attribution
   :no-members:

.. autofunction:: tatc.tradespace.attribute_variance

.. autoclass:: tatc.tradespace.Attribution
   :members: results, shares, summary, summary_of, differences, designs

Optimizers
==========

.. automodule:: tatc.tradespace.optimization
   :no-members:

.. autofunction:: tatc.tradespace.optimize_optuna

.. autofunction:: tatc.tradespace.optimize_pymoo

.. autoclass:: tatc.tradespace.Search
   :members: results, summary, designs, evaluated, criteria, state

.. autofunction:: tatc.tradespace.optimization.make_pymoo_problem

.. autoclass:: tatc.tradespace.optimization.Evaluations
   :members: evaluate, results, worst, designs, search

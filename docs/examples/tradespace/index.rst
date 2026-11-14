====================
Trade Space Tutorial
====================

These notebooks follow one trade study with the ``tatc.tradespace`` module
from start to finish: how many satellites, in how many orbit planes, an
imager needs to observe the Earth often enough. Each notebook adds one idea
to those before it, so read them in order:

#. :doc:`EvaluateTradeSpace`: declare a design space and a study, evaluate
   every design, and find the Pareto-efficient designs.
#. :doc:`ExplainTradeSpace`: measure sampling by its error, explain why
   designs differ, and find the decisions that drive the objectives.
#. :doc:`SearchTradeSpace`: find the efficient designs of larger spaces by
   screening with a cheaper study and by searching with optimizers.
#. :doc:`RankUncertainDesigns`: rank designs with the uncertainty of their
   replicates, refine and confirm the front, and require reliable
   performance.
#. :doc:`CompareScenarios`: compare designs under conditions that are not
   known in advance, and find designs robust to them.
#. :doc:`TradeCostAndUtility`: estimate life-cycle costs with uncertain
   growth, trade them within a budget, and weigh objectives by utility.
#. :doc:`PlanDeployment`: account for the phases of each design's
   deployment, and stage it.
#. :doc:`PlanForFailures`: simulate failures, choose restoration policies
   with the designs, and find what spreads the results.

The notebooks import the design space and the study of the first notebook,
and the parts of later notebooks that the notebooks after them use, from
``trade_space_example.py`` beside them. Each study keeps its results in a
store, from which a notebook takes the evaluations of the notebooks before it
where their studies are the same. The notebooks require::

  pip install tatc[tradespace]

.. toctree::
  :hidden:

  EvaluateTradeSpace.ipynb
  ExplainTradeSpace.ipynb
  SearchTradeSpace.ipynb
  RankUncertainDesigns.ipynb
  CompareScenarios.ipynb
  TradeCostAndUtility.ipynb
  PlanDeployment.ipynb
  PlanForFailures.ipynb

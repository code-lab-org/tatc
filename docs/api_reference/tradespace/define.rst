==============
Design Spaces
==============

.. automodule:: tatc.tradespace.space
   :no-members:

Design Spaces and Designs
=========================

.. autoclass:: tatc.tradespace.DesignSpace
   :members: enumerate, size, sample, design, designs_of, extend, suggest, unit_paths, uncertain_paths, encode, decode

.. autoclass:: tatc.tradespace.Design
   :members: flat, realize, values

Variables and Values
====================

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

.. autoclass:: tatc.tradespace.nodes.Variable
   :members: domain

Building Blocks of Constellations
=================================

Fragments of design spaces of common layouts of constellations, and of the
policies that restore their failures, whose values a study builds into
satellites (see :doc:`missions` for the layouts themselves).

.. autofunction:: tatc.tradespace.walker_space

.. autofunction:: tatc.tradespace.sso_planes_space

.. autofunction:: tatc.tradespace.planes_space

.. autofunction:: tatc.tradespace.policy_space

.. autofunction:: tatc.tradespace.build_constellation

.. autofunction:: tatc.tradespace.count_satellites

.. autofunction:: tatc.tradespace.count_planes

Views of Values
===============

The callables of a design space (bounds, derived values, and constraints)
and the values of designs are views of the values decided so far.

.. autoclass:: tatc.tradespace.values.Values

.. autoclass:: tatc.tradespace.values.ChoiceValues
   :members: choice, value

.. autoclass:: tatc.tradespace.values.RepeatValues
   :members: count

.. autoclass:: tatc.tradespace.values.Scope

.. autofunction:: tatc.tradespace.values.configure

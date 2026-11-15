==============
Design Spaces
==============

.. automodule:: tatc.tradespace.space
   :no-members:

Design Spaces and Designs
=========================

.. autoclass:: tatc.tradespace.DesignSpace
   :members: enumerate, size, sample, design, designs_of, extend, uncertain_paths

.. autoclass:: tatc.tradespace.Design
   :members: flat, uncertain_values, realize, values, decisions, scenario

.. autodata:: tatc.tradespace.space.RESERVED_NAMES
   :annotation:

.. autodata:: tatc.tradespace.space.STATUS_COLUMNS
   :annotation:

Variables and Values
====================

.. autoclass:: tatc.tradespace.Real

.. autoclass:: tatc.tradespace.Integer

.. autoclass:: tatc.tradespace.Categorical

.. autoclass:: tatc.tradespace.Fixed

.. autoclass:: tatc.tradespace.Derived

.. autoclass:: tatc.tradespace.Uncertain
   :members: uniform, normal, triangular, choice, seed

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

.. autofunction:: tatc.tradespace.soc_space

.. autofunction:: tatc.tradespace.mog_space

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

.. autofunction:: tatc.tradespace.values.generic_path

Reading Designs
===============

.. automodule:: tatc.tradespace.references
   :no-members:

.. autodata:: tatc.tradespace.references.Reference
   :annotation:

.. autofunction:: tatc.tradespace.value_of

.. autofunction:: tatc.tradespace.configure

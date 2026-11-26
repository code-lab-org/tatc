================
Sizing Functions
================

Sizing functions estimate what a spacecraft needs to keep its orbit for its
life: the velocity to raise its orbit from a launch's injection, to hold it
against drag, and to dispose of it, and the propellant these take. The drag
of the thermosphere and the decay of orbits are utility functions (see
:func:`tatc.utils.compute_atmospheric_density` and
:func:`tatc.utils.compute_orbit_lifetime`).

Velocity Budgets
----------------

.. autofunction:: tatc.sizing.compute_delta_v_budget

.. autopydantic_model:: tatc.sizing.DeltaVBudget

.. autofunction:: tatc.sizing.compute_disposal_delta_v

.. autopydantic_model:: tatc.sizing.DisposalRule

.. autodata:: tatc.sizing.budget.MAX_ECCENTRICITY

Propellant
----------

.. autofunction:: tatc.sizing.compute_propellant_mass

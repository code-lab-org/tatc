"""
Sizing of spacecraft for pre-Phase A studies: the velocity budgets of their
orbits (raising, holding against drag, and disposal), the propellant of those
budgets, and the mass and power of spacecraft sized from their payloads.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from .budget import (
    DeltaVBudget,
    DisposalRule,
    compute_delta_v_budget,
    compute_disposal_delta_v,
    compute_propellant_mass,
)
from .spacecraft import (
    LIFETIME_HORIZON,
    SizedSpacecraft,
    SpacecraftSizing,
    size_spacecraft,
)

__all__ = [
    "DeltaVBudget",
    "DisposalRule",
    "compute_delta_v_budget",
    "compute_disposal_delta_v",
    "compute_propellant_mass",
    "LIFETIME_HORIZON",
    "SizedSpacecraft",
    "SpacecraftSizing",
    "size_spacecraft",
]

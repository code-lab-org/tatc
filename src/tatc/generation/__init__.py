"""
Defines generation functions.
"""

from .cells import generate_equally_spaced_cells
from .points import (
    generate_equally_spaced_points,
    generate_fibonacci_lattice_points,
)

__all__ = [
    "generate_equally_spaced_cells",
    "generate_equally_spaced_points",
    "generate_fibonacci_lattice_points",
]
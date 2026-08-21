"""
Object schema for orbits defined by two line elements (TLE).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from typing import Literal

from pydantic import Field, model_validator

from .base import OrbitBase
from .gp import GeneralPerturbationsOrbit


class TwoLineElements(OrbitBase):
    """
    Orbit defined by two line elements (TLE).
    """

    type: Literal["tle"] = Field(default="tle", description="Orbit type discriminator.")
    tle: tuple[str, str] = Field(
        ...,
        description="The two TLE lines.",
        examples=[
            (
                "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
                "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
            )
        ],
    )

    def _compute_gp_orbit(self) -> GeneralPerturbationsOrbit:
        """
        Computes a general perturbations orbit representation of this
        orbit by parsing its TLE lines via SGP4.

        Returns:
            GeneralPerturbationsOrbit: the general perturbations orbit
        """
        return GeneralPerturbationsOrbit.from_tle(list(self.tle))

    @model_validator(mode="after")
    def _validate_tle(self) -> TwoLineElements:
        """
        Validates the TLE lines by parsing them, raising if they cannot
        be converted to a general perturbations orbit. Reuses
        `to_gp_orbit()` (rather than calling `_compute_gp_orbit()`
        directly) so the parsed result is cached, avoiding a redundant
        SGP4 fit the first time this orbit's gp orbit is requested.
        """
        self.to_gp_orbit()
        return self

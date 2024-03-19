# -*- coding: utf-8 -*-
"""
Object schemas for architectures.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from typing import List, Union

from pydantic import BaseModel, Field

from .satellite import Satellite, TrainConstellation, WalkerConstellation
from .point import GroundStation


class Architecture(BaseModel):
    """
    Mission architecture.
    """

    name: str = Field(..., description="Name of this mission.")
    satellites: List[Union[Satellite, TrainConstellation, WalkerConstellation]] = Field(
        [], description="List of member space systems."
    )
    stations: List[GroundStation] = Field(
        [], description="List of member ground stations."
    )

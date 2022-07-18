# -*- coding: utf-8 -*-
"""
Object schemas for architectures.

@author: Paul T. Grogan <pgrogan@stevens.edu>
"""

from typing import Optional, List, Union
from pydantic import BaseModel, Field
from datetime import timedelta

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

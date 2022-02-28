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
    Representation of a mission architecture.

    :param name: Name of this mission
    :type name: str
    :param satellites: List of member satellites
    :type satellites: list[:class:`tatc.schemas.satellite.Satellite`],
        list[:class:`tatc.schemas.satellite.TrainConstellation`], or
        list[:class:`tatc.schemas.satellite.WalkerConstellation`]
    :param stations: List of member ground stations
    :type stations: list[:class:`tatc.schemas.point.GroundStation`]
    """

    name: str = Field(..., description="Name of this mission.")
    satellites: List[Union[Satellite, TrainConstellation, WalkerConstellation]] = Field(
        [], description="List of member satellites."
    )
    stations: List[GroundStation] = Field(
        [], description="List of member ground stations."
    )

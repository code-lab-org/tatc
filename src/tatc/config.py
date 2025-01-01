# -*- coding: utf-8 -*-
"""
Configuration Settings.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import os
import pathlib

import yaml
from yaml.parser import ParserError
from pydantic import BaseModel, Field, ValidationError


class ConfigError(Exception):
    """Configuration error"""


class RuntimeConfiguration(BaseModel):
    """
    Runtime configuration settings.
    """

    footprint_points_elliptical: int = Field(
        32, description="Number of points for a SPICE elliptical footprint.", ge=4
    )
    footprint_points_rectangular_side: int = Field(
        8, description="Number of points for a SPICE rectangular footprint side.", ge=1
    )
    repeat_cycle_delta_position_m: float = Field(
        10000,
        description="Maximum difference in position (meters) for a valid repeat.",
        gt=0,
    )
    repeat_cycle_delta_velocity_m_per_s: float = Field(
        10,
        description="Maximum difference in velocity (meters/second) for a valid repeat.",
        gt=0,
    )
    repeat_cycle_search_elevation_deg: float = Field(
        88, description="Minimum elevation angle (degrees) for .", gt=0
    )
    repeat_cycle_search_duration_days: float = Field(
        30, description="Maximum duration for which to search for repeat cycles."
    )
    repeat_cycle_lazy_load: bool = Field(
        True, description="True, if a previously-computed repeat cycle should be used."
    )
    repeat_cycle_for_orbit_track: bool = Field(
        True,
        description="True, if a repeat cycle should be used to generate orbit tracks.",
    )
    repeat_cycle_for_observation_events: bool = Field(
        True,
        description="True, if a repeat cycle should be used to generate observation events.",
    )
    orbit_tle_lazy_load: bool = Field(
        True, description="True, if a previously-computed tle should be used."
    )


def load_yaml_config(path: pathlib.Path):
    """
    Load configuration settings from a YAML file.

    Args:
      path (Path): The file path to load.

    Returns:
      RuntimeConfiguration: The loaded configuration settings.
    """
    if not os.path.exists(path):
        raise ConfigError("Couldn't load config file (not found)")
    with open(path, "r", encoding="utf-8") as f:
        try:
            config = yaml.safe_load(f)
        except ParserError as err:
            raise ConfigError(f"Couldn't parse config file - {err}") from err
    try:
        return RuntimeConfiguration(**config)
    except ValidationError as err:
        raise ConfigError(f"Couldn't validate config file - {err}") from err


rc = load_yaml_config(
    os.path.join(os.path.dirname(__file__), "resources", "defaults.yml")
)

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

    footprint_points: int = Field(
        33, description="Number of points for a SPICE footprint.", gt=4
    )
    footprint_points_rect_side: int = Field(
        4, description="Number of points for a SPICE rectangular footprint side."
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

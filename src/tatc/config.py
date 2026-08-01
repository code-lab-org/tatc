"""
Configuration Settings.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import functools
import logging
import os
import pathlib

import yaml
from pydantic import BaseModel, Field, ValidationError
from yaml.parser import ParserError

logger = logging.getLogger(__name__)


class ConfigError(Exception):
    """Configuration error"""


class RuntimeConfiguration(BaseModel):
    """
    Runtime configuration settings.
    """

    footprint_points_elliptical: int = Field(
        default=32,
        description="Number of points for a SPICE elliptical footprint.",
        ge=4,
    )
    footprint_points_rectangular_side: int = Field(
        default=8,
        description="Number of points for a SPICE rectangular footprint side.",
        ge=1,
    )
    repeat_cycle_delta_position_m: float = Field(
        default=10000,
        description="Maximum difference in position (meters) for a valid repeat.",
        gt=0,
    )
    repeat_cycle_delta_velocity_m_per_s: float = Field(
        default=10,
        description="Maximum difference in velocity (meters/second) for a valid repeat.",
        gt=0,
    )
    repeat_cycle_search_elevation_deg: float = Field(
        default=88,
        description="Minimum elevation angle (degrees) for screening repeats.",
        gt=0,
    )
    repeat_cycle_search_duration_days: float = Field(
        default=30,
        description="Maximum duration for which to search for repeat cycles.",
    )
    repeat_cycle_lazy_load: bool = Field(
        default=True,
        description="True, if a previously-computed repeat cycle should be used.",
    )
    repeat_cycle_for_orbit_track: bool = Field(
        default=True,
        description="True, if a repeat cycle should be used to generate orbit tracks.",
    )
    repeat_cycle_for_observation_events: bool = Field(
        default=True,
        description="True, if a repeat cycle should be used to generate observation events.",
    )
    gp_orbit_lazy_load: bool = Field(
        default=True,
        description="True, if a previously-computed general perturbations orbit should be used.",
    )


def load_yaml_config(path: pathlib.Path) -> RuntimeConfiguration:
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


@functools.lru_cache(maxsize=1)
def get_rc() -> RuntimeConfiguration:
    """
    Return the process-wide runtime configuration, loading it from the
    packaged defaults file on first access and caching it thereafter.

    Returns:
      RuntimeConfiguration: The runtime configuration settings.
    """
    try:
        return load_yaml_config(pathlib.Path(__file__).parent / "resources" / "defaults.yml")
    except ConfigError as err:
        # fall back to default constructor, but warn since this masks a broken install
        logger.warning("Falling back to hard-coded runtime configuration defaults: %s", err)
        return RuntimeConfiguration()


def reset_rc() -> None:
    """
    Clear the cached runtime configuration so the next call to get_rc()
    reloads it from disk. Intended for tests that need to exercise
    load-time behavior (e.g. a missing packaged defaults file); not
    needed in normal use.
    """
    get_rc.cache_clear()


def __getattr__(name: str):
    # Preserve `tatc.config.rc` as a read path onto the lazy singleton, for
    # backward compatibility with code written against the old eager global.
    if name == "rc":
        return get_rc()
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

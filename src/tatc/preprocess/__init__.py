"""
Preprocessing utilities that derive TAT-C model inputs from external
datasets and tools: raster weights and terrain from geospatial datasets (e.g.
a digital elevation model), which require `pip install tatc[preprocess]`, and
the trajectories and instruments of aircraft from flight plans, recordings,
and sensor models (e.g. those of HyPlan, installed with
`pip install tatc[airborne]`). Each module is imported when first used, so
that it requires only its own dependencies.
"""

from typing import Any

_MODULES = {
    "aircraft_from_hyplan": "airborne",
    "estimate_attitude": "airborne",
    "instrument_from_sensor": "airborne",
    "sample_flight_plan": "airborne",
    "states_from_iwg1": "airborne",
    "trajectory_from_waypoints": "airborne",
    "read_raster_weights": "raster",
    "add_local_incidence_angles": "terrain",
    "compute_terrain_mask": "terrain",
    "compute_terrain_mask_for_platform": "terrain",
    "compute_terrain_mask_for_station": "terrain",
    "compute_terrain_slope_aspect": "terrain",
    "get_copernicus_dem_tile_urls": "terrain",
    "sample_dem_elevation": "terrain",
}


def __getattr__(name: str) -> Any:
    """
    Gets a utility from its module, imported when first used.
    """
    # pylint: disable=import-outside-toplevel
    module = _MODULES.get(name)
    if module == "airborne":
        from . import airborne as source
    elif module == "raster":
        from . import raster as source
    elif module == "terrain":
        from . import terrain as source
    else:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    return getattr(source, name)


def __dir__() -> list[str]:
    return sorted(list(globals()) + list(_MODULES))


__all__ = sorted(_MODULES)

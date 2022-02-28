import json
import geopandas as gpd
from pkg_resources import resource_stream
from shapely.geometry import shape
from tatc.generation.points import (
    generate_cubed_sphere_points,
    generate_fibonacci_lattice_points,
)
from tatc.generation.cells import generate_cubed_sphere_cells

from .schemas import KnownShape
from ..worker import app


def load_country_mask(iso_3166_1_alpha3: str):
    world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
    return world.query(f'iso_a3 == "{iso_3166_1_alpha3}"').to_crs("EPSG:4326").geometry


def load_known_shape(shape: KnownShape):
    if shape == KnownShape.conus:
        usa = gpd.read_file(
            resource_stream(__name__, "../../resources/cb_2020_us_state_20m.zip").name
        )
        return (
            usa[(usa.STUSPS != "AK") & (usa.STUSPS != "HI") & (usa.STUSPS != "PR")]
            .dissolve()
            .to_crs("EPSG:4326")
            .geometry
        )


@app.task
def generate_cubed_sphere_points_task(distance, mask):
    return generate_cubed_sphere_points(
        distance,
        load_country_mask(mask)
        if isinstance(mask, str) and len(mask) == 3
        else load_known_shape(mask)
        if isinstance(mask, str) and mask in KnownShape.__members__
        else shape(mask)
        if mask is not None
        else None,
    ).to_json(show_bbox=False, drop_id=True)


@app.task
def generate_fibonacci_lattice_points_task(distance, mask):
    return generate_fibonacci_lattice_points(
        distance,
        load_country_mask(mask)
        if isinstance(mask, str) and len(mask) == 3
        else load_known_shape(mask)
        if isinstance(mask, str) and mask in KnownShape.__members__
        else shape(mask)
        if mask is not None
        else None,
    ).to_json(show_bbox=False, drop_id=True)


@app.task
def generate_cubed_sphere_cells_task(distance, mask, strips):
    return generate_cubed_sphere_cells(
        distance,
        load_country_mask(mask)
        if isinstance(mask, str) and len(mask) == 3
        else load_known_shape(mask)
        if isinstance(mask, str) and mask in KnownShape.__members__
        else shape(mask)
        if mask is not None
        else None,
        strips,
    ).to_json(show_bbox=False, drop_id=True)

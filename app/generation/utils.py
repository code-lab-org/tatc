import json
import geopandas as gpd

from .schemas import (
    Point,
    PointGenerator,
    PointGeneratorMethod,
    CellGenerator,
    CellGeneratorMethod,
)
from .tasks import (
    generate_cubed_sphere_points_task,
    generate_fibonacci_lattice_points_task,
    generate_cubed_sphere_cells_task,
)


def generate_points(points):
    if isinstance(points, PointGenerator):
        if points.method == PointGeneratorMethod.cubed_square:
            return gpd.GeoDataFrame.from_features(
                json.loads(
                    generate_cubed_sphere_points_task.delay(
                        points.distance, points.dict().get("mask", None)
                    ).get()
                )
            ).apply(
                lambda r: Point(
                    id=r.point_id, latitude=r.geometry.y, longitude=r.geometry.x
                ),
                axis=1,
            )
        elif points.method == PointGeneratorMethod.fibonacci_lattice:
            return gpd.GeoDataFrame.from_features(
                json.loads(
                    generate_fibonacci_lattice_points_task.delay(
                        points.distance, points.dict().get("mask", None)
                    ).get()
                )
            ).apply(
                lambda r: Point(
                    id=r.point_id, latitude=r.geometry.y, longitude=r.geometry.x
                ),
                axis=1,
            )
    return points


def generate_cells(cells):
    if isinstance(cells, CellGenerator):
        if cells.method == CellGeneratorMethod.cubed_square:
            return gpd.GeoDataFrame.from_features(
                json.loads(
                    generate_cubed_sphere_cells_task.delay(
                        cells.distance,
                        cells.dict().get("mask", None),
                        cells.dict().get("strips", None),
                    ).get()
                )
            )
    return cells

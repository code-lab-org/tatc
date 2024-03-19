from datetime import datetime, timezone
import json
import time

import geopandas as gpd
import pandas as pd

from tatc.analysis import (
    collect_observations,
    collect_downlinks,
    compute_latencies,
    reduce_latencies,
    grid_latencies,
)
from tatc.generation import generate_cubed_sphere_points, generate_cubed_sphere_cells
from tatc.schemas import (
    Point,
    GroundStation,
    WalkerConstellation,
    TwoLineElements,
    Instrument,
)

from .base import TatcTestCase
from ..latency.schemas import LatencyAnalysisRequest
from ..generation.schemas import PointGenerator, CellGenerator


class LatencyAnalysisTestCase(TatcTestCase):
    def test_analyze_latency(self):
        points = PointGenerator(
            method="cubed_square",
            distance=5000e3,
        )
        cells = CellGenerator(
            method="cubed_square",
            distance=5000e3,
        )
        instrument = Instrument(name="Test", field_of_regard=180.0)
        orbit = TwoLineElements(
            tle=[
                "1 25544U 98067A   22171.11255782  .00008307  00000+0  15444-3 0  9992",
                "2 25544  51.6448 322.0970 0003980 282.3738 231.6559 15.49798078345636",
            ]
        )
        satellite = WalkerConstellation(
            name="Test",
            orbit=orbit,
            instruments=[instrument],
            number_satellites=4,
            number_planes=2,
        )
        station = GroundStation(
            name="Station 1", latitude=0, longitude=180, min_elevation_angle=10
        )
        start = datetime(2022, 6, 1, tzinfo=timezone.utc)
        end = datetime(2022, 6, 2, tzinfo=timezone.utc)
        response = self.client.post(
            "/analyze/latency",
            json=LatencyAnalysisRequest(
                satellites=[satellite],
                stations=[station],
                start=start,
                end=end,
                points=points,
                cells=cells,
            ).model_dump_json(),
        )
        self.assertEqual(response.status_code, 200)
        task_id = response.json().get("task_id")
        while not self.client.get(f"tasks/{task_id}/status").json().get("ready"):
            time.sleep(0.1)
        response = self.client.get(f"/analyze/latency/{task_id}")
        gdf_points = gpd.GeoDataFrame.from_features(
            json.loads(json.dumps(response.json().get("points"))), crs="EPSG:4326"
        )
        gdf_points = gdf_points.reindex(columns=sorted(gdf_points.columns))
        gdf_points["latency"] = gdf_points["latency"].astype("timedelta64[ns]")
        tatc_results_downlinks = pd.concat(
            [
                collect_downlinks([station], sat, start, end)
                for sat in satellite.generate_members()
            ],
            ignore_index=True,
        )
        tatc_results_observations = pd.concat(
            [
                collect_observations(point, sat, instrument, start, end)
                for point in generate_cubed_sphere_points(5000e3).apply(
                    lambda r: Point(
                        id=r.point_id,
                        latitude=r.geometry.y,
                        longitude=r.geometry.x,
                    ),
                    axis=1,
                )
                for sat in satellite.generate_members()
            ],
            ignore_index=True,
        )
        tatc_results_points = reduce_latencies(
            compute_latencies(tatc_results_observations, tatc_results_downlinks)
        )
        tatc_results_points = tatc_results_points.reindex(
            columns=sorted(tatc_results_points.columns)
        )
        self.assertTrue(gdf_points.equals(tatc_results_points))
        gdf_cells = gpd.GeoDataFrame.from_features(
            json.loads(json.dumps(response.json().get("cells"))), crs="EPSG:4326"
        )
        gdf_cells = gdf_cells.reindex(columns=sorted(gdf_cells.columns))
        gdf_cells["latency"] = gdf_cells["latency"].astype("timedelta64[ns]")
        tatc_results_cells = grid_latencies(
            tatc_results_points, generate_cubed_sphere_cells(5000e3)
        )
        tatc_results_cells = tatc_results_cells.reindex(
            columns=sorted(tatc_results_cells.columns)
        )
        self.assertTrue(gdf_cells.equals(tatc_results_cells))

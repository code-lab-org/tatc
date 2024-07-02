from datetime import datetime, timezone
import json
import time

import geopandas as gpd
import pandas as pd

from tatc.analysis import (
    collect_multi_observations,
    aggregate_observations,
    reduce_observations,
    grid_observations,
)
from tatc.generation import (
    generate_equally_spaced_points,
    generate_equally_spaced_cells,
)
from tatc.schemas import Point, WalkerConstellation, TwoLineElements, Instrument

from .base import TatcTestCase
from ..coverage.schemas import CoverageAnalysisRequest
from ..generation.schemas import PointGenerator, CellGenerator


class CoverageAnalysisTestCase(TatcTestCase):
    def test_analyze_coverage(self):
        points = PointGenerator(
            method="equally_spaced",
            distance=5000e3,
        )
        cells = CellGenerator(
            method="equally_spaced",
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
        start = datetime(2022, 6, 1, tzinfo=timezone.utc)
        end = datetime(2022, 6, 2, tzinfo=timezone.utc)
        response = self.client.post(
            "/analyze/coverage",
            content=CoverageAnalysisRequest(
                satellites=[satellite],
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
        response = self.client.get(f"/analyze/coverage/{task_id}")
        gdf_points = gpd.GeoDataFrame.from_features(
            json.loads(json.dumps(response.json().get("points"))), crs="EPSG:4326"
        )
        gdf_points = gdf_points.reindex(columns=sorted(gdf_points.columns))
        gdf_points["revisit"] = gdf_points["revisit"].astype("timedelta64[ns]")
        gdf_points["access"] = gdf_points["access"].astype("timedelta64[ns]")
        tatc_results_points = reduce_observations(
            pd.concat(
                [
                    aggregate_observations(
                        collect_multi_observations(
                            point, satellite.generate_members(), start, end
                        )
                    )
                    for point in generate_equally_spaced_points(5000e3).apply(
                        lambda r: Point(
                            id=r.point_id, latitude=r.geometry.y, longitude=r.geometry.x
                        ),
                        axis=1,
                    )
                ],
                ignore_index=True,
            )
        )
        tatc_results_points = tatc_results_points.reindex(
            columns=sorted(tatc_results_points.columns)
        )
        self.assertTrue(gdf_points.equals(tatc_results_points))
        gdf_cells = gpd.GeoDataFrame.from_features(
            json.loads(json.dumps(response.json().get("cells"))), crs="EPSG:4326"
        )
        gdf_cells = gdf_cells.reindex(columns=sorted(gdf_cells.columns))
        gdf_cells["revisit"] = gdf_cells["revisit"].astype("timedelta64[ns]")
        gdf_cells["access"] = gdf_cells["access"].astype("timedelta64[ns]")
        tatc_results_cells = grid_observations(
            tatc_results_points, generate_equally_spaced_cells(5000e3)
        )
        tatc_results_cells = tatc_results_cells.reindex(
            columns=sorted(tatc_results_cells.columns)
        )
        self.assertTrue(gdf_cells.equals(tatc_results_cells))

import json
import time
from geojson_pydantic import FeatureCollection
import geopandas as gpd
import pandas as pd
from datetime import datetime, timedelta, timezone

from tatc.analysis import collect_orbit_track
from tatc.schemas import Point, Satellite, TwoLineElements, Instrument

from .base import TatcTestCase
from ..tracking.schemas import OrbitTrackAnalysisRequest


class AnalyzeOrbitTrackTestCase(TatcTestCase):
    def test_analyze_orbit_track(self):
        point = Point(id=0, latitude=0, longitude=0)
        instrument = Instrument(name="Test", field_of_regard=180.0)
        orbit = TwoLineElements(
            tle=[
                "1 25544U 98067A   22171.11255782  .00008307  00000+0  15444-3 0  9992",
                "2 25544  51.6448 322.0970 0003980 282.3738 231.6559 15.49798078345636",
            ]
        )
        satellite = Satellite(name="Test", orbit=orbit, instruments=[instrument])
        times = [
            datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
            for i in range(10)
        ]
        response = self.client.post(
            "/analyze/orbit-track",
            OrbitTrackAnalysisRequest(
                satellite=satellite,
                instrument=instrument,
                times=times,
            ).json(),
        )
        self.assertEqual(response.status_code, 200)
        task_id = response.json().get("task_id")
        while not self.client.get(f"tasks/{task_id}/status").json().get("ready"):
            time.sleep(0.1)
        response = self.client.get(f"/analyze/orbit-track/{task_id}")
        gdf = gpd.GeoDataFrame.from_features(
            json.loads(response.content), crs="EPSG:4326"
        )
        gdf = gdf.reindex(columns=sorted(gdf.columns))
        gdf["time"] = gdf["time"].astype("datetime64[ns, utc]")
        tatc_results = collect_orbit_track(satellite, instrument, times)
        tatc_results = tatc_results.reindex(columns=sorted(tatc_results.columns))
        self.assertTrue(gdf.equals(tatc_results))

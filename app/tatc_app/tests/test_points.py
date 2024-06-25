import json

import geopandas as gpd

from tatc.generation import (
    generate_equally_spaced_points,
    generate_fibonacci_lattice_points,
)

from .base import TatcTestCase


class GeneratePointsTestCase(TatcTestCase):
    def test_generate_points_equally_spaced(self):
        response = self.client.post(
            "/generate/points",
            json={"method": "equally_spaced", "distance": 5000e3},
        )
        self.assertEqual(response.status_code, 200)
        gdf = gpd.GeoDataFrame.from_features(
            json.loads(response.content), crs="EPSG:4326"
        )
        gdf = gdf.reindex(columns=sorted(gdf.columns))
        tatc_results = generate_equally_spaced_points(5000e3)
        tatc_results = tatc_results.reindex(columns=sorted(tatc_results.columns))
        print(tatc_results.columns)
        self.assertTrue(gdf.equals(tatc_results))

    def test_generate_points_fibonacci_lattice(self):
        response = self.client.post(
            "/generate/points",
            json={"method": "fibonacci_lattice", "distance": 5000e3},
        )
        self.assertEqual(response.status_code, 200)
        gdf = gpd.GeoDataFrame.from_features(
            json.loads(response.content), crs="EPSG:4326"
        )
        gdf = gdf.reindex(columns=sorted(gdf.columns))
        tatc_results = generate_fibonacci_lattice_points(5000e3)
        tatc_results = tatc_results.reindex(columns=sorted(tatc_results.columns))
        print(tatc_results.columns)
        self.assertTrue(gdf.equals(tatc_results))

import json

import geopandas as gpd

from tatc.generation import generate_cubed_sphere_cells

from .base import TatcTestCase


class GenerateCellsTestCase(TatcTestCase):
    def test_generate_cells_cubed_sphere(self):
        response = self.client.post(
            "/generate/cells",
            json={"method": "cubed_square", "distance": 5000e3},
        )
        self.assertEqual(response.status_code, 200)
        gdf = gpd.GeoDataFrame.from_features(
            json.loads(response.content), crs="EPSG:4326"
        )
        gdf = gdf.reindex(columns=sorted(gdf.columns))
        tatc_results = generate_cubed_sphere_cells(5000e3)
        tatc_results = tatc_results.reindex(columns=sorted(tatc_results.columns))
        self.assertTrue(gdf.equals(tatc_results))

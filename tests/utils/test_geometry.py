"""
Unit tests for the tatc.utils.geometry module.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""
import unittest

import geopandas as gpd
from shapely.geometry import MultiPolygon, Polygon

from tatc.utils import normalize_geometry, project_polygon_to_elevation, split_polygon


class TestGeometry(unittest.TestCase):
    """
    Unit tests for the tatc.utils.geometry module.
    """
    def test_project_polygon_to_elevation_polygon(self):
        """
        Test that all exterior coordinates of a polygon are assigned the
        specified elevation as a z-coordinate, leaving x/y unchanged.
        """
        polygon = Polygon([(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)])
        result = project_polygon_to_elevation(polygon, 500)
        self.assertIsInstance(result, Polygon)
        self.assertEqual(
            list(result.exterior.coords),
            [(x, y, 500) for x, y in polygon.exterior.coords],
        )

    def test_project_polygon_to_elevation_polygon_with_hole(self):
        """
        Test that both exterior and interior ring coordinates are assigned
        the specified elevation.
        """
        exterior = [(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)]
        interior = [(2, 2), (2, 4), (4, 4), (4, 2), (2, 2)]
        polygon = Polygon(exterior, [interior])
        result = project_polygon_to_elevation(polygon, 250)
        self.assertEqual(
            list(result.exterior.coords), [(x, y, 250) for x, y in exterior]
        )
        self.assertEqual(len(list(result.interiors)), 1)
        self.assertEqual(
            list(result.interiors[0].coords), [(x, y, 250) for x, y in interior]
        )

    def test_project_polygon_to_elevation_overwrites_existing_z(self):
        """
        Test that an existing z-coordinate is replaced (not offset) by the
        specified elevation.
        """
        polygon = Polygon(
            [(0, 0, 100), (10, 0, 100), (10, 10, 100), (0, 10, 100), (0, 0, 100)]
        )
        result = project_polygon_to_elevation(polygon, 50)
        self.assertEqual(
            list(result.exterior.coords), [(x, y, 50) for x, y, _ in polygon.exterior.coords]
        )

    def test_project_polygon_to_elevation_negative(self):
        """
        Test that a negative elevation (below the WGS 84 geoid) is applied as-is.
        """
        polygon = Polygon([(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)])
        result = project_polygon_to_elevation(polygon, -500)
        self.assertEqual(
            list(result.exterior.coords),
            [(x, y, -500) for x, y in polygon.exterior.coords],
        )

    def test_project_polygon_to_elevation_multipolygon(self):
        """
        Test that a multipolygon projects each constituent polygon to the
        specified elevation and preserves the number of geometries.
        """
        polygon_a = Polygon([(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)])
        polygon_b = Polygon([(20, 0), (30, 0), (30, 10), (20, 10), (20, 0)])
        multipolygon = MultiPolygon([polygon_a, polygon_b])
        result = project_polygon_to_elevation(multipolygon, 1000)
        self.assertIsInstance(result, MultiPolygon)
        self.assertEqual(len(result.geoms), 2)
        for original, projected in zip(multipolygon.geoms, result.geoms):
            self.assertEqual(
                list(projected.exterior.coords),
                [(x, y, 1000) for x, y in original.exterior.coords],
            )

    def test_split_polygon_nominal_small(self):
        """
        Test that a polygon that does not cross the antimeridian or poles is not split.
        """
        polygon = Polygon([(-10, 10), (10, 10), (10, -10), (-10, -10), (-10, 10)])
        self.assertEqual(split_polygon(polygon), polygon)

    def test_split_polygon_nominal_large(self):
        """
        Test that a polygon that does not cross the antimeridian or poles is not split.
        """
        polygon = Polygon([(-90, 10), (90, 10), (90, -10), (-90, -10), (-90, 10)])
        self.assertEqual(split_polygon(polygon), polygon)

    def test_normalize_geometry_polygon(self):
        """
        Test that a polygon is normalized to a GeoDataFrame with the correct CRS and geometry.
        """
        polygon = Polygon([(-50, 10), (50, 10), (50, -10), (-50, -10), (-50, 10)])
        result = normalize_geometry(polygon)
        self.assertIsInstance(result, gpd.GeoDataFrame)
        self.assertEqual(result.crs, "EPSG:4326")
        self.assertEqual(len(result.index), 1)
        self.assertEqual(result.iloc[0].geometry, split_polygon(polygon))

    def test_normalize_geometry_polygon_invalid(self):
        """
        Test that an invalid polygon raises a ValueError when normalized.
        """
        polygon = Polygon([(-150, 0), (-50, 0), (0, 0), (50, 0), (150, 0)])
        with self.assertRaises(ValueError):
            normalize_geometry(polygon)

    def test_normalize_geometry_multipolygon(self):
        """
        Test that a multipolygon is normalized to a GeoDataFrame with the correct CRS and geometry.
        """
        geometries = [
            [(-50, 10), (50, 10), (50, -10), (-50, -10), (-50, 10)],
            [(-50, 30), (50, 30), (50, 20), (-50, 20), (-50, 30)],
        ]
        multipolygon = MultiPolygon([[geometry, []] for geometry in geometries])
        result = normalize_geometry(multipolygon)
        self.assertIsInstance(result, gpd.GeoDataFrame)
        self.assertEqual(result.crs, "EPSG:4326")
        self.assertEqual(len(result.index), 1)
        self.assertEqual(
            result.iloc[0].geometry,
            split_polygon(multipolygon),
        )

    def test_normalize_geometry_geoseries(self):
        """
        Test that a GeoSeries is normalized to a GeoDataFrame with the correct CRS and geometry.
        """
        polygon = Polygon([(-50, 10), (50, 10), (50, -10), (-50, -10), (-50, 10)])
        gs = gpd.GeoSeries(polygon, crs="EPSG:4326")
        result = normalize_geometry(gs)
        self.assertIsInstance(result, gpd.GeoDataFrame)
        self.assertEqual(result.crs, "EPSG:4326")
        self.assertEqual(len(result.index), 1)
        self.assertEqual(
            result.iloc[0].geometry,
            split_polygon(polygon),
        )

    def test_normalize_geometry_geodataframe(self):
        """
        Test that a GeoDataFrame is normalized to a GeoDataFrame with the correct CRS and geometry.
        """
        polygon = Polygon([(-50, 10), (50, 10), (50, -10), (-50, -10), (-50, 10)])
        df = gpd.GeoDataFrame(geometry=[polygon], index=[0], crs="EPSG:4326")
        result = normalize_geometry(df)
        self.assertIsInstance(result, gpd.GeoDataFrame)
        self.assertEqual(result.crs, "EPSG:4326")
        self.assertEqual(len(result.index), 1)
        self.assertEqual(
            result.iloc[0].geometry,
            split_polygon(polygon),
        )

    def test_split_polygon_north_pole(self):
        """
        Test that a polygon that crosses the north pole is split into two polygons.
        """
        polygon = Polygon([(-50, 95), (-20, 95), (-20, 85), (-50, 85), (-50, 95)])
        result = MultiPolygon(
            [
                Polygon([(130, 85), (160, 85), (160, 90), (130, 90), (130, 85)]),
                Polygon([(-20, 90), (-20, 85), (-50, 85), (-50, 90), (-20, 90)]),
            ]
        )
        self.assertTrue(split_polygon(polygon).equals(result))

    def test_split_polygon_south_pole(self):
        """
        Test that a polygon that crosses the south pole is split into two polygons.
        """
        polygon = Polygon(
            [(-150, -95), (-120, -95), (-120, -85), (-150, -85), (-150, -95)]
        )
        result = MultiPolygon(
            [
                Polygon(
                    [(-150, -90), (-150, -85), (-120, -85), (-120, -90), (-150, -90)]
                ),
                Polygon([(30, -85), (30, -90), (60, -90), (60, -85), (30, -85)]),
            ]
        )
        self.assertTrue(split_polygon(polygon).equals(result))

    def test_split_polygon_antimeridian_short_cw(self):
        """
        Test that a polygon that crosses the antimeridian is split into two polygons.
        """
        polygon = Polygon([(170, 10), (-170, 10), (-170, -10), (170, -10), (170, 10)])
        result = MultiPolygon(
            [
                Polygon([(170, -10), (170, 10), (180, 10), (180, -10), (170, -10)]),
                Polygon(
                    [(-180, -10), (-180, 10), (-170, 10), (-170, -10), (-180, -10)]
                ),
            ]
        )
        self.assertTrue(split_polygon(polygon).equals(result))

    def test_split_polygon_antimeridian_short_ccw(self):
        """
        Test that a polygon that crosses the antimeridian in a counter-clockwise direction is split into two polygons.
        """
        polygon = Polygon([(-170, 10), (170, 10), (170, -10), (-170, -10), (-170, 10)])
        result = MultiPolygon(
            [
                Polygon([(-180, 10), (-170, 10), (-170, -10), (-180, -10), (-180, 10)]),
                Polygon([(180, -10), (170, -10), (170, 10), (180, 10), (180, -10)]),
            ]
        )
        self.assertTrue(split_polygon(polygon).equals(result))

    def test_split_polygon_antimeridian_long(self):
        """
        Test that a polygon that crosses the antimeridian multiple times is split into multiple polygons.
        """
        polygon = Polygon(
            [
                (170, 10),
                (-170, 10),
                (-70, 10),
                (30, 10),
                (30, -10),
                (-70, -10),
                (-170, -10),
                (170, -10),
                (170, 10),
            ]
        )
        result = MultiPolygon(
            [
                Polygon([(170, 10), (180, 10), (180, -10), (170, -10), (170, 10)]),
                Polygon(
                    [
                        (-180, 10),
                        (-170, 10),
                        (-70, 10),
                        (30, 10),
                        (30, -10),
                        (-70, -10),
                        (-170, -10),
                        (-180, -10),
                        (-180, 10),
                    ]
                ),
            ]
        )
        self.assertTrue(split_polygon(polygon).equals(result))

    def test_split_polygon_north_pole_multipolygon(self):
        """
        Test that a multipolygon that crosses the north pole is split into multiple polygons.
        """
        polygon = MultiPolygon(
            [
                Polygon([(-150, 95), (-120, 95), (-120, 85), (-150, 85), (-150, 95)]),
                Polygon([(150, 95), (120, 95), (120, 85), (150, 85), (150, 95)]),
            ]
        )
        result = MultiPolygon(
            [
                Polygon([(-150, 85), (-150, 90), (-120, 90), (-120, 85), (-150, 85)]),
                Polygon([(30, 90), (30, 85), (60, 85), (60, 90), (30, 90)]),
                Polygon([(120, 85), (120, 90), (150, 90), (150, 85), (120, 85)]),
                Polygon([(-60, 90), (-60, 85), (-30, 85), (-30, 90), (-60, 90)]),
            ]
        )
        self.assertTrue(split_polygon(polygon).equals(result))

    def test_split_polygon_south_pole_multipolygon(self):
        """
        Test that a multipolygon that crosses the south pole is split into multiple polygons.
        """
        polygon = MultiPolygon(
            [
                Polygon(
                    [(-150, -95), (-120, -95), (-120, -85), (-150, -85), (-150, -95)]
                ),
                Polygon([(150, -95), (120, -95), (120, -85), (150, -85), (150, -95)]),
            ]
        )
        result = MultiPolygon(
            [
                Polygon(
                    [(-150, -90), (-150, -85), (-120, -85), (-120, -90), (-150, -90)]
                ),
                Polygon([(30, -85), (30, -90), (60, -90), (60, -85), (30, -85)]),
                Polygon([(120, -90), (120, -85), (150, -85), (150, -90), (120, -90)]),
                Polygon([(-60, -85), (-60, -90), (-30, -90), (-30, -85), (-60, -85)]),
            ]
        )
        self.assertTrue(split_polygon(polygon).equals(result))

    def test_split_polygon_antimeridian_multipolygon(self):
        """
        Test that a multipolygon that crosses the antimeridian is split into multiple polygons.
        """
        polygon = MultiPolygon(
            [
                Polygon([(170, 10), (-170, 10), (-170, -10), (170, -10), (170, 10)]),
                Polygon([(170, 50), (-170, 50), (-170, 30), (170, 30), (170, 50)]),
            ]
        )
        result = MultiPolygon(
            [
                Polygon([(170, -10), (170, 10), (180, 10), (180, -10), (170, -10)]),
                Polygon(
                    [(-180, -10), (-180, 10), (-170, 10), (-170, -10), (-180, -10)]
                ),
                Polygon([(170, 30), (170, 50), (180, 50), (180, 30), (170, 30)]),
                Polygon([(-180, 30), (-180, 50), (-170, 50), (-170, 30), (-180, 30)]),
            ]
        )
        self.assertTrue(split_polygon(polygon).equals(result))

    def test_split_polygon_north_pole_top_multipolygon(self):
        """
        Test that a polygon that crosses the north pole is split into two polygons.
        """
        polygon = Polygon(
            [
                (-150, 95),
                (-140, 95),
                (-140, 85),
                (-130, 85),
                (-130, 95),
                (-120, 95),
                (-120, 85),
                (-150, 85),
                (-150, 95),
            ]
        )
        result = MultiPolygon(
            [
                Polygon([(50, 90), (60, 90), (60, 85), (50, 85), (50, 90)]),
                Polygon([(30, 90), (40, 90), (40, 85), (30, 85), (30, 90)]),
                Polygon([(-150, 85), (-150, 90), (-140, 90), (-140, 85), (-150, 85)]),
                Polygon([(-130, 85), (-130, 90), (-120, 90), (-120, 85), (-130, 85)]),
            ]
        )
        self.assertTrue(split_polygon(polygon).equals(result))

    def test_split_polygon_south_pole_bottom_multipolygon(self):
        """
        Test that a polygon that crosses the south pole is split into two polygons.
        """
        polygon = Polygon(
            [
                (-150, -95),
                (-140, -95),
                (-140, -85),
                (-130, -85),
                (-130, -95),
                (-120, -95),
                (-120, -85),
                (-150, -85),
                (-150, -95),
            ]
        )
        result = MultiPolygon(
            [
                Polygon([(50, -85), (60, -85), (60, -90), (50, -90), (50, -85)]),
                Polygon([(30, -85), (40, -85), (40, -90), (30, -90), (30, -85)]),
                Polygon(
                    [(-150, -90), (-150, -85), (-140, -85), (-140, -90), (-150, -90)]
                ),
                Polygon(
                    [(-130, -85), (-120, -85), (-120, -90), (-130, -90), (-130, -85)]
                ),
            ]
        )
        self.assertTrue(split_polygon(polygon).equals(result))

    def test_split_polygon_antimeridian_right_multipolygon(self):
        """
        Test that a polygon that crosses the antimeridian is split into two polygons.
        """
        polygon = Polygon(
            [
                (170, 10),
                (-170, 10),
                (-170, 5),
                (170, 5),
                (170, -5),
                (-170, -5),
                (-170, -10),
                (170, -10),
                (170, 10),
            ]
        )
        result = MultiPolygon(
            [
                Polygon([(170, -10), (170, -5), (180, -5), (180, -10), (170, -10)]),
                Polygon([(170, 5), (170, 10), (180, 10), (180, 5), (170, 5)]),
                Polygon(
                    [(-180, -10), (-180, -5), (-170, -5), (-170, -10), (-180, -10)]
                ),
                Polygon([(-180, 5), (-180, 10), (-170, 10), (-170, 5), (-180, 5)]),
            ]
        )
        self.assertTrue(split_polygon(polygon).equals(result))

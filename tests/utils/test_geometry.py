"""
Unit tests for the tatc.utils.geometry module.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""
import unittest

import geopandas as gpd
from shapely.geometry import GeometryCollection, LineString, MultiPolygon, Point, Polygon

from tatc.utils import normalize_geometry, project_polygon_to_elevation, split_polygon
from tatc.utils.geometry import (
    _convert_collection_to_polygon,
    _split_polygon_antimeridian,
    _split_polygon_north_pole,
    _split_polygon_south_pole,
    _wrap_polygon_over_antimeridian,
    _wrap_polygon_over_north_pole,
    _wrap_polygon_over_south_pole,
)


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

    def test_wrap_polygon_over_north_pole_no_wrap_needed(self):
        """
        Test that a polygon entirely within [-90, 90] latitude is returned unchanged.
        """
        polygon = Polygon([(-10, 80), (10, 80), (10, 85), (-10, 85), (-10, 80)])
        result = _wrap_polygon_over_north_pole(polygon)
        self.assertIs(result, polygon)

    def test_wrap_polygon_over_north_pole_exactly_90(self):
        """
        Test that a polygon with latitude exactly at (but not exceeding) 90
        degrees is treated as not needing a wrap.
        """
        polygon = Polygon([(-10, 90), (10, 90), (10, 85), (-10, 85), (-10, 90)])
        result = _wrap_polygon_over_north_pole(polygon)
        self.assertIs(result, polygon)

    def test_wrap_polygon_over_north_pole_west_side(self):
        """
        Test that a polygon entirely above the pole, on the west side of the
        prime meridian, is wrapped to the opposite (east) longitude with
        latitude reflected below 90 degrees.
        """
        polygon = Polygon([(-50, 90), (-20, 90), (-20, 95), (-50, 95), (-50, 90)])
        result = _wrap_polygon_over_north_pole(polygon)
        self.assertTrue(result.is_valid)
        self.assertEqual(
            list(result.exterior.coords),
            [(130, 90), (160, 90), (160, 85), (130, 85), (130, 90)],
        )

    def test_wrap_polygon_over_north_pole_east_side(self):
        """
        Test that a polygon entirely above the pole, on the east side of the
        prime meridian, is wrapped to the opposite (west) longitude with
        latitude reflected below 90 degrees.
        """
        polygon = Polygon([(20, 90), (50, 90), (50, 95), (20, 95), (20, 90)])
        result = _wrap_polygon_over_north_pole(polygon)
        self.assertTrue(result.is_valid)
        self.assertEqual(
            list(result.exterior.coords),
            [(-160, 90), (-130, 90), (-130, 85), (-160, 85), (-160, 90)],
        )

    def test_wrap_polygon_over_north_pole_with_hole(self):
        """
        Test that interior ring (hole) coordinates are wrapped along with
        the exterior ring.
        """
        polygon = Polygon(
            [(-50, 90), (-20, 90), (-20, 95), (-50, 95), (-50, 90)],
            [[(-45, 91), (-45, 92), (-40, 92), (-40, 91), (-45, 91)]],
        )
        result = _wrap_polygon_over_north_pole(polygon)
        self.assertEqual(
            list(result.exterior.coords),
            [(130, 90), (160, 90), (160, 85), (130, 85), (130, 90)],
        )
        self.assertEqual(
            list(result.interiors[0].coords),
            [(135, 89), (135, 88), (140, 88), (140, 89), (135, 89)],
        )

    def test_wrap_polygon_over_north_pole_multipolygon(self):
        """
        Test that each polygon in a multipolygon is wrapped independently.
        """
        polygon_west = Polygon([(-50, 90), (-20, 90), (-20, 95), (-50, 95), (-50, 90)])
        polygon_east = Polygon([(20, 90), (50, 90), (50, 95), (20, 95), (20, 90)])
        result = _wrap_polygon_over_north_pole(
            MultiPolygon([polygon_west, polygon_east])
        )
        self.assertIsInstance(result, MultiPolygon)
        self.assertEqual(len(result.geoms), 2)
        self.assertEqual(
            list(result.geoms[0].exterior.coords),
            [(130, 90), (160, 90), (160, 85), (130, 85), (130, 90)],
        )
        self.assertEqual(
            list(result.geoms[1].exterior.coords),
            [(-160, 90), (-130, 90), (-130, 85), (-160, 85), (-160, 90)],
        )

    def test_wrap_polygon_over_north_pole_invalid_result_returns_original(self):
        """
        Test the documented limitation: if a polygon straddling 90 degrees
        latitude mixes vertices on both sides of the wrap (i.e. it is not
        confined to one side of the prime meridian above 90 degrees), the
        naive per-coordinate wrap can self-intersect. In that case the
        original, unwrapped polygon is returned instead of an invalid one.
        """
        polygon = Polygon([(-50, 95), (-20, 95), (-20, 85), (-50, 85), (-50, 95)])
        result = _wrap_polygon_over_north_pole(polygon)
        self.assertIs(result, polygon)

    def test_wrap_polygon_over_north_pole_unknown_geometry(self):
        """
        Test that an unsupported geometry type raises a ValueError.
        """
        with self.assertRaises(ValueError):
            _wrap_polygon_over_north_pole(Point(0, 0))

    def test_split_polygon_north_pole_helper_no_split_needed(self):
        """
        Test that a polygon entirely within [-90, 90] latitude is returned unchanged.
        """
        polygon = Polygon([(-10, 80), (10, 80), (10, 85), (-10, 85), (-10, 80)])
        result = _split_polygon_north_pole(polygon)
        self.assertIs(result, polygon)

    def test_split_polygon_north_pole_helper_single_side_of_meridian(self):
        """
        Test that a polygon crossing the north pole, confined to one side of
        the prime meridian, is split along the pole into two valid pieces
        wrapped to their correct latitude/longitude.
        """
        polygon = Polygon([(-50, 95), (-20, 95), (-20, 85), (-50, 85), (-50, 95)])
        result = _split_polygon_north_pole(polygon)
        self.assertIsInstance(result, MultiPolygon)
        self.assertEqual(len(result.geoms), 2)
        self.assertTrue(result.is_valid)
        self.assertTrue(
            result.equals(
                MultiPolygon(
                    [
                        Polygon(
                            [(130, 85), (160, 85), (160, 90), (130, 90), (130, 85)]
                        ),
                        Polygon(
                            [(-20, 90), (-20, 85), (-50, 85), (-50, 90), (-20, 90)]
                        ),
                    ]
                )
            )
        )

    def test_split_polygon_north_pole_helper_crosses_prime_meridian(self):
        """
        Test that a polygon crossing the north pole AND straddling the prime
        meridian is additionally split along the prime meridian before
        wrapping, since a single-side wrap is not possible for the part
        that spans both sides.
        """
        polygon = Polygon([(-10, 95), (10, 95), (10, 85), (-10, 85), (-10, 95)])
        result = _split_polygon_north_pole(polygon)
        self.assertIsInstance(result, MultiPolygon)
        self.assertEqual(len(result.geoms), 3)
        self.assertTrue(result.is_valid)
        self.assertTrue(
            result.equals(
                MultiPolygon(
                    [
                        Polygon(
                            [(10, 90), (10, 85), (-10, 85), (-10, 90), (10, 90)]
                        ),
                        Polygon(
                            [(170, 85), (180, 85), (180, 90), (170, 90), (170, 85)]
                        ),
                        Polygon(
                            [
                                (-180, 85),
                                (-170, 85),
                                (-170, 90),
                                (-180, 90),
                                (-180, 85),
                            ]
                        ),
                    ]
                )
            )
        )

    def test_split_polygon_north_pole_helper_multipolygon(self):
        """
        Test that each polygon of a multipolygon is split independently and
        the results are flattened into a single multipolygon.
        """
        no_split = Polygon([(-10, 80), (10, 80), (10, 85), (-10, 85), (-10, 80)])
        needs_split = Polygon(
            [(-50, 95), (-20, 95), (-20, 85), (-50, 85), (-50, 95)]
        )
        result = _split_polygon_north_pole(MultiPolygon([no_split, needs_split]))
        self.assertIsInstance(result, MultiPolygon)
        # one unchanged part plus two parts from the split polygon
        self.assertEqual(len(result.geoms), 3)

    def test_split_polygon_north_pole_helper_unknown_geometry(self):
        """
        Test that an unsupported geometry type raises a ValueError.
        """
        with self.assertRaises(ValueError):
            _split_polygon_north_pole(Point(0, 0))

    def test_wrap_polygon_over_south_pole_no_wrap_needed(self):
        """
        Test that a polygon entirely within [-90, 90] latitude is returned unchanged.
        """
        polygon = Polygon([(-10, -80), (10, -80), (10, -85), (-10, -85), (-10, -80)])
        result = _wrap_polygon_over_south_pole(polygon)
        self.assertIs(result, polygon)

    def test_wrap_polygon_over_south_pole_exactly_minus_90(self):
        """
        Test that a polygon with latitude exactly at (but not below) -90
        degrees is treated as not needing a wrap.
        """
        polygon = Polygon([(-10, -90), (10, -90), (10, -85), (-10, -85), (-10, -90)])
        result = _wrap_polygon_over_south_pole(polygon)
        self.assertIs(result, polygon)

    def test_wrap_polygon_over_south_pole_west_side(self):
        """
        Test that a polygon entirely below the pole, on the west side of the
        prime meridian, is wrapped to the opposite (east) longitude with
        latitude reflected above -90 degrees.
        """
        polygon = Polygon([(-50, -90), (-20, -90), (-20, -95), (-50, -95), (-50, -90)])
        result = _wrap_polygon_over_south_pole(polygon)
        self.assertTrue(result.is_valid)
        self.assertEqual(
            list(result.exterior.coords),
            [(130, -90), (160, -90), (160, -85), (130, -85), (130, -90)],
        )

    def test_wrap_polygon_over_south_pole_east_side(self):
        """
        Test that a polygon entirely below the pole, on the east side of the
        prime meridian, is wrapped to the opposite (west) longitude with
        latitude reflected above -90 degrees.
        """
        polygon = Polygon([(20, -90), (50, -90), (50, -95), (20, -95), (20, -90)])
        result = _wrap_polygon_over_south_pole(polygon)
        self.assertTrue(result.is_valid)
        self.assertEqual(
            list(result.exterior.coords),
            [(-160, -90), (-130, -90), (-130, -85), (-160, -85), (-160, -90)],
        )

    def test_wrap_polygon_over_south_pole_with_hole(self):
        """
        Test that interior ring (hole) coordinates are wrapped along with
        the exterior ring.
        """
        polygon = Polygon(
            [(-50, -90), (-20, -90), (-20, -95), (-50, -95), (-50, -90)],
            [[(-45, -91), (-45, -92), (-40, -92), (-40, -91), (-45, -91)]],
        )
        result = _wrap_polygon_over_south_pole(polygon)
        self.assertEqual(
            list(result.exterior.coords),
            [(130, -90), (160, -90), (160, -85), (130, -85), (130, -90)],
        )
        self.assertEqual(
            list(result.interiors[0].coords),
            [(135, -89), (135, -88), (140, -88), (140, -89), (135, -89)],
        )

    def test_wrap_polygon_over_south_pole_multipolygon(self):
        """
        Test that each polygon in a multipolygon is wrapped independently.
        """
        polygon_west = Polygon(
            [(-50, -90), (-20, -90), (-20, -95), (-50, -95), (-50, -90)]
        )
        polygon_east = Polygon(
            [(20, -90), (50, -90), (50, -95), (20, -95), (20, -90)]
        )
        result = _wrap_polygon_over_south_pole(
            MultiPolygon([polygon_west, polygon_east])
        )
        self.assertIsInstance(result, MultiPolygon)
        self.assertEqual(len(result.geoms), 2)
        self.assertEqual(
            list(result.geoms[0].exterior.coords),
            [(130, -90), (160, -90), (160, -85), (130, -85), (130, -90)],
        )
        self.assertEqual(
            list(result.geoms[1].exterior.coords),
            [(-160, -90), (-130, -90), (-130, -85), (-160, -85), (-160, -90)],
        )

    def test_wrap_polygon_over_south_pole_invalid_result_returns_original(self):
        """
        Test the documented limitation: if a polygon straddling -90 degrees
        latitude mixes vertices on both sides of the wrap (i.e. it is not
        confined to one side of the prime meridian below -90 degrees), the
        naive per-coordinate wrap can self-intersect. In that case the
        original, unwrapped polygon is returned instead of an invalid one.
        """
        polygon = Polygon([(-50, -95), (-20, -95), (-20, -85), (-50, -85), (-50, -95)])
        result = _wrap_polygon_over_south_pole(polygon)
        self.assertIs(result, polygon)

    def test_wrap_polygon_over_south_pole_unknown_geometry(self):
        """
        Test that an unsupported geometry type raises a ValueError.
        """
        with self.assertRaises(ValueError):
            _wrap_polygon_over_south_pole(Point(0, 0))

    def test_split_polygon_south_pole_helper_no_split_needed(self):
        """
        Test that a polygon entirely within [-90, 90] latitude is returned unchanged.
        """
        polygon = Polygon([(-10, -80), (10, -80), (10, -85), (-10, -85), (-10, -80)])
        result = _split_polygon_south_pole(polygon)
        self.assertIs(result, polygon)

    def test_split_polygon_south_pole_helper_single_side_of_meridian(self):
        """
        Test that a polygon crossing the south pole, confined to one side of
        the prime meridian, is split along the pole into two valid pieces
        wrapped to their correct latitude/longitude.
        """
        polygon = Polygon(
            [(-150, -95), (-120, -95), (-120, -85), (-150, -85), (-150, -95)]
        )
        result = _split_polygon_south_pole(polygon)
        self.assertIsInstance(result, MultiPolygon)
        self.assertEqual(len(result.geoms), 2)
        self.assertTrue(result.is_valid)
        self.assertTrue(
            result.equals(
                MultiPolygon(
                    [
                        Polygon(
                            [(-150, -90), (-150, -85), (-120, -85), (-120, -90), (-150, -90)]
                        ),
                        Polygon([(30, -85), (30, -90), (60, -90), (60, -85), (30, -85)]),
                    ]
                )
            )
        )

    def test_split_polygon_south_pole_helper_crosses_prime_meridian(self):
        """
        Test that a polygon crossing the south pole AND straddling the prime
        meridian is additionally split along the prime meridian before
        wrapping, since a single-side wrap is not possible for the part
        that spans both sides.
        """
        polygon = Polygon([(-10, -95), (10, -95), (10, -85), (-10, -85), (-10, -95)])
        result = _split_polygon_south_pole(polygon)
        self.assertIsInstance(result, MultiPolygon)
        self.assertEqual(len(result.geoms), 3)
        self.assertTrue(result.is_valid)
        self.assertTrue(
            result.equals(
                MultiPolygon(
                    [
                        Polygon(
                            [(-10, -90), (-10, -85), (10, -85), (10, -90), (-10, -90)]
                        ),
                        Polygon(
                            [
                                (-170, -90),
                                (-170, -85),
                                (-180, -85),
                                (-180, -90),
                                (-170, -90),
                            ]
                        ),
                        Polygon(
                            [(180, -85), (170, -85), (170, -90), (180, -90), (180, -85)]
                        ),
                    ]
                )
            )
        )

    def test_split_polygon_south_pole_helper_multipolygon(self):
        """
        Test that each polygon of a multipolygon is split independently and
        the results are flattened into a single multipolygon.
        """
        no_split = Polygon([(-10, -80), (10, -80), (10, -85), (-10, -85), (-10, -80)])
        needs_split = Polygon(
            [(-150, -95), (-120, -95), (-120, -85), (-150, -85), (-150, -95)]
        )
        result = _split_polygon_south_pole(MultiPolygon([no_split, needs_split]))
        self.assertIsInstance(result, MultiPolygon)
        # one unchanged part plus two parts from the split polygon
        self.assertEqual(len(result.geoms), 3)

    def test_split_polygon_south_pole_helper_unknown_geometry(self):
        """
        Test that an unsupported geometry type raises a ValueError.
        """
        with self.assertRaises(ValueError):
            _split_polygon_south_pole(Point(0, 0))

    def test_wrap_polygon_over_antimeridian_no_wrap_needed(self):
        """
        Test that a polygon entirely within [-180, 180] longitude is returned unchanged.
        """
        polygon = Polygon([(-170, 10), (170, 10), (170, -10), (-170, -10), (-170, 10)])
        result = _wrap_polygon_over_antimeridian(polygon)
        self.assertIs(result, polygon)

    def test_wrap_polygon_over_antimeridian_exactly_boundary(self):
        """
        Test that a polygon with longitude exactly at (but not exceeding)
        -180/180 degrees is treated as not needing a wrap.
        """
        polygon = Polygon([(-180, 10), (180, 10), (180, -10), (-180, -10), (-180, 10)])
        result = _wrap_polygon_over_antimeridian(polygon)
        self.assertIs(result, polygon)

    def test_wrap_polygon_over_antimeridian_negative_side(self):
        """
        Test that a polygon entirely at or below -180 degrees longitude is
        wrapped by adding 360 degrees.
        """
        polygon = Polygon(
            [(-200, 10), (-190, 10), (-190, -10), (-200, -10), (-200, 10)]
        )
        result = _wrap_polygon_over_antimeridian(polygon)
        self.assertTrue(result.is_valid)
        self.assertEqual(
            list(result.exterior.coords),
            [(160, 10), (170, 10), (170, -10), (160, -10), (160, 10)],
        )

    def test_wrap_polygon_over_antimeridian_positive_side(self):
        """
        Test that a polygon entirely at or above 180 degrees longitude is
        wrapped by subtracting 360 degrees.
        """
        polygon = Polygon([(190, 10), (200, 10), (200, -10), (190, -10), (190, 10)])
        result = _wrap_polygon_over_antimeridian(polygon)
        self.assertTrue(result.is_valid)
        self.assertEqual(
            list(result.exterior.coords),
            [(-170, 10), (-160, 10), (-160, -10), (-170, -10), (-170, 10)],
        )

    def test_wrap_polygon_over_antimeridian_with_hole(self):
        """
        Test that interior ring (hole) coordinates are wrapped along with
        the exterior ring.
        """
        polygon = Polygon(
            [(-200, 10), (-190, 10), (-190, -10), (-200, -10), (-200, 10)],
            [[(-198, 2), (-198, 4), (-196, 4), (-196, 2), (-198, 2)]],
        )
        result = _wrap_polygon_over_antimeridian(polygon)
        self.assertEqual(
            list(result.exterior.coords),
            [(160, 10), (170, 10), (170, -10), (160, -10), (160, 10)],
        )
        self.assertEqual(
            list(result.interiors[0].coords),
            [(162, 2), (162, 4), (164, 4), (164, 2), (162, 2)],
        )

    def test_wrap_polygon_over_antimeridian_multipolygon(self):
        """
        Test that each polygon in a multipolygon is wrapped independently.
        """
        polygon_neg = Polygon(
            [(-200, 10), (-190, 10), (-190, -10), (-200, -10), (-200, 10)]
        )
        polygon_pos = Polygon(
            [(190, 10), (200, 10), (200, -10), (190, -10), (190, 10)]
        )
        result = _wrap_polygon_over_antimeridian(
            MultiPolygon([polygon_neg, polygon_pos])
        )
        self.assertIsInstance(result, MultiPolygon)
        self.assertEqual(len(result.geoms), 2)
        self.assertEqual(
            list(result.geoms[0].exterior.coords),
            [(160, 10), (170, 10), (170, -10), (160, -10), (160, 10)],
        )
        self.assertEqual(
            list(result.geoms[1].exterior.coords),
            [(-170, 10), (-160, 10), (-160, -10), (-170, -10), (-170, 10)],
        )

    def test_wrap_polygon_over_antimeridian_straddles_both_sides_returns_original(
        self,
    ):
        """
        Test the documented limitation: if a polygon has coordinates both at
        or below -180 degrees and at or above 180 degrees, neither wrap
        direction applies uniformly, so the original polygon is returned
        unchanged.
        """
        polygon = Polygon(
            [(-190, 10), (190, 10), (190, -10), (-190, -10), (-190, 10)]
        )
        result = _wrap_polygon_over_antimeridian(polygon)
        self.assertIs(result, polygon)

    def test_wrap_polygon_over_antimeridian_unknown_geometry(self):
        """
        Test that an unsupported geometry type raises a ValueError.
        """
        with self.assertRaises(ValueError):
            _wrap_polygon_over_antimeridian(Point(0, 0))

    def test_convert_collection_to_polygon_single_polygon(self):
        """
        Test that a collection containing a single polygon returns that
        polygon directly, not wrapped in a MultiPolygon.
        """
        polygon = Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])
        result = _convert_collection_to_polygon(GeometryCollection([polygon]))
        self.assertIsInstance(result, Polygon)
        self.assertTrue(result.equals(polygon))

    def test_convert_collection_to_polygon_multiple_polygons(self):
        """
        Test that a collection containing multiple polygons returns a
        MultiPolygon preserving all of them.
        """
        polygon_a = Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])
        polygon_b = Polygon([(2, 0), (3, 0), (3, 1), (2, 1), (2, 0)])
        result = _convert_collection_to_polygon(
            GeometryCollection([polygon_a, polygon_b])
        )
        self.assertIsInstance(result, MultiPolygon)
        self.assertEqual(len(result.geoms), 2)
        self.assertTrue(result.geoms[0].equals(polygon_a))
        self.assertTrue(result.geoms[1].equals(polygon_b))

    def test_convert_collection_to_polygon_flattens_multipolygon_member(self):
        """
        Test that a MultiPolygon member of the collection is flattened so
        its constituent polygons are included individually.
        """
        polygon_a = Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])
        polygon_b = Polygon([(2, 0), (3, 0), (3, 1), (2, 1), (2, 0)])
        polygon_c = Polygon([(4, 0), (5, 0), (5, 1), (4, 1), (4, 0)])
        result = _convert_collection_to_polygon(
            GeometryCollection([polygon_a, MultiPolygon([polygon_b, polygon_c])])
        )
        self.assertIsInstance(result, MultiPolygon)
        self.assertEqual(len(result.geoms), 3)
        self.assertTrue(result.geoms[0].equals(polygon_a))
        self.assertTrue(result.geoms[1].equals(polygon_b))
        self.assertTrue(result.geoms[2].equals(polygon_c))

    def test_convert_collection_to_polygon_drops_points_and_lines(self):
        """
        Test that non-polygon geometries (points, lines) in the collection
        are dropped, leaving only the polygon.
        """
        polygon = Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])
        result = _convert_collection_to_polygon(
            GeometryCollection(
                [polygon, Point(0, 0), LineString([(0, 0), (1, 1)])]
            )
        )
        self.assertIsInstance(result, Polygon)
        self.assertTrue(result.equals(polygon))

    def test_convert_collection_to_polygon_no_polygons(self):
        """
        Test that a collection with no polygon members returns an empty MultiPolygon.
        """
        result = _convert_collection_to_polygon(
            GeometryCollection([Point(0, 0), LineString([(0, 0), (1, 1)])])
        )
        self.assertIsInstance(result, MultiPolygon)
        self.assertTrue(result.is_empty)

    def test_split_polygon_antimeridian_helper_no_split_needed(self):
        """
        Test that a polygon whose adjacent vertices never jump by 180
        degrees or more of longitude is returned unchanged.
        """
        polygon = Polygon([(-10, 10), (10, 10), (10, -10), (-10, -10), (-10, 10)])
        result = _split_polygon_antimeridian(polygon)
        self.assertIs(result, polygon)

    def test_split_polygon_antimeridian_helper_short_cw(self):
        """
        Test that a polygon crossing the antimeridian clockwise is split
        into two valid polygons on either side of it.
        """
        polygon = Polygon([(170, 10), (-170, 10), (-170, -10), (170, -10), (170, 10)])
        result = _split_polygon_antimeridian(polygon)
        self.assertTrue(result.is_valid)
        self.assertTrue(
            result.equals(
                MultiPolygon(
                    [
                        Polygon(
                            [(170, -10), (170, 10), (180, 10), (180, -10), (170, -10)]
                        ),
                        Polygon(
                            [
                                (-180, -10),
                                (-180, 10),
                                (-170, 10),
                                (-170, -10),
                                (-180, -10),
                            ]
                        ),
                    ]
                )
            )
        )

    def test_split_polygon_antimeridian_helper_short_ccw(self):
        """
        Test that a polygon crossing the antimeridian counter-clockwise is
        split into two valid polygons on either side of it.
        """
        polygon = Polygon([(-170, 10), (170, 10), (170, -10), (-170, -10), (-170, 10)])
        result = _split_polygon_antimeridian(polygon)
        self.assertTrue(result.is_valid)
        self.assertTrue(
            result.equals(
                MultiPolygon(
                    [
                        Polygon(
                            [(-180, 10), (-170, 10), (-170, -10), (-180, -10), (-180, 10)]
                        ),
                        Polygon(
                            [(180, -10), (170, -10), (170, 10), (180, 10), (180, -10)]
                        ),
                    ]
                )
            )
        )

    def test_split_polygon_antimeridian_helper_multiple_crossings(self):
        """
        Test that a polygon crossing the antimeridian multiple times is
        split into valid polygons on either side of it.
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
        result = _split_polygon_antimeridian(polygon)
        self.assertTrue(result.is_valid)
        self.assertTrue(
            result.equals(
                MultiPolygon(
                    [
                        Polygon(
                            [(170, 10), (180, 10), (180, -10), (170, -10), (170, 10)]
                        ),
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
            )
        )

    def test_split_polygon_antimeridian_helper_multipolygon(self):
        """
        Test that each polygon of a multipolygon is split independently and
        the results are flattened into a single multipolygon.
        """
        no_split = Polygon([(-10, 10), (10, 10), (10, -10), (-10, -10), (-10, 10)])
        needs_split = Polygon(
            [(170, 10), (-170, 10), (-170, -10), (170, -10), (170, 10)]
        )
        result = _split_polygon_antimeridian(MultiPolygon([no_split, needs_split]))
        self.assertIsInstance(result, MultiPolygon)
        # one unchanged part plus two parts from the split polygon
        self.assertEqual(len(result.geoms), 3)

    def test_split_polygon_antimeridian_helper_contains_pole(self):
        """
        Test the "contains a pole" special case: a polygon whose vertex
        longitudes wrap all the way around the globe (e.g. a polar cap) is
        reconstructed with a flattened pole edge and split along the prime
        meridian, rather than treated as a simple antimeridian crossing.

        Note: the raw result of this branch touches along the shared
        prime-meridian cut edge, so the combined MultiPolygon is reported
        as invalid even though each individual piece is valid; the public
        `split_polygon` function repairs this via `shapely.make_valid`.
        """
        polygon = Polygon([(45, 80), (135, 80), (-135, 80), (-45, 80), (45, 80)])
        result = _split_polygon_antimeridian(polygon)
        self.assertIsInstance(result, MultiPolygon)
        self.assertEqual(len(result.geoms), 2)
        for part in result.geoms:
            self.assertTrue(part.is_valid)
        self.assertFalse(result.is_valid)
        self.assertEqual(
            list(result.geoms[0].exterior.coords),
            [(0, 80), (-45, 80), (-135, 80), (-180, 80), (-180, 90), (0, 90), (0, 80)],
        )
        self.assertEqual(
            list(result.geoms[1].exterior.coords),
            [(0, 90), (180, 90), (180, 80), (135, 80), (45, 80), (0, 80), (0, 90)],
        )

    def test_split_polygon_antimeridian_helper_unknown_geometry(self):
        """
        Test that an unsupported geometry type raises a ValueError.
        """
        with self.assertRaises(ValueError):
            _split_polygon_antimeridian(Point(0, 0))

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

    def test_split_polygon_polar_cap_repairs_invalid_antimeridian_split(self):
        """
        Test that a polygon encircling a pole (all longitudes, without
        exceeding +/-90 degrees latitude) is split into a valid geometry.
        The antimeridian split's raw result touches itself along the
        prime-meridian cut, so this exercises the `make_valid` repair.
        """
        polygon = Polygon([(45, 80), (135, 80), (-135, 80), (-45, 80), (45, 80)])
        result = split_polygon(polygon)
        self.assertTrue(result.is_valid)
        self.assertTrue(result.contains(Point(100, 85)))
        self.assertFalse(result.contains(Point(100, 75)))

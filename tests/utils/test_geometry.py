"""
Unit tests for the tatc.utils.geometry module.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""
import unittest

import geopandas as gpd
from shapely.geometry import MultiPolygon, Point, Polygon

from tatc.utils import (
    geodesic_distance,
    get_planar_bounds,
    normalize_geometry,
    project_polygon_to_elevation,
    split_polygon,
)


class TestGeometry(unittest.TestCase):  # pylint: disable=too-many-public-methods
    """
    Unit tests for the tatc.utils.geometry module.
    """
    def test_geodesic_distance_same_point(self):
        """
        Test that the geodesic distance between a point and itself is zero.
        """
        self.assertAlmostEqual(geodesic_distance(10, 20, 10, 20), 0, delta=1e-6)

    def test_geodesic_distance_one_degree_at_equator(self):
        """
        Test the geodesic distance for one degree of longitude along the
        equator against the known WGS 84 equatorial circumference (the
        equator is itself a geodesic, so this distance is exact):
        2 * pi * EARTH_EQUATORIAL_RADIUS / 360 = 111319.4908 meters.
        """
        self.assertAlmostEqual(
            geodesic_distance(0, 0, 1, 0), 111319.4908, delta=0.01
        )

    def test_geodesic_distance_pole_to_equator(self):
        """
        Test the geodesic distance from the North pole to the equator
        against the published WGS 84 meridian quadrant length (a meridian
        is itself a geodesic, so this distance is exact): 10001965.7293
        meters.
        """
        self.assertAlmostEqual(
            geodesic_distance(0, 90, 0, 0), 10001965.7293, delta=0.01
        )

    def test_geodesic_distance_is_symmetric(self):
        """
        Test that the geodesic distance does not depend on point order.
        """
        self.assertAlmostEqual(
            geodesic_distance(-73.9857, 40.7484, -0.1278, 51.5074),
            geodesic_distance(-0.1278, 51.5074, -73.9857, 40.7484),
            delta=1e-6,
        )

    def test_geodesic_distance_across_antimeridian(self):
        """
        Test that the geodesic distance between points on either side of
        the antimeridian takes the short way across it, rather than the
        long way around through the prime meridian.
        """
        self.assertAlmostEqual(
            geodesic_distance(179, 0, -179, 0), 222638.9816, delta=0.01
        )

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

    def test_split_polygon_nominal_boundary(self):
        """
        Test that a polygon touching (but not exceeding) 90 degrees
        latitude is not split.
        """
        polygon = Polygon([(-10, 90), (10, 90), (10, 80), (-10, 80), (-10, 90)])
        self.assertEqual(split_polygon(polygon), polygon)

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

    def test_split_polygon_north_pole_crosses_prime_meridian(self):
        """
        Test that a polygon crossing the north pole AND straddling the
        prime meridian is split into three valid polygons: the part below
        the pole, and the two halves above it wrapped to either side.
        """
        polygon = Polygon([(-10, 95), (10, 95), (10, 85), (-10, 85), (-10, 95)])
        result = MultiPolygon(
            [
                Polygon([(10, 90), (10, 85), (-10, 85), (-10, 90), (10, 90)]),
                Polygon([(170, 85), (180, 85), (180, 90), (170, 90), (170, 85)]),
                Polygon(
                    [(-180, 85), (-170, 85), (-170, 90), (-180, 90), (-180, 85)]
                ),
            ]
        )
        self.assertTrue(split_polygon(polygon).equals(result))

    def test_split_polygon_north_pole_with_hole(self):
        """
        Test that an interior ring (hole) survives the split, wrapped
        along with the piece that contains it.
        """
        polygon = Polygon(
            [(-50, 95), (-20, 95), (-20, 85), (-50, 85), (-50, 95)],
            [[(-45, 91), (-45, 93), (-40, 93), (-40, 91), (-45, 91)]],
        )
        result = split_polygon(polygon)
        self.assertIsInstance(result, MultiPolygon)
        self.assertTrue(result.is_valid)
        wrapped_part = next(g for g in result.geoms if len(g.interiors) > 0)
        self.assertEqual(
            list(wrapped_part.interiors[0].coords),
            [(135, 89), (140, 89), (140, 87), (135, 87), (135, 89)],
        )

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

    def test_split_polygon_south_pole_crosses_prime_meridian(self):
        """
        Test that a polygon crossing the south pole AND straddling the
        prime meridian is split into three valid polygons: the part above
        the pole, and the two halves below it wrapped to either side.
        """
        polygon = Polygon([(-10, -95), (10, -95), (10, -85), (-10, -85), (-10, -95)])
        result = MultiPolygon(
            [
                Polygon([(-10, -90), (-10, -85), (10, -85), (10, -90), (-10, -90)]),
                Polygon(
                    [
                        (-170, -90),
                        (-170, -85),
                        (-180, -85),
                        (-180, -90),
                        (-170, -90),
                    ]
                ),
                Polygon([(180, -85), (170, -85), (170, -90), (180, -90), (180, -85)]),
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
        Test that a polygon that crosses the antimeridian in a
        counter-clockwise direction is split into two polygons.
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
        Test that a polygon that crosses the antimeridian multiple times is
        split into multiple polygons.
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

    def test_split_polygon_antimeridian_with_hole(self):
        """
        Test that an interior ring (hole) survives the antimeridian split,
        remaining in the piece that contains it.
        """
        polygon = Polygon(
            [(170, 10), (-170, 10), (-170, -10), (170, -10), (170, 10)],
            [[(175, 2), (175, 4), (177, 4), (177, 2), (175, 2)]],
        )
        result = split_polygon(polygon)
        self.assertIsInstance(result, MultiPolygon)
        self.assertTrue(result.is_valid)
        wrapped_part = next(g for g in result.geoms if len(g.interiors) > 0)
        self.assertEqual(
            list(wrapped_part.interiors[0].coords),
            [(175, 2), (177, 2), (177, 4), (175, 4), (175, 2)],
        )

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

    def test_split_polygon_polar_cap_drops_z_dimension(self):
        """
        Test that a polar cap polygon with a z-dimension does not raise an
        error: the antimeridian split's pole-containing case reconstructs
        the ring using synthetic 2D corner points, so any input z must be
        discarded rather than mixed in.
        """
        polygon = Polygon(
            [
                (45, 80, 100),
                (135, 80, 100),
                (-135, 80, 100),
                (-45, 80, 100),
                (45, 80, 100),
            ]
        )
        result = split_polygon(polygon)
        self.assertTrue(result.is_valid)
        self.assertFalse(result.has_z)

    def test_split_polygon_unknown_geometry(self):
        """
        Test that an unsupported geometry type raises a ValueError.
        """
        with self.assertRaises(ValueError):
            split_polygon(Point(0, 0))

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

    def test_get_planar_bounds_no_mask_returns_global_extent(self):
        """
        Test that omitting a mask returns the full global longitude/latitude
        extent.
        """
        self.assertEqual(get_planar_bounds(None), (-180, -90, 180, 90))

    def test_get_planar_bounds_polygon_mask_returns_its_bounds(self):
        """
        Test that a Polygon mask's bounds are returned directly.
        """
        mask = Polygon([[-100, 25], [-50, 25], [-50, -25], [-100, -25], [-100, 25]])
        self.assertEqual(get_planar_bounds(mask), (-100, -25, -50, 25))

    def test_get_planar_bounds_multipolygon_mask_returns_combined_bounds(self):
        """
        Test that a MultiPolygon mask's bounds span all its constituent
        polygons.
        """
        mask = MultiPolygon(
            [
                Polygon([[-100, 0], [-90, 0], [-90, 10], [-100, 10], [-100, 0]]),
                Polygon([[50, -20], [60, -20], [60, -10], [50, -10], [50, -20]]),
            ]
        )
        self.assertEqual(get_planar_bounds(mask), (-100, -20, 60, 10))

    def test_get_planar_bounds_mask_with_max_longitude_exactly_negative_180_expands_to_180(
        self,
    ):
        """
        Known-limitation regression test (not a full antimeridian fix, see
        get_planar_bounds's docstring): when a mask's raw bounding box has
        its maximum longitude exactly -180 (e.g. a mask expressed with
        unwrapped-negative longitude spanning -190 to -180), max_longitude
        is expanded to 180. This only patches max_longitude; min_longitude
        (-190 here) is left as-is, so the resulting bounds do not
        coherently describe the mask's actual extent.
        """
        mask = Polygon([[-190, -10], [-180, -10], [-180, 10], [-190, 10], [-190, -10]])
        min_longitude, _, max_longitude, _ = get_planar_bounds(mask)
        self.assertEqual(max_longitude, 180)
        self.assertEqual(min_longitude, -190)

    def test_get_planar_bounds_invalid_polygon_raises_value_error(self):
        """
        Test that a self-intersecting (invalid) Polygon mask raises a
        ValueError.
        """
        mask = Polygon([[-100, 25], [-50, 25], [-100, -25], [-50, -25], [-100, 25]])
        with self.assertRaises(ValueError):
            get_planar_bounds(mask)

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

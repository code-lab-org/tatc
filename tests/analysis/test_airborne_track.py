"""
Unit tests for airborne ground track analysis functions.

The expected values come from formulas (the geometry of a scan over a flat
surface, rotations of the scan rays, and the SPICE intersection of a ray with
the WGS 84 ellipsoid), not from the results of an earlier implementation.

@author: Naeun Kim <hkim439@asu.edu>
"""

import importlib
import inspect
import sys
import unittest
from unittest import mock

import numpy as np
import pandas as pd
from pyproj import Geod
from shapely.geometry import Polygon
from spiceypy.spiceypy import recgeo, surfpt

from tatc.analysis import collect_airborne_ground_track, compute_airborne_edge_points
from tatc.analysis.airborne_track import (
    _get_flat_edge_points,
    _get_flat_scan_offsets,
    _get_ned_to_earth_fixed,
)
from tatc.utils.ellipsoid import (
    compute_ellipsoid_intersection,
    geodetic_to_rectangular,
    rectangular_to_geodetic,
)
from tatc.constants import (
    EARTH_EQUATORIAL_RADIUS,
    EARTH_FLATTENING,
    EARTH_POLAR_RADIUS,
)

GEOD = Geod(ellps="WGS84")

LEVEL_COLUMNS = ["time", "aircraft", "instrument", "roll", "pitch", "valid_obs"]
STRIP_COLUMNS = [
    "time",
    "aircraft",
    "instrument",
    "roll",
    "pitch",
    "heading",
    "valid_obs",
]


def _trajectory(
    n=60,
    altitude=10000.0,
    speed=200.0,
    track=30.0,
    roll=0.0,
    pitch=0.0,
    heading=None,
):
    """
    Creates the trajectory of a straight flight at constant altitude, speed,
    and attitude, sampled every second from 33.4 N, 111.9 W.
    """
    longitude, latitude, _ = GEOD.fwd(
        np.full(n, -111.9),
        np.full(n, 33.4),
        np.full(n, track),
        np.arange(n) * speed,
    )
    trajectory = pd.DataFrame(
        {
            "time": pd.date_range(
                "2026-10-15 17:00:00", periods=n, freq="1s", tz="UTC"
            ),
            "latitude": latitude,
            "longitude": longitude,
            "altitude": altitude,
            "track": track,
            "ground_speed": speed,
            "roll": roll,
            "pitch": pitch,
        }
    )
    if heading is not None:
        trajectory["heading"] = heading
    return trajectory


def _areas(footprints):
    """
    Gets the areas (square meters) of footprints in a projected (UTM) CRS.
    """
    return footprints.to_crs(footprints.estimate_utm_crs()).geometry.area.to_numpy()


def _centroid_offsets(trajectory, footprints, index):
    """
    Gets the distances (meters) of a footprint's centroid from the aircraft
    along and to the right of the aircraft's track (flat-surface scale).
    """
    row = trajectory.iloc[index]
    centroid = footprints.geometry.iloc[index].centroid
    azimuth, _, distance = GEOD.inv(
        row["longitude"], row["latitude"], centroid.x, centroid.y
    )
    relative = np.radians(azimuth - row["track"])
    return distance * np.cos(relative), distance * np.sin(relative)


class TestAirborneTrack(unittest.TestCase):
    """
    Unit tests for `collect_airborne_ground_track`.
    """

    def test_empty_trajectory(self):
        """
        Tests that an empty trajectory gives an empty result.
        """
        footprints = collect_airborne_ground_track(pd.DataFrame(), 25.0, 1.0)
        self.assertTrue(footprints.empty)
        self.assertEqual(footprints.crs.to_epsg(), 4326)
        self.assertEqual(list(footprints.columns), LEVEL_COLUMNS + ["geometry"])

    def test_invalid_arguments(self):
        """
        Tests that invalid trajectories and arguments are rejected.
        """
        trajectory = _trajectory(n=5)
        with self.assertRaisesRegex(ValueError, "altitude"):
            collect_airborne_ground_track(trajectory.drop(columns="altitude"), 25, 1)
        for interval in (0, -1):
            with self.assertRaises(ValueError):
                collect_airborne_ground_track(trajectory, 25.0, interval)
        for angle in (0, -5, 90, 120):
            with self.assertRaises(ValueError):
                collect_airborne_ground_track(trajectory, angle, 1.0)
        with self.assertRaisesRegex(ValueError, "ellipsoid"):
            collect_airborne_ground_track(trajectory, 25.0, 1.0, geometry="globe")

    def test_spice_alias(self):
        """
        Tests that the `"spice"` geometry (the name in the earlier prototype)
        is an alias of `"ellipsoid"`: the same results with or without attitude.
        """
        trajectory = _trajectory(n=20, roll=15.0, pitch=3.0, heading=5.0)
        trajectory.loc[5:8, "altitude"] = -10.0
        for use_attitude in (False, True):
            alias = collect_airborne_ground_track(
                trajectory, 25.0, 1.0, use_attitude=use_attitude, geometry="spice"
            )
            expected = collect_airborne_ground_track(
                trajectory, 25.0, 1.0, use_attitude=use_attitude, geometry="ellipsoid"
            )
            self.assertEqual(len(alias), 20)
            self.assertFalse(alias["valid_obs"].all())
            pd.testing.assert_frame_equal(
                pd.DataFrame(alias.drop(columns="geometry")),
                pd.DataFrame(expected.drop(columns="geometry")),
            )
            for left, right in zip(alias.geometry, expected.geometry):
                self.assertTrue(left.equals_exact(right, 0.0))

    def test_level_rectangle_dimensions(self):
        """
        Tests that the level footprint is a rectangle as wide as the swath
        of a flat scan (2 h tan(a)) and as long as the distance flown.
        """
        footprints = collect_airborne_ground_track(
            _trajectory(altitude=10000.0, speed=200.0), 25.0, 1.0, "AC", "SCANNER"
        )
        self.assertEqual(list(footprints.columns), LEVEL_COLUMNS + ["geometry"])
        self.assertEqual(set(footprints["aircraft"]), {"AC"})
        self.assertEqual(set(footprints["instrument"]), {"SCANNER"})
        self.assertTrue(footprints["valid_obs"].all())
        swath = 2 * 10000.0 * np.tan(np.radians(25.0))
        for geometry in footprints.geometry:
            self.assertIsInstance(geometry, Polygon)
            corners = np.array(geometry.exterior.coords)
            sides = sorted(GEOD.inv(*corners[i], *corners[i + 1])[2] for i in range(4))
            self.assertAlmostEqual(sides[0], 200.0, delta=0.5)
            self.assertAlmostEqual(sides[1], 200.0, delta=0.5)
            self.assertAlmostEqual(sides[2], swath, delta=2.0)
            self.assertAlmostEqual(sides[3], swath, delta=2.0)
        # centered on the aircraft (the centroid is taken in longitude and
        # latitude, which is off by a few meters at this scale)
        for index in (0, 20, 59):
            along, across = _centroid_offsets(_trajectory(), footprints, index)
            self.assertAlmostEqual(along, 0.0, delta=3.0)
            self.assertAlmostEqual(across, 0.0, delta=3.0)

    def test_level_rectangle_does_not_apply_attitude(self):
        """
        Tests that the level rectangle reports but does not apply the roll
        and pitch angles of the trajectory.
        """
        level = collect_airborne_ground_track(_trajectory(), 25.0, 1.0)
        banked = collect_airborne_ground_track(_trajectory(roll=20.0, pitch=5.0), 25, 1)
        self.assertTrue((banked["roll"] == 20.0).all())
        self.assertTrue((banked["pitch"] == 5.0).all())
        for left, right in zip(level.geometry, banked.geometry):
            self.assertTrue(left.equals_exact(right, 1e-12))

    def test_flat_attitude_shifts_swath(self):
        """
        Tests that, over a flat surface, a roll shifts the center of the swath
        to the outside of a turn by h (tan(a - roll) + tan(-a - roll)) / 2 /
        cos(pitch) and a pitch shifts it forward by h tan(pitch).
        """
        altitude, half_angle, step = 10000.0, 25.0, 200.0
        for roll, pitch in ((20.0, 0.0), (-20.0, 0.0), (0.0, 5.0), (15.0, -4.0)):
            trajectory = _trajectory(altitude=altitude, roll=roll, pitch=pitch)
            footprints = collect_airborne_ground_track(
                trajectory, half_angle, 1.0, use_attitude=True
            )
            self.assertEqual(list(footprints.columns[:7]), STRIP_COLUMNS)
            self.assertTrue(footprints["valid_obs"].all())
            a, phi, theta = (np.radians(x) for x in (half_angle, roll, pitch))
            expected_across = (
                altitude * (np.tan(a - phi) + np.tan(-a - phi)) / 2 / np.cos(theta)
            )
            # a strip spans an observation interval, so its center is half a
            # step ahead of its edges
            expected_along = altitude * np.tan(theta) + step / 2
            along, across = _centroid_offsets(trajectory, footprints, 20)
            self.assertAlmostEqual(across, expected_across, delta=5.0)
            self.assertAlmostEqual(along, expected_along, delta=5.0)
            if roll != 0:
                # opposite to the roll: toward the outside of the turn
                self.assertLess(across * roll, 0)

    def test_heading_orients_scan(self):
        """
        Tests that the scan is perpendicular to the heading, which differs from
        the track with wind (a crab angle).
        """
        latitude, longitude, altitude = [33.4], [-111.9], [10000.0]
        for heading in (10.0, 100.0, 275.0):
            lat, lon, valid = compute_airborne_edge_points(
                latitude, longitude, altitude, [0.0], [0.0], [heading], 25.0
            )
            self.assertTrue(valid.all())
            azimuth = GEOD.inv(lon[0, 0], lat[0, 0], lon[0, 1], lat[0, 1])[0]
            difference = (azimuth - (heading + 90) + 180) % 360 - 180
            self.assertAlmostEqual(difference, 0.0, delta=0.1)
        trajectory = _trajectory(track=0.0, heading=10.0)
        crabbed = collect_airborne_ground_track(
            trajectory, 25.0, 1.0, use_attitude=True
        )
        straight = collect_airborne_ground_track(
            trajectory.drop(columns="heading"), 25.0, 1.0, use_attitude=True
        )
        self.assertTrue((crabbed["heading"] == 10.0).all())
        self.assertTrue((straight["heading"] == 0.0).all())
        self.assertFalse(
            crabbed.geometry.iloc[10].equals_exact(straight.geometry.iloc[10], 1e-9)
        )

    def test_attitude_columns(self):
        """
        Tests that the attitude is applied only where requested (zero otherwise).
        """
        trajectory = _trajectory(n=10, roll=12.0, pitch=3.0, heading=5.0)
        for geometry in ("flat", "ellipsoid"):
            level = collect_airborne_ground_track(
                trajectory, 25.0, 1.0, use_attitude=False, geometry=geometry
            )
            banked = collect_airborne_ground_track(
                trajectory, 25.0, 1.0, use_attitude=True, geometry=geometry
            )
            if geometry == "ellipsoid":
                self.assertTrue((level["roll"] == 0.0).all())
                self.assertTrue((level["pitch"] == 0.0).all())
                self.assertEqual(list(level.columns[:7]), STRIP_COLUMNS)
            self.assertTrue((banked["roll"] == 12.0).all())
            self.assertTrue((banked["pitch"] == 3.0).all())
            self.assertEqual(list(banked.columns[:7]), STRIP_COLUMNS)

    def test_sorts_by_time_and_carries_segment_name(self):
        """
        Tests that unsorted samples with naive times are sorted and returned in
        UTC, and that the segment name is carried to the results with attitude.
        """
        trajectory = _trajectory(n=20)
        trajectory["segment_name"] = ["leg"] * 10 + ["turn"] * 10
        trajectory["time"] = trajectory["time"].dt.tz_localize(None)
        shuffled = trajectory.sample(frac=1.0, random_state=1)
        for geometry in ("flat", "ellipsoid"):
            footprints = collect_airborne_ground_track(
                shuffled, 25.0, 1.0, use_attitude=True, geometry=geometry
            )
            self.assertTrue(footprints["time"].is_monotonic_increasing)
            self.assertEqual(str(footprints["time"].dt.tz), "UTC")
            self.assertEqual(
                list(footprints["segment_name"]), ["leg"] * 10 + ["turn"] * 10
            )

    def test_gap_is_extrapolated(self):
        """
        Tests that a strip before a gap in the samples (a gap longer than 1.5
        observation intervals, as between flight legs) spans one interval
        rather than the gap.
        """
        trajectory = _trajectory(n=60, roll=10.0)
        gapped = trajectory.drop(index=range(30, 40)).reset_index(drop=True)
        for geometry in ("flat", "ellipsoid"):
            footprints = collect_airborne_ground_track(
                gapped, 25.0, 1.0, use_attitude=True, geometry=geometry
            )
            areas = _areas(footprints)
            # the sample before the gap is index 29 and has the same area as
            # any other sample in steady flight (not that of a strip across
            # the ten-second gap, which is eleven times larger)
            self.assertAlmostEqual(areas[29] / areas[28], 1.0, delta=1e-3)
            self.assertTrue(footprints["valid_obs"].all())

    def test_surface_altitude_is_not_valid(self):
        """
        Tests that, in every geometry (including the default level rectangle),
        samples at or below the surface are not valid and have an empty
        footprint, while the samples above it are valid with valid shapes.
        """
        trajectory = _trajectory(n=12, roll=5.0)
        trajectory.loc[3:5, "altitude"] = 0.0
        trajectory.loc[8, "altitude"] = -100.0
        below = trajectory["altitude"] <= 0
        for use_attitude, geometry in (
            (False, "flat"),
            (True, "flat"),
            (False, "ellipsoid"),
            (True, "ellipsoid"),
        ):
            with self.subTest(use_attitude=use_attitude, geometry=geometry):
                footprints = collect_airborne_ground_track(
                    trajectory, 25.0, 1.0, use_attitude=use_attitude, geometry=geometry
                )
                self.assertEqual(len(footprints), 12)
                self.assertFalse(footprints.loc[below, "valid_obs"].any())
                self.assertTrue(footprints.loc[below].geometry.is_empty.all())
                # the sample before a sample at the surface has no next edge
                # in the strip geometries, so check the samples away from it
                away = ~below & ~below.shift(-1, fill_value=False)
                self.assertTrue(footprints.loc[away, "valid_obs"].all())
                self.assertTrue(footprints.loc[away].geometry.is_valid.all())
                self.assertFalse(footprints.loc[away].geometry.is_empty.any())

    def test_invalid_samples(self):
        """
        Tests that samples that cannot observe the ground are not valid.
        """
        for geometry in ("flat", "ellipsoid"):
            # a roll of 85 degrees points the left beam above the horizon
            footprints = collect_airborne_ground_track(
                _trajectory(n=10, roll=85.0),
                25.0,
                1.0,
                use_attitude=True,
                geometry=geometry,
            )
            self.assertFalse(footprints["valid_obs"].any())
            if geometry == "ellipsoid":
                self.assertTrue(footprints.geometry.is_empty.all())
            # at or below the surface
            for altitude in (0.0, -100.0):
                footprints = collect_airborne_ground_track(
                    _trajectory(n=10, altitude=altitude),
                    25.0,
                    1.0,
                    use_attitude=True,
                    geometry=geometry,
                )
                self.assertFalse(footprints["valid_obs"].any())
            # a stationary aircraft sweeps no area
            footprints = collect_airborne_ground_track(
                _trajectory(n=10, speed=0.0),
                25.0,
                1.0,
                use_attitude=True,
                geometry=geometry,
            )
            self.assertFalse(footprints["valid_obs"].any())
            self.assertEqual(len(footprints), 10)

    def test_flat_beam_limit(self):
        """
        Tests that the flat-surface formulas are valid up to 80 degrees from
        nadir: a roll of 50 degrees puts the left beam at 75 degrees (valid)
        and a roll of 60 degrees puts it at 85 degrees (not valid).
        """
        for roll, expected in ((50.0, True), (60.0, False)):
            footprints = collect_airborne_ground_track(
                _trajectory(n=10, roll=roll), 25.0, 1.0, use_attitude=True
            )
            self.assertEqual(footprints["valid_obs"].tolist(), [expected] * 10)
            self.assertFalse(footprints.geometry.is_empty.any())

    def test_hyplan_is_not_required(self):
        """
        Tests that the analysis does not use HyPlan (a trajectory from any
        source can be used): it runs with HyPlan blocked from being imported.
        """
        module = importlib.import_module("tatc.analysis.airborne_track")
        self.assertNotIn("hyplan", inspect.getsource(module).lower())
        with mock.patch.dict(sys.modules, {"hyplan": None}):
            with self.assertRaises(ImportError):
                importlib.import_module("hyplan")
            for use_attitude, geometry in (
                (False, "flat"),
                (True, "flat"),
                (True, "ellipsoid"),
            ):
                footprints = collect_airborne_ground_track(
                    _trajectory(n=10, roll=10.0),
                    25.0,
                    1.0,
                    use_attitude=use_attitude,
                    geometry=geometry,
                )
                self.assertEqual(len(footprints), 10)
                self.assertTrue(footprints["valid_obs"].all())


class TestAirborneEdgePoints(unittest.TestCase):
    """
    Unit tests for `compute_airborne_edge_points`.
    """

    def test_matches_spice_intersection(self):
        """
        Tests the edge points against the SPICE intersection of rays with the
        WGS 84 ellipsoid, for rays derived by hand at the equator and the
        prime meridian, where the local north, east, and down axes are the
        Earth-fixed z, y, and -x axes. With the nose north, a beam at scan
        angle s is (-cos(s) cos(pitch), sin(s), cos(s) sin(pitch)).
        """
        altitude, half_angle = 10000.0, 25.0
        for pitch in (0.0, 7.0, -12.0):
            lat, lon, valid = compute_airborne_edge_points(
                [0.0], [0.0], [altitude], [0.0], [pitch], [0.0], half_angle
            )
            self.assertTrue(valid.all())
            for k, scan in enumerate((-half_angle, half_angle)):
                s, theta = np.radians(scan), np.radians(pitch)
                ray = [-np.cos(s) * np.cos(theta), np.sin(s), np.cos(s) * np.sin(theta)]
                point = surfpt(
                    [EARTH_EQUATORIAL_RADIUS + altitude, 0.0, 0.0],
                    ray,
                    EARTH_EQUATORIAL_RADIUS,
                    EARTH_EQUATORIAL_RADIUS,
                    EARTH_POLAR_RADIUS,
                )
                expected_lon, expected_lat, _ = recgeo(
                    point, EARTH_EQUATORIAL_RADIUS, EARTH_FLATTENING
                )
                self.assertAlmostEqual(lon[0, k], np.degrees(expected_lon), delta=1e-9)
                self.assertAlmostEqual(lat[0, k], np.degrees(expected_lat), delta=1e-9)

    def test_roll_rotates_scan(self):
        """
        Tests that a roll rotates the scan angles: the edges at a roll are those
        of a level scan at (minus and plus) the half angle less the roll.
        """
        position = ([33.4], [-111.9], [10000.0])
        heading = [45.0]
        half_angle, roll = 25.0, 20.0
        banked = compute_airborne_edge_points(
            *position, [roll], [0.0], heading, half_angle
        )
        right = compute_airborne_edge_points(
            *position, [0.0], [0.0], heading, half_angle - roll
        )
        left = compute_airborne_edge_points(
            *position, [0.0], [0.0], heading, half_angle + roll
        )
        np.testing.assert_allclose(banked[0][0, 1], right[0][0, 1], atol=1e-10)
        np.testing.assert_allclose(banked[1][0, 1], right[1][0, 1], atol=1e-10)
        np.testing.assert_allclose(banked[0][0, 0], left[0][0, 0], atol=1e-10)
        np.testing.assert_allclose(banked[1][0, 0], left[1][0, 0], atol=1e-10)

    def test_invalid_rays(self):
        """
        Tests that rays that miss the Earth, and rays from at or below the
        surface, are not valid and have no point.
        """
        # a roll of 85 degrees points the left ray above the horizon
        lat, lon, valid = compute_airborne_edge_points(
            [33.4], [-111.9], [10000.0], [85.0], [0.0], [0.0], 25.0
        )
        self.assertEqual(valid.tolist(), [[False, True]])
        self.assertTrue(np.isnan(lat[0, 0]) and np.isnan(lon[0, 0]))
        self.assertTrue(np.isfinite(lat[0, 1]) and np.isfinite(lon[0, 1]))
        for altitude in (0.0, -50.0):
            lat, lon, valid = compute_airborne_edge_points(
                [33.4], [-111.9], [altitude], [0.0], [0.0], [0.0], 25.0
            )
            self.assertFalse(valid.any())
            self.assertTrue(np.isnan(lat).all() and np.isnan(lon).all())

    def test_elevation(self):
        """
        Tests that a higher surface brings the edges closer to the aircraft.
        """
        args = ([33.4], [-111.9], [10000.0], [0.0], [0.0], [0.0], 25.0)
        distances = []
        for elevation in (0.0, 2000.0):
            lat, lon, _ = compute_airborne_edge_points(*args, elevation=elevation)
            distances.append(GEOD.inv(-111.9, 33.4, lon[0, 1], lat[0, 1])[2])
        self.assertGreater(distances[0], distances[1])
        self.assertAlmostEqual(distances[0] / distances[1], 10000 / 8000, places=2)


class TestFlatAndEllipsoidAgreement(unittest.TestCase):
    """
    Compares the flat-surface and ellipsoid scan edges. The flat formulas
    assume a planar ground at the aircraft's altitude below it and ignore the
    Earth's curvature (the surface falls away from the beam, which lengthens
    it), so they agree with the ellipsoid within the curvature effect: about
    0.03 percent of the edge's ground range for a 25 degree beam and 0.22
    percent at 50 degrees from nadir (altitude 19.8 km), growing with the
    off-nadir angle. The tolerance is 0.4 percent of the ground range plus 1
    meter.
    """

    def test_edges_agree_within_curvature(self):
        """
        Tests the edges over a range of altitudes, rolls, and pitches.
        """
        half_angle, heading = 25.0, 60.0
        for altitude in (3000.0, 8000.0, 19800.0):
            for roll, pitch in ((0, 0), (25, 0), (-25, 0), (0, 5), (25, 5), (-15, -5)):
                with self.subTest(altitude=altitude, roll=roll, pitch=pitch):
                    lon = np.array([-111.9])
                    lat = np.array([33.4])
                    forward, left, right, valid = _get_flat_scan_offsets(
                        np.array([altitude]),
                        np.array([float(roll)]),
                        np.array([float(pitch)]),
                        half_angle,
                    )
                    self.assertTrue(valid.all())
                    flat = _get_flat_edge_points(
                        lon, lat, np.array([heading]), forward, left, right
                    )
                    ellipsoid_lat, ellipsoid_lon, found = compute_airborne_edge_points(
                        lat, lon, [altitude], [roll], [pitch], [heading], half_angle
                    )
                    self.assertTrue(found.all())
                    for k in range(2):
                        difference = GEOD.inv(
                            flat[2 * k][0],
                            flat[2 * k + 1][0],
                            ellipsoid_lon[0, k],
                            ellipsoid_lat[0, k],
                        )[2]
                        ground_range = GEOD.inv(
                            lon[0], lat[0], ellipsoid_lon[0, k], ellipsoid_lat[0, k]
                        )[2]
                        self.assertLess(difference, 0.004 * ground_range + 1.0)

    def test_footprints_agree_in_validity_and_size(self):
        """
        Tests that the footprints of the two geometries have the same
        validity and areas that agree within one percent.
        """
        trajectory = _trajectory(n=30, altitude=10000.0, roll=15.0, pitch=2.0)
        flat = collect_airborne_ground_track(
            trajectory, 25.0, 1.0, use_attitude=True, geometry="flat"
        )
        ellipsoid = collect_airborne_ground_track(
            trajectory, 25.0, 1.0, use_attitude=True, geometry="ellipsoid"
        )
        self.assertTrue(flat["valid_obs"].equals(ellipsoid["valid_obs"]))
        np.testing.assert_allclose(_areas(flat), _areas(ellipsoid), rtol=0.01)


class TestObservedElevation(unittest.TestCase):
    """
    Unit tests for the `elevation` argument of `collect_airborne_ground_track`:
    footprints where the scan meets a level above the ground (for example a
    layer of a storm), with the aircraft kept at its altitude.
    """

    MODES = (
        {"geometry": "flat", "use_attitude": False},
        {"geometry": "flat", "use_attitude": True},
        {"geometry": "ellipsoid", "use_attitude": True},
    )

    def test_zero_elevation_keeps_results(self):
        """
        Tests that an elevation of zero gives exactly the results without the
        argument, in every geometry, level and banked.
        """
        for roll in (0.0, 20.0):
            trajectory = _trajectory(n=20, altitude=5800.0, roll=roll, pitch=2.0)
            for mode in self.MODES:
                with self.subTest(roll=roll, **mode):
                    default = collect_airborne_ground_track(trajectory, 25.0, 1.0, **mode)
                    zero = collect_airborne_ground_track(
                        trajectory, 25.0, 1.0, elevation=0.0, **mode
                    )
                    self.assertTrue(
                        default.drop(columns="geometry").equals(
                            zero.drop(columns="geometry")
                        )
                    )
                    self.assertEqual(
                        list(default.geometry.to_wkb()), list(zero.geometry.to_wkb())
                    )

    def test_flat_width_matches_formula(self):
        """
        Tests that over a flat ground a straight, level flight at height h
        observes a level at elevation z across a swath 2 (h - z) tan(half
        angle) wide: each footprint's area is that width times the distance
        flown in an observation interval.
        """
        altitude, half_angle, speed = 5800.0, 25.0, 200.0
        trajectory = _trajectory(n=10, altitude=altitude, speed=speed)
        for elevation in (0.0, 3000.0, 5000.0):
            expected = 2 * (altitude - elevation) * np.tan(np.radians(half_angle))
            for use_attitude in (False, True):
                with self.subTest(elevation=elevation, use_attitude=use_attitude):
                    footprints = collect_airborne_ground_track(
                        trajectory,
                        half_angle,
                        1.0,
                        use_attitude=use_attitude,
                        geometry="flat",
                        elevation=elevation,
                    )
                    self.assertTrue(footprints["valid_obs"].all())
                    np.testing.assert_allclose(
                        _areas(footprints) / speed, expected, rtol=2e-3
                    )

    def test_ellipsoid_width_near_formula(self):
        """
        Tests that the ellipsoid swath at an elevation agrees with the flat
        formula 2 (h - z) tan(half angle) within the Earth's curvature (a few
        hundredths of a percent at these heights) and the map projection.
        """
        altitude, half_angle, speed = 5800.0, 25.0, 200.0
        trajectory = _trajectory(n=10, altitude=altitude, speed=speed)
        for elevation in (0.0, 3000.0, 5000.0):
            with self.subTest(elevation=elevation):
                footprints = collect_airborne_ground_track(
                    trajectory,
                    half_angle,
                    1.0,
                    use_attitude=True,
                    geometry="ellipsoid",
                    elevation=elevation,
                )
                self.assertTrue(footprints["valid_obs"].all())
                expected = 2 * (altitude - elevation) * np.tan(np.radians(half_angle))
                np.testing.assert_allclose(_areas(footprints) / speed, expected, rtol=3e-3)

    def test_level_at_or_above_aircraft_is_not_valid(self):
        """
        Tests that a downward-looking instrument cannot observe a level at or
        above the aircraft: those samples are not valid and have empty
        footprints in every geometry, while a level just below it is observed.
        """
        altitude = 5800.0
        trajectory = _trajectory(n=5, altitude=altitude)
        for mode in self.MODES:
            for elevation in (altitude, altitude + 500.0):
                with self.subTest(elevation=elevation, **mode):
                    footprints = collect_airborne_ground_track(
                        trajectory, 25.0, 1.0, elevation=elevation, **mode
                    )
                    self.assertFalse(footprints["valid_obs"].any())
                    self.assertTrue(footprints.geometry.is_empty.all())
            with self.subTest(elevation=altitude - 100.0, **mode):
                footprints = collect_airborne_ground_track(
                    trajectory, 25.0, 1.0, elevation=altitude - 100.0, **mode
                )
                self.assertTrue(footprints["valid_obs"].all())

    def test_elevation_must_be_finite(self):
        """
        Tests that a non-finite elevation is rejected.
        """
        trajectory = _trajectory(n=3)
        for elevation in (np.nan, np.inf):
            with self.assertRaisesRegex(ValueError, "elevation"):
                collect_airborne_ground_track(trajectory, 25.0, 1.0, elevation=elevation)

    def test_raised_surface_height(self):
        """
        Tests the approximation documented in `compute_airborne_edge_points`:
        the surface at an elevation is the ellipsoid with both semi-axes
        extended by it, which a vertical ray meets within 1 cm of that
        geodetic height up to 5 km (exactly at the equator and the poles).
        """
        for latitude in (0.0, 30.0, 45.0, 60.0, 90.0):
            for elevation in (3000.0, 5000.0):
                with self.subTest(latitude=latitude, elevation=elevation):
                    origin = geodetic_to_rectangular(
                        np.array([0.0]), np.array([latitude]), np.array([elevation + 15000.0])
                    )
                    down = _get_ned_to_earth_fixed(np.array([0.0]), np.array([latitude]))[
                        :, :, 2
                    ].T
                    point, found = compute_ellipsoid_intersection(
                        origin, down, elevation=elevation
                    )
                    self.assertTrue(found.all())
                    _, _, height = rectangular_to_geodetic(point)
                    self.assertLess(abs(height[0] - elevation), 0.01)


if __name__ == "__main__":
    unittest.main()

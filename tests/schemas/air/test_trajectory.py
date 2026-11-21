"""
Unit tests for the trajectory schema.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import timedelta

import numpy as np
import pandas as pd
from pydantic import ValidationError
from shapely.geometry import LineString
from skyfield.api import wgs84
from skyfield.framelib import itrs

from tatc.schemas import Trajectory
from tatc.utils.time import _to_time

from .common import GEOD, START, straight_trajectory, surface_distance


def _geodetic(track):
    """
    Gets the longitudes, latitudes (degrees), and altitudes (meters) of a track.
    """
    position = wgs84.geographic_position_of(track)
    return (
        position.longitude.degrees,
        position.latitude.degrees,
        position.elevation.m,
    )


def _ground_velocity(track):
    """
    Gets the azimuths (degrees) and speeds (meters per second) of the
    Earth-fixed velocities of a track.
    """
    position, velocity = track.frame_xyz_and_velocity(itrs)
    lon, lat = np.radians(_geodetic(track)[:2])
    v = np.array(velocity.m_per_s)
    east = -np.sin(lon) * v[0] + np.cos(lon) * v[1]
    north = (
        -np.sin(lat) * np.cos(lon) * v[0]
        - np.sin(lat) * np.sin(lon) * v[1]
        + np.cos(lat) * v[2]
    )
    del position
    return np.mod(np.degrees(np.arctan2(east, north)), 360), np.hypot(east, north)


class TestTrajectory(unittest.TestCase):
    """
    Unit tests for the trajectory schema.
    """

    def test_from_dataframe(self):
        """
        Tests that a data frame's rows are sorted by time, times without a
        time zone are taken as UTC, and other columns are ignored.
        """
        data = straight_trajectory(n=5, roll=2.0)
        data["segment_name"] = "Line 1"
        shuffled = data.iloc[[3, 0, 4, 1, 2]].copy()
        shuffled["time"] = shuffled["time"].dt.tz_localize(None)
        trajectory = Trajectory.from_dataframe(shuffled)
        self.assertEqual(trajectory.time, list(data["time"].dt.to_pydatetime()))
        self.assertEqual(trajectory.latitude, list(data["latitude"]))
        self.assertEqual(trajectory.roll, [2.0] * 5)
        self.assertIsNone(trajectory.pitch)
        self.assertTrue(trajectory.has_attitude())
        self.assertFalse(
            Trajectory.from_dataframe(data.drop(columns="roll")).has_attitude()
        )

    def test_from_dataframe_missing_columns(self):
        """
        Tests that a data frame without a required column is not accepted.
        """
        with self.assertRaisesRegex(ValueError, "altitude"):
            Trajectory.from_dataframe(straight_trajectory(n=3).drop(columns="altitude"))

    def test_from_dataframe_missing_values(self):
        """
        Tests that rows missing a value are omitted and that an optional column
        without any value is ignored.
        """
        data = straight_trajectory(n=5)
        data.loc[2, "altitude"] = np.nan
        data["roll"] = np.nan
        trajectory = Trajectory.from_dataframe(data)
        self.assertEqual(len(trajectory.time), 4)
        self.assertNotIn(data["time"][2], trajectory.time)
        self.assertIsNone(trajectory.roll)

    def test_from_dataframe_repeated_times(self):
        """
        Tests that a data frame with a repeated time (as in a recording) is not
        accepted.
        """
        data = straight_trajectory(n=3)
        with self.assertRaisesRegex(ValidationError, "strictly increasing"):
            Trajectory.from_dataframe(pd.concat([data, data.iloc[[1]]]))

    def test_non_finite_values(self):
        """
        Tests that values that are not finite are not accepted.
        """
        values = Trajectory.from_dataframe(
            straight_trajectory(n=3, roll=0.0)
        ).model_dump()
        for name in ("latitude", "altitude", "track", "roll"):
            for value in (np.nan, np.inf):
                with self.subTest(name=name, value=value):
                    with self.assertRaisesRegex(ValidationError, "finite"):
                        Trajectory(**{**values, name: [0.0, value, 0.0]})

    def test_validation(self):
        """
        Tests that the states must have every given quantity, increasing times,
        and latitudes in range.
        """
        values = Trajectory.from_dataframe(straight_trajectory(n=3)).model_dump()
        for name, value in (
            ("altitude", [1.0, 2.0]),
            ("roll", [0.0]),
            ("time", [START, START, START + timedelta(seconds=1)]),
            ("latitude", [91.0, 0.0, 0.0]),
            ("ground_speed", [-1.0, 0.0, 0.0]),
            ("max_gap", timedelta(0)),
        ):
            with self.subTest(name=name), self.assertRaises(ValidationError):
                Trajectory(**{**values, name: value})

    def test_serialization(self):
        """
        Tests that a trajectory is serialized to JSON and back.
        """
        trajectory = Trajectory.from_dataframe(straight_trajectory(n=4, heading=31.0))
        self.assertEqual(
            Trajectory.model_validate_json(trajectory.model_dump_json()), trajectory
        )

    def test_track_at_states(self):
        """
        Tests that the track at the times of the states has their positions and
        velocities (along the track at the ground speed).
        """
        data = straight_trajectory(n=10, altitude=8000.0, speed=150.0, track=300.0)
        trajectory = Trajectory.from_dataframe(data)
        track = trajectory.get_orbit_track_at_time(_to_time(trajectory.time))
        lon, lat, alt = _geodetic(track)
        np.testing.assert_allclose(lon, data["longitude"], atol=1e-8)
        np.testing.assert_allclose(lat, data["latitude"], atol=1e-8)
        np.testing.assert_allclose(alt, data["altitude"], atol=1e-3)
        azimuth, speed = _ground_velocity(track)
        np.testing.assert_allclose(azimuth, data["track"], atol=1e-6)
        np.testing.assert_allclose(speed, 150.0, rtol=1e-6)

    def test_track_between_states(self):
        """
        Tests that the track between the states of a flight along a geodesic
        stays on the geodesic at the distance flown, with the altitude
        interpolated linearly.
        """
        data = straight_trajectory(n=5, altitude=5000.0, speed=200.0, step=10.0)
        seconds = np.array([2.5, 13.0, 27.5, 39.0])
        times = _to_time([START + timedelta(seconds=s) for s in seconds])
        lon, lat, _ = _geodetic(
            Trajectory.from_dataframe(data).get_orbit_track_at_time(times)
        )
        azimuth, distance = GEOD.inv(
            np.full(4, data["longitude"][0]), np.full(4, data["latitude"][0]), lon, lat
        )[::2]
        np.testing.assert_allclose(
            distance,
            surface_distance(200.0 * seconds, 35.0, data["track"][0], 5000.0),
            atol=0.01,
        )
        np.testing.assert_allclose(azimuth, data["track"][0], atol=1e-3)
        data["altitude"] = [5000.0, 5100.0, 5200.0, 5300.0, 5400.0]
        _, _, alt = _geodetic(
            Trajectory.from_dataframe(data).get_orbit_track_at_time(times)
        )
        np.testing.assert_allclose(alt, 5000.0 + 10.0 * seconds, atol=1e-3)

    def test_track_beyond_states(self):
        """
        Tests that the track continues along the track at the ground speed
        before the first state, after the last state, and near a gap between
        states.
        """
        data = straight_trajectory(n=4, speed=100.0, track=90.0)
        # a gap of 10 minutes after the second state, before a leg to the south
        later = straight_trajectory(n=2, speed=100.0, track=180.0, latitude=34.0)
        later["time"] += timedelta(minutes=10)
        trajectory = Trajectory.from_dataframe(
            pd.concat([data.iloc[:2], later], ignore_index=True)
        )
        for state, offset in ((0, -2.0), (1, 0.5), (3, 3.0)):
            with self.subTest(state=state, offset=offset):
                time = trajectory.time[state] + timedelta(seconds=offset)
                track = trajectory.get_orbit_track_at_time(_to_time(time))
                lon, lat, _ = _geodetic(track)
                azimuth, _, distance = GEOD.inv(
                    trajectory.longitude[state], trajectory.latitude[state], lon, lat
                )
                # within the gap, the curve to the next leg bends the path by
                # about 0.2 meters
                self.assertAlmostEqual(
                    distance,
                    surface_distance(
                        100.0 * abs(offset),
                        trajectory.latitude[state],
                        trajectory.track[state],
                        10000.0,
                    ),
                    delta=0.2,
                )
                expected = trajectory.track[state] + (180 if offset < 0 else 0)
                self.assertAlmostEqual(
                    np.mod(azimuth - expected + 180, 360) - 180, 0, delta=0.5
                )

    def test_track_of_one_state(self):
        """
        Tests that a trajectory of one state continues along its track.
        """
        trajectory = Trajectory.from_dataframe(straight_trajectory(n=1, speed=50.0))
        track = trajectory.get_orbit_track_at_time(
            _to_time([START, START + timedelta(seconds=4)])
        )
        lon, lat, _ = _geodetic(track)
        self.assertAlmostEqual(
            GEOD.inv(lon[0], lat[0], lon[1], lat[1])[2],
            surface_distance(200.0, 35.0, 30.0, 10000.0),
            delta=0.01,
        )

    def test_stationary_track(self):
        """
        Tests that the velocity of a stationary aircraft still points along its
        track (so that it orients the views of its instruments).
        """
        data = straight_trajectory(n=3, speed=0.0, track=45.0)
        track = Trajectory.from_dataframe(data).get_orbit_track_at_time(
            _to_time(list(data["time"]) + [START + timedelta(seconds=0.5)])
        )
        azimuth, speed = _ground_velocity(track)
        np.testing.assert_allclose(azimuth, 45.0, atol=1e-6)
        self.assertTrue(np.all(speed > 0))

    def test_angles(self):
        """
        Tests that the attitude angles are interpolated linearly, with the
        heading relative to the track continuous across north.
        """
        data = straight_trajectory(n=3, track=0.0, step=10.0)
        data["track"] = [1.0, 1.0, 1.0]
        data["heading"] = [359.0, 3.0, 5.0]
        data["roll"] = [0.0, 10.0, 10.0]
        trajectory = Trajectory.from_dataframe(data)
        crab, roll, pitch = trajectory.get_angles(
            _to_time([START + timedelta(seconds=s) for s in (-5.0, 5.0, 15.0, 30.0)])
        )
        np.testing.assert_allclose(crab, [-2.0, 0.0, 3.0, 4.0], atol=1e-6)
        np.testing.assert_allclose(roll, [0.0, 5.0, 10.0, 10.0], atol=1e-6)
        np.testing.assert_allclose(pitch, 0.0)

    def test_contains(self):
        """
        Tests that a trajectory contains the times from its first to its last
        state, except within a gap longer than the maximum gap.
        """
        data = straight_trajectory(n=4, step=30.0)
        later = straight_trajectory(n=2, latitude=34.0)
        later["time"] += timedelta(minutes=10)
        trajectory = Trajectory.from_dataframe(pd.concat([data, later]))
        offsets = [-1.0, 0.0, 45.0, 90.0, 91.0, 300.0, 599.0, 600.0, 601.0, 602.0]
        times = [START + timedelta(seconds=s) for s in offsets]
        np.testing.assert_array_equal(
            trajectory.contains(times),
            [False, True, True, True, False, False, False, True, True, False],
        )
        self.assertTrue(trajectory.contains(START))
        np.testing.assert_array_equal(
            Trajectory.from_dataframe(pd.concat([data, later]), max_gap=None).contains(
                times
            ),
            [False] + [True] * 8 + [False],
        )

    def test_estimated_motion(self):
        """
        Tests that the track and ground speed, if not given, are estimated from
        the positions, within each part of the trajectory between gaps, and
        that a state without a neighbor within the maximum gap is not accepted.
        """
        data = straight_trajectory(n=10, speed=150.0, track=300.0, step=5.0)
        given = Trajectory.from_dataframe(data)
        estimated = Trajectory.from_dataframe(
            data.drop(columns=["track", "ground_speed"])
        )
        times = _to_time([START + timedelta(seconds=s) for s in (0.0, 12.5, 45.0)])
        np.testing.assert_allclose(
            _ground_velocity(estimated.get_orbit_track_at_time(times)),
            _ground_velocity(given.get_orbit_track_at_time(times)),
            rtol=1e-4,
        )
        later = straight_trajectory(n=2, speed=100.0, track=180.0, latitude=34.0)
        later["time"] += timedelta(minutes=10)
        parts = Trajectory.from_dataframe(
            pd.concat([data, later]).drop(columns=["track", "ground_speed"])
        )
        azimuth, speed = _ground_velocity(
            parts.get_orbit_track_at_time(_to_time(parts.time))
        )
        np.testing.assert_allclose(
            azimuth, list(data["track"]) + [180.0] * 2, atol=1e-3
        )
        np.testing.assert_allclose(speed, [150.0] * 10 + [100.0] * 2, rtol=1e-4)
        with self.assertRaisesRegex(ValidationError, "neighboring state"):
            Trajectory.from_dataframe(
                pd.concat([data, later.iloc[:1]]).drop(columns="track")
            )

    def test_estimated_motion_across_antimeridian_and_pole(self):
        """
        Tests that the track and ground speed estimated from the positions are
        those given across the anti-meridian and over the pole, where the
        longitude jumps.
        """
        for track, latitude, longitude in ((90.0, 10.0, 179.995), (0.0, 89.995, 0.0)):
            with self.subTest(latitude=latitude):
                data = straight_trajectory(
                    n=10, track=track, latitude=latitude, longitude=longitude
                )
                self.assertGreater(np.ptp(data["longitude"]), 90)
                given = Trajectory.from_dataframe(data)
                estimated = Trajectory.from_dataframe(
                    data.drop(columns=["track", "ground_speed"])
                )
                times = _to_time(
                    [START + timedelta(seconds=s) for s in (0.0, 4.5, 9.0)]
                )
                given_track = given.get_orbit_track_at_time(times)
                estimated_track = estimated.get_orbit_track_at_time(times)
                np.testing.assert_allclose(
                    estimated_track.velocity.m_per_s,
                    given_track.velocity.m_per_s,
                    atol=0.05,
                )
                np.testing.assert_allclose(
                    estimated_track.position.m, given_track.position.m, atol=0.01
                )

    def test_assignment(self):
        """
        Tests that assigned fields are validated and change the track.
        """
        trajectory = Trajectory.from_dataframe(straight_trajectory(n=3))
        before = trajectory.get_angles(_to_time(START))
        trajectory.roll = [5.0, 5.0, 5.0]
        np.testing.assert_allclose(trajectory.get_angles(_to_time(START))[1], 5.0)
        self.assertEqual(before[1], 0.0)
        with self.assertRaises(ValidationError):
            trajectory.pitch = [1.0]

    def test_from_waypoints(self):
        """
        Tests that a flight between two waypoints matches the states of an
        aircraft flying along the geodesic at its ground speed and altitude,
        ending at the last waypoint.
        """
        data = straight_trajectory(n=61, track=30.0)
        trajectory = Trajectory.from_waypoints(
            [
                (data.longitude.iloc[0], data.latitude.iloc[0]),
                (data.longitude.iloc[-1], data.latitude.iloc[-1]),
            ],
            START.to_pydatetime(),
            speed=200.0,
            altitude=10000.0,
        )
        self.assertEqual(len(trajectory.time), 61)
        self.assertIsNone(trajectory.max_gap)
        self.assertFalse(trajectory.has_attitude())
        np.testing.assert_allclose(trajectory.latitude, data.latitude, atol=1e-6)
        np.testing.assert_allclose(trajectory.longitude, data.longitude, atol=1e-6)
        np.testing.assert_allclose(trajectory.track, data.track, atol=1e-6)
        np.testing.assert_allclose(trajectory.ground_speed, 200.0)
        self.assertAlmostEqual(
            (trajectory.time[-1] - START).total_seconds(), 60.0, places=3
        )

    def test_from_waypoints_legs(self):
        """
        Tests that a flight through several waypoints (a line, with a
        repeated waypoint) visits each, at the altitude of each, with a final
        state at the last waypoint between time steps.
        """
        corner = GEOD.fwd(-110.0, 35.0, 90.0, 30e3)[:2]
        end = GEOD.fwd(*corner, 0.0, 20e3)[:2]
        trajectory = Trajectory.from_waypoints(
            LineString([(-110.0, 35.0), corner, corner, end]),
            START.to_pydatetime(),
            speed=150.0,
            altitude=[1000.0, 3000.0, 3000.0, 2000.0],
            time_step=timedelta(seconds=7),
        )
        latitude, longitude = np.array(trajectory.latitude), np.array(
            trajectory.longitude
        )
        _, _, to_corner = GEOD.inv(
            longitude,
            latitude,
            np.full(len(latitude), corner[0]),
            np.full(len(latitude), corner[1]),
        )
        self.assertLess(np.min(to_corner), 150.0 * 7)
        self.assertAlmostEqual(trajectory.longitude[-1], end[0], places=6)
        self.assertAlmostEqual(trajectory.latitude[-1], end[1], places=6)
        self.assertAlmostEqual(trajectory.altitude[0], 1000.0)
        self.assertAlmostEqual(trajectory.altitude[-1], 2000.0)
        self.assertLessEqual(max(trajectory.altitude), 3000.0)
        steps = np.diff([time.timestamp() for time in trajectory.time])
        np.testing.assert_allclose(steps[:-1], 7.0)
        self.assertLess(steps[-1], 7.0)
        # the flight takes about the distance over the speed
        duration = (trajectory.time[-1] - trajectory.time[0]).total_seconds()
        self.assertAlmostEqual(duration, 50e3 / 150.0, delta=1.0)

    def test_from_waypoints_invalid(self):
        """
        Tests that waypoints without two distinct points, a naive start, or
        a speed or time step that is not positive raise a ValueError.
        """
        start = START.to_pydatetime()
        waypoints = [(-110.0, 35.0), (-109.0, 35.0)]
        cases = {
            "one point": lambda: Trajectory.from_waypoints(
                [(-110.0, 35.0), (-110.0, 35.0)], start, 100.0, 1000.0
            ),
            "naive start": lambda: Trajectory.from_waypoints(
                waypoints, start.replace(tzinfo=None), 100.0, 1000.0
            ),
            "speed": lambda: Trajectory.from_waypoints(waypoints, start, 0.0, 1000.0),
            "time step": lambda: Trajectory.from_waypoints(
                waypoints, start, 100.0, 1000.0, timedelta(0)
            ),
        }
        for name, call in cases.items():
            with self.subTest(name):
                with self.assertRaises(ValueError):
                    call()


if __name__ == "__main__":
    unittest.main()

"""
Unit tests for the tatc.preprocess.airborne module.

The tests of flight plans and sensor models use HyPlan, if installed.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import importlib.util
import unittest
from datetime import timedelta
from types import SimpleNamespace

import geopandas as gpd
import numpy as np
import pandas as pd
from pyproj import Geod
from shapely.geometry import LineString, Point

from tatc.analysis import collect_airborne_ground_track
from tatc.preprocess import (
    estimate_attitude,
    instrument_from_sensor,
    sample_flight_plan,
    states_from_iwg1,
)
from tatc.schemas import Aircraft, Trajectory
from tatc.utils import compute_geoid_height

GEOD = Geod(ellps="WGS84")
TAKEOFF = pd.Timestamp("2026-10-15T17:00:00Z")
HAS_HYPLAN = importlib.util.find_spec("hyplan") is not None


class _Quantity:
    """
    A stand-in for a quantity with units (as of pint), in the units in which it
    is read.
    """

    def __init__(self, value):
        self.value = value

    def __truediv__(self, other):
        return _Quantity(self.value / other.value)

    def __rmul__(self, factor):
        return _Quantity(factor * self.value)

    def m_as(self, unit):  # pylint: disable=unused-argument
        """
        Gets the magnitude.
        """
        return self.value


def _line(longitude, latitude, azimuth, distance, n=50):
    """
    Gets a line string along a geodesic.
    """
    lon, lat, _ = GEOD.fwd(
        np.full(n, longitude),
        np.full(n, latitude),
        np.full(n, azimuth),
        np.linspace(0, distance, n),
    )
    return LineString(zip(lon, lat))


def _plan():
    """
    Gets a flight plan of two segments: a 60-second line north at 10,000 feet
    (with a crab angle) and a 30-second climb east to 11,000 feet.
    """
    first = _line(-112.0, 33.0, 0.0, 12000.0)
    end = first.coords[-1]
    return gpd.GeoDataFrame(
        {
            "segment_type": ["flight_line", "climb"],
            "segment_name": ["Line 1", "Climb"],
            "time_to_segment": [1.0, 0.5],
            "start_altitude": [10000.0, 10000.0],
            "end_altitude": [10000.0, 11000.0],
            "crab_angle_deg": [-5.0, np.nan],
        },
        geometry=[first, _line(end[0], end[1], 90.0, 6000.0)],
        crs="EPSG:4326",
    )


class TestSampleFlightPlan(unittest.TestCase):
    """
    Unit tests for `sample_flight_plan`.
    """

    def test_states(self):
        """
        Tests the states sampled along the segments of a plan.
        """
        states = sample_flight_plan(_plan(), TAKEOFF, geoid_height=-30.0)
        self.assertEqual(len(states), 91)
        self.assertEqual(states["time"].iloc[0], TAKEOFF)
        self.assertEqual(states["time"].iloc[-1], TAKEOFF + timedelta(seconds=90))
        line = states["segment_name"] == "Line 1"
        self.assertEqual(line.sum(), 60)
        np.testing.assert_allclose(states.loc[line, "ground_speed"], 200.0, rtol=1e-6)
        np.testing.assert_allclose(states.loc[~line, "ground_speed"], 200.0, rtol=1e-6)
        np.testing.assert_allclose(
            np.mod(states.loc[line, "track"] + 180, 360) - 180, 0.0, atol=1e-6
        )
        np.testing.assert_allclose(states.loc[line, "heading"], 355.0, atol=1e-6)
        np.testing.assert_allclose(
            states.loc[~line, "heading"], states.loc[~line, "track"]
        )
        np.testing.assert_allclose(
            states.loc[line, "altitude"], 10000 * 0.3048 - 30.0, atol=1e-9
        )
        self.assertAlmostEqual(states["altitude"].iloc[-1], 11000 * 0.3048 - 30.0)
        # by default, the EGM2008 geoid at each state
        egm2008 = sample_flight_plan(_plan(), TAKEOFF)
        np.testing.assert_allclose(
            egm2008["altitude"],
            states["altitude"]
            + 30.0
            + compute_geoid_height(egm2008["longitude"], egm2008["latitude"]),
        )
        # the states lie along the path at the distance flown
        distance = GEOD.inv(
            -112.0, 33.0, states["longitude"].iloc[30], states["latitude"].iloc[30]
        )[2]
        self.assertAlmostEqual(distance, 6000.0, delta=0.01)
        Trajectory.from_dataframe(states)

    def test_time_step(self):
        """
        Tests the time step between states and that it must be positive.
        """
        states = sample_flight_plan(_plan(), TAKEOFF, time_step=timedelta(seconds=10))
        self.assertEqual(len(states), 10)
        with self.assertRaises(ValueError):
            sample_flight_plan(_plan(), TAKEOFF, time_step=timedelta(0))

    def test_skipped_segments(self):
        """
        Tests that a segment without altitudes is skipped while its time
        elapses, that a segment without duration is skipped, and that an empty
        plan has no states.
        """
        line = _line(-112.0, 33.0, 0.0, 12000.0)
        plan = gpd.GeoDataFrame(
            {
                "time_to_segment": [1.0, 2.0, 0.0, 1.0],
                "start_altitude": [10000.0, np.nan, 10000.0, 10000.0],
                "end_altitude": [10000.0, np.nan, 10000.0, 10000.0],
            },
            geometry=[line] * 4,
            crs="EPSG:4326",
        )
        states = sample_flight_plan(plan, TAKEOFF, geoid_height=0.0)
        seconds = (states["time"] - TAKEOFF).dt.total_seconds()
        np.testing.assert_array_equal(seconds, list(range(60)) + list(range(180, 241)))
        # without the optional columns, the heading is the track
        self.assertTrue(states["segment_name"].isna().all())
        np.testing.assert_allclose(states["heading"], states["track"])
        empty = sample_flight_plan(plan.iloc[:0], TAKEOFF)
        self.assertTrue(empty.empty)
        self.assertEqual(
            list(empty.columns),
            [
                "time",
                "latitude",
                "longitude",
                "altitude",
                "track",
                "ground_speed",
                "heading",
                "segment_type",
                "segment_name",
            ],
        )

    def test_hold(self):
        """
        Tests that a path without length (as a hold) is stationary.
        """
        plan = gpd.GeoDataFrame(
            {"time_to_segment": [0.1], "start_altitude": [1e4], "end_altitude": [1e4]},
            geometry=[LineString([(-112.0, 33.0), (-112.0, 33.0)])],
            crs="EPSG:4326",
        )
        states = sample_flight_plan(plan, TAKEOFF)
        self.assertEqual(len(states), 7)
        np.testing.assert_allclose(states["ground_speed"], 0.0)
        np.testing.assert_allclose(
            states[["longitude", "latitude"]], [[-112.0, 33.0]] * 7
        )
        Trajectory.from_dataframe(states)

    @unittest.skipUnless(HAS_HYPLAN, "requires hyplan")
    def test_hyplan(self):
        """
        Tests the states of a flight plan computed by HyPlan.
        """
        import hyplan  # pylint: disable=import-outside-toplevel

        lines = [
            hyplan.FlightLine.center_length_azimuth(
                lat=33.4,
                lon=longitude,
                length=20 * hyplan.ureg.km,
                az=azimuth,
                altitude_msl=5000 * hyplan.ureg.meter,
                site_name=name,
            )
            for longitude, azimuth, name in ((-111.9, 0.0, "A"), (-111.7, 180.0, "B"))
        ]
        plan = hyplan.compute_flight_plan(
            hyplan.NASA_P3(),
            lines,
            takeoff_time=TAKEOFF.to_pydatetime(),
            wind_speed=20 * hyplan.ureg.knot,
            wind_direction=270.0,
        )
        states = sample_flight_plan(plan, TAKEOFF)
        duration = plan["time_to_segment"].sum() * 60
        self.assertAlmostEqual(len(states), duration, delta=1.0)
        self.assertEqual(set(states["segment_name"]), set(plan["segment_name"]))
        for name, track in (("A", 0.0), ("B", 180.0)):
            line = states["segment_name"] == name
            np.testing.assert_allclose(
                np.mod(states.loc[line, "track"] - track + 180, 360) - 180, 0, atol=0.2
            )
            crab = plan.loc[plan["segment_name"] == name, "crab_angle_deg"].iloc[0]
            self.assertNotEqual(crab, 0.0)
            np.testing.assert_allclose(
                np.mod(
                    states.loc[line, "heading"]
                    - states.loc[line, "track"]
                    - crab
                    + 180,
                    360,
                )
                - 180,
                0.0,
                atol=1e-6,
            )


class TestEstimateAttitude(unittest.TestCase):
    """
    Unit tests for `estimate_attitude`.
    """

    def _turn(self, n=121, speed=150.0, rate=3.0, climb=5.0):
        """
        Gets the states of a turn at a constant rate (degrees per second),
        climbing at a constant rate (meters per second).
        """
        t = np.arange(n, dtype=float)
        return pd.DataFrame(
            {
                "time": TAKEOFF + pd.to_timedelta(t, unit="s"),
                "altitude": 5000.0 + climb * t,
                "track": np.mod(350.0 + rate * t, 360.0),
                "ground_speed": speed,
            }
        )

    def test_turn(self):
        """
        Tests the roll and pitch of a steady, climbing turn.
        """
        states = self._turn()
        result = estimate_attitude(states)
        expected_roll = np.degrees(np.arctan(150.0 * np.radians(3.0) / 9.80665))
        np.testing.assert_allclose(result["roll"], expected_roll, rtol=1e-9)
        np.testing.assert_allclose(
            result["pitch"], np.degrees(np.arctan2(5.0, 150.0)), rtol=1e-9
        )
        pd.testing.assert_frame_equal(result[states.columns], states)
        left = estimate_attitude(self._turn(rate=-3.0))
        np.testing.assert_allclose(left["roll"], -expected_roll, rtol=1e-9)

    def test_smoothing_and_gaps(self):
        """
        Tests that a change of rate is smoothed over the smoothing duration,
        and that rates are not computed across a gap.
        """
        straight = self._turn(n=60, rate=0.0, climb=0.0)
        turn = self._turn(n=60, rate=3.0, climb=0.0)
        turn["time"] += timedelta(seconds=60)
        states = pd.concat([straight, turn], ignore_index=True)
        smoothed = estimate_attitude(states)["roll"]
        unsmoothed = estimate_attitude(states, smoothing=None)["roll"]
        self.assertTrue(smoothed.iloc[0] == 0 and smoothed.iloc[-1] > 30)
        self.assertGreater(smoothed.iloc[55], 0)
        self.assertEqual(unsmoothed.iloc[55], 0)
        turn["time"] += timedelta(minutes=10)
        gapped = estimate_attitude(pd.concat([straight, turn], ignore_index=True))
        np.testing.assert_allclose(gapped["roll"].iloc[:60], 0.0, atol=1e-12)
        with self.assertRaisesRegex(ValueError, "ground_speed"):
            estimate_attitude(states.drop(columns="ground_speed"))

    def test_invalid_states(self):
        """
        Tests that states missing a value, or not in increasing order of time,
        are not accepted.
        """
        states = self._turn(n=10)
        for name in ("time", "altitude", "track", "ground_speed"):
            with self.subTest(name=name):
                missing = states.copy()
                missing.loc[4, name] = None
                with self.assertRaisesRegex(ValueError, f"missing values.*{name}"):
                    estimate_attitude(missing)
        for order in ([1, 0, 2], [0, 1, 1, 2]):
            with self.subTest(order=order):
                with self.assertRaisesRegex(ValueError, "strictly increasing"):
                    estimate_attitude(states.iloc[order])

    def test_isolated_states(self):
        """
        Tests that states without a neighbor within the maximum gap are level.
        """
        states = self._turn(n=5)
        states["time"] = TAKEOFF + pd.to_timedelta(120 * np.arange(5), unit="s")
        result = estimate_attitude(states)
        np.testing.assert_allclose(result[["roll", "pitch"]], 0.0)
        self.assertTrue((estimate_attitude(states, max_gap=None)["roll"] > 0).all())

    def test_flight_path_pitch(self):
        """
        Tests that the pitch is that of the flight path, toward vertical for a
        climb without ground speed.
        """
        for speed, pitch in ((0.0, 90.0), (5.0, 45.0)):
            with self.subTest(speed=speed):
                result = estimate_attitude(self._turn(n=10, speed=speed, rate=0.0))
                np.testing.assert_allclose(result["pitch"], pitch, atol=1e-4)
                np.testing.assert_allclose(result["roll"], 0.0)


class TestStatesFromIwg1(unittest.TestCase):
    """
    Unit tests for `states_from_iwg1`.
    """

    def test_states(self):
        """
        Tests the units and columns of the states of a recording, and that
        rows missing a value are omitted.
        """
        data = pd.DataFrame(
            {
                "timestamp": pd.date_range("2026-10-15 17:00", periods=4, freq="s"),
                "latitude": [33.0, 33.001, np.nan, 33.003],
                "longitude": -112.0,
                "altitude": 19000.0,
                "altitude_gps_ft": 20000.0,
                "groundspeed": 300.0,
                "track": 1.0,
                "true_heading": 2.0,
                "roll_deg": 3.0,
                "pitch_deg": np.nan,
            }
        )
        states = states_from_iwg1(data, geoid_height=-30.0)
        self.assertEqual(
            list(states.columns),
            [
                "time",
                "latitude",
                "longitude",
                "altitude",
                "track",
                "ground_speed",
                "heading",
                "roll",
            ],
        )
        self.assertEqual(len(states), 3)
        self.assertEqual(str(states["time"].dt.tz), "UTC")
        np.testing.assert_allclose(states["altitude"], 20000 * 0.3048 - 30.0)
        np.testing.assert_allclose(states["ground_speed"], 300 * 1852 / 3600)
        trajectory = Trajectory.from_dataframe(states)
        self.assertEqual(trajectory.roll, [3.0] * 3)
        np.testing.assert_allclose(
            states_from_iwg1(data)["altitude"],
            20000 * 0.3048 + compute_geoid_height(-112.0, states["latitude"]),
        )
        with self.assertRaisesRegex(ValueError, "altitude_gps_ft"):
            states_from_iwg1(data.drop(columns="altitude_gps_ft"))


class TestInstrumentFromSensor(unittest.TestCase):
    """
    Unit tests for `instrument_from_sensor`.
    """

    def test_scanner(self):
        """
        Tests that a scanning sensor tilted to the right is a swept scan
        rolled to the right (a negative roll angle).
        """
        sensor = SimpleNamespace(
            name="scanner", ifov=0.1, swath_offset_angles=lambda: (-20.0, 30.0)
        )
        instrument = instrument_from_sensor(sensor)
        self.assertEqual(instrument.name, "scanner")
        self.assertEqual(instrument.cross_track_field_of_view, 50.0)
        self.assertEqual(instrument.along_track_field_of_view, 0.1)
        self.assertEqual(instrument.roll_angle, -5.0)
        self.assertTrue(instrument.is_swept)
        self.assertEqual(instrument.view_geometry, "scan")
        self.assertEqual(
            instrument_from_sensor(
                sensor, along_track_field_of_view=1.0
            ).along_track_field_of_view,
            1.0,
        )
        with self.assertRaisesRegex(ValueError, "along_track_field_of_view"):
            instrument_from_sensor(
                SimpleNamespace(name="x", swath_offset_angles=lambda: (-1.0, 1.0))
            )
        with self.assertRaisesRegex(TypeError, "not supported"):
            instrument_from_sensor(SimpleNamespace(name="profiler"))

    def test_along_track_field_of_view(self):
        """
        Tests the along-track field of view of a scanning lidar (its beam
        divergence) and radar (its beam width, the wavelength over the
        antenna's length or twice the azimuth resolution).
        """
        width = np.degrees(0.24 / 1.6)
        for attributes, expected in (
            ({"beam_divergence": _Quantity(0.05)}, 0.05),
            ({"lens": SimpleNamespace(divergence_mrad=1.0)}, np.degrees(0.001)),
            ({"wavelength": _Quantity(0.24), "antenna_length": _Quantity(1.6)}, width),
            (
                {"wavelength": _Quantity(0.24), "azimuth_resolution": _Quantity(0.8)},
                width,
            ),
        ):
            with self.subTest(attributes=sorted(attributes)):
                sensor = SimpleNamespace(
                    name="sensor",
                    swath_offset_angles=lambda: (-10.0, 10.0),
                    **attributes,
                )
                self.assertAlmostEqual(
                    instrument_from_sensor(sensor).along_track_field_of_view, expected
                )

    def test_frame_camera(self):
        """
        Tests that a frame camera tilted toward a direction (0 degrees forward,
        90 degrees to the left) views the ground at the tilt angle in that
        direction, and that its footprints are collected.
        """
        height, tilt = 1000.0, 20.0
        data = pd.DataFrame(
            {
                "time": TAKEOFF + pd.to_timedelta([0, 1], unit="s"),
                "latitude": [33.0, 33.0018],
                "longitude": -112.0,
                "altitude": height,
                "track": 0.0,
                "ground_speed": 200.0,
            }
        )
        for direction in (None, 0.0, 45.0, 90.0, 200.0):
            with self.subTest(direction=direction):
                camera = SimpleNamespace(
                    name="camera", fov_x=10.0, fov_y=8.0, tilt_angle=tilt
                )
                if direction is not None:
                    camera.tilt_direction = direction
                instrument = instrument_from_sensor(camera)
                self.assertTrue(instrument.is_rectangular)
                self.assertEqual(
                    (
                        instrument.cross_track_field_of_view,
                        instrument.along_track_field_of_view,
                    ),
                    (10.0, 8.0),
                )
                aircraft = Aircraft(
                    name="aircraft",
                    trajectory=Trajectory.from_dataframe(data),
                    instruments=[instrument],
                )
                center = aircraft.get_instrument().compute_footprint_center(
                    aircraft.get_track(), 0
                )
                azimuth, _, distance = GEOD.inv(
                    -112.0,
                    33.0,
                    center.longitude.degrees[0],
                    center.latitude.degrees[0],
                )
                # toward the tilt direction (counterclockwise from the track)
                angle = np.radians(0.0 if direction is None else direction)
                offset = height * np.tan(np.radians(tilt))
                left = -distance * np.sin(np.radians(azimuth))
                forward = distance * np.cos(np.radians(azimuth))
                self.assertAlmostEqual(left, offset * np.sin(angle), delta=0.5)
                self.assertAlmostEqual(forward, offset * np.cos(angle), delta=0.5)
                footprints = collect_airborne_ground_track(aircraft)
                self.assertTrue(footprints["valid_obs"].all())
                self.assertTrue(
                    footprints.geometry.iloc[0].contains(
                        Point(center.longitude.degrees[0], center.latitude.degrees[0])
                    )
                )

    @unittest.skipUnless(HAS_HYPLAN, "requires hyplan")
    def test_hyplan_sensors(self):
        """
        Tests the instruments of HyPlan's scanning sensors: their swaths span
        the cross-track angles of their edges.
        """
        import hyplan  # pylint: disable=import-outside-toplevel

        for name in ("HyTES", "AVIRIS3", "LVIS", "UAVSAR_Lband"):
            with self.subTest(name=name):
                sensor = hyplan.create_sensor(name)
                instrument = instrument_from_sensor(sensor)
                port, starboard = sensor.swath_offset_angles()
                # the left edge (roll + half field of view) is at -port
                self.assertAlmostEqual(
                    instrument.roll_angle + instrument.cross_track_field_of_view / 2,
                    -port,
                )
                self.assertAlmostEqual(
                    instrument.roll_angle - instrument.cross_track_field_of_view / 2,
                    -starboard,
                )
                self.assertGreater(instrument.along_track_field_of_view, 0)
        # a left-looking radar looks to the left
        self.assertGreater(
            instrument_from_sensor(hyplan.create_sensor("UAVSAR_Lband")).roll_angle, 0
        )

    @unittest.skipUnless(HAS_HYPLAN, "requires hyplan")
    def test_hyplan_frame_camera(self):
        """
        Tests that the view center of a tilted frame camera is where HyPlan
        projects it (the centroid of the footprint of a camera with a narrow
        field of view).
        """
        import hyplan  # pylint: disable=import-outside-toplevel

        ureg = hyplan.ureg
        times = TAKEOFF + pd.to_timedelta([0, 1], unit="s")
        data = pd.DataFrame(
            {
                "time": times,
                "latitude": [33.0, 33.0018],
                "longitude": -112.0,
                "altitude": 1000.0,
                "track": 0.0,
                "ground_speed": 200.0,
            }
        )
        for direction in (0.0, 45.0, 90.0, 200.0):
            with self.subTest(direction=direction):
                camera = hyplan.FrameCamera(
                    "camera",
                    sensor_width=1 * ureg.mm,
                    sensor_height=1 * ureg.mm,
                    focal_length=500 * ureg.mm,
                    resolution_x=6000,
                    resolution_y=4000,
                    frame_rate=1 * ureg.Hz,
                    f_speed=4.0,
                    tilt_angle=20.0,
                    tilt_direction=direction,
                )
                expected = camera.ground_footprint(1000 * ureg.m).centroid
                aircraft = Aircraft(
                    name="aircraft",
                    trajectory=Trajectory.from_dataframe(data),
                    instruments=[instrument_from_sensor(camera)],
                )
                center = aircraft.get_instrument().compute_footprint_center(
                    aircraft.get_track(), 0
                )
                azimuth, _, distance = GEOD.inv(
                    -112.0,
                    33.0,
                    center.longitude.degrees[0],
                    center.latitude.degrees[0],
                )
                right = distance * np.sin(np.radians(azimuth))
                forward = distance * np.cos(np.radians(azimuth))
                # within the curvature of the Earth over a flat projection
                self.assertAlmostEqual(right, expected.x, delta=1.0)
                self.assertAlmostEqual(forward, expected.y, delta=1.0)
                footprints = collect_airborne_ground_track(aircraft)
                self.assertTrue(footprints["valid_obs"].all())


if __name__ == "__main__":
    unittest.main()

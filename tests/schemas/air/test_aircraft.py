"""
Unit tests for the aircraft schema.

The expected attitudes come from the aerospace rotation matrices (heading, then
pitch, then roll) and the SPICE intersection of a ray with the WGS 84
ellipsoid.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import timedelta

import numpy as np
import pandas as pd
from pydantic import ValidationError
from spiceypy.spiceypy import recgeo, surfpt

from tatc.constants import (
    EARTH_EQUATORIAL_RADIUS,
    EARTH_FLATTENING,
    EARTH_POLAR_RADIUS,
)
from tatc.schemas import (
    Aircraft,
    ConicalInstrument,
    Instrument,
    PointedInstrument,
    Trajectory,
)
from tatc.schemas.air.aircraft import _to_body_angles
from tatc.utils.ellipsoid import geodetic_to_rectangular
from tatc.utils.projection import _rotate_view_axes

from .common import GEOD, straight_trajectory


def _body_to_ned(heading, pitch, roll):
    """
    Gets the rotation matrix from an aircraft's body frame (x forward, y right
    wing, z down) to the local north-east-down frame.
    """
    y, p, r = np.radians([heading, pitch, roll])
    yaw = np.array([[np.cos(y), -np.sin(y), 0], [np.sin(y), np.cos(y), 0], [0, 0, 1]])
    pit = np.array([[np.cos(p), 0, np.sin(p)], [0, 1, 0], [-np.sin(p), 0, np.cos(p)]])
    rol = np.array([[1, 0, 0], [0, np.cos(r), -np.sin(r)], [0, np.sin(r), np.cos(r)]])
    return yaw @ pit @ rol


def _ned_to_earth_fixed(longitude, latitude):
    """
    Gets the rotation matrix from the local north-east-down frame to the
    Earth-fixed frame.
    """
    lon, lat = np.radians([longitude, latitude])
    north = [-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)]
    east = [-np.sin(lon), np.cos(lon), 0]
    down = [-np.cos(lat) * np.cos(lon), -np.cos(lat) * np.sin(lon), -np.sin(lat)]
    return np.array([north, east, down]).T


class TestAircraft(unittest.TestCase):
    """
    Unit tests for the aircraft schema.
    """

    def test_body_angles(self):
        """
        Tests that the body angles orient the axes of the aircraft's body.
        """
        rng = np.random.default_rng(0)
        crab, roll, pitch = rng.uniform(-60, 60, (3, 50))
        yaw, body_roll, body_pitch = _to_body_angles(crab, roll, pitch)
        for i in range(50):
            # the axes in a frame of the track: forward, right, and down
            matrix = _body_to_ned(crab[i], pitch[i], roll[i])
            forward, left = np.array([1, 0, 0]), np.array([0, -1, 0])
            down = np.array([0, 0, 1])
            # the yawed nadir frame of TAT-C (yaw turns forward to the left)
            psi = np.radians(yaw[i])
            along = np.cos(psi) * forward + np.sin(psi) * left
            cross = np.cos(psi) * left - np.sin(psi) * forward
            boresight, along, cross = _rotate_view_axes(
                down, along, cross, body_roll[i], body_pitch[i]
            )
            np.testing.assert_allclose(boresight, matrix[:, 2], atol=1e-12)
            np.testing.assert_allclose(along, matrix[:, 0], atol=1e-12)
            np.testing.assert_allclose(cross, -matrix[:, 1], atol=1e-12)

    def test_views_match_spice_intersection(self):
        """
        Tests that the view of an instrument rolled from the body of an
        aircraft in any attitude meets the ellipsoid where the ray rotated by
        the aircraft's attitude does.
        """
        rng = np.random.default_rng(1)
        for i in range(20):
            latitude, longitude = rng.uniform(-70, 70), rng.uniform(-180, 180)
            heading, roll, pitch = rng.uniform(-20, 20), *rng.uniform(-40, 40, 2)
            mounting = rng.uniform(-30, 30)
            data = straight_trajectory(
                n=2,
                altitude=rng.uniform(1000, 15000),
                track=rng.uniform(0, 360),
                latitude=latitude,
                longitude=longitude,
                roll=roll,
                pitch=pitch,
            )
            data["heading"] = data["track"] + heading
            aircraft = Aircraft(
                name="aircraft",
                trajectory=Trajectory.from_dataframe(data),
                instruments=[
                    PointedInstrument(
                        cross_track_field_of_view=1,
                        along_track_field_of_view=1,
                        roll_angle=mounting,
                    )
                ],
            )
            center = aircraft.get_instrument().compute_footprint_center(
                aircraft.get_track(aircraft.trajectory.time[0])
            )
            # a positive roll angle of the instrument looks to the left wing
            ray = _body_to_ned(data["heading"][0], pitch, roll) @ np.array(
                [0, -np.sin(np.radians(mounting)), np.cos(np.radians(mounting))]
            )
            position = geodetic_to_rectangular(
                data["longitude"][0], data["latitude"][0], data["altitude"][0]
            )
            point = surfpt(
                position,
                _ned_to_earth_fixed(data["longitude"][0], data["latitude"][0]) @ ray,
                EARTH_EQUATORIAL_RADIUS,
                EARTH_EQUATORIAL_RADIUS,
                EARTH_POLAR_RADIUS,
            )
            lon, lat, _ = recgeo(point, EARTH_EQUATORIAL_RADIUS, EARTH_FLATTENING)
            with self.subTest(i=i):
                self.assertAlmostEqual(
                    center.latitude.degrees, np.degrees(lat), places=7
                )
                self.assertAlmostEqual(
                    center.longitude.degrees, np.degrees(lon), places=7
                )

    def test_get_instrument(self):
        """
        Tests that instruments are bound to the aircraft's attitude (level if
        its trajectory gives none, as for an unbound instrument) and observe at
        the times that the trajectory contains.
        """
        instrument = PointedInstrument(
            cross_track_field_of_view=40, along_track_field_of_view=1
        )
        data = straight_trajectory(n=3)
        later = straight_trajectory(n=3, latitude=36.0)
        later["time"] += timedelta(minutes=10)
        level = Aircraft(
            name="a",
            trajectory=Trajectory.from_dataframe(pd.concat([data, later])),
            instruments=[instrument],
        )
        bound = level.get_instrument()
        self.assertIsNotNone(bound.platform_attitude)
        self.assertIsNone(level.instruments[0].platform_attitude)
        track = level.get_track(level.trajectory.time)
        np.testing.assert_allclose(
            bound.compute_footprint_center(track).latitude.degrees,
            instrument.compute_footprint_center(track).latitude.degrees,
            atol=1e-12,
        )
        offsets = [-1.0, 0.0, 2.0, 3.0, 300.0, 600.0, 602.0, 603.0]
        times = [data["time"][0] + timedelta(seconds=s) for s in offsets]
        np.testing.assert_array_equal(
            bound.is_observing(level.get_track(times)),
            [False, True, True, False, False, True, True, False],
        )

    def test_check_instruments(self):
        """
        Tests that instruments with angle profiles or an inertial velocity
        frame are not accepted.
        """
        trajectory = Trajectory.from_dataframe(straight_trajectory(n=3))
        for kwargs in (
            {"roll_angle_profile": [(0, 10)]},
            {"pitch_angle_profile": [(0, 10)]},
            {"velocity_frame": "inertial"},
        ):
            with self.subTest(**{k: str(v) for k, v in kwargs.items()}):
                with self.assertRaises(ValidationError):
                    Aircraft(
                        name="a",
                        trajectory=trajectory,
                        instruments=[
                            PointedInstrument(
                                cross_track_field_of_view=10,
                                along_track_field_of_view=10,
                                **kwargs,
                            )
                        ],
                    )
        with self.assertRaises(ValidationError):
            Aircraft(
                name="a",
                trajectory=trajectory,
                instruments=[
                    ConicalInstrument(
                        cone_angle=45,
                        along_track_field_of_view=2,
                        velocity_frame="inertial",
                    )
                ],
            )

    def test_serialization(self):
        """
        Tests that an aircraft is serialized to JSON and back.
        """
        aircraft = Aircraft(
            name="a",
            trajectory=Trajectory.from_dataframe(straight_trajectory(n=3, pitch=2.0)),
            instruments=[
                Instrument(name="cone", field_of_regard=30),
                PointedInstrument(
                    name="scanner",
                    cross_track_field_of_view=40,
                    along_track_field_of_view=1,
                    view_geometry="scan",
                    is_swept=True,
                ),
            ],
        )
        self.assertEqual(
            Aircraft.model_validate_json(aircraft.model_dump_json()), aircraft
        )

    def test_roll_and_pitch_conventions(self):
        """
        Tests that the roll and pitch of an aircraft (right wing down and nose
        up positive) shift a nadir view as the roll and pitch angles of an
        instrument do: to the left and forward.
        """
        data = straight_trajectory(n=2, track=0.0, altitude=5000.0)
        for name, angles in (("roll", (10.0, 0.0)), ("pitch", (0.0, 10.0))):
            with self.subTest(name=name):
                rotated = Aircraft(
                    name="rotated",
                    trajectory=Trajectory.from_dataframe(data.assign(**{name: 10.0})),
                    instruments=[
                        PointedInstrument(
                            cross_track_field_of_view=1, along_track_field_of_view=1
                        )
                    ],
                )
                pointed = Aircraft(
                    name="pointed",
                    trajectory=Trajectory.from_dataframe(data),
                    instruments=[
                        PointedInstrument(
                            cross_track_field_of_view=1,
                            along_track_field_of_view=1,
                            roll_angle=angles[0],
                            pitch_angle=angles[1],
                        )
                    ],
                )
                centers = [
                    aircraft.get_instrument().compute_footprint_center(
                        aircraft.get_track(), 0
                    )
                    for aircraft in (rotated, pointed)
                ]
                np.testing.assert_allclose(
                    centers[0].longitude.degrees,
                    centers[1].longitude.degrees,
                    atol=1e-9,
                )
                np.testing.assert_allclose(
                    centers[0].latitude.degrees, centers[1].latitude.degrees, atol=1e-9
                )
                azimuth = GEOD.inv(
                    data["longitude"][0],
                    data["latitude"][0],
                    centers[0].longitude.degrees[0],
                    centers[0].latitude.degrees[0],
                )[0]
                # to the left (west) of a northward track, or forward (north)
                self.assertAlmostEqual(
                    azimuth, -90.0 if name == "roll" else 0.0, delta=0.01
                )


if __name__ == "__main__":
    unittest.main()

"""
Unit tests for the TrainConstellation schema.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""
import unittest
from datetime import timedelta

from pydantic import ValidationError

from tatc import constants
from tatc.schemas import CircularOrbit, Instrument, TrainConstellation


class TestTrainConstellation(unittest.TestCase):
    """
    Unit tests for the TrainConstellation schema.
    """
    def setUp(self):
        self.test_data = {
            "name": "Test Constellation",
            "orbit": {
                "mean_altitude": "400000",
                "inclination": 51.6,
                "right_ascension_ascending_node": 180,
                "true_anomaly": 180,
            },
            "instruments": [{"name": "Test Instrument", "field_of_regard": 25.0}],
            "number_satellites": 4,
            "interval": timedelta(minutes=10),
        }
        self.test_orbit = CircularOrbit(**self.test_data.get("orbit"))
        self.test_con_rgt = TrainConstellation(**self.test_data, repeat_ground_track=True)
        self.test_con_nrgt = TrainConstellation(**self.test_data, repeat_ground_track=False)

    def test_good_data(self):
        """
        Test that the TrainConstellation object can be created from valid data.
        """
        self.assertEqual(self.test_con_rgt.name, self.test_data.get("name"))
        self.assertEqual(
            self.test_con_rgt.orbit, self.test_orbit
        )
        self.assertEqual(len(self.test_con_rgt.instruments), 1)
        self.assertEqual(
            self.test_con_rgt.instruments[0],
            Instrument(**self.test_data.get("instruments")[0]),
        )
        self.assertEqual(
            self.test_con_rgt.number_satellites, self.test_data.get("number_satellites")
        )
        self.assertEqual(self.test_con_rgt.interval, self.test_data.get("interval"))
        self.assertEqual(
            self.test_con_rgt.repeat_ground_track,
            True,
        )

    def test_get_delta_mean_anomaly_repeat_ground_track_tle(self):
        """
        Test that the delta mean anomaly can be retrieved from the TrainConstellation 
        object for a repeat ground track TLE orbit.
        """
        self.assertAlmostEqual(
            self.test_con_rgt.get_delta_mean_anomaly(),
            -360 * self.test_con_rgt.interval / self.test_orbit.get_orbit_period(),
            delta=0.001,
        )

    def test_get_delta_mean_anomaly_no_repeat_ground_track_tle(self):
        """
        Test that the delta mean anomaly can be retrieved from the TrainConstellation
        object for a non-repeat ground track TLE orbit.
        """
        self.assertAlmostEqual(
            self.test_con_nrgt.get_delta_mean_anomaly(),
            -360 * self.test_con_nrgt.interval / self.test_orbit.get_orbit_period(),
            delta=0.001,
        )

    def test_get_delta_raan_repeat_ground_track_tle(self):
        """
        Test that the delta right ascension of ascending node can be retrieved from the
        TrainConstellation object for a repeat ground track TLE orbit.
        """
        self.assertAlmostEqual(
            self.test_con_rgt.get_delta_raan(),
            360
            * self.test_con_rgt.interval.total_seconds()
            / constants.EARTH_SIDEREAL_DAY_S,
            delta=1e-9,
        )

    def test_get_delta_raan_repeat_ground_track_uses_sidereal_day(self):
        """
        Test that the repeat ground track RAAN correction is based on the
        sidereal day (Earth's rotation relative to inertial space), not
        the 24-hour solar day: the two differ by about 0.27%, which would
        otherwise leave each trailing satellite's ascending node at a
        slightly different Earth-fixed longitude than the one ahead of it
        rather than truly repeating the ground track. This is a regression
        test for a bug where `timedelta(days=1)` (the solar day) was used
        instead of `constants.EARTH_SIDEREAL_DAY_S`.
        """
        con = TrainConstellation(
            **{**self.test_data, "interval": timedelta(hours=1)},
            repeat_ground_track=True,
        )
        solar_day_raan = 360 * (con.interval / timedelta(days=1))
        self.assertNotAlmostEqual(con.get_delta_raan(), solar_day_raan, delta=1e-4)

    def test_get_delta_raan_no_repeat_ground_track_tle(self):
        """
        Test that the delta right ascension of ascending node can be retrieved from the
        TrainConstellation object for a non-repeat ground track TLE orbit.
        """
        self.assertEqual(self.test_con_nrgt.get_delta_raan(), 0.0)

    def helper_test_generate_members(self, constellation):
        """
        Helper function to test that the members of a TrainConstellation object can be
        generated correctly.
        """
        members = constellation.generate_members()
        self.assertEqual(len(members), constellation.number_satellites)
        for i in range(len(members) - 1):
            self.assertAlmostEqual(
                (
                    members[i + 1].orbit.get_mean_anomaly()
                    - members[i].orbit.get_mean_anomaly()
                )
                % 360,
                constellation.get_delta_mean_anomaly() % 360,
                delta=0.001,
            )
            self.assertAlmostEqual(
                (
                    (
                        members[i + 1].orbit.get_right_ascension_ascending_node()
                        - members[i].orbit.get_right_ascension_ascending_node()
                    )
                    % 360
                    if constellation.orbit.type == "tle"
                    else (
                        members[i + 1].orbit.right_ascension_ascending_node
                        - members[i].orbit.right_ascension_ascending_node
                    )
                    % 360
                ),
                constellation.get_delta_raan() % 360,
                delta=0.001,
            )

    def test_generate_members_repeat_ground_track_tle(self):
        """
        Test that the members of a TrainConstellation object can be generated correctly
        for a repeat ground track TLE orbit.
        """
        self.helper_test_generate_members(self.test_con_rgt)

    def test_generate_members_non_repeat_ground_track_tle(self):
        """
        Test that the members of a TrainConstellation object can be generated correctly
        for a non-repeat ground track TLE orbit.
        """
        self.helper_test_generate_members(self.test_con_nrgt)

    def test_generate_members_first_satellite_matches_lead_orbit(self):
        """
        Test that the first generated member (index 0, zero mean anomaly
        and RAAN offset) matches the constellation's own lead orbit.
        """
        members = self.test_con_rgt.generate_members()
        self.assertEqual(members[0].orbit, self.test_orbit)

    def test_generate_members_instruments_independent_per_satellite(self):
        """
        Test that each generated member gets its own deep-copied
        instruments list, not one shared (aliased) across members or with
        the constellation itself: mutating one member's instrument must
        not affect another member's or the constellation's.
        """
        members = self.test_con_rgt.generate_members()
        self.assertIsNot(members[0].instruments, members[1].instruments)
        self.assertIsNot(members[0].instruments, self.test_con_rgt.instruments)
        members[0].instruments[0].name = "Renamed"
        self.assertEqual(members[1].instruments[0].name, "Test Instrument")
        self.assertEqual(self.test_con_rgt.instruments[0].name, "Test Instrument")

    def test_type_defaults_to_train(self):
        """
        Test that omitting type defaults to the "train" discriminator.
        """
        self.assertEqual(self.test_con_rgt.type, "train")

    def test_type_rejects_other_values(self):
        """
        Test that type only accepts the "train" literal, rejecting other
        space system type discriminators.
        """
        with self.assertRaises(ValidationError):
            TrainConstellation(**self.test_data, type="walker")

    def test_number_satellites_default(self):
        """
        Test that omitting number_satellites defaults to a single satellite.
        """
        con = TrainConstellation(
            name="Test Constellation",
            orbit=self.test_data["orbit"],
            interval=timedelta(minutes=10),
        )
        self.assertEqual(con.number_satellites, 1)

    def test_number_satellites_must_be_at_least_one(self):
        """
        Test that number_satellites must be at least 1.
        """
        with self.assertRaises(ValidationError):
            TrainConstellation(**{**self.test_data, "number_satellites": 0})

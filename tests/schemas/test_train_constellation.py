"""
Unit tests for the TrainConstellation schema.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""
import unittest
from datetime import timedelta

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
        self.assertEqual(
            self.test_con_rgt.get_delta_raan(),
            360 * self.test_con_rgt.interval / timedelta(days=1),
        )

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

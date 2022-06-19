import unittest

from datetime import timedelta
import numpy as np

from tatc.schemas import WalkerConstellation, TwoLineElements, Instrument


class TestWalkerConstellation(unittest.TestCase):
    def setUp(self):
        self.d420_data = {
            "name": "Test Constellation",
            "configuration": "delta",
            "orbit": {
                "tle": [
                    "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
                    "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
                ]
            },
            "instruments": [{"name": "Test Instrument", "field_of_regard": 25.0}],
            "number_satellites": 4,
            "number_planes": 2,
            "relative_spacing": 0,
        }
        self.d420_con = WalkerConstellation(**self.d420_data)
        self.s420_data = {
            "name": "Test Constellation",
            "configuration": "star",
            "orbit": {
                "tle": [
                    "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
                    "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
                ]
            },
            "instruments": [{"name": "Test Instrument", "field_of_regard": 25.0}],
            "number_satellites": 4,
            "number_planes": 2,
            "relative_spacing": 0,
        }
        self.s420_con = WalkerConstellation(**self.s420_data)
        self.d520_data = {
            "name": "Test Constellation",
            "configuration": "delta",
            "orbit": {
                "altitude": "400000",
                "inclination": 51.6,
                "right_ascension_ascending_node": 180,
                "true_anomaly": 180,
            },
            "instruments": [{"name": "Test Instrument", "field_of_regard": 25.0}],
            "number_satellites": 5,
            "number_planes": 2,
            "relative_spacing": 1,
        }
        self.d520_con = WalkerConstellation(**self.d520_data)

    def normalize_angle(self, angle):
        return np.mod(360 + angle, 360)

    def test_constructor(self):
        self.assertEqual(self.d420_con.name, self.d420_data.get("name"))
        self.assertEqual(
            self.d420_con.orbit, TwoLineElements(**self.d420_data.get("orbit"))
        )
        self.assertEqual(len(self.d420_con.instruments), 1)
        self.assertEqual(
            self.d420_con.instruments[0],
            Instrument(**self.d420_data.get("instruments")[0]),
        )
        self.assertEqual(
            self.d420_con.configuration, self.d420_data.get("configuration")
        )
        self.assertEqual(
            self.d420_con.number_satellites, self.d420_data.get("number_satellites")
        )
        self.assertEqual(
            self.d420_con.number_planes, self.d420_data.get("number_planes")
        )
        self.assertEqual(
            self.d420_con.relative_spacing, self.d420_data.get("relative_spacing")
        )

    def test_get_satellites_per_plane(self):
        self.assertEqual(
            self.d420_con.get_satellites_per_plane(),
            np.ceil(self.d420_con.number_satellites / self.d420_con.number_planes),
        )
        self.assertEqual(
            self.d520_con.get_satellites_per_plane(),
            np.ceil(self.d520_con.number_satellites / self.d520_con.number_planes),
        )
        self.assertEqual(
            self.s420_con.get_satellites_per_plane(),
            np.ceil(self.s420_con.number_satellites / self.s420_con.number_planes),
        )

    def test_get_delta_mean_anomaly_within_planes(self):
        self.assertEqual(
            self.d420_con.get_delta_mean_anomaly_within_planes(),
            360 / self.d420_con.get_satellites_per_plane(),
        )
        self.assertEqual(
            self.d520_con.get_delta_mean_anomaly_within_planes(),
            360 / self.d520_con.get_satellites_per_plane(),
        )
        self.assertEqual(
            self.s420_con.get_delta_mean_anomaly_within_planes(),
            360 / self.s420_con.get_satellites_per_plane(),
        )

    def test_get_delta_mean_anomaly_between_planes(self):
        self.assertEqual(
            self.d420_con.get_delta_mean_anomaly_between_planes(),
            self.d420_con.relative_spacing * 360 / self.d420_con.number_satellites,
        )
        self.assertEqual(
            self.d520_con.get_delta_mean_anomaly_between_planes(),
            self.d520_con.relative_spacing * 360 / self.d520_con.number_satellites,
        )
        self.assertEqual(
            self.s420_con.get_delta_mean_anomaly_between_planes(),
            self.s420_con.relative_spacing * 360 / self.s420_con.number_satellites,
        )

    def test_get_delta_raan_between_planes_delta(self):
        self.assertEqual(
            self.d420_con.get_delta_raan_between_planes(),
            360 / self.d420_con.number_planes,
        )
        self.assertEqual(
            self.d520_con.get_delta_raan_between_planes(),
            360 / self.d520_con.number_planes,
        )

    def test_get_delta_raan_between_planes_star(self):
        self.assertEqual(
            self.s420_con.get_delta_raan_between_planes(),
            180 / self.s420_con.number_planes,
        )

    def helper_test_generate_members(self, constellation):
        members = constellation.generate_members()
        self.assertEqual(len(members), constellation.number_satellites)
        for i in range(len(members) - 1):
            sat_in_plane = np.mod(i, constellation.get_satellites_per_plane())
            next_sat_in_plane = np.mod(i + 1, constellation.get_satellites_per_plane())
            plane = i // constellation.get_satellites_per_plane()
            next_plane = (i + 1) // constellation.get_satellites_per_plane()
            if plane == next_plane:
                self.assertAlmostEqual(
                    self.normalize_angle(
                        members[i + 1].orbit.get_mean_anomaly()
                        - members[i].orbit.get_mean_anomaly()
                    ),
                    self.normalize_angle(
                        constellation.get_delta_mean_anomaly_within_planes()
                    ),
                    delta=0.001,
                )
                self.assertAlmostEqual(
                    self.normalize_angle(
                        members[i + 1].orbit.get_right_ascension_ascending_node()
                        - members[i].orbit.get_right_ascension_ascending_node()
                        if constellation.orbit.type == "tle"
                        else members[i + 1].orbit.right_ascension_ascending_node
                        - members[i].orbit.right_ascension_ascending_node
                    ),
                    0.0,
                    delta=0.001,
                )
            else:
                self.assertAlmostEqual(
                    self.normalize_angle(
                        members[i + 1].orbit.get_mean_anomaly()
                        - members[i].orbit.get_mean_anomaly()
                    ),
                    self.normalize_angle(
                        constellation.get_delta_mean_anomaly_between_planes()
                        + constellation.get_delta_mean_anomaly_within_planes()
                        * (next_sat_in_plane - sat_in_plane)
                    ),
                    delta=0.001,
                )
                self.assertAlmostEqual(
                    self.normalize_angle(
                        members[i + 1].orbit.get_right_ascension_ascending_node()
                        - members[i].orbit.get_right_ascension_ascending_node()
                        if constellation.orbit.type == "tle"
                        else members[i + 1].orbit.right_ascension_ascending_node
                        - members[i].orbit.right_ascension_ascending_node
                    ),
                    self.normalize_angle(constellation.get_delta_raan_between_planes()),
                    delta=0.001,
                )

    def test_generate_members_delta_tle(self):
        self.helper_test_generate_members(self.d420_con)

    def test_generate_members_delta_circular(self):
        self.helper_test_generate_members(self.d520_con)

    def test_generate_members_star_tle(self):
        self.helper_test_generate_members(self.s420_con)

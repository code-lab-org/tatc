import unittest

from datetime import timedelta

from tatc.schemas import TrainConstellation, TwoLineElements, Instrument


class TestTrainConstellation(unittest.TestCase):
    def setUp(self):
        self.test_data_1 = {
            "name": "Test Constellation",
            "orbit": {
                "tle": [
                    "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
                    "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
                ]
            },
            "instruments": [{"name": "Test Instrument", "field_of_regard": 25.0}],
            "number_satellites": 4,
            "interval": timedelta(minutes=10),
            "repeat_ground_track": True,
        }
        self.test_con_1 = TrainConstellation(**self.test_data_1)
        self.test_data_2 = {
            "name": "Test Constellation",
            "orbit": {
                "tle": [
                    "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
                    "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
                ]
            },
            "instruments": [{"name": "Test Instrument", "field_of_regard": 25.0}],
            "number_satellites": 4,
            "interval": timedelta(minutes=10),
            "repeat_ground_track": False,
        }
        self.test_con_2 = TrainConstellation(**self.test_data_2)

    def test_good_data(self):
        self.assertEqual(self.test_con_1.name, self.test_data_1.get("name"))
        self.assertEqual(
            self.test_con_1.orbit, TwoLineElements(**self.test_data_1.get("orbit"))
        )
        self.assertEqual(len(self.test_con_1.instruments), 1)
        self.assertEqual(
            self.test_con_1.instruments[0],
            Instrument(**self.test_data_1.get("instruments")[0]),
        )
        self.assertEqual(
            self.test_con_1.number_satellites, self.test_data_1.get("number_satellites")
        )
        self.assertEqual(self.test_con_1.interval, self.test_data_1.get("interval"))
        self.assertEqual(
            self.test_con_1.repeat_ground_track,
            self.test_data_1.get("repeat_ground_track"),
        )

    def test_get_delta_mean_anomaly(self):
        self.assertAlmostEqual(
            self.test_con_1.get_delta_mean_anomaly(),
            360 * self.test_con_1.interval / self.test_con_1.orbit.get_orbit_period(),
            delta=0.001,
        )
        self.assertAlmostEqual(
            self.test_con_2.get_delta_mean_anomaly(),
            360 * self.test_con_2.interval / self.test_con_2.orbit.get_orbit_period(),
            delta=0.001,
        )

    def test_get_delta_raan_repeat_ground_track(self):
        self.assertEqual(
            self.test_con_1.get_delta_raan(),
            -1 * 360 * self.test_con_1.interval / timedelta(days=1),
        )

    def test_get_delta_raan_no_repeat_ground_track(self):
        self.assertEqual(self.test_con_2.get_delta_raan(), 0.0)

    def test_generate_members_repeat_ground_track(self):
        members = self.test_con_1.generate_members()
        self.assertEqual(len(members), 4)
        for i in range(3):
            self.assertAlmostEqual(
                members[i + 1].orbit.get_mean_anomaly()
                - members[i].orbit.get_mean_anomaly(),
                self.test_con_1.get_delta_mean_anomaly(),
                delta=0.001,
            )
            self.assertAlmostEqual(
                members[i + 1].orbit.get_right_ascension_ascending_node()
                - members[i].orbit.get_right_ascension_ascending_node(),
                self.test_con_1.get_delta_raan(),
                delta=0.001,
            )

    def test_generate_members_non_repeat_ground_track(self):
        members = self.test_con_2.generate_members()
        self.assertEqual(len(members), 4)
        for i in range(3):
            self.assertAlmostEqual(
                members[i + 1].orbit.get_mean_anomaly()
                - members[i].orbit.get_mean_anomaly(),
                self.test_con_1.get_delta_mean_anomaly(),
                delta=0.001,
            )
            self.assertEqual(
                members[i + 1].orbit.get_right_ascension_ascending_node()
                - members[i].orbit.get_right_ascension_ascending_node(),
                0.0,
            )

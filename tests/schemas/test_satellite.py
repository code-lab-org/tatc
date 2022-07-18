import unittest

from tatc.schemas import Satellite, TwoLineElements, Instrument


class TestSatellite(unittest.TestCase):
    def setUp(self):
        self.test_data = {
            "name": "Test Satellite",
            "orbit": {
                "tle": [
                    "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
                    "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
                ]
            },
            "instruments": [{"name": "Test Instrument", "field_of_regard": 25.0}],
        }
        self.test_sat = Satellite(**self.test_data)

    def test_good_data(self):
        self.assertEqual(self.test_sat.name, self.test_data.get("name"))
        self.assertEqual(
            self.test_sat.orbit, TwoLineElements(**self.test_data.get("orbit"))
        )
        self.assertEqual(len(self.test_sat.instruments), 1)
        self.assertEqual(
            self.test_sat.instruments[0],
            Instrument(**self.test_data.get("instruments")[0]),
        )

    def test_generate_members(self):
        members = self.test_sat.generate_members()
        self.assertEqual(len(members), 1)
        self.assertEqual(members[0], self.test_sat)

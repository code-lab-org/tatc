import unittest

from tatc.schemas.orbit import TwoLineElements
from pydantic import ValidationError


class TestTLE(unittest.TestCase):
    def test_good_data(self):
        good_data = {
            "tle": [
                "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
                "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
            ]
        }
        orbit = TwoLineElements(**good_data)
        self.assertEqual(orbit.tle[0], good_data.get("tle")[0])
        self.assertEqual(orbit.tle[1], good_data.get("tle")[1])

    def test_bad_data_invalid(self):
        bad_data = {"tle": ["not valid", "not valid"]}
        with self.assertRaises(ValidationError):
            TwoLineElements(**bad_data)

    def test_bad_data_checksums(self):
        bad_data = {
            "tle": [
                "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9994",
                "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
            ]
        }
        with self.assertRaises(ValidationError):
            TwoLineElements(**bad_data)

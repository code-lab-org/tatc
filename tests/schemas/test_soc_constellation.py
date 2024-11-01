import unittest

from datetime import timedelta
import numpy as np
from pydantic import ValidationError

from tatc.utils import field_of_regard_to_swath_width
from tatc.schemas import SOCConstellation, CircularOrbit, Instrument

class TestSOCConstellation(unittest.TestCase):
    def setUp(self):
        self.d420_data = {
            "name": "Test Constellation",
            "orbit": {
                "tle": [
                    "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
                    "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
                ],
                "altitude": 780000,
                "inclination": 86.4,
            },
            "instruments": [{"name": "Test Instrument", "field_of_regard": 180.0}],
            "swath_width": field_of_regard_to_swath_width(altitude=780000, field_of_regard=150),
            "packing_distance": 1
        }
        self.d420_con = SOCConstellation(**self.d420_data)
    
    def test_constructor(self):
        self.assertEqual(self.d420_con.name, self.d420_data.get("name"))
        self.assertEqual(
            self.d420_con.orbit, CircularOrbit(**self.d420_data.get("orbit"))
        )
        self.assertEqual(len(self.d420_con.instruments), 1)
        self.assertEqual(
            self.d420_con.instruments[0],
            Instrument(**self.d420_data.get("instruments")[0]),
        )
        self.assertEqual(
            self.d420_con.swath_width, self.d420_data.get("swath_width")
        )
        self.assertEqual(
            self.d420_con.packing_distance, self.d420_data.get("packing_distance")
        )

    def test_get_num_satellites(self):
        self.assertEqual(
            len(self.d420_con.generate_members()), self.d420_con.generate_walker().number_satellites
        )
    
    def test_get_satellites_per_plane(self):
        self.assertEqual(
            self.d420_con.generate_walker().number_satellites / self.d420_con.generate_walker().number_planes,
            7
        )
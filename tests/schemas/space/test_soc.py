"""
Unit tests for the SOCConstellation schema.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""
import unittest

from tatc.schemas import CircularOrbit, Instrument, SOCConstellation
from tatc.utils import field_of_regard_to_swath_width


class TestSOCConstellation(unittest.TestCase):
    """
    Unit tests for the SOCConstellation schema.
    """
    def setUp(self):
        self.d420_data = {
            "name": "Test Constellation",
            "orbit": {
                "type": "circular",
                "mean_altitude": 780000,
                "inclination": 86.4,
                "epoch": "2000-01-01T00:00:00Z",
            },
            "instruments": [{"name": "Test Instrument", "field_of_regard": 150.0}],
            "swath_width": field_of_regard_to_swath_width(
                altitude=780000, field_of_regard=150
            ),
            "packing_distance": 1,
        }
        self.d420_con = SOCConstellation(**self.d420_data)

    def test_constructor(self):
        """
        Test that the SOCConstellation schema correctly initializes with valid data.
        """
        self.assertEqual(self.d420_con.name, self.d420_data.get("name"))
        self.assertEqual(
            self.d420_con.orbit, CircularOrbit(**self.d420_data.get("orbit"))
        )
        self.assertEqual(len(self.d420_con.instruments), 1)
        self.assertEqual(
            self.d420_con.instruments[0],
            Instrument(**self.d420_data.get("instruments")[0]),
        )
        self.assertEqual(self.d420_con.swath_width, self.d420_data.get("swath_width"))
        self.assertEqual(
            self.d420_con.packing_distance, self.d420_data.get("packing_distance")
        )

    def test_get_num_satellites(self):
        """
        Test that the SOCConstellation schema correctly calculates the number of satellites.
        """
        self.assertEqual(
            len(self.d420_con.generate_members()),
            self.d420_con.generate_walker().number_satellites,
        )

    def test_get_satellites_per_plane(self):
        """
        Test that the SOCConstellation schema correctly calculates the number of satellites per plane.
        """
        self.assertEqual(
            self.d420_con.generate_walker().number_satellites
            / self.d420_con.generate_walker().number_planes,
            7,
        )

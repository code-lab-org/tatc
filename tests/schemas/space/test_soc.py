"""
Unit tests for the SOCConstellation schema.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

from pydantic import ValidationError

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

    def test_type_defaults_to_soc(self):
        """
        Test that omitting type defaults to the "soc" discriminator.
        """
        self.assertEqual(self.d420_con.type, "soc")

    def test_type_rejects_other_values(self):
        """
        Test that type only accepts the "soc" literal, rejecting other
        space system type discriminators.
        """
        with self.assertRaises(ValidationError):
            SOCConstellation(**self.d420_data, type="walker")

    def test_swath_width_must_be_positive(self):
        """
        Test that swath_width must be positive.
        """
        with self.assertRaises(ValidationError):
            SOCConstellation(**{**self.d420_data, "swath_width": 0})

    def test_packing_distance_bounds(self):
        """
        Test that packing_distance must be in (0, 1]: values above 1 would
        space footprint centers farther apart than the footprint diameter,
        leaving gaps and violating the continuous "streets of coverage"
        design goal.
        """
        SOCConstellation(**{**self.d420_data, "packing_distance": 1.0})
        with self.assertRaises(ValidationError):
            SOCConstellation(**{**self.d420_data, "packing_distance": 0})
        with self.assertRaises(ValidationError):
            SOCConstellation(**{**self.d420_data, "packing_distance": 1.1})

    def test_generate_walker_planes_are_hex_offset(self):
        """
        Test that adjacent planes are offset by half a within-plane
        satellite spacing, realizing the staggered hexagonal packing
        implied by the sqrt(3) row spacing (Eq. 24 in Anderson et al.
        2022) rather than a plain rectangular grid of planes. This is a
        regression test for a bug where the generated WalkerConstellation
        left relative_spacing at its default of 0 (no offset).
        """
        walker = self.d420_con.generate_walker()
        self.assertEqual(
            walker.get_delta_mean_anomaly_between_planes(),
            walker.get_delta_mean_anomaly_within_planes() / 2,
        )

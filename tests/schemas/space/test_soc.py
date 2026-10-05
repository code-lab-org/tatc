"""
Unit tests for the SOCConstellation schema.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import math
import unittest

from pydantic import ValidationError

from tatc.constants import EARTH_MEAN_RADIUS
from tatc.schemas import CircularOrbit, Instrument, SOCConstellation
from tatc.utils import field_of_regard_to_swath_width


def iridium_swath_width(altitude=780000, min_elevation_angle=8.2):
    """
    Footprint diameter (m) above a minimum elevation angle from an altitude.
    """
    eta = math.degrees(
        math.asin(
            EARTH_MEAN_RADIUS
            / (EARTH_MEAN_RADIUS + altitude)
            * math.cos(math.radians(min_elevation_angle))
        )
    )
    return 2 * math.radians(90 - min_elevation_angle - eta) * EARTH_MEAN_RADIUS


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
        # the same constellation with the Walker delta pattern
        self.d420_delta = SOCConstellation(**self.d420_data, polar=False)
        self.iridium = SOCConstellation(
            name="Iridium",
            orbit=CircularOrbit(mean_altitude=780000, inclination=86.4),
            swath_width=iridium_swath_width(),
            packing_distance=1,
        )

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
        for con in (self.d420_con, self.d420_delta):
            self.assertEqual(
                len(con.generate_members()),
                con.get_satellites_per_plane() * con.get_number_planes(),
            )

    def test_get_satellites_per_plane(self):
        """
        Test that the SOCConstellation schema correctly calculates the number of satellites per plane.
        """
        self.assertEqual(
            self.d420_delta.generate_walker().number_satellites
            / self.d420_delta.generate_walker().number_planes,
            7,
        )
        self.assertEqual(self.d420_delta.get_satellites_per_plane(), 7)

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
        walker = self.d420_delta.generate_walker()
        self.assertEqual(
            walker.get_delta_mean_anomaly_between_planes(),
            walker.get_delta_mean_anomaly_within_planes() / 2,
        )

    def test_is_polar(self):
        """
        Test that the polar pattern is used by default for near-polar
        orbits (inclination within 10 degrees of 90 degrees), and can be
        forced on or off.
        """
        self.assertTrue(self.d420_con.is_polar())
        self.assertFalse(self.d420_delta.is_polar())
        for inclination, polar in [
            (53, False),
            (79.9, False),
            (80, True),
            (98.6, True),
            (100.1, False),
        ]:
            orbit = {**self.d420_data["orbit"], "inclination": inclination}
            con = SOCConstellation(**{**self.d420_data, "orbit": orbit})
            self.assertEqual(con.is_polar(), polar, inclination)
        orbit = {**self.d420_data["orbit"], "inclination": 53}
        self.assertTrue(
            SOCConstellation(
                **{**self.d420_data, "orbit": orbit, "polar": True}
            ).is_polar()
        )

    def test_polar_design_reproduces_iridium(self):
        """
        Test that the polar design for Iridium's altitude (780 km) and
        minimum elevation angle (8.2 degrees) has Iridium's 6 planes of 11
        satellites, with co-rotating planes about 31.4 degrees apart and a
        seam of about 23 degrees (Iridium: 31.6 and 22 degrees), spanning
        180 degrees.
        """
        satellites_per_plane, number_planes, spacing, seam = (
            self.iridium.get_polar_design()
        )
        self.assertEqual(satellites_per_plane, 11)
        self.assertEqual(number_planes, 6)
        self.assertAlmostEqual(spacing, 31.4, delta=0.1)
        self.assertAlmostEqual(seam, 23.0, delta=0.1)
        self.assertAlmostEqual((number_planes - 1) * spacing + seam, 180)
        self.assertEqual(self.iridium.get_satellites_per_plane(), 11)
        self.assertEqual(self.iridium.get_number_planes(), 6)

    def test_polar_design_satisfies_coverage_conditions(self):
        """
        Test that the polar design's plane spacings do not exceed the
        maxima for continuous coverage (gamma + c for co-rotating planes and
        2c across the seam) and that no design with fewer satellites
        satisfies them.
        """
        gamma = self.iridium.get_footprint_angle()
        s, p, spacing, seam = self.iridium.get_polar_design()
        c = math.degrees(
            math.acos(math.cos(math.radians(gamma)) / math.cos(math.pi / s))
        )
        self.assertLessEqual(spacing, gamma + c)
        self.assertLessEqual(seam, 2 * c)
        for s_other in range(math.floor(180 / gamma) + 1, 40):
            c_other = math.degrees(
                math.acos(math.cos(math.radians(gamma)) / math.cos(math.pi / s_other))
            )
            p_other = 1 + math.ceil((180 - 2 * c_other) / (gamma + c_other) - 1e-9)
            self.assertGreaterEqual(s_other * p_other, s * p)

    def test_polar_packing_distance_adds_margin(self):
        """
        Test that a packing distance below 1 reduces the effective footprint
        of the polar design, requiring more satellites.
        """
        margin = self.iridium.model_copy(update={"packing_distance": 0.9})
        self.assertGreater(
            margin.get_satellites_per_plane() * margin.get_number_planes(), 66
        )

    def test_polar_generate_members(self):
        """
        Test that polar members occupy planes spaced by the co-rotating
        spacing, with adjacent planes offset by half the in-plane spacing.
        """
        members = self.iridium.generate_members()
        self.assertEqual(len(members), 66)
        _, _, spacing, _ = self.iridium.get_polar_design()
        raan = [m.orbit.right_ascension_ascending_node for m in members]
        anomaly = [m.orbit.true_anomaly for m in members]
        for plane in range(6):
            self.assertAlmostEqual(raan[11 * plane], plane * spacing, places=6)
            self.assertAlmostEqual(
                (anomaly[11 * plane] - plane * 180 / 11) % 360, 0, places=6
            )
        self.assertAlmostEqual(anomaly[1] - anomaly[0], 360 / 11, places=6)

    def test_polar_generate_walker_raises(self):
        """
        Test that a polar design, whose planes are unequally spaced, cannot
        be generated as a Walker constellation.
        """
        with self.assertRaises(ValueError):
            self.iridium.generate_walker()

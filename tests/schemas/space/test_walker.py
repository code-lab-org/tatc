"""
Unit tests for the WalkerConstellation schema.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

import numpy as np
from pydantic import ValidationError

from tatc.schemas import CircularOrbit, Instrument, WalkerConstellation


class TestWalkerConstellation(unittest.TestCase):
    """
    Unit tests for the WalkerConstellation schema.
    """

    def setUp(self):
        self.d420_data = {
            "name": "Test Constellation",
            "configuration": "delta",
            "orbit": {
                "mean_altitude": "400000",
                "inclination": 51.6,
                "right_ascension_ascending_node": 180,
                "true_anomaly": 180,
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
                "mean_altitude": "400000",
                "inclination": 51.6,
                "right_ascension_ascending_node": 180,
                "true_anomaly": 180,
            },
            "instruments": [{"name": "Test Instrument", "field_of_regard": 25.0}],
            "number_satellites": 4,
            "number_planes": 2,
            "relative_spacing": 0,
        }
        self.s420_con = WalkerConstellation(**self.s420_data)

    def normalize_angle(self, angle):
        """
        Normalize an angle to the range [0, 360) degrees.
        """
        return np.mod(360 + angle, 360)

    def test_invalid_planes(self):
        """
        Test that the WalkerConstellation schema raises a ValidationError when the number of planes is invalid.
        """
        bad_data = {
            "name": "Test Constellation",
            "configuration": "delta",
            "orbit": {
                "mean_altitude": "400000",
                "inclination": 51.6,
                "right_ascension_ascending_node": 180,
                "true_anomaly": 180,
            },
            "instruments": [{"name": "Test Instrument", "field_of_regard": 25.0}],
            "number_satellites": 5,
            "number_planes": 6,
            "relative_spacing": 1,
        }
        with self.assertRaises(ValidationError):
            WalkerConstellation(**bad_data)

    def test_invalid_relative_spacing(self):
        """
        Test that the WalkerConstellation schema raises a ValidationError when the relative spacing is invalid.
        """
        bad_data = {
            "name": "Test Constellation",
            "configuration": "delta",
            "orbit": {
                "mean_altitude": "400000",
                "inclination": 51.6,
                "right_ascension_ascending_node": 180,
                "true_anomaly": 180,
            },
            "instruments": [{"name": "Test Instrument", "field_of_regard": 25.0}],
            "number_satellites": 5,
            "number_planes": 2,
            "relative_spacing": 5,
        }
        with self.assertRaises(ValidationError):
            WalkerConstellation(**bad_data)

    def test_number_planes_equal_to_number_satellites_is_valid(self):
        """
        Test that number_planes equal to number_satellites (the boundary
        of number_planes_le_number_satellites) is valid, not rejected.
        """
        con = WalkerConstellation(**{**self.d420_data, "number_planes": 4})
        self.assertEqual(con.number_planes, 4)

    def test_relative_spacing_equal_to_number_planes_minus_one_is_valid(self):
        """
        Test that relative_spacing equal to number_planes - 1 (the
        boundary of relative_spacing_lt_number_planes) is valid, not
        rejected.
        """
        con = WalkerConstellation(**{**self.d420_data, "relative_spacing": 1})
        self.assertEqual(con.relative_spacing, 1)

    def test_defaults(self):
        """
        Test that optional fields default to a single-satellite,
        single-plane Walker Delta constellation.
        """
        con = WalkerConstellation(
            name="Test Constellation", orbit=self.d420_data["orbit"]
        )
        self.assertEqual(con.configuration, "delta")
        self.assertEqual(con.number_satellites, 1)
        self.assertEqual(con.number_planes, 1)
        self.assertEqual(con.relative_spacing, 0)

    def test_type_defaults_to_walker(self):
        """
        Test that omitting type defaults to the "walker" discriminator.
        """
        self.assertEqual(self.d420_con.type, "walker")

    def test_type_rejects_other_values(self):
        """
        Test that type only accepts the "walker" literal, rejecting other
        space system type discriminators.
        """
        with self.assertRaises(ValidationError):
            WalkerConstellation(**self.d420_data, type="train")

    def test_constructor(self):
        """
        Test that the WalkerConstellation schema correctly initializes with valid data.
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
        """
        Test that the number of satellites per plane can be calculated correctly.
        """
        self.assertEqual(
            self.d420_con.get_satellites_per_plane(),
            np.ceil(self.d420_con.number_satellites / self.d420_con.number_planes),
        )
        self.assertEqual(
            self.s420_con.get_satellites_per_plane(),
            np.ceil(self.s420_con.number_satellites / self.s420_con.number_planes),
        )

    def test_get_delta_mean_anomaly_within_planes(self):
        """
        Test that the delta mean anomaly within planes can be calculated correctly.
        """
        self.assertEqual(
            self.d420_con.get_delta_mean_anomaly_within_planes(),
            360 / self.d420_con.get_satellites_per_plane(),
        )
        self.assertEqual(
            self.s420_con.get_delta_mean_anomaly_within_planes(),
            360 / self.s420_con.get_satellites_per_plane(),
        )

    def test_get_delta_mean_anomaly_between_planes(self):
        """
        Test that the delta mean anomaly between planes can be calculated correctly.
        """
        self.assertEqual(
            self.d420_con.get_delta_mean_anomaly_between_planes(),
            self.d420_con.relative_spacing * 360 / self.d420_con.number_satellites,
        )
        self.assertEqual(
            self.s420_con.get_delta_mean_anomaly_between_planes(),
            self.s420_con.relative_spacing * 360 / self.s420_con.number_satellites,
        )

    def test_get_delta_raan_between_planes_delta(self):
        """
        Test that the delta right ascension of ascending node between planes can be calculated correctly for delta configuration.
        """
        self.assertEqual(
            self.d420_con.get_delta_raan_between_planes(),
            360 / self.d420_con.number_planes,
        )

    def test_get_delta_raan_between_planes_star(self):
        """
        Test that the delta right ascension of ascending node between planes can be calculated correctly for star configuration.
        """
        self.assertEqual(
            self.s420_con.get_delta_raan_between_planes(),
            180 / self.s420_con.number_planes,
        )

    def helper_test_generate_members(self, constellation):
        """
        Helper function to test that the WalkerConstellation schema correctly
        generates constellation members with the specified parameters
        """
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

    def test_generate_members_delta(self):
        """
        Test that the WalkerConstellation schema correctly generates
        constellation members for delta configuration with TLE orbit.
        """
        self.helper_test_generate_members(self.d420_con)

    def test_generate_members_star(self):
        """
        Test that the WalkerConstellation schema correctly generates
        constellation members for star configuration with TLE orbit.
        """
        self.helper_test_generate_members(self.s420_con)

    def test_generate_members_first_satellite_matches_lead_orbit(self):
        """
        Test that the first generated member (index 0, zero mean anomaly
        and RAAN offset) matches the constellation's own lead orbit.
        """
        members = self.d420_con.generate_members()
        self.assertEqual(members[0].orbit, CircularOrbit(**self.d420_data["orbit"]))

    def test_generate_members_instruments_independent_per_satellite(self):
        """
        Test that each generated member gets its own deep-copied
        instruments list, not one shared (aliased) across members or with
        the constellation itself: mutating one member's instrument must
        not affect another member's or the constellation's.
        """
        members = self.d420_con.generate_members()
        self.assertIsNot(members[0].instruments, members[1].instruments)
        self.assertIsNot(members[0].instruments, self.d420_con.instruments)
        members[0].instruments[0].name = "Renamed"
        self.assertEqual(members[1].instruments[0].name, "Test Instrument")
        self.assertEqual(self.d420_con.instruments[0].name, "Test Instrument")

    def test_generate_members_uneven_plane_distribution(self):
        """
        Test the documented "(max) number of satellites per plane"
        generalization beyond the strict Walker definition: when
        number_satellites does not divide evenly by number_planes,
        get_satellites_per_plane() ceils, so some planes get fewer
        satellites than others (5 satellites / 2 planes -> 3 + 2, not a
        ValidationError), and every plane shares the same mean-anomaly
        slot spacing (360 / 3 degrees here) rather than each plane
        re-spreading its own satellites evenly across the full 360
        degrees. The under-full plane's satellites are therefore
        clustered within a portion of the plane, not evenly spaced
        around it.
        """
        con = WalkerConstellation(**{**self.d420_data, "number_satellites": 5})
        self.assertEqual(con.get_satellites_per_plane(), 3)
        members = con.generate_members()
        self.assertEqual(len(members), 5)
        plane_0 = [m.orbit.right_ascension_ascending_node for m in members[0:3]]
        plane_1 = [m.orbit.right_ascension_ascending_node for m in members[3:5]]
        self.assertEqual(len(set(plane_0)), 1)
        self.assertEqual(len(set(plane_1)), 1)
        self.assertNotEqual(plane_0[0], plane_1[0])
        # plane 1's 2 satellites are spaced by a full plane's slot spacing
        # (120 degrees, the same as within plane 0), not by 180 degrees
        # (360 / 2, which would evenly spread only 2 satellites)
        self.assertAlmostEqual(
            self.normalize_angle(
                members[4].orbit.get_mean_anomaly()
                - members[3].orbit.get_mean_anomaly()
            ),
            con.get_delta_mean_anomaly_within_planes(),
            delta=0.001,
        )
        self.assertNotAlmostEqual(
            self.normalize_angle(
                members[4].orbit.get_mean_anomaly()
                - members[3].orbit.get_mean_anomaly()
            ),
            180.0,
            delta=0.001,
        )

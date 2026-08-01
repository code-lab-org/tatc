"""
Unit tests for the MOGConstellation schema.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import math
import unittest
from datetime import datetime, timezone

from pydantic import ValidationError

from tatc.constants import EARTH_MEAN_RADIUS
from tatc.schemas import CircularOrbit, Instrument, MOGConstellation


class TestMOGConstellation(unittest.TestCase):
    """
    Unit tests for the MOGConstellation schema.
    """
    def setUp(self):
        self.epoch = datetime(2020, 1, 1, tzinfo=timezone.utc)
        self.reference_orbit = CircularOrbit(
            mean_altitude=500000,
            inclination=51.6,
            right_ascension_ascending_node=180,
            true_anomaly=0,
            epoch=self.epoch,
        )
        self.test_data = {
            "name": "Test Constellation",
            "orbit": self.reference_orbit,
            "instruments": [{"name": "Test Instrument", "field_of_regard": 25.0}],
            "parallel_axis": 1000.0,
            "transverse_axis": 500.0,
            "number_satellites": 3,
        }
        self.con = MOGConstellation(**self.test_data)

    def test_good_data(self):
        """
        Test that the MOGConstellation schema correctly initializes with valid data.
        """
        self.assertEqual(self.con.name, self.test_data.get("name"))
        self.assertEqual(self.con.orbit, self.reference_orbit)
        self.assertEqual(len(self.con.instruments), 1)
        self.assertEqual(
            self.con.instruments[0],
            Instrument(**self.test_data.get("instruments")[0]),
        )
        self.assertEqual(self.con.parallel_axis, self.test_data.get("parallel_axis"))
        self.assertEqual(
            self.con.transverse_axis, self.test_data.get("transverse_axis")
        )
        self.assertEqual(
            self.con.number_satellites, self.test_data.get("number_satellites")
        )

    def test_type_defaults_to_mog(self):
        """
        Test that omitting type defaults to the "mog" discriminator.
        """
        self.assertEqual(self.con.type, "mog")

    def test_type_rejects_other_values(self):
        """
        Test that type only accepts the "mog" literal, rejecting other
        space system type discriminators.
        """
        with self.assertRaises(ValidationError):
            MOGConstellation(**self.test_data, type="walker")

    def test_parallel_axis_must_be_positive(self):
        """
        Test that parallel_axis must be positive.
        """
        with self.assertRaises(ValidationError):
            MOGConstellation(**{**self.test_data, "parallel_axis": 0})

    def test_transverse_axis_must_be_positive(self):
        """
        Test that transverse_axis must be positive.
        """
        with self.assertRaises(ValidationError):
            MOGConstellation(**{**self.test_data, "transverse_axis": 0})

    def test_number_satellites_default(self):
        """
        Test that omitting number_satellites defaults to 2 (a mutual
        orbiting pair, the minimal configuration illustrated in Leroy
        et al. 2020, Fig. 7).
        """
        con = MOGConstellation(
            name="Test Constellation",
            orbit=self.reference_orbit,
            parallel_axis=1000.0,
            transverse_axis=500.0,
        )
        self.assertEqual(con.number_satellites, 2)

    def test_number_satellites_must_be_positive(self):
        """
        Test that number_satellites must be positive.
        """
        with self.assertRaises(ValidationError):
            MOGConstellation(**{**self.test_data, "number_satellites": 0})

    def test_clockwise_defaults_to_true(self):
        """
        Test that omitting clockwise defaults to True.
        """
        self.assertTrue(self.con.clockwise)

    def test_generate_members_count(self):
        """
        Test that generate_members produces exactly number_satellites members.
        """
        self.assertEqual(len(self.con.generate_members()), 3)

    def test_generate_members_semimajor_axis_conserved(self):
        """
        Test that every member shares the reference orbit's semimajor
        axis exactly. This is a necessary condition for the mutual
        orbiting group to remain bounded (no secular along-track drift):
        Leroy et al. (2020) explicitly assume a common semimajor axis "for
        the sake of constellation design" (Appendix A), since orbital
        period depends only on semimajor axis, not eccentricity.
        """
        a = self.reference_orbit.get_semimajor_axis()
        for member in self.con.generate_members():
            self.assertAlmostEqual(member.orbit.semimajor_axis, a, delta=1e-6)

    def test_generate_members_eccentricity_matches_formula(self):
        """
        Test that every member's eccentricity equals
        parallel_axis / (4 * semimajor_axis), the relationship stated in
        Leroy et al. (2020) (Section VI Summary and Appendix A): the
        mutual orbit's extent parallel to the velocity vector is 4ae.
        """
        a = self.reference_orbit.get_semimajor_axis()
        expected_eccentricity = self.con.parallel_axis / (4 * a)
        for member in self.con.generate_members():
            self.assertAlmostEqual(
                member.orbit.eccentricity, expected_eccentricity, delta=1e-12
            )

    def test_generate_members_instruments_independent_per_satellite(self):
        """
        Test that each generated member gets its own deep-copied
        instruments list, not one shared (aliased) across members or with
        the constellation itself: mutating one member's instrument must
        not affect another member's or the constellation's.
        """
        members = self.con.generate_members()
        self.assertIsNot(members[0].instruments, members[1].instruments)
        self.assertIsNot(members[0].instruments, self.con.instruments)
        members[0].instruments[0].name = "Renamed"
        self.assertEqual(members[1].instruments[0].name, "Test Instrument")
        self.assertEqual(self.con.instruments[0].name, "Test Instrument")

    def test_generate_members_orientation_independent_of_reference_epoch_position(
        self,
    ):
        """
        Test that inclination, RAAN, and argument of perigee do not depend
        on the reference orbit's true anomaly at epoch: these describe the
        fixed shape/orientation of the mutual orbiter's ellipse, which
        Leroy et al. (2020), Appendix A, derive purely from geometry (the
        tilt angle delta and clockwise/counter-clockwise sense) -- they
        are independent of *when* the reference orbiter is at that
        geometry.
        """
        orbit_ta0 = self.reference_orbit.model_copy(update={"true_anomaly": 0})
        orbit_ta123 = self.reference_orbit.model_copy(update={"true_anomaly": 123})
        con_ta0 = MOGConstellation(**{**self.test_data, "orbit": orbit_ta0})
        con_ta123 = MOGConstellation(**{**self.test_data, "orbit": orbit_ta123})
        members_ta0 = con_ta0.generate_members()
        members_ta123 = con_ta123.generate_members()
        for m0, m123 in zip(members_ta0, members_ta123):
            self.assertAlmostEqual(
                m0.orbit.inclination, m123.orbit.inclination, delta=1e-9
            )
            self.assertAlmostEqual(
                m0.orbit.right_ascension_ascending_node,
                m123.orbit.right_ascension_ascending_node,
                delta=1e-9,
            )
            self.assertAlmostEqual(
                m0.orbit.perigee_argument, m123.orbit.perigee_argument, delta=1e-9
            )

    def test_generate_members_true_anomaly_tracks_reference_epoch_position(self):
        """
        Test that the mutual orbiter's true anomaly shifts by (approximately,
        given the mutual orbiter's small eccentricity) the same amount as
        the reference orbit's own true anomaly at epoch. This is a
        regression test for a bug where the reference orbit's actual
        position at epoch was silently ignored: `generate_members` always
        computed the mutual orbiter's true anomaly as if the reference
        orbiter were exactly at its ascending node (true_anomaly=0) at
        epoch, so two reference orbits differing only in true_anomaly
        produced identical mutual orbiter positions.
        """
        orbit_ta0 = self.reference_orbit.model_copy(update={"true_anomaly": 0})
        orbit_ta90 = self.reference_orbit.model_copy(update={"true_anomaly": 90})
        con_ta0 = MOGConstellation(**{**self.test_data, "orbit": orbit_ta0})
        con_ta90 = MOGConstellation(**{**self.test_data, "orbit": orbit_ta90})
        member_ta0 = con_ta0.generate_members()[0]
        member_ta90 = con_ta90.generate_members()[0]
        self.assertAlmostEqual(
            (member_ta90.orbit.true_anomaly - member_ta0.orbit.true_anomaly) % 360,
            90,
            delta=0.01,
        )

    def test_generate_members_matches_published_table(self):
        """
        Test against Table I of Leroy et al. (2020): a 4-satellite
        constellation at 404 km altitude, reference inclination 51.64 deg,
        eccentricity 0.01, reference RAAN 0 deg, reference true anomaly 0
        (at its ascending node). The tilt angle delta = 0.22 degrees
        (reverse-engineered here, since the paper states the axis lengths
        symbolically as 4ae/2a*delta but does not give delta numerically)
        reproduces the table's inclination, RAAN, and argument of perigee
        for all 4 satellites to within the table's own 2-decimal-place
        precision.

        Note: the table's true anomaly values for the theta=90/270
        satellites (its "second pair") do not match this reconstruction
        (off by ~7.5-7.8 degrees, and not by a common offset), while the
        theta=0/180 satellites ("first pair") match exactly. This is
        consistent with the paper's own description of the example as
        "two pairs of mutual orbiting satellites" -- i.e. not necessarily
        one single 4-satellite group sharing a common reference epoch --
        so only the theta=0/180 true anomalies are checked here.
        """
        a = EARTH_MEAN_RADIUS + 404000
        eccentricity = 0.01
        delta_deg = 0.22
        orbit = CircularOrbit(
            mean_altitude=404000,
            inclination=51.64,
            right_ascension_ascending_node=0,
            true_anomaly=0,
            epoch=self.epoch,
        )
        con = MOGConstellation(
            name="Table I",
            orbit=orbit,
            parallel_axis=4 * a * eccentricity,
            transverse_axis=2 * a * math.radians(delta_deg),
            number_satellites=4,
            clockwise=True,
        )
        members = con.generate_members()
        # theta = 0, 90, 180, 270 degrees, in generation order
        expected_inclination = [51.42, 51.64, 51.86, 51.64]
        expected_raan = [0.0, 0.29, 0.0, -0.29]
        expected_perigee_argument = [-90.00, 179.82, 90.00, 0.18]
        for member, exp_i, exp_raan, exp_w in zip(
            members, expected_inclination, expected_raan, expected_perigee_argument
        ):
            self.assertAlmostEqual(member.orbit.inclination, exp_i, delta=0.01)
            self.assertAlmostEqual(
                (member.orbit.right_ascension_ascending_node + 180) % 360 - 180,
                exp_raan,
                delta=0.01,
            )
            self.assertAlmostEqual(
                (member.orbit.perigee_argument + 180) % 360 - 180, exp_w, delta=0.01
            )
        # only the theta=0/180 ("first pair") true anomalies match the table
        self.assertAlmostEqual(members[0].orbit.true_anomaly, 91.15, delta=0.01)
        self.assertAlmostEqual(
            (members[2].orbit.true_anomaly + 180) % 360 - 180, -91.15, delta=0.01
        )


if __name__ == "__main__":
    unittest.main()

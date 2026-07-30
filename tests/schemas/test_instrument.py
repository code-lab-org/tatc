"""
Unit tests for the Instrument schema.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

from skyfield.api import EarthSatellite

from tatc.constants import timescale
from tatc.schemas import CircularOrbit, Instrument


class TestInstrument(unittest.TestCase):
    """
    Unit tests for the Instrument schema.
    """
    def setUp(self):
        noon_utc = datetime(2020, 3, 20, 12, tzinfo=timezone.utc)
        self.test_time = timescale.from_datetime(noon_utc)
        self.test_sat_1 = EarthSatellite.from_satrec(
            CircularOrbit(
                mean_altitude=400000,
                true_anomaly=0,
                epoch=noon_utc,
                inclination=0.0,
                right_ascension_ascending_node=0.0,
            ).to_gp_orbit().elements[0].to_satrec(),
            timescale
        )
        self.test_sat_2 = EarthSatellite.from_satrec(
            CircularOrbit(
                mean_altitude=400000,
                true_anomaly=0,
                epoch=noon_utc,
                inclination=0.0,
                right_ascension_ascending_node=80.0,
            ).to_gp_orbit().elements[0].to_satrec(),
            timescale,
        )
        self.test_sat_3 = EarthSatellite.from_satrec(
            CircularOrbit(
                mean_altitude=400000,
                true_anomaly=0,
                epoch=noon_utc,
                inclination=0.0,
                right_ascension_ascending_node=100.0,
            ).to_gp_orbit().elements[0].to_satrec(),
            timescale,
        )
        self.test_sat_4 = EarthSatellite.from_satrec(
            CircularOrbit(
                mean_altitude=400000,
                true_anomaly=0,
                epoch=noon_utc,
                inclination=0.0,
                right_ascension_ascending_node=180.0,
            ).to_gp_orbit().elements[0].to_satrec(),
            timescale,
        )
        self.test_sat_5 = EarthSatellite.from_satrec(
            CircularOrbit(
                mean_altitude=400000,
                true_anomaly=0,
                epoch=noon_utc,
                inclination=45.0,
                right_ascension_ascending_node=0.0,
            ).to_gp_orbit().elements[0].to_satrec(),
            timescale,
        )

    def test_good_data(self):
        """
        Test that an Instrument can be created with valid data.
        """
        good_data = {
            "name": "Test Instrument",
            "field_of_regard": 20.0,
            "min_access_time": timedelta(seconds=10),
            "req_self_sunlit": None,
            "req_target_sunlit": None,
        }
        o = Instrument(**good_data)
        self.assertEqual(o.name, good_data.get("name"))
        self.assertEqual(o.field_of_regard, good_data.get("field_of_regard"))
        self.assertEqual(o.min_access_time, good_data.get("min_access_time"))
        self.assertEqual(o.req_self_sunlit, good_data.get("req_self_sunlit"))
        self.assertEqual(o.req_target_sunlit, good_data.get("req_target_sunlit"))

    def test_get_swath_width(self):
        """
        Test that the swath width can be computed from the field of regard.
        """
        o = Instrument(name="GMI", field_of_regard=15.0)
        self.assertAlmostEqual(o.get_swath_width(705000), 185815, delta=1.0)

    def test_get_min_elevation_angle(self):
        """
        Test that the minimum elevation angle can be computed from the field of regard.
        """
        o = Instrument(name="GMI", field_of_regard=15.0)
        self.assertAlmostEqual(o.get_min_elevation_angle(705000), 81.66446, delta=0.01)

    def test_valid_observation_no_constraints(self):
        """
        Test that an observation is valid when there are no constraints 
        on sunlit conditions.
        """
        o = Instrument(name="Test Instrument")
        self.assertTrue(o.is_valid_observation(self.test_sat_1.at(self.test_time)).all()) # type: ignore
        self.assertTrue(o.is_valid_observation(self.test_sat_2.at(self.test_time)).all()) # type: ignore
        self.assertTrue(o.is_valid_observation(self.test_sat_3.at(self.test_time)).all()) # type: ignore
        self.assertTrue(o.is_valid_observation(self.test_sat_4.at(self.test_time)).all()) # type: ignore

    def test_valid_observation_self_sunlit(self):
        """
        Test that an observation is valid when the instrument requires 
        self-sunlit conditions.
        """
        o = Instrument(name="Test Instrument", req_self_sunlit=True)
        self.assertTrue(o.is_valid_observation(self.test_sat_1.at(self.test_time)).all()) # type: ignore
        self.assertTrue(o.is_valid_observation(self.test_sat_2.at(self.test_time)).all()) # type: ignore
        self.assertTrue(o.is_valid_observation(self.test_sat_3.at(self.test_time)).all()) # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_4.at(self.test_time)).all()) # type: ignore

    def test_valid_observation_self_not_sunlit(self):
        """
        Test that an observation is valid when the instrument requires 
        self-not-sunlit conditions.
        """
        o = Instrument(name="Test Instrument", req_self_sunlit=False)
        self.assertFalse(o.is_valid_observation(self.test_sat_1.at(self.test_time)).all()) # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_2.at(self.test_time)).all()) # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_3.at(self.test_time)).all()) # type: ignore
        self.assertTrue(o.is_valid_observation(self.test_sat_4.at(self.test_time)).all()) # type: ignore

    def test_valid_observation_target_sunlit(self):
        """
        Test that an observation is valid when the instrument requires 
        target-sunlit conditions.
        """
        o = Instrument(name="Test Instrument", req_target_sunlit=True)
        self.assertTrue(o.is_valid_observation(self.test_sat_1.at(self.test_time)).all()) # type: ignore
        self.assertTrue(o.is_valid_observation(self.test_sat_2.at(self.test_time)).all()) # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_3.at(self.test_time)).all()) # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_4.at(self.test_time)).all()) # type: ignore

    def test_valid_observation_target_not_sunlit(self):
        """
        Test that an observation is valid when the instrument requires 
        target-not-sunlit conditions.
        """
        o = Instrument(name="Test Instrument", req_target_sunlit=False)
        self.assertFalse(o.is_valid_observation(self.test_sat_1.at(self.test_time)).all()) # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_2.at(self.test_time)).all()) # type: ignore
        self.assertTrue(o.is_valid_observation(self.test_sat_3.at(self.test_time)).all()) # type: ignore
        self.assertTrue(o.is_valid_observation(self.test_sat_4.at(self.test_time)).all()) # type: ignore

    def test_valid_observation_self_sunlit_target_sunlit(self):
        """
        Test that an observation is valid when the instrument requires 
        both self-sunlit and target-sunlit conditions."""
        o = Instrument(
            name="Test Instrument", req_self_sunlit=True, req_target_sunlit=True
        )
        self.assertTrue(o.is_valid_observation(self.test_sat_1.at(self.test_time)).all()) # type: ignore
        self.assertTrue(o.is_valid_observation(self.test_sat_2.at(self.test_time)).all()) # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_3.at(self.test_time)).all()) # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_4.at(self.test_time)).all()) # type: ignore

    def test_valid_observation_self_not_sunlit_target_sunlit(self):
        """
        Test that an observation is valid when the instrument requires 
        self-not-sunlit and target-sunlit conditions.
        """
        o = Instrument(
            name="Test Instrument", req_self_sunlit=False, req_target_sunlit=True
        )
        self.assertFalse(o.is_valid_observation(self.test_sat_1.at(self.test_time)).all()) # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_2.at(self.test_time)).all()) # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_3.at(self.test_time)).all()) # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_4.at(self.test_time)).all()) # type: ignore

    def test_valid_observation_self_sunlit_target_not_sunlit(self):
        """
        Test that an observation is valid when the instrument requires 
        self-sunlit and target-not-sunlit conditions.
        """
        o = Instrument(
            name="Test Instrument", req_self_sunlit=True, req_target_sunlit=False
        )
        self.assertFalse(o.is_valid_observation(self.test_sat_1.at(self.test_time)).all()) # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_2.at(self.test_time)).all()) # type: ignore
        self.assertTrue(o.is_valid_observation(self.test_sat_3.at(self.test_time)).all()) # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_4.at(self.test_time)).all()) # type: ignore

    def test_valid_observation_self_not_sunlit_target_not_sunlit(self):
        """
        Test that an observation is valid when the instrument requires 
        both self-not-sunlit and target-not-sunlit conditions.
        """
        o = Instrument(
            name="Test Instrument", req_self_sunlit=False, req_target_sunlit=False
        )
        self.assertFalse(o.is_valid_observation(self.test_sat_1.at(self.test_time)).all()) # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_2.at(self.test_time)).all()) # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_3.at(self.test_time)).all()) # type: ignore
        self.assertTrue(o.is_valid_observation(self.test_sat_4.at(self.test_time)).all()) # type: ignore

    def test_valid_observation_self_sunlit_vector(self):
        """
        Test that an observation is valid when the instrument requires 
        self-sunlit conditions for a vector of times.
        """
        o = Instrument(name="Test Instrument", req_self_sunlit=True)
        times = timescale.utc(2020, 3, 20, [11, 12, 13]) # type: ignore
        results = o.is_valid_observation(self.test_sat_1.at(times)) # type: ignore
        self.assertEqual(len(results), 3)
        self.assertFalse(results[0])
        self.assertTrue(results[1])
        self.assertFalse(results[2])

    def test_valid_observation_target_sunlit_vector(self):
        """
        Test that an observation is valid when the instrument requires 
        target-sunlit conditions for a vector of times.
        """
        o = Instrument(name="Test Instrument", req_target_sunlit=True)
        times = timescale.utc(2020, 3, 20, [11, 12, 13]) # type: ignore
        results = o.is_valid_observation(self.test_sat_1.at(times)) # type: ignore
        self.assertEqual(len(results), 3)
        self.assertFalse(results[0])
        self.assertTrue(results[1])
        self.assertFalse(results[2])

    def test_valid_observation_self_sunlit_vector_inclined(self):
        """
        Test that an observation is valid when the instrument requires self-sunlit conditions for a vector of times with an inclined orbit.
        """
        o = Instrument(name="Test Instrument", req_self_sunlit=True)
        times = timescale.utc(2020, 3, 20, [11, 12, 13]) # type: ignore
        results = o.is_valid_observation(self.test_sat_5.at(times)) # type: ignore
        self.assertEqual(len(results), 3)
        self.assertFalse(results[0])
        self.assertTrue(results[1])
        self.assertFalse(results[2])

    def test_valid_observation_target_sunlit_vector_inclined(self):
        """
        Test that an observation is valid when the instrument requires target-sunlit conditions for a vector of times with an inclined orbit.
        """
        o = Instrument(name="Test Instrument", req_target_sunlit=True)
        times = timescale.utc(2020, 3, 20, [11, 12, 13]) # type: ignore
        results = o.is_valid_observation(self.test_sat_5.at(times)) # type: ignore
        self.assertEqual(len(results), 3)
        self.assertFalse(results[0])
        self.assertTrue(results[1])
        self.assertFalse(results[2])

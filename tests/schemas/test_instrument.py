import unittest

from tatc.schemas import Instrument, CircularOrbit
from pydantic import ValidationError
from datetime import datetime, timedelta, timezone
from tatc.constants import timescale
from skyfield.api import wgs84, EarthSatellite
from sgp4.api import Satrec, WGS72


class TestInstrument(unittest.TestCase):
    def setUp(self):
        noon_utc = datetime(2020, 3, 20, 12, tzinfo=timezone.utc)
        self.test_time = timescale.from_datetime(noon_utc)
        test_orbit_1 = CircularOrbit(
            altitude=400000,
            true_anomaly=0,
            epoch=noon_utc,
            inclination=0.0,
            right_ascension_ascending_node=140.0,
        ).to_tle()
        self.test_sat_1 = EarthSatellite.from_satrec(
            Satrec.twoline2rv(test_orbit_1.tle[0], test_orbit_1.tle[1], WGS72),
            timescale,
        )
        test_orbit_2 = CircularOrbit(
            altitude=400000,
            true_anomaly=0,
            epoch=noon_utc,
            inclination=0.0,
            right_ascension_ascending_node=225.0,
        ).to_tle()
        self.test_sat_2 = EarthSatellite.from_satrec(
            Satrec.twoline2rv(test_orbit_2.tle[0], test_orbit_2.tle[1], WGS72),
            timescale,
        )
        test_orbit_3 = CircularOrbit(
            altitude=400000,
            true_anomaly=0,
            epoch=noon_utc,
            inclination=0.0,
            right_ascension_ascending_node=235.0,
        ).to_tle()
        self.test_sat_3 = EarthSatellite.from_satrec(
            Satrec.twoline2rv(test_orbit_3.tle[0], test_orbit_3.tle[1], WGS72),
            timescale,
        )
        test_orbit_4 = CircularOrbit(
            altitude=400000,
            true_anomaly=0,
            epoch=noon_utc,
            inclination=0.0,
            right_ascension_ascending_node=320.0,
        ).to_tle()
        self.test_sat_4 = EarthSatellite.from_satrec(
            Satrec.twoline2rv(test_orbit_4.tle[0], test_orbit_4.tle[1], WGS72),
            timescale,
        )

    def test_good_data(self):
        good_data = {
            "name": "Test Instrument",
            "field_of_regard": 20.0,
            "min_access_time": timedelta(seconds=10),
            "req_self_sunlit": None,
            "req_target_sunlit": None,
            "duty_cycle": 1.0,
            "duty_cycle_scheme": "fixed",
        }
        o = Instrument(**good_data)
        self.assertEqual(o.name, good_data.get("name"))
        self.assertEqual(o.field_of_regard, good_data.get("field_of_regard"))
        self.assertEqual(o.min_access_time, good_data.get("min_access_time"))
        self.assertEqual(o.req_self_sunlit, good_data.get("req_self_sunlit"))
        self.assertEqual(o.req_target_sunlit, good_data.get("req_target_sunlit"))
        self.assertEqual(o.duty_cycle, good_data.get("duty_cycle"))
        self.assertEqual(o.duty_cycle_scheme, good_data.get("duty_cycle_scheme"))

    def test_valid_observation_no_constraints(self):
        o = Instrument(name="Test Instrument")
        self.assertTrue(o.is_valid_observation(self.test_sat_1, self.test_time))
        self.assertTrue(o.is_valid_observation(self.test_sat_2, self.test_time))
        self.assertTrue(o.is_valid_observation(self.test_sat_3, self.test_time))
        self.assertTrue(o.is_valid_observation(self.test_sat_4, self.test_time))

    def test_valid_observation_self_sunlit(self):
        o = Instrument(name="Test Instrument", req_self_sunlit=True)
        self.assertTrue(o.is_valid_observation(self.test_sat_1, self.test_time))
        self.assertTrue(o.is_valid_observation(self.test_sat_2, self.test_time))
        self.assertTrue(o.is_valid_observation(self.test_sat_3, self.test_time))
        self.assertFalse(o.is_valid_observation(self.test_sat_4, self.test_time))

    def test_valid_observation_self_not_sunlit(self):
        o = Instrument(name="Test Instrument", req_self_sunlit=False)
        self.assertFalse(o.is_valid_observation(self.test_sat_1, self.test_time))
        self.assertFalse(o.is_valid_observation(self.test_sat_2, self.test_time))
        self.assertFalse(o.is_valid_observation(self.test_sat_3, self.test_time))
        self.assertTrue(o.is_valid_observation(self.test_sat_4, self.test_time))

    def test_valid_observation_target_sunlit(self):
        o = Instrument(name="Test Instrument", req_target_sunlit=True)
        self.assertTrue(o.is_valid_observation(self.test_sat_1, self.test_time))
        self.assertTrue(o.is_valid_observation(self.test_sat_2, self.test_time))
        self.assertFalse(o.is_valid_observation(self.test_sat_3, self.test_time))
        self.assertFalse(o.is_valid_observation(self.test_sat_4, self.test_time))

    def test_valid_observation_target_not_sunlit(self):
        o = Instrument(name="Test Instrument", req_target_sunlit=False)
        self.assertFalse(o.is_valid_observation(self.test_sat_1, self.test_time))
        self.assertFalse(o.is_valid_observation(self.test_sat_2, self.test_time))
        self.assertTrue(o.is_valid_observation(self.test_sat_3, self.test_time))
        self.assertTrue(o.is_valid_observation(self.test_sat_4, self.test_time))

    def test_valid_observation_self_sunlit_target_sunlit(self):
        o = Instrument(
            name="Test Instrument", req_self_sunlit=True, req_target_sunlit=True
        )
        self.assertTrue(o.is_valid_observation(self.test_sat_1, self.test_time))
        self.assertTrue(o.is_valid_observation(self.test_sat_2, self.test_time))
        self.assertFalse(o.is_valid_observation(self.test_sat_3, self.test_time))
        self.assertFalse(o.is_valid_observation(self.test_sat_4, self.test_time))

    def test_valid_observation_self_not_sunlit_target_sunlit(self):
        o = Instrument(
            name="Test Instrument", req_self_sunlit=False, req_target_sunlit=True
        )
        self.assertFalse(o.is_valid_observation(self.test_sat_1, self.test_time))
        self.assertFalse(o.is_valid_observation(self.test_sat_2, self.test_time))
        self.assertFalse(o.is_valid_observation(self.test_sat_3, self.test_time))
        self.assertFalse(o.is_valid_observation(self.test_sat_4, self.test_time))

    def test_valid_observation_self_sunlit_target_not_sunlit(self):
        o = Instrument(
            name="Test Instrument", req_self_sunlit=True, req_target_sunlit=False
        )
        self.assertFalse(o.is_valid_observation(self.test_sat_1, self.test_time))
        self.assertFalse(o.is_valid_observation(self.test_sat_2, self.test_time))
        self.assertTrue(o.is_valid_observation(self.test_sat_3, self.test_time))
        self.assertFalse(o.is_valid_observation(self.test_sat_4, self.test_time))

    def test_valid_observation_self_not_sunlit_target_not_sunlit(self):
        o = Instrument(
            name="Test Instrument", req_self_sunlit=False, req_target_sunlit=False
        )
        self.assertFalse(o.is_valid_observation(self.test_sat_1, self.test_time))
        self.assertFalse(o.is_valid_observation(self.test_sat_2, self.test_time))
        self.assertFalse(o.is_valid_observation(self.test_sat_3, self.test_time))
        self.assertTrue(o.is_valid_observation(self.test_sat_4, self.test_time))

    def test_valid_observation_self_sunlit_vector(self):
        o = Instrument(name="Test Instrument", req_self_sunlit=True)
        times = timescale.utc(2020, 3, 20, [11, 12, 13])
        results = o.is_valid_observation(self.test_sat_1, times)
        self.assertEqual(len(results), 3)
        self.assertFalse(results[0])
        self.assertTrue(results[1])
        self.assertFalse(results[2])

    def test_valid_observation_target_sunlit_vector(self):
        o = Instrument(name="Test Instrument", req_target_sunlit=True)
        times = timescale.utc(2020, 3, 20, [11, 12, 13])
        results = o.is_valid_observation(self.test_sat_1, times)
        self.assertEqual(len(results), 3)
        self.assertFalse(results[0])
        self.assertTrue(results[1])
        self.assertFalse(results[2])

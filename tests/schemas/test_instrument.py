import unittest

from tatc.schemas import Instrument
from pydantic import ValidationError
from datetime import datetime, timedelta, timezone
from tatc.constants import timescale
from skyfield.api import wgs84


class TestInstrument(unittest.TestCase):
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
        noon_utc = timescale.utc(2020, 3, 20, 12)
        test_point_1 = wgs84.latlon(0, 0, 400000).at(noon_utc)
        test_point_2 = wgs84.latlon(0, 85, 400000).at(noon_utc)
        test_point_3 = wgs84.latlon(0, 95, 400000).at(noon_utc)
        test_point_4 = wgs84.latlon(0, 180, 400000).at(noon_utc)
        self.assertTrue(o.is_valid_observation(noon_utc, test_point_1))
        self.assertTrue(o.is_valid_observation(noon_utc, test_point_2))
        self.assertTrue(o.is_valid_observation(noon_utc, test_point_3))
        self.assertTrue(o.is_valid_observation(noon_utc, test_point_4))

    def test_valid_observation_self_sunlit(self):
        o = Instrument(name="Test Instrument", req_self_sunlit=True)
        noon_utc = timescale.utc(2020, 3, 20, 12)
        test_point_1 = wgs84.latlon(0, 0, 400000).at(noon_utc)
        test_point_2 = wgs84.latlon(0, 85, 400000).at(noon_utc)
        test_point_3 = wgs84.latlon(0, 95, 400000).at(noon_utc)
        test_point_4 = wgs84.latlon(0, 180, 400000).at(noon_utc)
        self.assertTrue(o.is_valid_observation(noon_utc, test_point_1))
        self.assertTrue(o.is_valid_observation(noon_utc, test_point_2))
        self.assertTrue(o.is_valid_observation(noon_utc, test_point_3))
        self.assertFalse(o.is_valid_observation(noon_utc, test_point_4))

    def test_valid_observation_self_not_sunlit(self):
        o = Instrument(name="Test Instrument", req_self_sunlit=False)
        noon_utc = timescale.utc(2020, 3, 20, 12)
        test_point_1 = wgs84.latlon(0, 0, 400000).at(noon_utc)
        test_point_2 = wgs84.latlon(0, 85, 400000).at(noon_utc)
        test_point_3 = wgs84.latlon(0, 95, 400000).at(noon_utc)
        test_point_4 = wgs84.latlon(0, 180, 400000).at(noon_utc)
        self.assertFalse(o.is_valid_observation(noon_utc, test_point_1))
        self.assertFalse(o.is_valid_observation(noon_utc, test_point_2))
        self.assertFalse(o.is_valid_observation(noon_utc, test_point_3))
        self.assertTrue(o.is_valid_observation(noon_utc, test_point_4))

    def test_valid_observation_target_sunlit(self):
        o = Instrument(name="Test Instrument", req_target_sunlit=True)
        noon_utc = timescale.utc(2022, 3, 20, 12)
        test_point_1 = wgs84.latlon(0, 0, 400000).at(noon_utc)
        test_point_2 = wgs84.latlon(0, 85, 400000).at(noon_utc)
        test_point_3 = wgs84.latlon(0, 95, 400000).at(noon_utc)
        test_point_4 = wgs84.latlon(0, 180, 400000).at(noon_utc)
        self.assertTrue(o.is_valid_observation(noon_utc, test_point_1))
        self.assertTrue(o.is_valid_observation(noon_utc, test_point_2))
        self.assertFalse(o.is_valid_observation(noon_utc, test_point_3))
        self.assertFalse(o.is_valid_observation(noon_utc, test_point_4))

    def test_valid_observation_target_not_sunlit(self):
        o = Instrument(name="Test Instrument", req_target_sunlit=False)
        noon_utc = timescale.utc(2022, 3, 20, 12)
        test_point_1 = wgs84.latlon(0, 0, 400000).at(noon_utc)
        test_point_2 = wgs84.latlon(0, 85, 400000).at(noon_utc)
        test_point_3 = wgs84.latlon(0, 95, 400000).at(noon_utc)
        test_point_4 = wgs84.latlon(0, 180, 400000).at(noon_utc)
        self.assertFalse(o.is_valid_observation(noon_utc, test_point_1))
        self.assertFalse(o.is_valid_observation(noon_utc, test_point_2))
        self.assertTrue(o.is_valid_observation(noon_utc, test_point_3))
        self.assertTrue(o.is_valid_observation(noon_utc, test_point_4))

    def test_valid_observation_self_sunlit_target_sunlit(self):
        o = Instrument(name="Test Instrument", req_self_sunlit=True, req_target_sunlit=True)
        noon_utc = timescale.utc(2022, 3, 20, 12)
        test_point_1 = wgs84.latlon(0, 0, 400000).at(noon_utc)
        test_point_2 = wgs84.latlon(0, 85, 400000).at(noon_utc)
        test_point_3 = wgs84.latlon(0, 95, 400000).at(noon_utc)
        test_point_4 = wgs84.latlon(0, 180, 400000).at(noon_utc)
        self.assertTrue(o.is_valid_observation(noon_utc, test_point_1))
        self.assertTrue(o.is_valid_observation(noon_utc, test_point_2))
        self.assertFalse(o.is_valid_observation(noon_utc, test_point_3))
        self.assertFalse(o.is_valid_observation(noon_utc, test_point_4))

    def test_valid_observation_self_not_sunlit_target_sunlit(self):
        o = Instrument(name="Test Instrument", req_self_sunlit=False, req_target_sunlit=True)
        noon_utc = timescale.utc(2022, 3, 20, 12)
        test_point_1 = wgs84.latlon(0, 0, 400000).at(noon_utc)
        test_point_2 = wgs84.latlon(0, 85, 400000).at(noon_utc)
        test_point_3 = wgs84.latlon(0, 95, 400000).at(noon_utc)
        test_point_4 = wgs84.latlon(0, 180, 400000).at(noon_utc)
        self.assertFalse(o.is_valid_observation(noon_utc, test_point_1))
        self.assertFalse(o.is_valid_observation(noon_utc, test_point_2))
        self.assertFalse(o.is_valid_observation(noon_utc, test_point_3))
        self.assertFalse(o.is_valid_observation(noon_utc, test_point_4))

    def test_valid_observation_self_sunlit_target_not_sunlit(self):
        o = Instrument(name="Test Instrument", req_self_sunlit=True, req_target_sunlit=False)
        noon_utc = timescale.utc(2022, 3, 20, 12)
        test_point_1 = wgs84.latlon(0, 0, 400000).at(noon_utc)
        test_point_2 = wgs84.latlon(0, 85, 400000).at(noon_utc)
        test_point_3 = wgs84.latlon(0, 95, 400000).at(noon_utc)
        test_point_4 = wgs84.latlon(0, 180, 400000).at(noon_utc)
        self.assertFalse(o.is_valid_observation(noon_utc, test_point_1))
        self.assertFalse(o.is_valid_observation(noon_utc, test_point_2))
        self.assertTrue(o.is_valid_observation(noon_utc, test_point_3))
        self.assertFalse(o.is_valid_observation(noon_utc, test_point_4))

    def test_valid_observation_self_not_sunlit_target_not_sunlit(self):
        o = Instrument(name="Test Instrument", req_self_sunlit=False, req_target_sunlit=False)
        noon_utc = timescale.utc(2022, 3, 20, 12)
        test_point_1 = wgs84.latlon(0, 0, 400000).at(noon_utc)
        test_point_2 = wgs84.latlon(0, 85, 400000).at(noon_utc)
        test_point_3 = wgs84.latlon(0, 95, 400000).at(noon_utc)
        test_point_4 = wgs84.latlon(0, 180, 400000).at(noon_utc)
        self.assertFalse(o.is_valid_observation(noon_utc, test_point_1))
        self.assertFalse(o.is_valid_observation(noon_utc, test_point_2))
        self.assertFalse(o.is_valid_observation(noon_utc, test_point_3))
        self.assertTrue(o.is_valid_observation(noon_utc, test_point_4))

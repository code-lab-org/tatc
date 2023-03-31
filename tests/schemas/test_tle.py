import unittest

from pydantic import ValidationError
from datetime import datetime, timezone

from tatc.schemas import TwoLineElements


class TestTLE(unittest.TestCase):
    def setUp(self):
        self.test_tle = TwoLineElements(
            tle=[
                "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
                "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
            ]
        )

    def test_good_data(self):
        good_data = {
            "tle": [
                "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
                "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
            ]
        }
        o = TwoLineElements(**good_data)
        self.assertEqual(o.tle[0], good_data.get("tle")[0])
        self.assertEqual(o.tle[1], good_data.get("tle")[1])

    def test_bad_data_line_1_length(self):
        bad_data = {
            "tle": [
                "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  999",
                "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
            ]
        }
        with self.assertRaises(ValidationError):
            TwoLineElements(**bad_data)

    def test_bad_data_line_2_length(self):
        bad_data = {
            "tle": [
                "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
                "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.4895753428675",
            ]
        }
        with self.assertRaises(ValidationError):
            TwoLineElements(**bad_data)

    def test_bad_data_line_1_format(self):
        bad_data = {
            "tle": [
                "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  99Q3",
                "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
            ]
        }
        with self.assertRaises(ValidationError):
            TwoLineElements(**bad_data)

    def test_bad_data_line_2_format(self):
        bad_data = {
            "tle": [
                "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
                "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.489575342867Q4",
            ]
        }
        with self.assertRaises(ValidationError):
            TwoLineElements(**bad_data)

    def test_bad_data_line_1_checksums(self):
        bad_data = {
            "tle": [
                "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9994",
                "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
            ]
        }
        with self.assertRaises(ValidationError):
            TwoLineElements(**bad_data)

    def test_bad_data_line_2_checksums(self):
        bad_data = {
            "tle": [
                "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
                "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286753",
            ]
        }
        with self.assertRaises(ValidationError):
            TwoLineElements(**bad_data)

    def test_get_catalog_number(self):
        self.assertEqual(self.test_tle.get_catalog_number(), 25544)

    def test_get_classification(self):
        self.assertEqual(self.test_tle.get_classification(), "U")

    def test_get_international_designator(self):
        self.assertEqual(self.test_tle.get_international_designator(), "1998-067A")

    def test_get_epoch(self):
        self.assertEqual(
            self.test_tle.get_epoch(),
            datetime(2021, 6, 5, 7, 19, 36, 128928, tzinfo=timezone.utc),
        )

    def test_get_first_derivative_mean_motion(self):
        self.assertEqual(self.test_tle.get_first_derivative_mean_motion(), 0.00003432)

    def test_get_second_derivative_mean_motion(self):
        self.assertEqual(self.test_tle.get_second_derivative_mean_motion(), 0.0)

    def test_get_b_star(self):
        self.assertEqual(self.test_tle.get_b_star(), 0.000070541)

    def test_get_ephemeris_type(self):
        self.assertEqual(self.test_tle.get_ephemeris_type(), 0)

    def test_get_element_set_number(self):
        self.assertEqual(self.test_tle.get_element_set_number(), 999)

    def test_get_inclination(self):
        self.assertEqual(self.test_tle.get_inclination(), 51.6455)

    def test_get_right_ascension_ascending_node(self):
        self.assertEqual(self.test_tle.get_right_ascension_ascending_node(), 41.4969)

    def test_get_eccentricity(self):
        self.assertEqual(self.test_tle.get_eccentricity(), 0.0003508)

    def test_get_perigee_argument(self):
        self.assertEqual(self.test_tle.get_perigee_argument(), 68.0432)

    def test_get_mean_anomaly(self):
        self.assertEqual(self.test_tle.get_mean_anomaly(), 78.3395)

    def test_get_mean_motion(self):
        self.assertAlmostEqual(self.test_tle.get_mean_motion(), 15.48957534)

    def test_get_revolution_number_at_epoch(self):
        self.assertEqual(self.test_tle.get_revolution_number_at_epoch(), 28675)

    def test_get_semimajor_axis(self):
        self.assertAlmostEqual(self.test_tle.get_semimajor_axis(), 6797911, delta=1.0)

    def test_get_altitude(self):
        self.assertAlmostEqual(self.test_tle.get_altitude(), 426902, delta=1.0)

    def test_get_true_anomaly(self):
        self.assertAlmostEqual(self.test_tle.get_true_anomaly(), 78.3788725993742)

    def test_get_derived_orbit(self):
        derived_tle = self.test_tle.get_derived_orbit(20, 10)
        self.assertAlmostEqual(
            derived_tle.get_mean_anomaly(),
            self.test_tle.get_mean_anomaly() + 20,
            delta=0.001,
        )
        self.assertAlmostEqual(
            derived_tle.get_right_ascension_ascending_node(),
            self.test_tle.get_right_ascension_ascending_node() + 10,
            delta=0.001,
        )

    def test_get_tle(self):
        self.assertAlmostEqual(self.test_tle.to_tle(), self.test_tle)

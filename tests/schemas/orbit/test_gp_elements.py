"""
Unit tests for the GeneralPerturbationsElements schema.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import csv
import io
import json
import unittest
from datetime import datetime, timezone

from pydantic import ValidationError
from skyfield.api import EarthSatellite

from tatc.schemas.orbit.gp import GeneralPerturbationsElements


class TestGeneralPerturbationsElements(unittest.TestCase):
    """
    Unit tests for the GeneralPerturbationsElements schema.
    """

    def setUp(self):
        self.test_tle = (
            "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
            "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
        )
        self.test_elements = GeneralPerturbationsElements.from_tle(self.test_tle)

    def test_good_data(self):
        """
        Test that the GeneralPerturbationsElements schema correctly
        initializes with valid data.
        """
        good_data = {
            "epoch": datetime(2022, 1, 1, tzinfo=timezone.utc),
            "mean_motion": 0.06,
            "eccentricity": 0.001,
            "inclination": 51.6,
            "ra_of_asc_node": 40,
            "arg_of_pericenter": 68,
            "mean_anomaly": 78,
        }
        e = GeneralPerturbationsElements(**good_data)
        for field, value in good_data.items():
            self.assertEqual(getattr(e, field), value)

    def test_defaults(self):
        """
        Test that all optional fields default correctly when omitted.
        """
        e = GeneralPerturbationsElements(
            epoch=datetime(2022, 1, 1, tzinfo=timezone.utc),
            mean_motion=0.06,
            eccentricity=0.001,
            inclination=51.6,
            ra_of_asc_node=40,
            arg_of_pericenter=68,
            mean_anomaly=78,
        )
        self.assertIsNone(e.object_name)
        self.assertEqual(e.norad_cat_id, 0)
        self.assertEqual(e.bstar, 0)
        self.assertEqual(e.mean_motion_dot, 0)
        self.assertEqual(e.mean_motion_ddot, 0)
        self.assertEqual(e.classification, "U")
        self.assertEqual(e.international_designator, "00000A")
        self.assertEqual(e.ephemeris_type, 0)
        self.assertEqual(e.element_set_num, 0)
        self.assertEqual(e.revolution_num, 0)

    def test_bad_required_field_missing(self):
        """
        Test that omitting any required field raises a ValidationError.
        """
        required = {
            "epoch": datetime(2022, 1, 1, tzinfo=timezone.utc),
            "mean_motion": 0.06,
            "eccentricity": 0.001,
            "inclination": 51.6,
            "ra_of_asc_node": 40,
            "arg_of_pericenter": 68,
            "mean_anomaly": 78,
        }
        for field in required:
            with self.subTest(field=field):
                data = {k: v for k, v in required.items() if k != field}
                with self.assertRaises(ValidationError):
                    GeneralPerturbationsElements(**data)

    def test_bad_mean_motion_non_positive(self):
        """
        Test that a zero or negative mean_motion is rejected.
        """
        with self.assertRaises(ValidationError):
            GeneralPerturbationsElements(
                epoch=datetime(2022, 1, 1, tzinfo=timezone.utc),
                mean_motion=0,
                eccentricity=0,
                inclination=0,
                ra_of_asc_node=0,
                arg_of_pericenter=0,
                mean_anomaly=0,
            )

    def test_bad_eccentricity_out_of_range(self):
        """
        Test that an eccentricity outside [0, 1] is rejected.
        """
        base = dict(
            epoch=datetime(2022, 1, 1, tzinfo=timezone.utc),
            mean_motion=1,
            inclination=0,
            ra_of_asc_node=0,
            arg_of_pericenter=0,
            mean_anomaly=0,
        )
        with self.assertRaises(ValidationError):
            GeneralPerturbationsElements(eccentricity=-0.1, **base)
        with self.assertRaises(ValidationError):
            GeneralPerturbationsElements(eccentricity=1.1, **base)

    def test_bad_angle_fields_out_of_range(self):
        """
        Test that inclination, ra_of_asc_node, arg_of_pericenter, and
        mean_anomaly each reject values outside their documented ranges.
        """
        base = dict(
            epoch=datetime(2022, 1, 1, tzinfo=timezone.utc),
            mean_motion=1,
            eccentricity=0,
            inclination=0,
            ra_of_asc_node=0,
            arg_of_pericenter=0,
            mean_anomaly=0,
        )
        with self.assertRaises(ValidationError):
            GeneralPerturbationsElements(**{**base, "inclination": 180.1})
        with self.assertRaises(ValidationError):
            GeneralPerturbationsElements(**{**base, "inclination": -0.1})
        with self.assertRaises(ValidationError):
            GeneralPerturbationsElements(**{**base, "ra_of_asc_node": 360})
        with self.assertRaises(ValidationError):
            GeneralPerturbationsElements(**{**base, "arg_of_pericenter": 360})
        with self.assertRaises(ValidationError):
            GeneralPerturbationsElements(**{**base, "mean_anomaly": 360})

    def test_get_orbit_period_matches_published_iss_period(self):
        """
        Test get_orbit_period against the ISS's well-known published
        orbital period of approximately 93 minutes.
        """
        self.assertAlmostEqual(
            self.test_elements.get_orbit_period().total_seconds() / 60,
            93,
            delta=1.0,
        )

    def test_get_semimajor_axis(self):
        """
        Test get_semimajor_axis against the ISS's well-known published
        semimajor axis of approximately 6,798 km.
        """
        self.assertAlmostEqual(
            self.test_elements.get_semimajor_axis(), 6797911, delta=1.0
        )

    def test_get_mean_altitude(self):
        """
        Test get_mean_altitude against the ISS's well-known published
        mean altitude of approximately 427 km.
        """
        self.assertAlmostEqual(
            self.test_elements.get_mean_altitude(), 426902, delta=1.0
        )

    def test_get_true_anomaly(self):
        """
        Test that get_true_anomaly converts mean anomaly to true anomaly
        using the class's own eccentricity.
        """
        self.assertAlmostEqual(
            self.test_elements.get_true_anomaly(), 78.3788725993742, delta=1e-6
        )

    def test_from_tle(self):
        """
        Test that from_tle correctly parses all fields from a real TLE.
        """
        e = self.test_elements
        self.assertEqual(e.norad_cat_id, 25544)
        self.assertEqual(
            e.epoch, datetime(2021, 6, 5, 7, 19, 36, 128928, tzinfo=timezone.utc)
        )
        self.assertEqual(e.inclination, 51.6455)
        self.assertEqual(e.ra_of_asc_node, 41.4969)
        self.assertEqual(e.eccentricity, 0.0003508)
        self.assertEqual(e.arg_of_pericenter, 68.0432)
        self.assertEqual(e.mean_anomaly, 78.3395)
        self.assertAlmostEqual(e.mean_motion * 86400 / 360, 15.48957534, delta=1e-9)
        self.assertEqual(e.bstar, 0.000070541)
        self.assertEqual(e.international_designator, "98067A")
        self.assertEqual(e.element_set_num, 999)
        self.assertEqual(e.revolution_num, 28675)

    def test_to_tle_round_trip(self):
        """
        Test that converting to TLE and back reproduces the original TLE
        lines exactly.
        """
        tle_out = self.test_elements.to_tle()
        self.assertEqual(tle_out[0], self.test_tle[0])
        self.assertEqual(tle_out[1], self.test_tle[1])

    def test_from_satrec_to_satrec_round_trip(self):
        """
        Test that converting to a Satrec and back preserves every field
        (to within floating-point precision from the radian/degree and
        per-minute/per-second unit conversions).
        """
        satrec = self.test_elements.to_satrec()
        round_tripped = GeneralPerturbationsElements.from_satrec(satrec)
        for field in type(self.test_elements).model_fields:
            with self.subTest(field=field):
                original = getattr(self.test_elements, field)
                restored = getattr(round_tripped, field)
                if isinstance(original, float):
                    self.assertAlmostEqual(original, restored, delta=1e-9)
                else:
                    self.assertEqual(original, restored)

    def test_to_omm_dict_contents(self):
        """
        Test that to_omm_dict produces the expected OMM field names and
        values, including the object name and mean motion in the OMM
        convention of revolutions/day (rather than this class's own
        degrees/second).
        """
        elements = self.test_elements.model_copy(update={"object_name": "ISS (ZARYA)"})
        omm_dict = elements.to_omm_dict()
        self.assertEqual(omm_dict["OBJECT_NAME"], "ISS (ZARYA)")
        self.assertEqual(omm_dict["NORAD_CAT_ID"], 25544)
        self.assertEqual(omm_dict["INCLINATION"], 51.6455)
        self.assertAlmostEqual(omm_dict["MEAN_MOTION"], 15.48957534286754, delta=1e-6)

    def test_from_omm_dict_round_trip_including_object_name(self):
        """
        Regression test: from_omm_dict must preserve object_name, not
        just the orbital elements. object_name has no equivalent on the
        Satrec object from_omm_dict constructs internally, so it must be
        restored directly from the OMM dictionary rather than being lost.
        """
        elements = self.test_elements.model_copy(update={"object_name": "ISS (ZARYA)"})
        omm_dict = elements.to_omm_dict()
        round_tripped = GeneralPerturbationsElements.from_omm_dict(omm_dict)
        self.assertEqual(round_tripped.object_name, "ISS (ZARYA)")
        self.assertEqual(round_tripped.norad_cat_id, elements.norad_cat_id)
        self.assertEqual(round_tripped.inclination, elements.inclination)

    def test_from_omm_dict_object_name_none_when_absent(self):
        """
        Test that from_omm_dict leaves object_name as None when the OMM
        dictionary has no OBJECT_NAME key, rather than raising an error.
        """
        omm_dict = self.test_elements.to_omm_dict()
        del omm_dict["OBJECT_NAME"]
        round_tripped = GeneralPerturbationsElements.from_omm_dict(omm_dict)
        self.assertIsNone(round_tripped.object_name)

    def test_from_omm_csv_round_trip_including_object_name(self):
        """
        Test that from_omm_csv preserves object_name through a CSV
        round trip, via the same from_omm_dict path.
        """
        elements = self.test_elements.model_copy(update={"object_name": "ISS (ZARYA)"})
        omm_dict = elements.to_omm_dict()
        buf = io.StringIO()
        writer = csv.DictWriter(buf, fieldnames=list(omm_dict.keys()))
        writer.writeheader()
        writer.writerow(omm_dict)
        round_tripped = GeneralPerturbationsElements.from_omm_csv(
            buf.getvalue().splitlines()
        )
        self.assertEqual(round_tripped.object_name, "ISS (ZARYA)")
        self.assertEqual(round_tripped.norad_cat_id, elements.norad_cat_id)

    def test_from_omm_csv_only_uses_first_row(self):
        """
        Test that from_omm_csv uses only the first data row when given
        multiple rows, since it returns a single GeneralPerturbationsElements
        rather than a list.
        """
        omm_dict = self.test_elements.to_omm_dict()
        other_dict = dict(omm_dict, NORAD_CAT_ID=99999)
        buf = io.StringIO()
        writer = csv.DictWriter(buf, fieldnames=list(omm_dict.keys()))
        writer.writeheader()
        writer.writerow(omm_dict)
        writer.writerow(other_dict)
        round_tripped = GeneralPerturbationsElements.from_omm_csv(
            buf.getvalue().splitlines()
        )
        self.assertEqual(round_tripped.norad_cat_id, 25544)

    def test_from_omm_csv_empty_raises_value_error(self):
        """
        Test that from_omm_csv raises a ValueError when given no lines.
        """
        with self.assertRaises(ValueError):
            GeneralPerturbationsElements.from_omm_csv([])

    def test_from_omm_json_round_trip_including_object_name(self):
        """
        Test that from_omm_json preserves object_name through a JSON
        round trip, via the same from_omm_dict path.
        """
        elements = self.test_elements.model_copy(update={"object_name": "ISS (ZARYA)"})
        omm_json = json.dumps([elements.to_omm_dict()])
        round_tripped = GeneralPerturbationsElements.from_omm_json(omm_json)
        self.assertEqual(round_tripped.object_name, "ISS (ZARYA)")
        self.assertEqual(round_tripped.norad_cat_id, elements.norad_cat_id)

    def test_from_omm_json_only_uses_first_entry(self):
        """
        Test that from_omm_json uses only the first entry when given
        multiple, since it returns a single GeneralPerturbationsElements
        rather than a list.
        """
        omm_dict = self.test_elements.to_omm_dict()
        other_dict = dict(omm_dict, NORAD_CAT_ID=99999)
        omm_json = json.dumps([omm_dict, other_dict])
        round_tripped = GeneralPerturbationsElements.from_omm_json(omm_json)
        self.assertEqual(round_tripped.norad_cat_id, 25544)

    def test_from_omm_json_empty_raises_value_error(self):
        """
        Test that from_omm_json raises a ValueError when given an empty
        JSON array.
        """
        with self.assertRaises(ValueError):
            GeneralPerturbationsElements.from_omm_json("[]")

    def test_to_skyfield_returns_configured_earth_satellite(self):
        """
        Test that to_skyfield returns a Skyfield EarthSatellite whose
        catalog number and epoch match this element set, confirming it
        is usable to propagate this orbital state via SGP4.
        """
        satellite = self.test_elements.to_skyfield()
        self.assertIsInstance(satellite, EarthSatellite)
        self.assertEqual(satellite.model.satnum, 25544)
        self.assertAlmostEqual(
            satellite.epoch.utc_datetime().timestamp(),
            self.test_elements.epoch.timestamp(),
            delta=1e-3,
        )


if __name__ == "__main__":
    unittest.main()

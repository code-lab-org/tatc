"""
Unit tests for the GeneralPerturbationsElements schema.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import csv
import io
import json
import unittest
from datetime import datetime, timedelta, timezone

from pydantic import ValidationError
from skyfield.api import EarthSatellite

from tatc import constants
from tatc.schemas.orbit.gp_elements import GeneralPerturbationsElements


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


class TestGetRepeatCycle(unittest.TestCase):
    """
    Unit tests for GeneralPerturbationsElements.get_repeat_cycle. This is
    the single-element repeat-cycle computation itself; see
    tests/schemas/orbit/test_gp.py's TestGetRepeatCycle for
    GeneralPerturbationsOrbit's multi-element consistency-checking
    aggregation on top of this.
    """

    def setUp(self):
        # a real Landsat-8 element set (NORAD 39084, fetched from
        # Celestrak), whose published repeat ground track is exactly 233
        # orbits every 16 days (USGS)
        self.landsat_8_tle = (
            "1 39084U 13008A   26213.27824675  .00000294  00000+0  75333-4 0  9990",
            "2 39084  98.2277 282.8718 0001275  92.4910 267.6434 14.57104473704466",
        )
        # a real Sentinel-2A element set (NORAD 40697), whose published
        # repeat ground track is exactly 143 orbits every 10 days (ESA)
        self.sentinel_2a_tle = (
            "1 40697U 15028A   26213.24967738  .00000092  00000+0  51774-4 0  9990",
            "2 40697  98.5671 287.4570 0001329  92.6662 267.4673 14.30818788580177",
        )
        # a real Sentinel-1A element set (NORAD 39634) from early 2026,
        # while the satellite was still actively operated (its mission
        # concluded 2026-06-29): published repeat ground track is exactly
        # 175 orbits every 12 days (ESA)
        self.sentinel_1a_tle = (
            "1 39634U 14016A   26001.19041520  .00000521  00000-0  12012-3 0  9995",
            "2 39634  98.1805  11.1453 0001276  85.2963 274.8383 14.59199668625660",
        )
        # a real GPS BIIR-5 element set (NORAD 26407): GPS orbits are
        # designed to repeat their ground track every sidereal day (two
        # ~12-hour orbits per day), a much shorter cycle than the
        # sun-synchronous imaging orbits above, at MEO altitude (~20,200 km)
        self.gps_tle = (
            "1 26407U 00040A   26213.32131914  .00000072  00000+0  00000+0 0  9998",
            "2 26407  54.8467 213.3697 0120005 302.9740 169.1529  2.00558010190856",
        )
        # a real Molniya 3-8 element set (NORAD 10455): a highly eccentric,
        # critical-inclination (~63.4 degree) orbit that, like GPS, repeats
        # every sidereal day by design, but is geometrically nothing like
        # the near-circular orbits above
        self.molniya_tle = (
            "1 10455U 77105A   26212.97315042  .00000728  00000+0  00000+0 0  9994",
            "2 10455  63.8024 172.4301 6701249 276.1687  17.0849  2.00778403357346",
        )
        # the ISS is not designed for a repeat ground track
        self.iss_tle = (
            "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
            "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
        )

    def test_landsat_8_repeat_cycle_matches_published_16_days(self):
        """
        Test against Landsat-8's published repeat ground track of 233
        orbits every 16 days (USGS). A small tolerance accounts for the
        real orbit's minor drift between station-keeping maneuvers.
        """
        elements = GeneralPerturbationsElements.from_tle(self.landsat_8_tle)
        repeat_cycle = elements.get_repeat_cycle()
        self.assertIsNotNone(repeat_cycle)
        self.assertAlmostEqual(repeat_cycle.total_seconds() / 86400, 16, delta=0.1)

    def test_sentinel_2a_repeat_cycle_matches_published_10_days(self):
        """
        Test against Sentinel-2A's published repeat ground track of 143
        orbits every 10 days (ESA), at a different altitude/inclination
        than Landsat-8, using the default tolerances.
        """
        elements = GeneralPerturbationsElements.from_tle(self.sentinel_2a_tle)
        repeat_cycle = elements.get_repeat_cycle()
        self.assertIsNotNone(repeat_cycle)
        self.assertAlmostEqual(repeat_cycle.total_seconds() / 86400, 10, delta=0.1)

    def test_sentinel_1a_repeat_cycle_matches_published_12_days(self):
        """
        Test against Sentinel-1A's published repeat ground track of 175
        orbits every 12 days (ESA), at a different altitude/inclination
        than Landsat-8, using the default tolerances.
        """
        elements = GeneralPerturbationsElements.from_tle(self.sentinel_1a_tle)
        repeat_cycle = elements.get_repeat_cycle()
        self.assertIsNotNone(repeat_cycle)
        self.assertAlmostEqual(repeat_cycle.total_seconds() / 86400, 12, delta=0.1)

    def test_gps_repeat_cycle_matches_sidereal_day(self):
        """
        Test against GPS's designed repeat ground track of one sidereal
        day (~23h56m), using the default tolerances. Confirms the search
        also correctly resolves a very short repeat cycle (day 1), not
        just the multi-week cycles above, and at a completely different
        (MEO) altitude regime.
        """
        elements = GeneralPerturbationsElements.from_tle(self.gps_tle)
        repeat_cycle = elements.get_repeat_cycle()
        self.assertIsNotNone(repeat_cycle)
        self.assertAlmostEqual(
            repeat_cycle.total_seconds() / 86400,
            constants.EARTH_SIDEREAL_DAY_S / 86400,
            delta=0.01,
        )

    def test_molniya_repeat_cycle_matches_sidereal_day(self):
        """
        Test against a Molniya orbit's designed repeat ground track of
        one sidereal day, at a highly eccentric, critical-inclination
        geometry very different from every other case here. Public
        Molniya TLEs are published with an epoch near perigee (where
        ground radar tracking is easiest), and perigee is where a highly
        eccentric orbit moves fastest (~10 km/s here, vs ~1.5 km/s at
        apogee) -- so a given timing imprecision in the analytic
        prediction translates into a much larger position/velocity error
        than for the near-circular orbits above. This is a property of
        the epoch's orbital phase, not of the orbit's true repeatability,
        so this test widens the tolerance rather than the global default.
        """
        elements = GeneralPerturbationsElements.from_tle(self.molniya_tle)
        repeat_cycle = elements.get_repeat_cycle(
            max_delta_position=80000, max_delta_velocity=25
        )
        self.assertIsNotNone(repeat_cycle)
        self.assertAlmostEqual(
            repeat_cycle.total_seconds() / 86400,
            constants.EARTH_SIDEREAL_DAY_S / 86400,
            delta=0.01,
        )

    def test_non_repeating_orbit_returns_none(self):
        """
        Test that the ISS's orbit, which is not designed for a repeat
        ground track, finds no repeat at the default tolerances (10 km,
        3 m/s). The ISS does have a coincidental ~4-day near-repeat
        (delta position 10.6 km, delta velocity 8.0 m/s) -- narrowly
        outside the position tolerance, and well outside the velocity
        tolerance, so both criteria independently reject it. A looser
        position tolerance alone (e.g. 30 km, to accommodate a less
        precisely maintained real repeat orbit) would accept this
        incidental match on position, which is exactly why velocity is
        checked too rather than relying on position alone.
        """
        elements = GeneralPerturbationsElements.from_tle(self.iss_tle)
        self.assertIsNone(elements.get_repeat_cycle())

    def test_lazy_load_reuses_cached_result(self):
        """
        Test that calling get_repeat_cycle() twice with lazy_load=True
        (the default) returns the identical cached timedelta rather than
        recomputing.
        """
        elements = GeneralPerturbationsElements.from_tle(self.landsat_8_tle)
        first = elements.get_repeat_cycle()
        second = elements.get_repeat_cycle()
        self.assertIs(first, second)

    def test_lazy_load_false_forces_recomputation(self):
        """
        Test that lazy_load=False recomputes rather than reusing the
        cached result (a fresh but equal timedelta).
        """
        elements = GeneralPerturbationsElements.from_tle(self.landsat_8_tle)
        first = elements.get_repeat_cycle()
        second = elements.get_repeat_cycle(lazy_load=False)
        self.assertEqual(first, second)
        self.assertIsNot(first, second)

    def test_too_short_search_duration_returns_none(self):
        """
        Test that a max_search_duration shorter than the true repeat
        cycle (16 days) cannot find it.
        """
        elements = GeneralPerturbationsElements.from_tle(self.landsat_8_tle)
        repeat_cycle = elements.get_repeat_cycle(
            max_search_duration=timedelta(days=10), lazy_load=False
        )
        self.assertIsNone(repeat_cycle)

    def test_too_tight_tolerance_returns_none(self):
        """
        Test that an unrealistically tight position/velocity tolerance
        rejects even the real repeat cycle.
        """
        elements = GeneralPerturbationsElements.from_tle(self.landsat_8_tle)
        repeat_cycle = elements.get_repeat_cycle(
            max_delta_position=1, max_delta_velocity=0.001, lazy_load=False
        )
        self.assertIsNone(repeat_cycle)


if __name__ == "__main__":
    unittest.main()

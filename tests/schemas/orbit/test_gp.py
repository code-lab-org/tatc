"""
Unit tests for the GeneralPerturbationsOrbit schema.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import csv
import io
import json
import unittest
from datetime import datetime, timedelta, timezone

import numpy as np
from pydantic import ValidationError
from skyfield.api import wgs84

from tatc import config, constants
from tatc.schemas import GeneralPerturbationsOrbit, Point


class TestGPOrbit(unittest.TestCase):
    """
    Unit tests for the GeneralPerturbationsOrbit schema.
    """

    def setUp(self):
        self.test_tle = [
            "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
            "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
        ]
        self.test_orbit = GeneralPerturbationsOrbit.from_tle(self.test_tle)
        # a clearly-distinguishable second element, for testing index-based
        # access on a multi-element orbit
        self.second_element = self.test_orbit.elements[0].model_copy(
            update={
                "norad_cat_id": 99999,
                "epoch": datetime(2021, 6, 6, 7, 19, 36, tzinfo=timezone.utc),
                "inclination": 45.0,
                "eccentricity": 0.001,
                "ra_of_asc_node": 100.0,
                "arg_of_pericenter": 10.0,
                "mean_anomaly": 20.0,
                "bstar": 0.0001,
                "mean_motion_dot": 0.00001,
                "mean_motion_ddot": 0.000001,
            }
        )
        self.multi_element_orbit = GeneralPerturbationsOrbit(
            elements=[self.test_orbit.elements[0], self.second_element]
        )

    def test_bad_elements_empty(self):
        """
        Test that an empty elements list is rejected, since every getter
        assumes at least one element exists.
        """
        with self.assertRaises(ValidationError):
            GeneralPerturbationsOrbit(elements=[])

    def test_bad_elements_missing(self):
        """
        Test that the required elements field must be provided.
        """
        with self.assertRaises(ValidationError):
            GeneralPerturbationsOrbit()

    def test_elements_sorted_by_epoch_on_construction(self):
        """
        Regression test: elements provided out of epoch order must be
        sorted ascending by epoch at construction time, since
        get_closest_element_index relies on np.searchsorted, which
        silently produces incorrect results if its input is not sorted.
        """
        base = self.test_orbit.elements[0]
        e_jan = base.model_copy(
            update={"epoch": datetime(2022, 1, 1, tzinfo=timezone.utc)}
        )
        e_mar = base.model_copy(
            update={"epoch": datetime(2022, 3, 1, tzinfo=timezone.utc)}
        )
        e_feb = base.model_copy(
            update={"epoch": datetime(2022, 2, 1, tzinfo=timezone.utc)}
        )
        o = GeneralPerturbationsOrbit(elements=[e_jan, e_mar, e_feb])
        self.assertEqual(
            [el.epoch for el in o.elements], [e_jan.epoch, e_feb.epoch, e_mar.epoch]
        )

    def test_get_catalog_number(self):
        """
        Test that the catalog number can be retrieved from the GeneralPerturbationsOrbit object.
        """
        self.assertEqual(self.test_orbit.get_catalog_number(), 25544)

    def test_get_epoch(self):
        """
        Test that the epoch can be retrieved from the GeneralPerturbationsOrbit object.
        """
        self.assertEqual(
            self.test_orbit.get_epoch(),
            datetime(2021, 6, 5, 7, 19, 36, 128928, tzinfo=timezone.utc),
        )

    def test_get_mean_motion_dot(self):
        """
        Test that the first derivative of the mean motion can be retrieved from the GeneralPerturbationsOrbit object.
        """
        self.assertEqual(
            self.test_orbit.get_mean_motion_dot() * (24 * 60 * 60) ** 2 / 360,
            0.00003432,
        )

    def test_get_mean_motion_ddot(self):
        """
        Test that the second derivative of the mean motion can be retrieved from the GeneralPerturbationsOrbit object.
        """
        self.assertEqual(
            self.test_orbit.get_mean_motion_ddot() * (24 * 60 * 60) ** 3 / 360, 0.0
        )

    def test_get_bstar(self):
        """
        Test that the B* drag term can be retrieved from the GeneralPerturbationsOrbit object.
        """
        self.assertEqual(self.test_orbit.get_bstar(), 0.000070541)

    def test_get_inclination(self):
        """
        Test that the inclination can be retrieved from the GeneralPerturbationsOrbit object.
        """
        self.assertEqual(self.test_orbit.get_inclination(), 51.6455)

    def test_get_right_ascension_ascending_node(self):
        """
        Test that the right ascension of the ascending node can be retrieved from the GeneralPerturbationsOrbit object.
        """
        self.assertEqual(self.test_orbit.get_right_ascension_ascending_node(), 41.4969)

    def test_get_eccentricity(self):
        """
        Test that the eccentricity can be retrieved from the GeneralPerturbationsOrbit object.
        """
        self.assertEqual(self.test_orbit.get_eccentricity(), 0.0003508)

    def test_get_perigee_argument(self):
        """
        Test that the argument of perigee can be retrieved from the GeneralPerturbationsOrbit object.
        """
        self.assertEqual(self.test_orbit.get_perigee_argument(), 68.0432)

    def test_get_mean_anomaly(self):
        """
        Test that the mean anomaly can be retrieved from the GeneralPerturbationsOrbit object.
        """
        self.assertEqual(self.test_orbit.get_mean_anomaly(), 78.3395)

    def test_get_mean_motion(self):
        """
        Test that the mean motion can be retrieved from the GeneralPerturbationsOrbit object.
        """
        self.assertAlmostEqual(
            self.test_orbit.get_mean_motion() * (24 * 60 * 60) / 360, 15.48957534
        )

    def test_get_orbit_period(self):
        """
        Test that the orbit period can be retrieved from the GeneralPerturbationsOrbit object
        using the ISS's well-known published orbital period of approximately 93 minutes.
        """
        self.assertAlmostEqual(
            self.test_orbit.get_orbit_period().total_seconds() / 60, 93, delta=1.0
        )

    def test_get_semimajor_axis(self):
        """
        Test that the semi-major axis can be retrieved from the GeneralPerturbationsOrbit object.
        """
        self.assertAlmostEqual(self.test_orbit.get_semimajor_axis(), 6797911, delta=1.0)

    def test_get_mean_altitude(self):
        """
        Test that the mean altitude can be retrieved from the GeneralPerturbationsOrbit object.
        """
        self.assertAlmostEqual(self.test_orbit.get_mean_altitude(), 426902, delta=1.0)

    def test_get_true_anomaly(self):
        """
        Test that the true anomaly can be retrieved from the GeneralPerturbationsOrbit object.
        """
        self.assertAlmostEqual(
            self.test_orbit.elements[0].get_true_anomaly(), 78.3788725993742
        )

    def test_getters_default_to_first_element(self):
        """
        Test that every index-aware getter defaults to index=0 (the first
        element), matching the design intent that a single-element orbit
        should be usable without ever passing an explicit index.
        """
        self.assertEqual(
            self.multi_element_orbit.get_catalog_number(),
            self.multi_element_orbit.get_catalog_number(0),
        )
        self.assertEqual(
            self.multi_element_orbit.get_inclination(),
            self.multi_element_orbit.get_inclination(0),
        )

    def test_getters_select_specified_element_by_index(self):
        """
        Test that passing index=1 on a multi-element orbit retrieves
        values from the second element, not the first -- the core
        multi-element design intent of this class.
        """
        o = self.multi_element_orbit
        self.assertEqual(o.get_catalog_number(1), 99999)
        self.assertEqual(o.get_inclination(1), 45.0)
        self.assertEqual(o.get_eccentricity(1), 0.001)
        self.assertEqual(o.get_right_ascension_ascending_node(1), 100.0)
        self.assertEqual(o.get_perigee_argument(1), 10.0)
        self.assertEqual(o.get_mean_anomaly(1), 20.0)
        self.assertEqual(o.get_bstar(1), 0.0001)
        self.assertEqual(o.get_mean_motion_dot(1), 0.00001)
        self.assertEqual(o.get_mean_motion_ddot(1), 0.000001)
        self.assertEqual(
            o.get_epoch(1), datetime(2021, 6, 6, 7, 19, 36, tzinfo=timezone.utc)
        )
        # first element (index 0) is unaffected
        self.assertEqual(o.get_catalog_number(0), 25544)

    def test_from_tle_multiple_elements(self):
        """
        Test that from_tle builds one element per TLE pair when given a
        flat list of multiple concatenated TLEs.
        """
        combined_tle_lines = self.test_tle + self.test_tle
        o = GeneralPerturbationsOrbit.from_tle(combined_tle_lines)
        self.assertEqual(len(o.elements), 2)

    def test_from_tle_odd_number_of_lines_raises_value_error(self):
        """
        Test that from_tle raises a clear ValueError (rather than an
        unhelpful IndexError) when given an odd number of TLE lines.
        """
        with self.assertRaises(ValueError):
            GeneralPerturbationsOrbit.from_tle(self.test_tle + [self.test_tle[0]])

    def test_from_omm_csv_multiple_elements(self):
        """
        Test that from_omm_csv builds one element per CSV row, unlike
        GeneralPerturbationsElements.from_omm_csv, which only uses the
        first row.
        """
        omm_dict_1 = self.test_orbit.elements[0].to_omm_dict()
        omm_dict_2 = self.second_element.to_omm_dict()
        buf = io.StringIO()
        writer = csv.DictWriter(buf, fieldnames=list(omm_dict_1.keys()))
        writer.writeheader()
        writer.writerow(omm_dict_1)
        writer.writerow(omm_dict_2)
        o = GeneralPerturbationsOrbit.from_omm_csv(buf.getvalue().splitlines())
        self.assertEqual(len(o.elements), 2)
        self.assertEqual(o.get_catalog_number(0), 25544)
        self.assertEqual(o.get_catalog_number(1), 99999)

    def test_from_omm_csv_empty_rejected(self):
        """
        Test that from_omm_csv with no data rows is rejected via the
        elements field's min_length constraint (rather than silently
        building a zero-element orbit).
        """
        with self.assertRaises(ValidationError):
            GeneralPerturbationsOrbit.from_omm_csv([])

    def test_from_omm_json_multiple_elements(self):
        """
        Test that from_omm_json builds one element per JSON entry, unlike
        GeneralPerturbationsElements.from_omm_json, which only uses the
        first entry.
        """
        omm_json = json.dumps(
            [
                self.test_orbit.elements[0].to_omm_dict(),
                self.second_element.to_omm_dict(),
            ]
        )
        o = GeneralPerturbationsOrbit.from_omm_json(omm_json)
        self.assertEqual(len(o.elements), 2)
        self.assertEqual(o.get_catalog_number(0), 25544)
        self.assertEqual(o.get_catalog_number(1), 99999)

    def test_from_omm_json_empty_rejected(self):
        """
        Test that from_omm_json with an empty JSON array is rejected via
        the elements field's min_length constraint.
        """
        with self.assertRaises(ValidationError):
            GeneralPerturbationsOrbit.from_omm_json("[]")

    def test_get_element_epochs(self):
        """
        Test that get_element_epochs returns the epoch of every element,
        in element order.
        """
        epochs = self.multi_element_orbit.get_element_epochs()
        self.assertEqual(len(epochs), 2)
        self.assertEqual(epochs[0], self.test_orbit.elements[0].epoch)
        self.assertEqual(epochs[1], self.second_element.epoch)

    def test_get_derived_orbit(self):
        """
        Test that a derived orbit can be created from the GeneralPerturbationsOrbit object.
        """
        derived_orbit = self.test_orbit.get_derived_orbit(20, 10)
        self.assertAlmostEqual(
            derived_orbit.get_mean_anomaly(),
            self.test_orbit.get_mean_anomaly() + 20,
            delta=0.001,
        )
        self.assertAlmostEqual(
            derived_orbit.get_right_ascension_ascending_node(),
            self.test_orbit.get_right_ascension_ascending_node() + 10,
            delta=0.001,
        )

    def test_get_derived_orbit_shifts_every_element(self):
        """
        Test that get_derived_orbit applies the same mean anomaly and
        RAAN perturbation to every element in a multi-element orbit, not
        just the first.
        """
        derived_orbit = self.multi_element_orbit.get_derived_orbit(20, 10)
        self.assertEqual(len(derived_orbit.elements), 2)
        for i in range(2):
            with self.subTest(index=i):
                self.assertAlmostEqual(
                    derived_orbit.get_mean_anomaly(i),
                    self.multi_element_orbit.get_mean_anomaly(i) + 20,
                    delta=0.001,
                )
                self.assertAlmostEqual(
                    derived_orbit.get_right_ascension_ascending_node(i),
                    self.multi_element_orbit.get_right_ascension_ascending_node(i) + 10,
                    delta=0.001,
                )

    def test_get_derived_orbit_does_not_mutate_original(self):
        """
        Test that get_derived_orbit does not mutate the original orbit's
        elements (since it deep-copies each element before modifying it).
        """
        original_mean_anomaly = self.test_orbit.get_mean_anomaly()
        self.test_orbit.get_derived_orbit(20, 10)
        self.assertEqual(self.test_orbit.get_mean_anomaly(), original_mean_anomaly)


class TestGetClosestElementIndex(unittest.TestCase):
    """
    Unit tests for GeneralPerturbationsOrbit.get_closest_element_index.
    """

    def setUp(self):
        tle = [
            "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
            "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
        ]
        base = GeneralPerturbationsOrbit.from_tle(tle).elements[0]
        self.epoch_0 = datetime(2022, 1, 1, tzinfo=timezone.utc)
        self.epoch_1 = datetime(2022, 1, 3, tzinfo=timezone.utc)
        self.epoch_2 = datetime(2022, 1, 5, tzinfo=timezone.utc)
        self.multi_element_orbit = GeneralPerturbationsOrbit(
            elements=[
                base.model_copy(update={"epoch": self.epoch_0}),
                base.model_copy(update={"epoch": self.epoch_1}),
                base.model_copy(update={"epoch": self.epoch_2}),
            ]
        )
        self.single_element_orbit = GeneralPerturbationsOrbit(
            elements=[base.model_copy(update={"epoch": self.epoch_0})]
        )

    def test_none_always_returns_first_index(self):
        """
        Test that at_times=None always selects index 0, regardless of
        how many elements exist.
        """
        self.assertEqual(self.single_element_orbit.get_closest_element_index(None), 0)
        self.assertEqual(self.multi_element_orbit.get_closest_element_index(None), 0)

    def test_single_element_orbit_always_returns_zero(self):
        """
        Test that a single-element orbit always returns index 0 for any
        query time, since there is only one element to choose from.
        """
        far_future = self.epoch_0 + timedelta(days=3650)
        self.assertEqual(
            self.single_element_orbit.get_closest_element_index(far_future), 0
        )

    def test_query_exactly_at_an_epoch_returns_that_index(self):
        """
        Test that querying exactly at an element's epoch returns that
        element's index.
        """
        o = self.multi_element_orbit
        self.assertEqual(o.get_closest_element_index(self.epoch_0), 0)
        self.assertEqual(o.get_closest_element_index(self.epoch_1), 1)
        self.assertEqual(o.get_closest_element_index(self.epoch_2), 2)

    def test_query_before_first_epoch_clamps_to_first_index(self):
        """
        Test that a query time before every element's epoch returns the
        first (nearest) index.
        """
        query = self.epoch_0 - timedelta(days=365)
        self.assertEqual(self.multi_element_orbit.get_closest_element_index(query), 0)

    def test_query_after_last_epoch_clamps_to_last_index(self):
        """
        Test that a query time after every element's epoch returns the
        last (nearest) index.
        """
        query = self.epoch_2 + timedelta(days=365)
        self.assertEqual(self.multi_element_orbit.get_closest_element_index(query), 2)

    def test_query_closer_to_earlier_epoch(self):
        """
        Test that a query time closer to an earlier epoch than the next
        one returns the earlier index.
        """
        query = self.epoch_0 + timedelta(hours=1)  # much closer to epoch_0 than epoch_1
        self.assertEqual(self.multi_element_orbit.get_closest_element_index(query), 0)

    def test_query_closer_to_later_epoch(self):
        """
        Test that a query time closer to a later epoch than the previous
        one returns the later index.
        """
        query = self.epoch_1 - timedelta(hours=1)  # much closer to epoch_1 than epoch_0
        self.assertEqual(self.multi_element_orbit.get_closest_element_index(query), 1)

    def test_query_at_exact_midpoint_breaks_tie_toward_later_index(self):
        """
        Test the documented tie-breaking convention: a query exactly
        equidistant between two epochs resolves to the later index
        (since the implementation only prefers the earlier neighbor when
        it is strictly closer).
        """
        midpoint = self.epoch_0 + (self.epoch_1 - self.epoch_0) / 2
        self.assertEqual(
            self.multi_element_orbit.get_closest_element_index(midpoint), 1
        )

    def test_list_input_returns_list_of_indices(self):
        """
        Test that a list of query times returns a list of indices, one
        per query, matching each query's closest element independently.
        """
        queries = [
            self.epoch_0 - timedelta(days=365),
            self.epoch_0 + (self.epoch_1 - self.epoch_0) / 2,
            self.epoch_2 + timedelta(days=365),
        ]
        self.assertEqual(
            self.multi_element_orbit.get_closest_element_index(queries), [0, 1, 2]
        )

    def test_ndarray_input_returns_list_of_indices(self):
        """
        Test that a numpy datetime64 array of query times is accepted
        and returns the same result as an equivalent list of datetimes.
        """
        # strip tzinfo (all test epochs are UTC) before building the
        # array, since numpy warns when directly casting timezone-aware
        # datetimes to datetime64
        queries = np.array(
            [
                self.epoch_0.replace(tzinfo=None),
                self.epoch_1.replace(tzinfo=None),
                self.epoch_2.replace(tzinfo=None),
            ],
            dtype="datetime64[ns]",
        )
        self.assertEqual(
            self.multi_element_orbit.get_closest_element_index(queries), [0, 1, 2]
        )

    def test_empty_list_input_returns_empty_list(self):
        """
        Test that an empty list of query times returns an empty list of
        indices, rather than raising an error.
        """
        self.assertEqual(self.multi_element_orbit.get_closest_element_index([]), [])


class TestGetClosestElement(unittest.TestCase):
    """
    Unit tests for GeneralPerturbationsOrbit.get_closest_element. Since
    this is a thin wrapper mapping get_closest_element_index's result
    onto self.elements, index-selection edge cases (ties, boundaries,
    sorting) are covered by TestGetClosestElementIndex; these tests focus
    on the element-mapping behavior itself.
    """

    def setUp(self):
        tle = [
            "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
            "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
        ]
        base = GeneralPerturbationsOrbit.from_tle(tle).elements[0]
        self.epoch_0 = datetime(2022, 1, 1, tzinfo=timezone.utc)
        self.epoch_1 = datetime(2022, 1, 3, tzinfo=timezone.utc)
        self.element_0 = base.model_copy(
            update={"norad_cat_id": 11111, "epoch": self.epoch_0}
        )
        self.element_1 = base.model_copy(
            update={"norad_cat_id": 22222, "epoch": self.epoch_1}
        )
        self.multi_element_orbit = GeneralPerturbationsOrbit(
            elements=[self.element_0, self.element_1]
        )

    def test_none_returns_first_element(self):
        """
        Test that at_times=None returns the first element.
        """
        result = self.multi_element_orbit.get_closest_element(None)
        self.assertIs(result, self.element_0)

    def test_scalar_returns_matching_element(self):
        """
        Test that a scalar query time returns the actual closest element
        object (by identity), not merely its index.
        """
        near_epoch_1 = self.epoch_1 - timedelta(hours=1)
        result = self.multi_element_orbit.get_closest_element(near_epoch_1)
        self.assertIs(result, self.element_1)
        self.assertEqual(result.norad_cat_id, 22222)

    def test_list_input_returns_list_of_matching_elements(self):
        """
        Test that a list of query times returns a list of the
        corresponding closest element objects, in query order.
        """
        queries = [
            self.epoch_0 + timedelta(hours=1),
            self.epoch_1 - timedelta(hours=1),
        ]
        result = self.multi_element_orbit.get_closest_element(queries)
        self.assertEqual(len(result), 2)
        self.assertIs(result[0], self.element_0)
        self.assertIs(result[1], self.element_1)

    def test_empty_list_input_returns_empty_list(self):
        """
        Test that an empty list of query times returns an empty list of
        elements, rather than raising an error.
        """
        self.assertEqual(self.multi_element_orbit.get_closest_element([]), [])

    def test_single_element_orbit_always_returns_that_element(self):
        """
        Test that a single-element orbit always returns its one element,
        regardless of the query time.
        """
        o = GeneralPerturbationsOrbit(elements=[self.element_0])
        far_future = self.epoch_0 + timedelta(days=3650)
        self.assertIs(o.get_closest_element(far_future), self.element_0)


class TestGetOrbitTrackAtTime(unittest.TestCase):
    """
    Unit tests for GeneralPerturbationsOrbit.get_orbit_track_at_time. This
    method is mostly a thin wrapper around Skyfield's own SGP4 propagation
    (EarthSatellite.at); the specific behavior worth validating here is
    that a multi-element orbit selects and propagates the correct element
    per query time, both for a single scalar time and for a vectorized
    batch spanning multiple elements' regions.
    """

    def setUp(self):
        tle = [
            "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
            "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
        ]
        base = GeneralPerturbationsOrbit.from_tle(tle).elements[0]
        self.epoch_0 = datetime(2022, 1, 1, tzinfo=timezone.utc)
        self.epoch_1 = datetime(2022, 1, 3, tzinfo=timezone.utc)
        self.epoch_2 = datetime(2022, 1, 5, tzinfo=timezone.utc)
        self.element_0 = base.model_copy(update={"epoch": self.epoch_0})
        self.element_1 = base.model_copy(update={"epoch": self.epoch_1})
        self.element_2 = base.model_copy(update={"epoch": self.epoch_2})
        self.multi_element_orbit = GeneralPerturbationsOrbit(
            elements=[self.element_0, self.element_1, self.element_2]
        )
        self.single_element_orbit = GeneralPerturbationsOrbit(
            elements=[self.element_0]
        )

    def test_single_element_scalar_time_matches_direct_propagation(self):
        """
        Test that a single-element orbit's result at a scalar time is
        identical to propagating that element directly.
        """
        t = constants.timescale.from_datetime(self.epoch_0 + timedelta(hours=5))
        expected = self.element_0.to_skyfield().at(t)
        actual = self.single_element_orbit.get_orbit_track_at_time(t)
        self.assertTrue(np.array_equal(actual.position.km, expected.position.km))
        self.assertTrue(
            np.array_equal(actual.velocity.km_per_s, expected.velocity.km_per_s)
        )

    def test_single_element_vector_time_matches_direct_propagation(self):
        """
        Test that a single-element orbit's result at a vector of times is
        identical to propagating that element directly.
        """
        times = [self.epoch_0 + timedelta(hours=h) for h in (1, 5, 10)]
        t = constants.timescale.from_datetimes(times)
        expected = self.element_0.to_skyfield().at(t)
        actual = self.single_element_orbit.get_orbit_track_at_time(t)
        self.assertTrue(np.array_equal(actual.position.km, expected.position.km))

    def test_multi_element_scalar_time_uses_nearest_element(self):
        """
        Test that a scalar query time near epoch_1 is propagated using
        element_1, not element_0 -- verified both by matching element_1's
        own direct propagation, and by confirming element_0 would have
        given a different (wrong) answer, since the two elements carry
        different epochs and so a different elapsed-time offset from the
        same absolute query time.
        """
        t = constants.timescale.from_datetime(self.epoch_1 - timedelta(hours=1))
        actual = self.multi_element_orbit.get_orbit_track_at_time(t)
        expected = self.element_1.to_skyfield().at(t)
        self.assertTrue(np.array_equal(actual.position.km, expected.position.km))
        wrong = self.element_0.to_skyfield().at(t)
        self.assertFalse(np.array_equal(actual.position.km, wrong.position.km))

    def test_multi_element_vector_time_selects_nearest_element_per_time(self):
        """
        Test that a vectorized query spanning all three elements' regions
        propagates each time with its own nearest element, rather than
        applying a single element to the whole batch -- exercising the
        implementation's grouped-by-unique-index vectorization.
        """
        times = [
            self.epoch_0 + timedelta(hours=1),  # nearest element_0
            self.epoch_1 - timedelta(hours=1),  # nearest element_1
            self.epoch_1 + timedelta(hours=1),  # nearest element_1
            self.epoch_2 + timedelta(hours=1),  # nearest element_2
        ]
        expected_elements = [
            self.element_0,
            self.element_1,
            self.element_1,
            self.element_2,
        ]
        t = constants.timescale.from_datetimes(times)
        actual = self.multi_element_orbit.get_orbit_track_at_time(t)
        for i, (time, element) in enumerate(zip(times, expected_elements)):
            expected = element.to_skyfield().at(
                constants.timescale.from_datetime(time)
            )
            self.assertTrue(
                np.array_equal(actual.position.km[:, i], expected.position.km),
                f"time index {i} did not match its nearest element",
            )

    def test_multi_element_vector_time_all_same_nearest_element(self):
        """
        Test the edge case where every time in a vectorized query shares
        the same nearest element, so the grouped-by-unique-index loop
        runs exactly once, still returning correct per-time results.
        """
        times = [self.epoch_1 + timedelta(minutes=m) for m in (10, 20, 30)]
        t = constants.timescale.from_datetimes(times)
        expected = self.element_1.to_skyfield().at(t)
        actual = self.multi_element_orbit.get_orbit_track_at_time(t)
        self.assertTrue(np.array_equal(actual.position.km, expected.position.km))


class TestGetOrbitTrack(unittest.TestCase):
    """
    Unit tests for GeneralPerturbationsOrbit.get_orbit_track. This method
    only builds a Skyfield Time from the given datetime(s) and delegates
    to get_orbit_track_at_time (see TestGetOrbitTrackAtTime for
    propagation and multi-element selection coverage), so these tests
    focus on the datetime-to-Time dispatch itself: a scalar datetime vs.
    a list, including the list-of-one edge case.
    """

    def setUp(self):
        tle = [
            "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
            "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
        ]
        self.orbit = GeneralPerturbationsOrbit.from_tle(tle)
        self.epoch = self.orbit.elements[0].epoch

    def test_scalar_datetime_matches_get_orbit_track_at_time(self):
        """
        Test that a single datetime produces the same result as calling
        get_orbit_track_at_time with the equivalent scalar Skyfield Time.
        """
        time = self.epoch + timedelta(hours=2)
        expected = self.orbit.get_orbit_track_at_time(
            constants.timescale.from_datetime(time)
        )
        actual = self.orbit.get_orbit_track(time)
        self.assertTrue(np.array_equal(actual.position.km, expected.position.km))

    def test_list_of_datetimes_matches_get_orbit_track_at_time(self):
        """
        Test that a list of datetimes produces the same result as calling
        get_orbit_track_at_time with the equivalent vector Skyfield Time.
        """
        times = [self.epoch + timedelta(hours=h) for h in (1, 2, 3)]
        expected = self.orbit.get_orbit_track_at_time(
            constants.timescale.from_datetimes(times)
        )
        actual = self.orbit.get_orbit_track(times)
        self.assertTrue(np.array_equal(actual.position.km, expected.position.km))

    def test_scalar_datetime_returns_scalar_shaped_result(self):
        """
        Test that a scalar datetime returns a scalar (non-vectorized)
        position, distinguishing it from an equal-valued single-item list
        (see test_list_of_one_datetime_returns_vector_shaped_result).
        """
        time = self.epoch + timedelta(hours=1)
        result = self.orbit.get_orbit_track(time)
        self.assertEqual(result.position.km.shape, (3,))

    def test_list_of_one_datetime_returns_vector_shaped_result(self):
        """
        Test that a single-item list is still treated as a vector query
        (shape (3, 1)), not collapsed to the scalar shape a bare datetime
        would produce -- the dispatch is based on the input's type
        (datetime vs. list), not its length.
        """
        time = self.epoch + timedelta(hours=1)
        result = self.orbit.get_orbit_track([time])
        self.assertEqual(result.position.km.shape, (3, 1))
        # and the single entry's value matches the scalar-input case
        scalar_result = self.orbit.get_orbit_track(time)
        self.assertTrue(
            np.array_equal(result.position.km[:, 0], scalar_result.position.km)
        )


class TestGetGeographicPositionAtTime(unittest.TestCase):
    """
    Unit tests for GeneralPerturbationsOrbit.get_geographic_position_at_time.
    This method is almost the same as get_orbit_track_at_time (already
    covered elsewhere), but optionally substitutes a nearby, epoch-relative
    time for a distant one when a repeat cycle is known, trading a little
    ground-track drift for much less accumulated SGP4 propagation error.
    These tests focus on that substitution: when it applies, whether it
    preserves the sign of the time offset from epoch, and when it falls
    back to direct propagation.
    """

    def setUp(self):
        # Landsat-8: real, single-element, actively repeating orbit
        landsat_8_tle = [
            "1 39084U 13008A   26213.27824675  .00000294  00000+0  75333-4 0  9990",
            "2 39084  98.2277 282.8718 0001275  92.4910 267.6434 14.57104473704466",
        ]
        self.repeat_orbit = GeneralPerturbationsOrbit.from_tle(landsat_8_tle)
        self.epoch = self.repeat_orbit.get_epoch()
        self.repeat_cycle = self.repeat_orbit.get_repeat_cycle()
        self.assertIsNotNone(self.repeat_cycle)

        # ISS: single-element, not designed for a repeat ground track
        iss_tle = [
            "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
            "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
        ]
        self.non_repeat_orbit = GeneralPerturbationsOrbit.from_tle(iss_tle)

        # a multi-element orbit (elements otherwise identical to Landsat-8,
        # just re-epoched), to confirm the substitution never applies
        base = self.repeat_orbit.elements[0]
        self.multi_element_orbit = GeneralPerturbationsOrbit(
            elements=[
                base.model_copy(update={"epoch": self.epoch}),
                base.model_copy(update={"epoch": self.epoch + timedelta(days=1)}),
            ]
        )

    def test_try_repeat_false_matches_direct_propagation(self):
        """
        Test that try_repeat=False always returns the true, directly
        propagated position, even for an orbit with a known repeat cycle.
        """
        t = constants.timescale.from_datetime(self.epoch + self.repeat_cycle * 3)
        actual = self.repeat_orbit.get_geographic_position_at_time(
            t, try_repeat=False
        )
        expected = wgs84.geographic_position_of(
            self.repeat_orbit.get_orbit_track_at_time(t)
        )
        self.assertEqual(actual.latitude.degrees, expected.latitude.degrees)
        self.assertEqual(actual.longitude.degrees, expected.longitude.degrees)

    def test_try_repeat_true_within_first_cycle_matches_direct_propagation(self):
        """
        Test that, for a time less than one repeat cycle after epoch, the
        substitution is a no-op (the wrapped offset equals the original
        offset), so try_repeat=True and try_repeat=False agree exactly.
        """
        t = constants.timescale.from_datetime(self.epoch + timedelta(hours=5))
        with_repeat = self.repeat_orbit.get_geographic_position_at_time(
            t, try_repeat=True
        )
        direct = self.repeat_orbit.get_geographic_position_at_time(
            t, try_repeat=False
        )
        self.assertEqual(with_repeat.latitude.degrees, direct.latitude.degrees)
        self.assertEqual(with_repeat.longitude.degrees, direct.longitude.degrees)

    def test_try_repeat_true_substitutes_epoch_relative_time_for_far_future(self):
        """
        Test that a time several repeat cycles in the future is
        substituted with the equivalent epoch-relative offset (t's offset
        from epoch, wrapped modulo the repeat cycle) rather than
        propagated directly -- verified against a manual replica of that
        formula, and confirmed to differ from true direct propagation
        (which accumulates more SGP4 error over the longer elapsed time).
        """
        far_future = self.epoch + self.repeat_cycle * 3 + timedelta(hours=5)
        t = constants.timescale.from_datetime(far_future)
        actual = self.repeat_orbit.get_geographic_position_at_time(
            t, try_repeat=True
        )

        offset_days = (far_future - self.epoch) / timedelta(days=1)
        cycle_days = self.repeat_cycle / timedelta(days=1)
        wrapped_offset = timedelta(days=float(np.mod(offset_days, cycle_days)))
        expected = wgs84.geographic_position_of(
            self.repeat_orbit.elements[0]
            .to_skyfield()
            .at(constants.timescale.from_datetime(self.epoch + wrapped_offset))
        )
        self.assertAlmostEqual(
            actual.latitude.degrees, expected.latitude.degrees, places=9
        )
        self.assertAlmostEqual(
            actual.longitude.degrees, expected.longitude.degrees, places=9
        )

        direct = wgs84.geographic_position_of(
            self.repeat_orbit.get_orbit_track_at_time(t)
        )
        self.assertNotEqual(actual.latitude.degrees, direct.latitude.degrees)

    def test_try_repeat_true_preserves_sign_for_time_before_epoch(self):
        """
        Regression test: a query 2.5 repeat cycles *before* epoch must
        wrap to -0.5 cycles (epoch minus half a cycle), not +0.5 cycles.
        A naive numpy np.mod() on the raw (negative) offset always
        returns a non-negative result, which would silently wrap to the
        wrong side of the repeat cycle -- a real bug this implementation
        avoids by explicitly reapplying the offset's original sign.
        """
        query_time = self.epoch - self.repeat_cycle * 2.5
        t = constants.timescale.from_datetime(query_time)
        actual = self.repeat_orbit.get_geographic_position_at_time(
            t, try_repeat=True
        )

        correct = wgs84.geographic_position_of(
            self.repeat_orbit.elements[0]
            .to_skyfield()
            .at(
                constants.timescale.from_datetime(
                    self.epoch - self.repeat_cycle * 0.5
                )
            )
        )
        wrong = wgs84.geographic_position_of(
            self.repeat_orbit.elements[0]
            .to_skyfield()
            .at(
                constants.timescale.from_datetime(
                    self.epoch + self.repeat_cycle * 0.5
                )
            )
        )
        self.assertAlmostEqual(
            actual.latitude.degrees, correct.latitude.degrees, places=9
        )
        self.assertNotAlmostEqual(
            actual.latitude.degrees, wrong.latitude.degrees, places=2
        )

    def test_try_repeat_true_ignored_for_multi_element_orbit(self):
        """
        Test that the repeat-cycle substitution never applies to a
        multi-element orbit, even with try_repeat=True and a detectable
        single-element repeat cycle -- it always falls back to direct,
        per-time nearest-element propagation (get_orbit_track_at_time).
        """
        t = constants.timescale.from_datetime(self.epoch + timedelta(days=40))
        actual = self.multi_element_orbit.get_geographic_position_at_time(
            t, try_repeat=True
        )
        expected = wgs84.geographic_position_of(
            self.multi_element_orbit.get_orbit_track_at_time(t)
        )
        self.assertEqual(actual.latitude.degrees, expected.latitude.degrees)
        self.assertEqual(actual.longitude.degrees, expected.longitude.degrees)

    def test_try_repeat_true_falls_back_when_no_repeat_cycle_found(self):
        """
        Test that try_repeat=True falls back to direct propagation for an
        orbit (the ISS) with no detectable repeat cycle.
        """
        t = constants.timescale.from_datetime(
            self.non_repeat_orbit.get_epoch() + timedelta(days=10)
        )
        actual = self.non_repeat_orbit.get_geographic_position_at_time(
            t, try_repeat=True
        )
        expected = wgs84.geographic_position_of(
            self.non_repeat_orbit.get_orbit_track_at_time(t)
        )
        self.assertEqual(actual.latitude.degrees, expected.latitude.degrees)
        self.assertEqual(actual.longitude.degrees, expected.longitude.degrees)

    def test_try_repeat_none_uses_config_default(self):
        """
        Test that omitting try_repeat follows
        config.get_rc().repeat_cycle_for_orbit_track.
        """
        t = constants.timescale.from_datetime(
            self.epoch + self.repeat_cycle * 3 + timedelta(hours=5)
        )
        original = config.get_rc().repeat_cycle_for_orbit_track
        try:
            config.get_rc().repeat_cycle_for_orbit_track = True
            with_default_true = self.repeat_orbit.get_geographic_position_at_time(t)
            config.get_rc().repeat_cycle_for_orbit_track = False
            with_default_false = self.repeat_orbit.get_geographic_position_at_time(t)
        finally:
            config.get_rc().repeat_cycle_for_orbit_track = original
        expected_true = self.repeat_orbit.get_geographic_position_at_time(
            t, try_repeat=True
        )
        expected_false = self.repeat_orbit.get_geographic_position_at_time(
            t, try_repeat=False
        )
        self.assertEqual(
            with_default_true.latitude.degrees, expected_true.latitude.degrees
        )
        self.assertEqual(
            with_default_false.latitude.degrees, expected_false.latitude.degrees
        )

    def test_vectorized_time_substitution_matches_per_time_computation(self):
        """
        Test that a vectorized query spanning multiple repeat cycles, on
        both sides of epoch, produces the same result as computing each
        time individually -- exercising the array branch of the epoch-
        relative substitution (the tests above use the scalar branch).
        """
        query_times = [
            self.epoch + timedelta(hours=5),
            self.epoch + self.repeat_cycle * 3 + timedelta(hours=5),
            self.epoch - self.repeat_cycle * 2.5,
        ]
        t_vector = constants.timescale.from_datetimes(query_times)
        actual = self.repeat_orbit.get_geographic_position_at_time(
            t_vector, try_repeat=True
        )
        for i, query_time in enumerate(query_times):
            expected = self.repeat_orbit.get_geographic_position_at_time(
                constants.timescale.from_datetime(query_time), try_repeat=True
            )
            self.assertAlmostEqual(
                actual.latitude.degrees[i], expected.latitude.degrees, places=9
            )
            self.assertAlmostEqual(
                actual.longitude.degrees[i], expected.longitude.degrees, places=9
            )


class TestPartitionByElementIndex(unittest.TestCase):
    """
    Unit tests for GeneralPerturbationsOrbit.partition_by_element_index.
    """

    def setUp(self):
        tle = [
            "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
            "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
        ]
        base = GeneralPerturbationsOrbit.from_tle(tle).elements[0]
        self.epoch_0 = datetime(2022, 1, 1, tzinfo=timezone.utc)
        self.epoch_1 = datetime(2022, 1, 3, tzinfo=timezone.utc)
        self.epoch_2 = datetime(2022, 1, 5, tzinfo=timezone.utc)
        self.multi_element_orbit = GeneralPerturbationsOrbit(
            elements=[
                base.model_copy(update={"epoch": self.epoch_0}),
                base.model_copy(update={"epoch": self.epoch_1}),
                base.model_copy(update={"epoch": self.epoch_2}),
            ]
        )

    def test_uses_get_element_epochs_cache(self):
        """
        Test that partition_by_element_index reads epochs through
        get_element_epochs (populating its lazy-load cache) rather than
        re-deriving them independently, so a subsequent call reuses the
        same cached list.
        """
        self.assertIsNone(self.multi_element_orbit.__dict__.get("element_epochs"))
        start = self.epoch_0 - timedelta(days=365)
        end = self.epoch_2 + timedelta(days=365)
        self.multi_element_orbit.partition_by_element_index(start, end)
        cached = self.multi_element_orbit.__dict__.get("element_epochs")
        self.assertIsNotNone(cached)
        self.assertEqual(cached, [self.epoch_0, self.epoch_1, self.epoch_2])

    def test_single_element_orbit(self):
        """
        Test that a single-element orbit returns exactly one segment
        covering the whole window, assigned to element 0.
        """
        o = GeneralPerturbationsOrbit(
            elements=[
                self.multi_element_orbit.elements[0].model_copy(
                    update={"epoch": self.epoch_0}
                )
            ]
        )
        start = self.epoch_0 - timedelta(days=1)
        end = self.epoch_0 + timedelta(days=1)
        boundary_times, segment_indices = o.partition_by_element_index(start, end)
        self.assertEqual(boundary_times, [start, end])
        self.assertEqual(segment_indices, [0])

    def test_window_spanning_all_elements(self):
        """
        Regression test for a bug where this method previously crashed
        with a TypeError comparing Python datetime against numpy
        datetime64 whenever called with more than one element -- the
        exact scenario this method exists for. A window spanning well
        before the first and well after the last epoch should produce
        one segment per element, in order.
        """
        start = self.epoch_0 - timedelta(days=365)
        end = self.epoch_2 + timedelta(days=365)
        boundary_times, segment_indices = (
            self.multi_element_orbit.partition_by_element_index(start, end)
        )
        self.assertEqual(len(boundary_times), 4)
        self.assertEqual(boundary_times[0], start)
        self.assertEqual(boundary_times[-1], end)
        self.assertEqual(segment_indices, [0, 1, 2])

    def test_boundary_times_are_plain_datetimes(self):
        """
        Regression test: every boundary time (including the internally
        computed midpoints, not just start/end) must be a plain Python
        datetime, not a numpy.datetime64 -- the downstream caller
        (get_observation_events) passes each one directly to
        skyfield's Timescale.from_datetime(), which requires a real
        datetime object.
        """
        start = self.epoch_0 - timedelta(days=365)
        end = self.epoch_2 + timedelta(days=365)
        boundary_times, _ = self.multi_element_orbit.partition_by_element_index(
            start, end
        )
        for t in boundary_times:
            with self.subTest(t=t):
                self.assertIsInstance(t, datetime)

    def test_window_within_single_elements_region(self):
        """
        Test that a narrow window falling entirely within one element's
        closest-region (no epoch midpoint inside the window) returns a
        single segment for that element, not one per element.
        """
        start = self.epoch_1 - timedelta(hours=1)
        end = self.epoch_1 + timedelta(hours=1)
        boundary_times, segment_indices = (
            self.multi_element_orbit.partition_by_element_index(start, end)
        )
        self.assertEqual(boundary_times, [start, end])
        self.assertEqual(segment_indices, [1])

    def test_window_spanning_only_first_midpoint(self):
        """
        Test a window that includes the first epoch midpoint but not the
        second, producing two segments (elements 0 and 1).
        """
        start = self.epoch_0
        end = self.epoch_1 + timedelta(hours=12)
        boundary_times, segment_indices = (
            self.multi_element_orbit.partition_by_element_index(start, end)
        )
        self.assertEqual(len(boundary_times), 3)
        self.assertEqual(boundary_times[0], start)
        self.assertEqual(boundary_times[-1], end)
        self.assertEqual(segment_indices, [0, 1])

    def test_number_of_indices_is_one_less_than_boundary_times(self):
        """
        Test the length contract: N+1 boundary times always yield
        exactly N segment indices, for a variety of window sizes.
        """
        cases = [
            (self.epoch_0 - timedelta(days=365), self.epoch_2 + timedelta(days=365)),
            (self.epoch_1 - timedelta(hours=1), self.epoch_1 + timedelta(hours=1)),
            (self.epoch_0, self.epoch_1 + timedelta(hours=12)),
        ]
        for start, end in cases:
            with self.subTest(start=start, end=end):
                boundary_times, segment_indices = (
                    self.multi_element_orbit.partition_by_element_index(start, end)
                )
                self.assertEqual(len(segment_indices), len(boundary_times) - 1)

    def test_get_observation_events_does_not_crash_for_multi_element_orbit(self):
        """
        Smoke test for the underlying bug: get_observation_events, the
        sole caller of partition_by_element_index, must not crash when
        given a multi-element orbit. Full coverage of
        get_observation_events itself is a separate, later piece of work.
        """
        point = Point(id=0, latitude=40.0, longitude=-74.0)
        start = self.epoch_0
        end = self.epoch_0 + timedelta(hours=12)
        _, events = self.multi_element_orbit.get_observation_events(
            point, start, end, min_elevation_angle=10, try_repeat=False
        )
        self.assertGreater(len(events), 0)


class TestGetRepeatCycle(unittest.TestCase):
    """
    Unit tests for GeneralPerturbationsOrbit.get_repeat_cycle.
    """

    def setUp(self):
        # a real Landsat-8 element set (NORAD 39084, fetched from
        # Celestrak), whose published repeat ground track is exactly 233
        # orbits every 16 days (USGS)
        self.landsat_8_tle = [
            "1 39084U 13008A   26213.27824675  .00000294  00000+0  75333-4 0  9990",
            "2 39084  98.2277 282.8718 0001275  92.4910 267.6434 14.57104473704466",
        ]
        # a real Sentinel-2A element set (NORAD 40697), whose published
        # repeat ground track is exactly 143 orbits every 10 days (ESA)
        self.sentinel_2a_tle = [
            "1 40697U 15028A   26213.24967738  .00000092  00000+0  51774-4 0  9990",
            "2 40697  98.5671 287.4570 0001329  92.6662 267.4673 14.30818788580177",
        ]
        # a real Sentinel-1A element set (NORAD 39634) from early 2026,
        # while the satellite was still actively operated (its mission
        # concluded 2026-06-29): published repeat ground track is exactly
        # 175 orbits every 12 days (ESA)
        self.sentinel_1a_tle = [
            "1 39634U 14016A   26001.19041520  .00000521  00000-0  12012-3 0  9995",
            "2 39634  98.1805  11.1453 0001276  85.2963 274.8383 14.59199668625660",
        ]
        # a real GPS BIIR-5 element set (NORAD 26407): GPS orbits are
        # designed to repeat their ground track every sidereal day (two
        # ~12-hour orbits per day), a much shorter cycle than the
        # sun-synchronous imaging orbits above, at MEO altitude (~20,200 km)
        self.gps_tle = [
            "1 26407U 00040A   26213.32131914  .00000072  00000+0  00000+0 0  9998",
            "2 26407  54.8467 213.3697 0120005 302.9740 169.1529  2.00558010190856",
        ]
        # a real Molniya 3-8 element set (NORAD 10455): a highly eccentric,
        # critical-inclination (~63.4 degree) orbit that, like GPS, repeats
        # every sidereal day by design, but is geometrically nothing like
        # the near-circular orbits above
        self.molniya_tle = [
            "1 10455U 77105A   26212.97315042  .00000728  00000+0  00000+0 0  9994",
            "2 10455  63.8024 172.4301 6701249 276.1687  17.0849  2.00778403357346",
        ]
        # the ISS is not designed for a repeat ground track
        self.iss_tle = [
            "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
            "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
        ]

    def test_landsat_8_repeat_cycle_matches_published_16_days(self):
        """
        Test against Landsat-8's published repeat ground track of 233
        orbits every 16 days (USGS). A small tolerance accounts for the
        real orbit's minor drift between station-keeping maneuvers.
        """
        orbit = GeneralPerturbationsOrbit.from_tle(self.landsat_8_tle)
        repeat_cycle = orbit.get_repeat_cycle()
        self.assertIsNotNone(repeat_cycle)
        self.assertAlmostEqual(repeat_cycle.total_seconds() / 86400, 16, delta=0.1)

    def test_sentinel_2a_repeat_cycle_matches_published_10_days(self):
        """
        Test against Sentinel-2A's published repeat ground track of 143
        orbits every 10 days (ESA), at a different altitude/inclination
        than Landsat-8, using the default tolerances.
        """
        orbit = GeneralPerturbationsOrbit.from_tle(self.sentinel_2a_tle)
        repeat_cycle = orbit.get_repeat_cycle()
        self.assertIsNotNone(repeat_cycle)
        self.assertAlmostEqual(repeat_cycle.total_seconds() / 86400, 10, delta=0.1)

    def test_sentinel_1a_repeat_cycle_matches_published_12_days(self):
        """
        Test against Sentinel-1A's published repeat ground track of 175
        orbits every 12 days (ESA), at a different altitude/inclination
        than Landsat-8, using the default tolerances.
        """
        orbit = GeneralPerturbationsOrbit.from_tle(self.sentinel_1a_tle)
        repeat_cycle = orbit.get_repeat_cycle()
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
        orbit = GeneralPerturbationsOrbit.from_tle(self.gps_tle)
        repeat_cycle = orbit.get_repeat_cycle()
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
        orbit = GeneralPerturbationsOrbit.from_tle(self.molniya_tle)
        repeat_cycle = orbit.get_repeat_cycle(
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
        orbit = GeneralPerturbationsOrbit.from_tle(self.iss_tle)
        self.assertIsNone(orbit.get_repeat_cycle())

    def test_lazy_load_reuses_cached_result(self):
        """
        Test that calling get_repeat_cycle() twice with lazy_load=True
        (the default) returns the identical cached timedelta rather than
        recomputing.
        """
        orbit = GeneralPerturbationsOrbit.from_tle(self.landsat_8_tle)
        first = orbit.get_repeat_cycle()
        second = orbit.get_repeat_cycle()
        self.assertIs(first, second)

    def test_lazy_load_false_forces_recomputation(self):
        """
        Test that lazy_load=False recomputes rather than reusing the
        cached result (a fresh but equal timedelta).
        """
        orbit = GeneralPerturbationsOrbit.from_tle(self.landsat_8_tle)
        first = orbit.get_repeat_cycle()
        second = orbit.get_repeat_cycle(lazy_load=False)
        self.assertEqual(first, second)
        self.assertIsNot(first, second)

    def test_too_short_search_duration_returns_none(self):
        """
        Test that a max_search_duration shorter than the true repeat
        cycle (16 days) cannot find it.
        """
        orbit = GeneralPerturbationsOrbit.from_tle(self.landsat_8_tle)
        repeat_cycle = orbit.get_repeat_cycle(
            max_search_duration=timedelta(days=10), lazy_load=False
        )
        self.assertIsNone(repeat_cycle)

    def test_too_tight_tolerance_returns_none(self):
        """
        Test that an unrealistically tight position/velocity tolerance
        rejects even the real repeat cycle.
        """
        orbit = GeneralPerturbationsOrbit.from_tle(self.landsat_8_tle)
        repeat_cycle = orbit.get_repeat_cycle(
            max_delta_position=1, max_delta_velocity=0.001, lazy_load=False
        )
        self.assertIsNone(repeat_cycle)

    def test_multi_element_consistent_repeat_cycles_are_combined(self):
        """
        Test that a multi-element orbit reports a repeat cycle when every
        element's own repeat cycle agrees within the consistency
        threshold. The second element here is a near-identical copy (mean
        motion nudged by 1e-6, a negligible fitting-noise-scale change)
        of the real Landsat-8 element, so both independently resolve to
        ~16 days, well within the default 1-hour threshold.
        """
        base = GeneralPerturbationsOrbit.from_tle(self.landsat_8_tle).elements[0]
        nearly_identical = base.model_copy(
            update={"mean_motion": base.mean_motion * (1 + 1e-6)}
        )
        orbit = GeneralPerturbationsOrbit(elements=[base, nearly_identical])
        repeat_cycle = orbit.get_repeat_cycle()
        self.assertIsNotNone(repeat_cycle)
        self.assertAlmostEqual(repeat_cycle.total_seconds() / 86400, 16, delta=0.1)

    def test_multi_element_valid_but_inconsistent_repeat_cycles_return_none(self):
        """
        Test that a multi-element orbit returns None when its elements
        each have a valid repeat cycle, but they disagree well beyond the
        consistency threshold -- simulating a maneuver partway through
        the orbit's history that changed its fundamental repeat behavior.
        A 0.05% mean motion change shifts the second element's best
        (smallest-drift) candidate from 16 days to 7 days, at a loosened
        but still fairly tight tolerance (40 km) chosen so that the
        original element still resolves cleanly to 16 days too (its own
        other candidate days all drift well over 100 km, so 40 km can't
        accidentally admit any of them) -- isolating the consistency
        check from the tolerance itself.
        """
        base = GeneralPerturbationsOrbit.from_tle(self.landsat_8_tle).elements[0]
        maneuvered = base.model_copy(
            update={"mean_motion": base.mean_motion * 1.0005}
        )
        orbit = GeneralPerturbationsOrbit(elements=[base, maneuvered])
        max_delta_position = 40000
        max_delta_velocity = 5
        base_cycle = base.get_repeat_cycle(max_delta_position, max_delta_velocity)
        maneuvered_cycle = maneuvered.get_repeat_cycle(
            max_delta_position, max_delta_velocity
        )
        self.assertIsNotNone(base_cycle)
        self.assertIsNotNone(maneuvered_cycle)
        self.assertNotEqual(base_cycle, maneuvered_cycle)
        self.assertIsNone(
            orbit.get_repeat_cycle(max_delta_position, max_delta_velocity)
        )

    def test_multi_element_one_element_without_repeat_cycle_returns_none(self):
        """
        Test that a multi-element orbit returns None if any element has
        no repeat cycle at all (here, a 0.1% mean motion change from the
        real Landsat-8 element, which at default tolerances finds no
        commensurate day within the search duration), even though the
        other element's own repeat cycle is perfectly valid.
        """
        base = GeneralPerturbationsOrbit.from_tle(self.landsat_8_tle).elements[0]
        maneuvered = base.model_copy(update={"mean_motion": base.mean_motion * 1.001})
        orbit = GeneralPerturbationsOrbit(elements=[base, maneuvered])
        self.assertIsNotNone(base.get_repeat_cycle())
        self.assertIsNone(maneuvered.get_repeat_cycle())
        self.assertIsNone(orbit.get_repeat_cycle())

    def test_consistency_threshold_override_changes_accept_reject_boundary(self):
        """
        Test that consistency_threshold controls the accept/reject
        boundary directly: the same pair of elements (differing by
        ~1.59 seconds in their independently-computed repeat cycles) is
        rejected under a 1-second threshold and accepted under a
        2-second threshold.
        """
        base = GeneralPerturbationsOrbit.from_tle(self.landsat_8_tle).elements[0]
        nearly_identical = base.model_copy(
            update={"mean_motion": base.mean_motion * (1 + 1e-6)}
        )
        orbit = GeneralPerturbationsOrbit(elements=[base, nearly_identical])
        self.assertIsNone(
            orbit.get_repeat_cycle(
                consistency_threshold=timedelta(seconds=1), lazy_load=False
            )
        )
        self.assertIsNotNone(
            orbit.get_repeat_cycle(
                consistency_threshold=timedelta(seconds=2), lazy_load=False
            )
        )


class TestGetObservationEvents(unittest.TestCase):
    """
    Unit tests for GeneralPerturbationsOrbit.get_observation_events.
    Focuses on which of the three strategies (repeat-cycle tiling,
    multi-element partitioning, direct propagation) gets used and when,
    since each one is separately exercised elsewhere (find_events itself
    is Skyfield's; partition_by_element_index has its own tests in
    TestPartitionByElementIndex; get_repeat_cycle has its own tests in
    TestGetRepeatCycle).
    """

    def setUp(self):
        landsat_8_tle = [
            "1 39084U 13008A   26213.27824675  .00000294  00000+0  75333-4 0  9990",
            "2 39084  98.2277 282.8718 0001275  92.4910 267.6434 14.57104473704466",
        ]
        self.base = GeneralPerturbationsOrbit.from_tle(landsat_8_tle).elements[0]
        self.point = Point(id=0, latitude=40.0, longitude=-105.0)

    def test_repeat_cycle_tiles_events_across_cycles(self):
        """
        Test that, when a repeat cycle applies, consecutive cycles'
        events are separated by exactly the repeat cycle duration --
        the structural signature of copy-pasting one cycle's events
        forward, rather than directly propagating the whole period (which
        would not produce exactly period-spaced events, since real
        orbital dynamics drift slightly cycle to cycle). Uses a
        two-element orbit (a near-identical copy, mean motion nudged by
        1e-6) to also confirm this applies even with multiple elements,
        as long as they agree on one repeat cycle.
        """
        nearly_identical = self.base.model_copy(
            update={"mean_motion": self.base.mean_motion * (1 + 1e-6)}
        )
        orbit = GeneralPerturbationsOrbit(elements=[self.base, nearly_identical])
        repeat_cycle = orbit.get_repeat_cycle()
        self.assertIsNotNone(repeat_cycle)
        start = self.base.epoch
        end = start + repeat_cycle * 2 + timedelta(hours=1)
        times, _ = orbit.get_observation_events(
            self.point, start, end, min_elevation_angle=10, try_repeat=True
        )
        utc_times = times.utc_datetime()
        events_per_cycle = int(np.sum(utc_times < start + repeat_cycle))
        self.assertGreater(events_per_cycle, 0)
        self.assertGreater(len(utc_times), events_per_cycle)
        for i in range(len(utc_times) - events_per_cycle):
            self.assertEqual(
                utc_times[i + events_per_cycle] - utc_times[i], repeat_cycle
            )

    def test_repeat_cycle_uses_element_closest_to_start(self):
        """
        Test that the repeat-cycle strategy propagates from whichever
        element is closest to `start`, not always the first (earliest)
        element -- important once a multi-element orbit's repeat cycle is
        trusted regardless of element count. The second element here
        (epoch 20 days after the first, but otherwise the same orbit) is
        closest to a `start` 25 days after the first element's epoch.
        """
        second = self.base.model_copy(
            update={
                "epoch": self.base.epoch + timedelta(days=20),
                "mean_motion": self.base.mean_motion * (1 + 1e-6),
            }
        )
        orbit = GeneralPerturbationsOrbit(elements=[self.base, second])
        self.assertIs(
            orbit.get_closest_element(self.base.epoch + timedelta(days=25)), second
        )
        repeat_cycle = orbit.get_repeat_cycle()
        self.assertIsNotNone(repeat_cycle)
        start = self.base.epoch + timedelta(days=25)
        end = start + repeat_cycle + timedelta(hours=1)

        actual_times, _ = orbit.get_observation_events(
            self.point, start, end, min_elevation_angle=10, try_repeat=True
        )
        topos = wgs84.latlon(
            self.point.latitude, self.point.longitude, self.point.elevation
        )
        t_0 = constants.timescale.from_datetime(start)
        repeat_t_1 = constants.timescale.from_datetime(start + repeat_cycle)
        expected_times, _ = second.to_skyfield().find_events(
            topos, t_0, repeat_t_1, 10
        )
        wrong_times, _ = self.base.to_skyfield().find_events(
            topos, t_0, repeat_t_1, 10
        )
        # sanity check that using the wrong (first) element would actually
        # have given a different answer, so this test is a meaningful
        # discriminator, not a coincidence
        self.assertFalse(
            np.array_equal(expected_times.utc_datetime(), wrong_times.utc_datetime())
        )
        self.assertTrue(
            np.array_equal(
                actual_times.utc_datetime()[: len(expected_times)],
                expected_times.utc_datetime(),
            )
        )

    def test_repeat_cycle_ignored_when_shorter_period_requested(self):
        """
        Test that a period shorter than the repeat cycle skips the
        tiling strategy entirely and falls through to direct propagation
        -- verified by exact equality with try_repeat=False, since both
        end up on the same code path in that case.
        """
        orbit = GeneralPerturbationsOrbit(elements=[self.base])
        repeat_cycle = orbit.get_repeat_cycle()
        self.assertIsNotNone(repeat_cycle)
        start = self.base.epoch
        end = start + timedelta(days=5)
        self.assertLess(end - start, repeat_cycle)
        with_repeat = orbit.get_observation_events(
            self.point, start, end, min_elevation_angle=10, try_repeat=True
        )
        without_repeat = orbit.get_observation_events(
            self.point, start, end, min_elevation_angle=10, try_repeat=False
        )
        self.assertTrue(
            np.array_equal(
                with_repeat[0].utc_datetime(), without_repeat[0].utc_datetime()
            )
        )
        self.assertTrue(np.array_equal(with_repeat[1], without_repeat[1]))

    def test_multi_element_partition_used_when_no_consistent_repeat_cycle(self):
        """
        Test that a multi-element orbit whose elements do not agree on a
        repeat cycle (simulating a maneuver) falls through to per-time
        nearest-element partitioning even with try_repeat=True --
        verified by exact equality with try_repeat=False, since both end
        up on the same partition_by_element_index-based code path when
        get_repeat_cycle() returns None.
        """
        maneuvered = self.base.model_copy(
            update={"mean_motion": self.base.mean_motion * 1.001}
        )
        orbit = GeneralPerturbationsOrbit(elements=[self.base, maneuvered])
        self.assertIsNone(orbit.get_repeat_cycle())
        start = self.base.epoch
        end = start + timedelta(days=5)
        with_repeat = orbit.get_observation_events(
            self.point, start, end, min_elevation_angle=10, try_repeat=True
        )
        without_repeat = orbit.get_observation_events(
            self.point, start, end, min_elevation_angle=10, try_repeat=False
        )
        self.assertGreater(len(with_repeat[1]), 0)
        self.assertTrue(
            np.array_equal(
                with_repeat[0].utc_datetime(), without_repeat[0].utc_datetime()
            )
        )
        self.assertTrue(np.array_equal(with_repeat[1], without_repeat[1]))


class TestToGpOrbit(unittest.TestCase):
    """
    Unit tests for GeneralPerturbationsOrbit.to_gp_orbit.
    """

    def setUp(self):
        tle = [
            "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
            "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
        ]
        self.orbit = GeneralPerturbationsOrbit.from_tle(tle)

    def test_returns_self(self):
        """
        Test that to_gp_orbit() returns this exact instance (identity,
        not merely an equal copy), since a GeneralPerturbationsOrbit
        already is its own general perturbations representation.
        """
        self.assertIs(self.orbit.to_gp_orbit(), self.orbit)

    def test_accepts_lazy_load_without_effect(self):
        """
        Test that to_gp_orbit() accepts the lazy_load parameter (for
        interface parity with OrbitBase.to_gp_orbit, used polymorphically
        across the AllOrbits union) without error, and that it has no
        effect on the result -- still just this instance, regardless of
        the value passed.
        """
        self.assertIs(self.orbit.to_gp_orbit(lazy_load=True), self.orbit)
        self.assertIs(self.orbit.to_gp_orbit(lazy_load=False), self.orbit)


if __name__ == "__main__":
    unittest.main()

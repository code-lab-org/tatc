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


if __name__ == "__main__":
    unittest.main()

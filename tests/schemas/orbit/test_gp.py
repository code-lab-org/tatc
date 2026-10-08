"""
Unit tests for the GeneralPerturbationsOrbit schema.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import copy
import csv
import io
import json
import pickle
import unittest
import warnings
from datetime import datetime, timedelta, timezone
from unittest.mock import patch

import numpy as np
from pydantic import ValidationError
from shapely.geometry import Point as ShapelyPoint
from skyfield.api import wgs84
from skyfield.framelib import itrs

from tatc import config, constants
from tatc.schemas import GeneralPerturbationsOrbit, Point
from tatc.utils.computation import _run
from tatc.utils.propagation import _complete_passes, _find_events

REPEATING = {"remove_drag": True, "repeat_cycle": "auto"}
"""Options to propagate an orbit maintained on its (found) repeat ground track."""


def direct(orbit: GeneralPerturbationsOrbit) -> GeneralPerturbationsOrbit:
    """Gets a copy of an orbit propagated directly (without a repeat cycle)."""
    return orbit.model_copy(update={"repeat_cycle": None})


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

    def test_skyfield_time_input_matches_datetime_input(self):
        """
        Test that a Skyfield time (array or scalar) is accepted and
        selects the same indices as the equivalent datetimes, including
        the tie-breaking convention at a midpoint.
        """
        queries = [
            self.epoch_0 - timedelta(days=365),
            self.epoch_0 + timedelta(hours=1),
            self.epoch_0 + (self.epoch_1 - self.epoch_0) / 2,
            self.epoch_1 - timedelta(hours=1),
            self.epoch_2 + timedelta(days=365),
        ]
        o = self.multi_element_orbit
        self.assertEqual(
            o.get_closest_element_index(constants.timescale.from_datetimes(queries)),
            o.get_closest_element_index(queries),
        )
        self.assertEqual(
            o.get_closest_element_index(constants.timescale.from_datetime(queries[3])),
            1,
        )


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
        self.single_element_orbit = GeneralPerturbationsOrbit(elements=[self.element_0])

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
            expected = element.to_skyfield().at(constants.timescale.from_datetime(time))
            # within a millimeter (nutation angles are interpolated in
            # propagation, see _interpolate_nutation), whereas another element
            # differs by kilometers
            self.assertTrue(
                np.allclose(
                    actual.position.km[:, i], expected.position.km, rtol=0, atol=1e-6
                ),
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
        self.repeat_orbit = GeneralPerturbationsOrbit.from_tle(
            landsat_8_tle, **REPEATING
        )
        self.epoch = self.repeat_orbit.get_epoch()
        self.repeat_cycle = self.repeat_orbit.get_repeat_cycle()
        self.assertIsNotNone(self.repeat_cycle)

        # ISS: single-element, not designed for a repeat ground track
        iss_tle = [
            "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
            "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
        ]
        self.non_repeat_orbit = GeneralPerturbationsOrbit.from_tle(iss_tle, **REPEATING)

        # a multi-element orbit (elements otherwise identical to Landsat-8,
        # just re-epoched), to confirm the substitution never applies
        base = self.repeat_orbit.elements[0]
        self.multi_element_orbit = GeneralPerturbationsOrbit(
            elements=[
                base.model_copy(update={"epoch": self.epoch}),
                base.model_copy(update={"epoch": self.epoch + timedelta(days=1)}),
            ],
            **REPEATING,
        )

    def test_default_orbit_matches_direct_propagation(self):
        """
        Test that an orbit with the default options (`repeat_cycle=None`,
        with drag) always returns the true, directly propagated position,
        even if its element has a repeat cycle.
        """
        t = constants.timescale.from_datetime(self.epoch + self.repeat_cycle * 3)
        default = GeneralPerturbationsOrbit(elements=self.repeat_orbit.elements)
        self.assertIsNone(default.repeat_cycle)
        self.assertFalse(default.remove_drag)
        # propagated with drag more than 30 days from the epoch
        with self.assertWarns(UserWarning):
            actual = default.get_geographic_position_at_time(t)
        expected = wgs84.geographic_position_of(
            self.repeat_orbit.elements[0].to_skyfield().at(t)
        )
        self.assertEqual(actual.latitude.degrees, expected.latitude.degrees)
        self.assertEqual(actual.longitude.degrees, expected.longitude.degrees)

    def test_repeat_within_first_cycle_uses_maintained_element(self):
        """
        Test that, for a time less than one repeat cycle after epoch, the
        time is not shifted, but propagated with the element maintained on
        its repeat ground track (see get_repeat_element), which stays close
        to the element's own (direct) propagation near the epoch.
        """
        t = constants.timescale.from_datetime(self.epoch + timedelta(hours=5))
        with_repeat = self.repeat_orbit.get_geographic_position_at_time(t)
        maintained = wgs84.geographic_position_of(
            self.repeat_orbit.get_repeat_element().to_skyfield().at(t)
        )
        unrepeated = direct(self.repeat_orbit).get_geographic_position_at_time(t)
        self.assertEqual(with_repeat.latitude.degrees, maintained.latitude.degrees)
        self.assertEqual(with_repeat.longitude.degrees, maintained.longitude.degrees)
        self.assertLess(
            np.linalg.norm(np.array(with_repeat.itrs_xyz.m) - unrepeated.itrs_xyz.m),
            5e3,
        )

    def test_repeat_continuous_across_cycles(self):
        """
        Test that the repeated orbit track joins at the ends of each repeat
        cycle: the element maintained on its repeat ground track returns
        close to its initial position (within about 1 km, from the eccentricity
        terms as the argument of perigee precesses), whereas the element's
        own mean motion, slightly off the exact repeat, would leave a jump
        of tens of kilometers (for this Landsat 8 element set, 34 km).
        """
        d = timedelta(seconds=0.5)

        def position(time):
            return np.array(
                self.repeat_orbit.get_geographic_position_at_time(
                    constants.timescale.from_datetime(time)
                ).itrs_xyz.m
            )

        for boundary in (
            self.epoch + self.repeat_cycle,
            self.epoch - self.repeat_cycle,
        ):
            step = np.linalg.norm(position(boundary + d) - position(boundary - d))
            later = boundary + timedelta(seconds=10)
            normal = np.linalg.norm(position(later + d) - position(later - d))
            self.assertLess(step - normal, 1.5e3)
        own = self.repeat_orbit.elements[0].without_drag()
        jump = np.linalg.norm(
            own.to_skyfield()
            .at(constants.timescale.from_datetime(self.epoch + self.repeat_cycle))
            .frame_xyz(itrs)
            .m
            - own.to_skyfield()
            .at(constants.timescale.from_datetime(self.epoch))
            .frame_xyz(itrs)
            .m
        )
        self.assertGreater(jump, 20e3)

    def test_get_repeat_element(self):
        """
        Test that the maintained element is the element adjusted to the
        exact repeat (cached), for each element of a multi-element orbit,
        and None without a repeat cycle.
        """
        element = self.repeat_orbit.get_repeat_element()
        self.assertEqual(
            element,
            self.repeat_orbit.elements[0].get_repeat_element(self.repeat_cycle),
        )
        self.assertIs(self.repeat_orbit.get_repeat_element(), element)
        self.assertEqual(element.bstar, 0)
        self.assertEqual(
            self.multi_element_orbit.get_repeat_element(1),
            self.multi_element_orbit.elements[1].get_repeat_element(self.repeat_cycle),
        )
        no_repeat = self.repeat_orbit.model_copy(
            update={
                "elements": [
                    # a realistic orbit without a repeat cycle within the
                    # search duration (mean motion is in radians per minute)
                    self.repeat_orbit.elements[0].model_copy(
                        update={
                            "mean_motion": self.repeat_orbit.elements[0].mean_motion
                            * 1.002
                        }
                    )
                ]
            }
        )
        self.assertIsNone(no_repeat.get_repeat_element())

    def test_repeat_substitutes_epoch_relative_time_for_far_future(self):
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
        actual = self.repeat_orbit.get_geographic_position_at_time(t)

        offset_days = (far_future - self.epoch) / timedelta(days=1)
        cycle_days = self.repeat_cycle / timedelta(days=1)
        wrapped_offset = timedelta(days=float(np.mod(offset_days, cycle_days)))
        expected = wgs84.geographic_position_of(
            self.repeat_orbit.get_repeat_element()
            .to_skyfield()
            .at(constants.timescale.from_datetime(self.epoch + wrapped_offset))
        )
        self.assertAlmostEqual(
            actual.latitude.degrees, expected.latitude.degrees, places=9
        )
        self.assertAlmostEqual(
            actual.longitude.degrees, expected.longitude.degrees, places=9
        )

        unrepeated = wgs84.geographic_position_of(
            direct(self.repeat_orbit).get_orbit_track_at_time(t)
        )
        self.assertNotEqual(actual.latitude.degrees, unrepeated.latitude.degrees)

    def test_repeat_preserves_sign_for_time_before_epoch(self):
        """
        Regression test: a query 2.5 repeat cycles *before* epoch must
        wrap to -0.5 cycles (epoch minus half a cycle), not +0.5 cycles.
        A naive numpy np.mod() on the raw (negative) offset always
        returns a non-negative result, which would silently wrap to the
        wrong side of the repeat cycle. Since the maintained element repeats
        (nearly) exactly, the two sides give (nearly) the same position
        (within about 1 km here).
        """
        query_time = self.epoch - self.repeat_cycle * 2.5
        t = constants.timescale.from_datetime(query_time)
        actual = self.repeat_orbit.get_geographic_position_at_time(t)

        correct = wgs84.geographic_position_of(
            self.repeat_orbit.get_repeat_element()
            .to_skyfield()
            .at(constants.timescale.from_datetime(self.epoch - self.repeat_cycle * 0.5))
        )
        wrong = wgs84.geographic_position_of(
            self.repeat_orbit.get_repeat_element()
            .to_skyfield()
            .at(constants.timescale.from_datetime(self.epoch + self.repeat_cycle * 0.5))
        )
        self.assertAlmostEqual(
            actual.latitude.degrees, correct.latitude.degrees, places=9
        )
        self.assertLess(
            np.linalg.norm(np.array(correct.itrs_xyz.m) - wrong.itrs_xyz.m), 5e3
        )

    def test_repeat_multi_element_orbit(self):
        """
        Test that a multi-element orbit is repeated with its first element
        before the first epoch and with its last element after the last
        epoch (as single-element orbits of those elements are), and
        propagated directly with the closest element between them.
        """
        first, last = self.multi_element_orbit.elements
        times = [
            self.epoch - timedelta(days=40),
            self.epoch + timedelta(hours=6),
            self.epoch + timedelta(days=41),
        ]
        expected_orbits = [
            GeneralPerturbationsOrbit(elements=[first], **REPEATING),
            None,
            GeneralPerturbationsOrbit(elements=[last], **REPEATING),
        ]
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            actual = self.multi_element_orbit.get_geographic_position(times)
        for i, (time, orbit) in enumerate(zip(times, expected_orbits)):
            t = constants.timescale.from_datetime(time)
            expected = (
                wgs84.geographic_position_of(first.to_skyfield(remove_drag=True).at(t))
                if orbit is None
                else orbit.get_geographic_position_at_time(t)
            )
            np.testing.assert_allclose(
                np.array(actual.itrs_xyz.m)[:, i], expected.itrs_xyz.m, atol=1e-3
            )

    def test_falls_back_when_no_repeat_cycle_found(self):
        """
        Test that an orbit falls back to direct propagation (without drag)
        if no repeat cycle is found (the ISS).
        """
        t = constants.timescale.from_datetime(
            self.non_repeat_orbit.get_epoch() + timedelta(days=10)
        )
        actual = self.non_repeat_orbit.get_geographic_position_at_time(t)
        expected = wgs84.geographic_position_of(
            self.non_repeat_orbit.elements[0].to_skyfield(remove_drag=True).at(t)
        )
        self.assertEqual(actual.latitude.degrees, expected.latitude.degrees)
        self.assertEqual(actual.longitude.degrees, expected.longitude.degrees)

    def test_repeat_cycle_field_selects_propagation(self):
        """
        Test that the `repeat_cycle` field selects how the orbit is
        propagated: "auto" and a declared repeat cycle repeat the orbit
        track, while None (the default) propagates the element directly.
        """
        t = constants.timescale.from_datetime(
            self.epoch + self.repeat_cycle * 3 + timedelta(hours=5)
        )
        self.assertEqual(self.repeat_orbit.repeat_cycle, "auto")
        repeated = self.repeat_orbit.get_geographic_position_at_time(t)
        declared = self.repeat_orbit.model_copy(
            update={"repeat_cycle": timedelta(days=16)}
        ).get_geographic_position_at_time(t)
        direct_position = direct(self.repeat_orbit).get_geographic_position_at_time(t)
        expected = wgs84.geographic_position_of(
            self.repeat_orbit.elements[0].to_skyfield(remove_drag=True).at(t)
        )
        np.testing.assert_allclose(repeated.itrs_xyz.m, declared.itrs_xyz.m, atol=1e-3)
        np.testing.assert_allclose(direct_position.itrs_xyz.m, expected.itrs_xyz.m)
        self.assertGreater(
            np.linalg.norm(np.array(repeated.itrs_xyz.m) - direct_position.itrs_xyz.m),
            10e3,
        )
        self.assertIsNone(direct(self.repeat_orbit).get_repeat_cycle())
        self.assertIsNone(direct(self.repeat_orbit).get_repeat_element())

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
        actual = self.repeat_orbit.get_geographic_position_at_time(t_vector)
        for i, query_time in enumerate(query_times):
            expected = self.repeat_orbit.get_geographic_position_at_time(
                constants.timescale.from_datetime(query_time)
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
        start = self.epoch_0 - timedelta(days=365)
        end = self.epoch_2 + timedelta(days=365)
        with patch.object(
            GeneralPerturbationsOrbit,
            "get_element_epochs",
            autospec=True,
            side_effect=GeneralPerturbationsOrbit.get_element_epochs,
        ) as get_element_epochs:
            self.multi_element_orbit.partition_by_element_index(start, end)
        get_element_epochs.assert_called()
        cached = self.multi_element_orbit.get_element_epochs()
        self.assertEqual(cached, [self.epoch_0, self.epoch_1, self.epoch_2])
        self.assertIs(self.multi_element_orbit.get_element_epochs(), cached)

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
        _, events = direct(self.multi_element_orbit).get_observation_events(
            point, start, end, min_elevation_angle=10
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
        orbit = GeneralPerturbationsOrbit.from_tle(self.landsat_8_tle, **REPEATING)
        repeat_cycle = orbit.get_repeat_cycle()
        self.assertIsNotNone(repeat_cycle)
        self.assertAlmostEqual(repeat_cycle.total_seconds() / 86400, 16, delta=0.1)

    def test_sentinel_2a_repeat_cycle_matches_published_10_days(self):
        """
        Test against Sentinel-2A's published repeat ground track of 143
        orbits every 10 days (ESA), at a different altitude/inclination
        than Landsat-8, using the default tolerances.
        """
        orbit = GeneralPerturbationsOrbit.from_tle(self.sentinel_2a_tle, **REPEATING)
        repeat_cycle = orbit.get_repeat_cycle()
        self.assertIsNotNone(repeat_cycle)
        self.assertAlmostEqual(repeat_cycle.total_seconds() / 86400, 10, delta=0.1)

    def test_sentinel_1a_repeat_cycle_matches_published_12_days(self):
        """
        Test against Sentinel-1A's published repeat ground track of 175
        orbits every 12 days (ESA), at a different altitude/inclination
        than Landsat-8, using the default tolerances.
        """
        orbit = GeneralPerturbationsOrbit.from_tle(self.sentinel_1a_tle, **REPEATING)
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
        orbit = GeneralPerturbationsOrbit.from_tle(self.gps_tle, **REPEATING)
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
        orbit = GeneralPerturbationsOrbit.from_tle(self.molniya_tle, **REPEATING)
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

    def test_semimajor_axis_tolerance_confirms_element_in_maintenance_band(self):
        """
        Test that an element whose semimajor axis is 90 m from Landsat-8's
        (within a maintenance band such as Landsat 9's 200 m), whose ground
        track drifts beyond the position tolerance over 16 days, is
        confirmed by the semimajor axis tolerance, but not without it.
        """
        base = GeneralPerturbationsOrbit.from_tle(self.landsat_8_tle).elements[0]
        a = base.get_semimajor_axis()
        offset = base.model_copy(
            update={"mean_motion": base.mean_motion * ((a + 90) / a) ** -1.5}
        )
        repeat_cycle = offset.get_repeat_cycle()
        self.assertIsNotNone(repeat_cycle)
        self.assertAlmostEqual(repeat_cycle.total_seconds() / 86400, 16, delta=0.1)
        self.assertIsNone(offset.get_repeat_cycle(max_delta_semimajor_axis=0))

    def test_long_search_does_not_admit_chance_repeat(self):
        """
        Test that a search of 100 days finds no repeat cycle for an ICESat-2
        element set, whose semimajor axis is 67 m from its 91-day repeat but
        within 14 m of a chance repeat after 62 days: the semimajor axis
        tolerance does not apply to such long repeat cycles.
        """
        orbit = GeneralPerturbationsOrbit.from_tle(
            [
                "1 43613U 18070A   26277.57046165  .00001587  00000+0  57409-4 0  9990",
                "2 43613  91.9994 320.5065 0005609  71.5392 288.6467 15.28297940449160",
            ],
            **REPEATING,
        )
        self.assertIsNone(
            orbit.get_repeat_cycle(max_search_duration=timedelta(days=100))
        )

    def test_cache_follows_runtime_configuration(self):
        """
        Test that the orbit's cached repeat cycle is recomputed when the
        runtime configuration of the search options changes (the defaults
        are resolved before the cache key is built).
        """
        orbit = GeneralPerturbationsOrbit.from_tle(
            [
                "1 39084U 13008A   26213.27824675  .00000294  00000+0  75333-4 0  9990",
                "2 39084  98.2277 282.8718 0001275  92.4910 267.6434 14.57104473704466",
            ],
            **REPEATING,
        )
        self.assertIsNotNone(orbit.get_repeat_cycle())
        rc = config.get_rc()
        original = (
            rc.repeat_cycle_delta_position_m,
            rc.repeat_cycle_delta_semimajor_axis_m,
        )
        try:
            rc.repeat_cycle_delta_position_m = 1
            rc.repeat_cycle_delta_semimajor_axis_m = 0
            self.assertIsNone(orbit.get_repeat_cycle())
        finally:
            rc.repeat_cycle_delta_position_m, rc.repeat_cycle_delta_semimajor_axis_m = (
                original
            )
        self.assertIsNotNone(orbit.get_repeat_cycle())

    def test_cache_follows_search_options(self):
        """
        Test that a cached repeat cycle is not reused for a search with
        other options.
        """
        orbit = GeneralPerturbationsOrbit.from_tle(self.landsat_8_tle, **REPEATING)
        self.assertIsNotNone(orbit.get_repeat_cycle())
        self.assertIsNone(
            orbit.get_repeat_cycle(max_search_duration=timedelta(days=10))
        )
        self.assertIsNone(
            orbit.elements[0].get_repeat_cycle(max_search_duration=timedelta(days=10))
        )

    def test_reuses_cached_result(self):
        """
        Test that calling get_repeat_cycle() twice returns the identical
        cached timedelta rather than recomputing.
        """
        orbit = GeneralPerturbationsOrbit.from_tle(self.landsat_8_tle, **REPEATING)
        first = orbit.get_repeat_cycle()
        second = orbit.get_repeat_cycle()
        self.assertIs(first, second)

    def test_too_short_search_duration_returns_none(self):
        """
        Test that a max_search_duration shorter than the true repeat
        cycle (16 days) cannot find it.
        """
        orbit = GeneralPerturbationsOrbit.from_tle(self.landsat_8_tle, **REPEATING)
        repeat_cycle = orbit.get_repeat_cycle(max_search_duration=timedelta(days=10))
        self.assertIsNone(repeat_cycle)

    def test_too_tight_tolerance_returns_none(self):
        """
        Test that an unrealistically tight position/velocity tolerance
        rejects even the real repeat cycle.
        """
        orbit = GeneralPerturbationsOrbit.from_tle(self.landsat_8_tle, **REPEATING)
        repeat_cycle = orbit.get_repeat_cycle(
            max_delta_position=1, max_delta_velocity=0.001
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
        orbit = GeneralPerturbationsOrbit(
            elements=[base, nearly_identical], **REPEATING
        )
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
        maneuvered = base.model_copy(update={"mean_motion": base.mean_motion * 1.0005})
        orbit = GeneralPerturbationsOrbit(elements=[base, maneuvered], **REPEATING)
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
        no repeat cycle at all (here, a 0.2% mean motion change from the
        real Landsat-8 element, which at default tolerances finds no
        commensurate day within the search duration), even though the
        other element's own repeat cycle is perfectly valid.
        """
        base = GeneralPerturbationsOrbit.from_tle(self.landsat_8_tle).elements[0]
        maneuvered = base.model_copy(update={"mean_motion": base.mean_motion * 1.002})
        orbit = GeneralPerturbationsOrbit(elements=[base, maneuvered], **REPEATING)
        self.assertIsNotNone(base.get_repeat_cycle())
        self.assertIsNone(maneuvered.get_repeat_cycle())
        self.assertIsNone(orbit.get_repeat_cycle())

    def test_consistency_threshold_override_changes_accept_reject_boundary(self):
        """
        Test that consistency_threshold controls the accept/reject
        boundary directly: the same pair of elements of a
        non-sun-synchronous orbit (SWOT), whose inclinations differ by
        0.002 degrees so that their repeat cycles (whole numbers of nodal
        days) differ by about 1.7 seconds, is rejected under a threshold
        below that difference and accepted under one above it.
        """
        base = GeneralPerturbationsOrbit.from_tle(
            [
                "1 54754U 22173A   26277.59902338  .00000094  00000+0  65386-4 0  9997",
                "2 54754  77.6084 247.6064 0000264 143.3265 216.7905 14.00173063194463",
            ]
        ).elements[0]
        nearly_identical = base.model_copy(
            update={"inclination": base.inclination + 0.002}
        )
        difference = abs(base.get_repeat_cycle() - nearly_identical.get_repeat_cycle())
        self.assertGreater(difference, timedelta(seconds=1))
        orbit = GeneralPerturbationsOrbit(
            elements=[base, nearly_identical], **REPEATING
        )
        self.assertIsNone(orbit.get_repeat_cycle(consistency_threshold=difference / 2))
        self.assertIsNotNone(
            orbit.get_repeat_cycle(consistency_threshold=difference * 2)
        )


class TestGetObservationEvents(unittest.TestCase):
    """
    Unit tests for GeneralPerturbationsOrbit.get_observation_events.
    Focuses on which strategy (repeat tracks before the first and after the
    last element's epoch, or direct propagation partitioned by the closest
    element) is used and when, since each one is separately exercised
    elsewhere (find_events itself is Skyfield's; partition_by_element_index
    has its own tests in TestPartitionByElementIndex; get_repeat_cycle has
    its own tests in TestGetRepeatCycle).
    """

    def setUp(self):
        landsat_8_tle = [
            "1 39084U 13008A   26213.27824675  .00000294  00000+0  75333-4 0  9990",
            "2 39084  98.2277 282.8718 0001275  92.4910 267.6434 14.57104473704466",
        ]
        self.base = GeneralPerturbationsOrbit.from_tle(landsat_8_tle).elements[0]
        self.point = Point(id=0, latitude=40.0, longitude=-105.0)

    def test_shapely_point_matches_point(self):
        """
        Test that the events for a shapely point (longitude, latitude, and
        elevation) are those for the equivalent TAT-C point.
        """
        orbit = GeneralPerturbationsOrbit(elements=[self.base])
        start = self.base.epoch
        end = start + timedelta(days=1)
        for point, shapely_point in [
            (self.point, ShapelyPoint(-105.0, 40.0)),
            (
                Point(latitude=40.0, longitude=-105.0, elevation=1600),
                ShapelyPoint(-105.0, 40.0, 1600),
            ),
        ]:
            expected_times, expected_codes = orbit.get_observation_events(
                point, start, end, min_elevation_angle=10
            )
            times, codes = orbit.get_observation_events(
                shapely_point, start, end, min_elevation_angle=10
            )
            self.assertGreater(len(codes), 0)
            np.testing.assert_array_equal(codes, expected_codes)
            np.testing.assert_array_equal(times.tt, expected_times.tt)

    def test_repeat_cycle_tiles_events_across_cycles(self):
        """
        Test that, when a repeat cycle applies, consecutive cycles'
        events are separated by exactly the repeat cycle duration --
        the structural signature of copy-pasting one cycle's events
        forward, rather than directly propagating the whole period (which
        would not produce exactly period-spaced events, since real
        orbital dynamics drift slightly cycle to cycle).
        """
        orbit = GeneralPerturbationsOrbit(elements=[self.base], **REPEATING)
        repeat_cycle = orbit.get_repeat_cycle()
        self.assertIsNotNone(repeat_cycle)
        start = self.base.epoch
        end = start + repeat_cycle * 2 + timedelta(hours=1)
        times, _ = orbit.get_observation_events(
            self.point, start, end, min_elevation_angle=10
        )
        utc_times = times.utc_datetime()
        events_per_cycle = int(np.sum(utc_times < start + repeat_cycle))
        self.assertGreater(events_per_cycle, 0)
        self.assertGreater(len(utc_times), events_per_cycle)
        for i in range(len(utc_times) - events_per_cycle):
            self.assertEqual(
                utc_times[i + events_per_cycle] - utc_times[i], repeat_cycle
            )

    def test_multi_element_orbit_repeats_first_and_last_elements(self):
        """
        Test that a multi-element orbit (such as a history of element sets)
        is repeated with its first element before the first epoch and with
        its last element after the last epoch, as single-element orbits of
        those elements are, and propagated directly with the closest element
        between the epochs, as without a repeat cycle.
        """
        second = self.base.model_copy(
            update={
                "epoch": self.base.epoch + timedelta(days=20),
                "mean_motion": self.base.mean_motion * (1 + 1e-6),
            }
        )
        orbit = GeneralPerturbationsOrbit(elements=[self.base, second], **REPEATING)
        first_epoch, last_epoch = self.base.epoch, second.epoch
        start, end = first_epoch - timedelta(days=40), last_epoch + timedelta(days=40)
        times, codes = orbit.get_observation_events(
            self.point, start, end, min_elevation_angle=10
        )
        times = np.array(times.utc_datetime())
        margin = timedelta(hours=1)
        for reference, lower, upper in [
            (
                GeneralPerturbationsOrbit(elements=[self.base], **REPEATING),
                start,
                first_epoch - margin,
            ),
            (direct(orbit), first_epoch + margin, last_epoch - margin),
            (
                GeneralPerturbationsOrbit(elements=[second], **REPEATING),
                last_epoch + margin,
                end,
            ),
        ]:
            expected_times, expected_codes = reference.get_observation_events(
                self.point, lower, upper, min_elevation_angle=10
            )
            selected = (times >= lower) & (times <= upper)
            self.assertGreater(len(expected_codes), 0)
            self.assertEqual(codes[selected].tolist(), expected_codes.tolist())
            differences = np.array(
                [
                    (actual - expected).total_seconds()
                    for actual, expected in zip(
                        times[selected], expected_times.utc_datetime()
                    )
                ]
            )
            np.testing.assert_allclose(differences[expected_codes != 1], 0, atol=1e-2)
            np.testing.assert_allclose(differences[expected_codes == 1], 0, atol=1)

    def _expected_repeated_events(self, orbit, start, end):
        """
        Expected repeated events: the events of the repeat cycle just after
        (or before) the epoch, directly propagated with the element
        maintained on its repeat ground track, shifted by each whole number
        of repeat cycles that maps them into the period on the same side of
        the epoch.
        """
        epoch, repeat_cycle = orbit.get_epoch(), orbit.get_repeat_cycle()
        maintained = GeneralPerturbationsOrbit(
            elements=[orbit.get_repeat_element()], **REPEATING
        )
        expected = []
        for after in (True, False):
            window = (
                (epoch, epoch + repeat_cycle)
                if after
                else (epoch - repeat_cycle, epoch)
            )
            times, codes = direct(maintained).get_observation_events(
                self.point, *window, min_elevation_angle=10
            )
            for time, code in zip(times.utc_datetime(), codes):
                for cycles in range(-50, 51):
                    shifted = time + cycles * repeat_cycle
                    if (
                        start <= shifted <= end
                        and (shifted >= epoch) == after
                        and int(np.trunc((shifted - epoch) / repeat_cycle)) == cycles
                    ):
                        expected.append((shifted, int(code)))
        return sorted(set(expected))

    def _assert_repeated_events(self, orbit, start, end, times, codes):
        """
        Asserts that events match the expected repeated events, with rise
        and set times within 10 ms (they are refined to a millisecond) and
        culmination times within 1 s (they are not refined), as the times
        depend slightly on the period over which events are found.
        """
        expected = self._expected_repeated_events(orbit, start, end)
        self.assertGreater(len(expected), 0)
        self.assertEqual(codes.tolist(), [code for _, code in expected])
        differences = np.array(
            [
                (actual - time).total_seconds()
                for actual, (time, _) in zip(times.utc_datetime(), expected)
            ]
        )
        np.testing.assert_allclose(differences[codes != 1], 0, atol=1e-2)
        np.testing.assert_allclose(differences[codes == 1], 0, atol=1)

    def test_repeat_cycle_anchored_at_epoch(self):
        """
        Test that repeated events are those of the repeat cycle just after
        the element's epoch (not of the cycle starting at `start`), shifted
        by whole repeat cycles, for a period that starts ten cycles after
        the epoch: the element is never propagated more than one cycle
        from its epoch.
        """
        orbit = GeneralPerturbationsOrbit(elements=[self.base], **REPEATING)
        repeat_cycle = orbit.get_repeat_cycle()
        start = self.base.epoch + 10 * repeat_cycle + timedelta(hours=3)
        end = start + repeat_cycle
        times, codes = orbit.get_observation_events(
            self.point, start, end, min_elevation_angle=10
        )
        self._assert_repeated_events(orbit, start, end, times, codes)

    def test_repeat_cycle_anchored_at_epoch_before_epoch(self):
        """
        Test that, before the element's epoch, repeated events are those of
        the repeat cycle just before the epoch, shifted by whole repeat
        cycles, for a period that spans the epoch.
        """
        orbit = GeneralPerturbationsOrbit(elements=[self.base], **REPEATING)
        repeat_cycle = orbit.get_repeat_cycle()
        start = self.base.epoch - 3 * repeat_cycle - timedelta(hours=5)
        end = self.base.epoch + timedelta(days=2)
        times, codes = orbit.get_observation_events(
            self.point, start, end, min_elevation_angle=10
        )
        self._assert_repeated_events(orbit, start, end, times, codes)

    def _get_overhead_point(self, orbit, time):
        """Gets the point beneath an orbit's (directly propagated) track at a time."""
        subpoint = wgs84.subpoint_of(direct(orbit).get_orbit_track(time))
        return Point(
            id=0,
            latitude=float(subpoint.latitude.degrees),
            longitude=float(subpoint.longitude.degrees),
        )

    def test_find_events_pass_culminating_outside_period(self):
        """
        Test that the set of a pass that culminates before the start of the
        period, and the rise of one that culminates after its end, are found
        (Skyfield's find_events finds neither).
        """
        orbit = GeneralPerturbationsOrbit(elements=[self.base], **REPEATING)
        culmination = self.base.epoch + timedelta(hours=1)
        point = self._get_overhead_point(orbit, culmination)
        topos = wgs84.latlon(point.latitude, point.longitude)
        satellite = self.base.to_skyfield()
        times, codes = _run(
            _find_events(
                satellite,
                topos,
                constants.timescale.from_datetime(culmination - timedelta(minutes=30)),
                constants.timescale.from_datetime(culmination + timedelta(minutes=30)),
                80,
            )
        )
        self.assertEqual(codes.tolist(), [0, 1, 2])
        rise, _, set_ = times.utc_datetime()
        after = _run(
            _find_events(
                satellite,
                topos,
                constants.timescale.from_datetime(culmination + timedelta(seconds=2)),
                constants.timescale.from_datetime(culmination + timedelta(minutes=30)),
                80,
            )
        )
        self.assertEqual(after[1].tolist(), [2])
        self.assertLess(abs((after[0].utc_datetime()[0] - set_).total_seconds()), 1e-2)
        before = _run(
            _find_events(
                satellite,
                topos,
                constants.timescale.from_datetime(culmination - timedelta(minutes=30)),
                constants.timescale.from_datetime(culmination - timedelta(seconds=2)),
                80,
            )
        )
        self.assertEqual(before[1].tolist(), [0])
        self.assertLess(abs((before[0].utc_datetime()[0] - rise).total_seconds()), 1e-2)

    def test_multi_element_pass_spanning_element_switch(self):
        """
        Test that a pass spanning the switch between two elements (here,
        identical, so that the orbit track is continuous) has the same
        events as with a single element, although each element's period
        contains only part of the pass.
        """
        single = GeneralPerturbationsOrbit(elements=[self.base], **REPEATING)
        multi = GeneralPerturbationsOrbit(
            elements=[self.base, self.base.model_copy()], **REPEATING
        )
        point = self._get_overhead_point(
            single, self.base.epoch + timedelta(seconds=10)
        )
        start = self.base.epoch - timedelta(minutes=30)
        end = self.base.epoch + timedelta(minutes=30)
        expected = direct(single).get_observation_events(point, start, end, 80)
        actual = direct(multi).get_observation_events(point, start, end, 80)
        self.assertEqual(expected[1].tolist(), [0, 1, 2])
        self.assertEqual(actual[1].tolist(), expected[1].tolist())
        np.testing.assert_allclose(actual[0].tt, expected[0].tt, atol=1e-2 / 86400)

    def test_multi_element_pass_ends_at_element_switch(self):
        """
        Test that a pass ends at the switch between two elements if the
        satellite is above the minimum elevation angle with the first
        element but not with the second (here, ahead along the orbit).
        """
        ahead = self.base.model_copy(
            update={"mean_anomaly": (self.base.mean_anomaly + 3) % 360}
        )
        orbit = GeneralPerturbationsOrbit(elements=[self.base, ahead], **REPEATING)
        point = self._get_overhead_point(
            GeneralPerturbationsOrbit(elements=[self.base], **REPEATING),
            self.base.epoch,
        )
        start = self.base.epoch - timedelta(minutes=30)
        end = self.base.epoch + timedelta(minutes=30)
        times, codes = direct(orbit).get_observation_events(point, start, end, 80)
        self.assertEqual(codes.tolist()[-1], 2)
        self.assertIn(0, codes.tolist()[:-1])
        self.assertLess(
            abs((times.utc_datetime()[-1] - self.base.epoch).total_seconds()), 1e-3
        )

    def test_warns_when_propagating_with_drag_far_from_epoch(self):
        """
        Test that propagating an element directly with drag far from its
        epoch (the ISS, with the default options) warns, unless it is
        propagated without drag (including when no repeat cycle is found) or
        near its epoch.
        """
        iss_tle = [
            "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
            "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
        ]
        orbit = GeneralPerturbationsOrbit.from_tle(iss_tle)
        start = orbit.get_epoch() + timedelta(days=60)
        end = start + timedelta(hours=6)
        with self.assertWarns(UserWarning):
            orbit.get_observation_events(self.point, start, end, 10)
        with self.assertWarns(UserWarning):
            orbit.get_geographic_position(start)
        epoch = orbit.get_epoch()
        for quiet_orbit, period in [
            (GeneralPerturbationsOrbit.from_tle(iss_tle, **REPEATING), (start, end)),
            (
                GeneralPerturbationsOrbit.from_tle(iss_tle, remove_drag=True),
                (start, end),
            ),
            (orbit, (epoch, epoch + timedelta(days=1))),
        ]:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                quiet_orbit.get_observation_events(self.point, *period, 10)
            self.assertFalse(
                any("with drag" in str(w.message) for w in caught),
                (period, quiet_orbit.repeat_cycle),
            )

    def test_repeat_cycle_used_within_first_cycle(self):
        """
        Test that a period within one repeat cycle of the epoch is also
        propagated with the element maintained on its repeat ground track,
        so that its events do not depend on the length of the period.
        """
        orbit = GeneralPerturbationsOrbit(elements=[self.base], **REPEATING)
        repeat_cycle = orbit.get_repeat_cycle()
        self.assertIsNotNone(repeat_cycle)
        start = self.base.epoch + timedelta(days=1)
        end = start + timedelta(days=5)
        times, codes = orbit.get_observation_events(
            self.point, start, end, min_elevation_angle=10
        )
        self._assert_repeated_events(orbit, start, end, times, codes)
        longer = orbit.get_observation_events(
            self.point, start, end + 2 * repeat_cycle, 10
        )
        self.assertEqual(longer[1][: len(codes)].tolist(), codes.tolist())
        differences = np.array(
            [
                (a - b).total_seconds()
                for a, b in zip(longer[0].utc_datetime(), times.utc_datetime())
            ]
        )
        np.testing.assert_allclose(differences[codes != 1], 0, atol=1e-2)
        np.testing.assert_allclose(differences[codes == 1], 0, atol=1)

    def test_multi_element_partition_used_between_epochs(self):
        """
        Test that a multi-element orbit whose elements do not agree on a
        repeat cycle (simulating a maneuver) uses per-time nearest-element
        partitioning between its first and last epochs -- verified by exact
        equality with direct propagation (`repeat_cycle=None`), since both
        use the same partition_by_element_index-based code path.
        """
        maneuvered = self.base.model_copy(
            update={
                "epoch": self.base.epoch + timedelta(days=5),
                "mean_motion": self.base.mean_motion * 1.002,
            }
        )
        orbit = GeneralPerturbationsOrbit(elements=[self.base, maneuvered], **REPEATING)
        self.assertIsNone(orbit.get_repeat_cycle())
        start = self.base.epoch
        end = start + timedelta(days=5)
        with_repeat = orbit.get_observation_events(
            self.point, start, end, min_elevation_angle=10
        )
        without_repeat = direct(orbit).get_observation_events(
            self.point, start, end, min_elevation_angle=10
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


if __name__ == "__main__":
    unittest.main()

    def test_rise_and_set_times_are_refined(self):
        """
        Test that rise and set times are where the elevation angle equals
        the minimum, to within a millisecond, and do not depend on the
        length of the requested period: Skyfield's `find_events` alone
        returns some of them seconds early or late over long periods,
        because it stops refining once the first of its unequal brackets
        converges.
        """
        orbit = GeneralPerturbationsOrbit(elements=[self.base])
        topos = wgs84.latlon(self.point.latitude, self.point.longitude)
        start = self.base.epoch
        times, events = direct(orbit).get_observation_events(
            self.point, start, start + timedelta(days=2), 35
        )
        self.assertGreater(np.sum(events == 0), 1)
        satellite = self.base.to_skyfield()
        for time, event in zip(times.utc_datetime(), events):
            if event == 1:
                continue
            ts = constants.timescale.from_datetimes(
                [time - timedelta(milliseconds=1), time + timedelta(milliseconds=1)]
            )
            elevation = (satellite - topos).at(ts).altaz()[0].degrees
            # the elevation angle crosses the minimum within the millisecond
            self.assertLess((elevation[0] - 35) * (elevation[1] - 35), 0)
            short_times, short_events = direct(orbit).get_observation_events(
                self.point,
                time - timedelta(minutes=10),
                time + timedelta(minutes=10),
                35,
            )
            matching = short_times.utc_datetime()[short_events == event]
            self.assertLess(
                min(abs((t - time).total_seconds()) for t in matching), 1e-3
            )


class TestRemoveDragAndRepeatCycle(unittest.TestCase):
    """
    Unit tests for removing the drag terms of a GP orbit and declaring the
    repeat cycle of a maintained orbit.
    """

    def setUp(self):
        # ICESat-2 (NORAD 43613) is maintained on a 91-day repeat ground
        # track (1,387 orbits), which is too long to be found by the
        # repeat-cycle search from a single element set
        self.icesat2_tle = [
            "1 43613U 18070A   26277.57046165  .00001587  00000+0  57409-4 0  9990",
            "2 43613  91.9994 320.5065 0005609  71.5392 288.6467 15.28297940449160",
        ]

    def test_defaults_keep_drag(self):
        """
        Test that, by default, the elements' drag terms are kept and they
        are propagated directly, without a repeat cycle.
        """
        orbit = GeneralPerturbationsOrbit.from_tle(self.icesat2_tle)
        self.assertFalse(orbit.remove_drag)
        self.assertIsNone(orbit.repeat_cycle)
        self.assertIsNone(orbit.get_repeat_cycle())
        self.assertNotEqual(orbit.get_bstar(), 0)

    def test_remove_drag_ignores_drag_terms(self):
        """
        Test that removing drag keeps the elements' drag terms but
        propagates them without drag, so the propagated orbit matches at the
        epoch but no longer decays (the two diverge over 30 days).
        """
        drag = GeneralPerturbationsOrbit.from_tle(self.icesat2_tle)
        no_drag = GeneralPerturbationsOrbit.from_tle(self.icesat2_tle, remove_drag=True)
        self.assertEqual(no_drag.elements, drag.elements)
        self.assertNotEqual(no_drag.get_bstar(), 0)
        epoch = drag.get_epoch()
        at_epoch = [
            direct(o).get_orbit_track(epoch).position.m for o in (drag, no_drag)
        ]
        np.testing.assert_allclose(at_epoch[0], at_epoch[1], atol=1e-3)
        earlier = epoch - timedelta(days=30)
        positions = [
            direct(o).get_orbit_track(earlier).position.m for o in (drag, no_drag)
        ]
        self.assertGreater(np.linalg.norm(positions[0] - positions[1]), 10e3)
        np.testing.assert_allclose(
            positions[1],
            drag.elements[0]
            .without_drag()
            .to_skyfield()
            .at(constants.timescale.from_datetime(earlier))
            .position.m,
        )

    def test_repeat_cycle_values(self):
        """
        Test that the repeat cycle accepts None, or "auto" or a positive
        duration (with `remove_drag`), and is preserved by serialization.
        """
        for remove_drag, repeat_cycle in [
            (False, None),
            (True, None),
            (True, "auto"),
            (True, timedelta(days=91)),
        ]:
            orbit = GeneralPerturbationsOrbit.from_tle(
                self.icesat2_tle, remove_drag=remove_drag, repeat_cycle=repeat_cycle
            )
            self.assertEqual(orbit.repeat_cycle, repeat_cycle)
            restored = GeneralPerturbationsOrbit.model_validate_json(
                orbit.model_dump_json()
            )
            self.assertEqual(restored.repeat_cycle, repeat_cycle)
        for repeat_cycle in ["never", timedelta(0), timedelta(days=-1)]:
            with self.assertRaises(ValidationError):
                GeneralPerturbationsOrbit.from_tle(
                    self.icesat2_tle, remove_drag=True, repeat_cycle=repeat_cycle
                )

    def test_repeat_cycle_requires_remove_drag(self):
        """
        Test that a repeat cycle cannot be declared (or found) for an orbit
        with drag.
        """
        with self.assertRaises(ValidationError):
            GeneralPerturbationsOrbit.from_tle(self.icesat2_tle, repeat_cycle="auto")
        with self.assertRaises(ValidationError):
            GeneralPerturbationsOrbit.from_tle(
                self.icesat2_tle, repeat_cycle=timedelta(days=91)
            )

    def test_declared_repeat_cycle_refined_to_nodal_days(self):
        """
        Test that a declared repeat cycle of 91 days is refined to 91
        nodal days (from the SGP4 secular rates), within a minute of
        ICESat-2's repeat cycle measured from the equator crossings of the
        same reference ground tracks in consecutive cycles (90 days,
        19:39:53 to 19:40:15), whereas the repeat-cycle search finds none.
        """
        orbit = GeneralPerturbationsOrbit.from_tle(
            self.icesat2_tle, remove_drag=True, repeat_cycle=timedelta(days=91)
        )
        model = orbit.get_repeat_element().to_satrec()
        nodal_day = (
            2
            * np.pi
            / (2 * np.pi / constants.EARTH_SIDEREAL_DAY_S * 60 - model.nodedot)
            * 60
        )
        repeat_cycle = orbit.get_repeat_cycle()
        self.assertAlmostEqual(repeat_cycle.total_seconds(), 91 * nodal_day, places=3)
        measured = timedelta(days=90, hours=19, minutes=40, seconds=5)
        self.assertLess(abs((repeat_cycle - measured).total_seconds()), 60)
        self.assertIsNone(
            GeneralPerturbationsOrbit.from_tle(
                self.icesat2_tle, remove_drag=True
            ).get_repeat_cycle()
        )

    def test_declared_repeat_cycle_adjusts_mean_motion(self):
        """
        Test that a declared repeat cycle keeps the orbit's elements, but
        propagates them maintained on the repeat ground track: with the mean
        motion adjusted (here, by 67 m of semimajor axis) so that 1,387
        nodal periods span the refined repeat cycle, so that the satellite
        returns to its initial Earth-fixed position after the repeat cycle
        (within 1 km), whereas with the element's own mean motion it misses
        by more than 800 km (about 113 s along track), a discontinuity at
        the end of each repeated cycle.
        """
        declared = GeneralPerturbationsOrbit.from_tle(
            self.icesat2_tle, remove_drag=True, repeat_cycle=timedelta(days=91)
        )
        no_drag = GeneralPerturbationsOrbit.from_tle(self.icesat2_tle, remove_drag=True)
        self.assertEqual(declared.elements, no_drag.elements)
        maintained = declared.get_repeat_element()
        repeat_cycle = declared.get_repeat_cycle()
        nodal_period, _ = maintained.get_nodal_period_and_day()
        self.assertAlmostEqual(
            1387 * nodal_period, repeat_cycle.total_seconds(), places=3
        )
        self.assertAlmostEqual(
            maintained.get_semimajor_axis() - no_drag.get_semimajor_axis(),
            67,
            delta=1,
        )
        self.assertEqual(maintained.bstar, 0)
        epoch = declared.get_epoch()
        misses = []
        for element in (maintained, no_drag.elements[0].without_drag()):
            positions = [
                np.array(
                    element.to_skyfield()
                    .at(constants.timescale.from_datetime(t))
                    .frame_xyz(itrs)
                    .m
                )
                for t in (epoch, epoch + repeat_cycle)
            ]
            misses.append(np.linalg.norm(positions[1] - positions[0]))
        self.assertLess(misses[0], 1e3)
        self.assertGreater(misses[1], 800e3)

    def test_get_repeat_element_idempotent(self):
        """
        Test that adjusting an element already adjusted to a repeat cycle
        leaves it unchanged.
        """
        element = GeneralPerturbationsOrbit.from_tle(
            self.icesat2_tle, remove_drag=True, repeat_cycle=timedelta(days=91)
        ).get_repeat_element()
        again = element.get_repeat_element(timedelta(days=91))
        self.assertEqual(again.mean_motion, element.mean_motion)
        self.assertEqual(again.inclination, element.inclination)

    def test_declared_repeat_cycle_sun_synchronous_solar_days(self):
        """
        Test that a declared repeat cycle of a sun-synchronous orbit
        (Landsat 9) is refined to whole mean solar days (16 days exactly),
        as the orbit is maintained at a constant local time of ascending
        node, rather than to whole nodal days from the elements' nodal
        precession, which differs from a solar day by a fraction of a second.
        """
        landsat_9_tle = [
            "1 49260U 21088A   26277.17558559  .00000158  00000+0  45237-4 0  9997",
            "2 49260  98.2176 345.9720 0001438  89.7555 270.3809 14.57106411266879",
        ]
        orbit = GeneralPerturbationsOrbit.from_tle(
            landsat_9_tle, remove_drag=True, repeat_cycle=timedelta(days=16)
        )
        self.assertEqual(orbit.get_repeat_cycle(), timedelta(days=16))
        self.assertEqual(
            orbit.elements[0].refine_repeat_cycle(timedelta(days=15, hours=20)),
            timedelta(days=16),
        )

    def test_repeat_element_sun_synchronous_inclination(self):
        """
        Test that the element of a sun-synchronous orbit maintained on its
        repeat ground track has its inclination adjusted (by a small
        fraction of a degree) so that its nodal day is exactly a mean solar
        day, so that the repeated ground track joins within about 1 km at the
        end of each repeat cycle (rather than about 4 km, from the elements'
        nodal precession).
        """
        landsat_8_tle = [
            "1 39084U 13008A   26213.27824675  .00000294  00000+0  75333-4 0  9990",
            "2 39084  98.2277 282.8718 0001275  92.4910 267.6434 14.57104473704466",
        ]
        orbit = GeneralPerturbationsOrbit.from_tle(
            landsat_8_tle, remove_drag=True, repeat_cycle=timedelta(days=16)
        )
        element = orbit.elements[0]
        maintained = orbit.get_repeat_element()
        _, nodal_day = maintained.get_nodal_period_and_day()
        self.assertAlmostEqual(nodal_day, 86400, places=5)
        self.assertNotEqual(maintained.inclination, element.inclination)
        self.assertLess(abs(maintained.inclination - element.inclination), 0.05)
        epoch = orbit.get_epoch()
        positions = [
            np.array(
                maintained.to_skyfield()
                .at(constants.timescale.from_datetime(t))
                .frame_xyz(itrs)
                .m
            )
            for t in (epoch, epoch + timedelta(days=16))
        ]
        self.assertLess(np.linalg.norm(positions[1] - positions[0]), 1.5e3)

    def test_declared_repeat_cycle_used_for_observation_events(self):
        """
        Test that the declared repeat cycle is used to repeat observation
        events: those two repeat cycles after the epoch are those of the
        maintained element in the first cycle, shifted by two cycles.
        """
        orbit = GeneralPerturbationsOrbit.from_tle(
            self.icesat2_tle, remove_drag=True, repeat_cycle=timedelta(days=91)
        )
        repeat_cycle = orbit.get_repeat_cycle()
        point = Point(id=0, latitude=70, longitude=-45)
        start = orbit.get_epoch() + 2 * repeat_cycle + timedelta(days=3)
        end = start + timedelta(days=2)
        times, codes = orbit.get_observation_events(point, start, end, 10)
        maintained = GeneralPerturbationsOrbit(elements=[orbit.get_repeat_element()])
        expected_times, expected_codes = direct(maintained).get_observation_events(
            point, start - 2 * repeat_cycle, end - 2 * repeat_cycle, 10
        )
        self.assertGreater(len(codes), 0)
        self.assertEqual(codes.tolist(), expected_codes.tolist())
        differences = np.array(
            [
                (a - b - 2 * repeat_cycle).total_seconds()
                for a, b in zip(times.utc_datetime(), expected_times.utc_datetime())
            ]
        )
        np.testing.assert_allclose(differences[codes != 1], 0, atol=1e-2)
        np.testing.assert_allclose(differences[codes == 1], 0, atol=1)

    def test_repeat_cycle_cache_follows_fields(self):
        """
        Test that the cached repeat cycle is recomputed for a copy with a
        different declared repeat cycle (which model_copy copies along with
        the cache).
        """
        orbit = GeneralPerturbationsOrbit.from_tle(
            self.icesat2_tle, remove_drag=True, repeat_cycle=timedelta(days=91)
        )
        orbit.get_repeat_cycle()
        copied = orbit.model_copy(update={"repeat_cycle": timedelta(days=30)})
        self.assertEqual(
            copied.get_repeat_cycle(),
            orbit.elements[0]
            .get_repeat_element(timedelta(days=30))
            .refine_repeat_cycle(timedelta(days=30)),
        )

    def test_pickle_and_copy_after_propagation(self):
        """
        Test that an orbit can be pickled (as to send it to another process),
        deep copied, and derived after it is propagated, which caches
        Skyfield satellites that cannot be pickled.
        """
        orbit = GeneralPerturbationsOrbit.from_tle(
            self.icesat2_tle, remove_drag=True, repeat_cycle=timedelta(days=91)
        )
        time = orbit.get_epoch() + timedelta(days=200)
        position = orbit.get_orbit_track(time).position.m
        for other in (
            pickle.loads(pickle.dumps(orbit)),
            copy.deepcopy(orbit),
            orbit.get_derived_orbit(0, 0),
        ):
            self.assertEqual(other, orbit)
            np.testing.assert_allclose(
                other.get_orbit_track(time).position.m, position, atol=1e-6
            )

    def test_options_preserved(self):
        """
        Test that the options are preserved by derived orbits, OMM
        constructors, and serialization.
        """
        orbit = GeneralPerturbationsOrbit.from_tle(
            self.icesat2_tle, remove_drag=True, repeat_cycle=timedelta(days=91)
        )
        derived = orbit.get_derived_orbit(10, 5)
        self.assertTrue(derived.remove_drag)
        self.assertEqual(derived.repeat_cycle, timedelta(days=91))
        restored = GeneralPerturbationsOrbit.model_validate_json(
            orbit.model_dump_json()
        )
        self.assertEqual(restored, orbit)
        omm = json.dumps([orbit.elements[0].to_omm_dict()])
        from_omm = GeneralPerturbationsOrbit.from_omm_json(
            omm, remove_drag=True, repeat_cycle=timedelta(days=91)
        )
        self.assertEqual(from_omm.get_repeat_cycle(), orbit.get_repeat_cycle())


class TestCompletePasses(unittest.TestCase):
    """
    Unit tests for `_complete_passes`, with a synthetic elevation angle (less
    the minimum) that is positive within passes of a half-width `w` about
    culminations at `c`: `1 - ((t - c) / w)^2`, at the nearest culmination.
    """

    def setUp(self):
        self.culminations = np.array([10.0, 20.0])
        self.width = 1.0

        def excess(x):
            x = np.asarray(x, dtype=float)
            yield constants.timescale.tt_jd(x)
            nearest = self.culminations[
                np.argmin(np.abs(x[..., None] - self.culminations), axis=-1)
            ]
            return 1 - ((x - nearest) / self.width) ** 2

        self.excess = excess

    def complete(self, jd, events, t_0=0.0, t_1=30.0):
        f_0, f_1 = _run(self.excess(np.array([t_0, t_1])))
        return _run(
            _complete_passes(
                self.excess,
                np.array(jd, dtype=float),
                np.array(events),
                t_0,
                t_1,
                f_0,
                f_1,
            )
        )

    def test_complete_passes_unchanged(self):
        """
        Test that complete passes are unchanged.
        """
        jd, events = self.complete([9, 10, 11, 19, 20, 21], [0, 1, 2, 0, 1, 2])
        np.testing.assert_array_equal(jd, [9, 10, 11, 19, 20, 21])
        np.testing.assert_array_equal(events, [0, 1, 2, 0, 1, 2])

    def test_missing_rise(self):
        """
        Test that a missing rise (a culmination and set without a rise) is
        found before the culmination.
        """
        jd, events = self.complete([10, 11, 19, 20, 21], [1, 2, 0, 1, 2])
        np.testing.assert_array_equal(events, [0, 1, 2, 0, 1, 2])
        self.assertAlmostEqual(jd[0], 9, delta=1e-6)

    def test_missing_set(self):
        """
        Test that a missing set (a rise and culmination followed by another
        rise, or by the end while below the minimum) is found after the
        culmination.
        """
        jd, events = self.complete([9, 10, 19, 20], [0, 1, 0, 1])
        np.testing.assert_array_equal(events, [0, 1, 2, 0, 1, 2])
        self.assertAlmostEqual(jd[2], 11, delta=1e-6)
        self.assertAlmostEqual(jd[5], 21, delta=1e-6)

    def test_lone_culmination(self):
        """
        Test that the rise and set of a lone culmination are found.
        """
        jd, events = self.complete([10], [1], t_1=15.0)
        np.testing.assert_array_equal(events, [0, 1, 2])
        np.testing.assert_allclose(jd, [9, 10, 11], atol=1e-6)

    def test_culmination_below_minimum(self):
        """
        Test that a culmination below the minimum, and a set without a
        culmination above it, are dropped.
        """
        self.width = 1e-9  # passes too narrow to reach the minimum at these times
        jd, events = self.complete([10.5, 11], [1, 2], t_1=15.0)
        self.assertEqual(len(jd), 0)
        self.assertEqual(len(events), 0)

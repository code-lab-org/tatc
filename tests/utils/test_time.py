"""
Unit tests for the tatc.utils.time module.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""
import unittest
from datetime import datetime, timedelta, timezone

import numpy as np

from tatc.utils import to_datetime64_ns


class TestTime(unittest.TestCase):
    """
    Unit tests for the tatc.utils.time module.
    """
    def test_to_datetime64_ns_scalar_utc(self):
        """
        Test that a UTC datetime is converted to a naive datetime64[ns].
        """
        value = datetime(2022, 6, 1, 12, tzinfo=timezone.utc)
        result = to_datetime64_ns(value)
        self.assertEqual(result, np.datetime64("2022-06-01T12:00:00", "ns"))

    def test_to_datetime64_ns_scalar_non_utc(self):
        """
        Test that a non-UTC timezone-aware datetime is normalized to UTC
        before conversion to datetime64[ns].
        """
        offset_tz = timezone(timedelta(hours=-5))
        value = datetime(2022, 6, 1, 7, tzinfo=offset_tz)
        result = to_datetime64_ns(value)
        self.assertEqual(result, np.datetime64("2022-06-01T12:00:00", "ns"))

    def test_to_datetime64_ns_list(self):
        """
        Test that a list of UTC datetimes is converted to an array of
        naive datetime64[ns] values.
        """
        values = [
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 2, tzinfo=timezone.utc),
        ]
        result = to_datetime64_ns(values)
        expected = np.array(
            ["2022-06-01T00:00:00", "2022-06-02T00:00:00"], dtype="datetime64[ns]"
        )
        np.testing.assert_array_equal(result, expected)

    def test_to_datetime64_ns_ndarray_passthrough(self):
        """
        Test that an existing datetime64 array is cast to datetime64[ns]
        without modification.
        """
        values = np.array(["2022-06-01", "2022-06-02"], dtype="datetime64[D]")
        result = to_datetime64_ns(values)
        expected = values.astype("datetime64[ns]")
        np.testing.assert_array_equal(result, expected)

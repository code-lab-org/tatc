"""
Unit tests for computations run together.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

# pylint: disable=protected-access

import unittest
from datetime import datetime, timedelta, timezone

import numpy as np
from skyfield.framelib import itrs

from tatc.constants import timescale
from tatc.schemas import SunSynchronousOrbit
from tatc.utils import computation


class TestRunTogether(unittest.TestCase):
    """
    Unit tests for computations run together.
    """

    def test_results_equal_those_run_alone(self):
        """
        Test that computations run together (sharing the Earth orientation
        quantities of their times) return the results of each run alone.
        """
        orbit = SunSynchronousOrbit(
            altitude=705e3, equator_crossing_time="13:30"
        ).to_gp_orbit()
        start = datetime(2026, 1, 1, tzinfo=timezone.utc)

        def positions(offset: float) -> computation.TimeRequest:
            t = timescale.from_datetimes(
                [start + timedelta(seconds=offset + 60 * i) for i in range(50)]
            )
            yield t
            return np.array(orbit.get_orbit_track_at_time(t).frame_xyz(itrs).m)

        together = computation._run_together([positions(k) for k in range(3)])
        for k, result in enumerate(together):
            np.testing.assert_allclose(
                result, computation._run(positions(k)), rtol=0, atol=1e-6
            )


if __name__ == "__main__":
    unittest.main()

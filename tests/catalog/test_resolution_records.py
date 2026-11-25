"""
Regression tests that recompute the records of the resolution configurations of SAR
entries from the published resolutions they are compared with (see
`validation/catalog/validate_catalog_csar.py`).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

import geopandas as gpd
import pandas as pd

from tatc.analysis import compute_sample_sizes
from tatc.catalog import get_entry
from tatc.catalog.standard import get_gates
from tatc.schemas import CircularOrbit
from tatc.utils import VelocityFrame, compute_projected_ray_position

PUBLISHED_EPOCH = datetime(2026, 10, 1, tzinfo=timezone.utc)
SIZE_TOLERANCE = 0.03


class TestCsarResolutions(unittest.TestCase):
    """
    Regression tests of the resolutions of Sentinel-1's IW sub-swaths, against the
    Sentinel-1 Product Definition's (Table 7-6, ground range detected high
    resolution): instrument, altitude (km), inclination (deg), roll (deg), and
    resolutions in azimuth and ground range (km).
    """

    resolutions = (
        ("C-SAR IW1", 693, 98.18, -29.3332, 0.0225, 0.0204),
        ("C-SAR IW2", 693, 98.18, -33.9851, 0.0226, 0.0203),
        ("C-SAR IW3", 693, 98.18, -38.0423, 0.0226, 0.0205),
    )

    def test_published_resolutions(self):
        """
        Test that the resolutions meet the gate as recorded.
        """
        entry = get_entry("csar")
        errors = []
        for name, altitude, inclination, angle, along, cross in self.resolutions:
            orbit = CircularOrbit(
                altitude=altitude * 1e3, inclination=inclination, epoch=PUBLISHED_EPOCH
            )
            satellite = entry.build(
                configuration="resolution_iw", orbit=orbit
            ).to_satellite(orbit, name=entry.name)
            time = PUBLISHED_EPOCH + timedelta(minutes=10)
            target = compute_projected_ray_position(
                orbit.to_gp_orbit().get_orbit_track(time),
                0,
                0,
                angle,
                0,
                velocity_frame=VelocityFrame.INERTIAL,
            )
            sizes = compute_sample_sizes(
                gpd.GeoDataFrame(
                    {
                        "satellite": [satellite.name],
                        "instrument": [name],
                        "epoch": pd.DatetimeIndex([time]),
                    },
                    geometry=gpd.points_from_xy(
                        [target.longitude.degrees], [target.latitude.degrees]
                    ),
                    crs="EPSG:4326",
                ),
                satellite,
            )
            errors += [
                abs(sizes.sample_along.iloc[0] / 1e3 / along - 1),
                abs(sizes.sample_cross.iloc[0] / 1e3 / cross - 1),
            ]
        within = 100 * sum(e <= SIZE_TOLERANCE for e in errors) / len(errors)
        metrics = {
            m.name: m.value
            for m in entry.record.configurations["resolution_iw"].metrics
        }
        self.assertEqual(within, metrics["sample_size_within_3pct"])
        self.assertAlmostEqual(
            100 * max(errors), metrics["sample_size_largest_error"], places=6
        )
        (gate,) = get_gates(entry.instrument_class, "sample_size")
        self.assertGreaterEqual(within, gate.limit)


if __name__ == "__main__":
    unittest.main()

"""
Unit tests for airborne ground track analysis functions.

The expected values come from the geometry of a scan over a flat surface,
which the footprints match to within the Earth's curvature.

@author: Naeun Kim <hkim439@asu.edu>
@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import timedelta

import numpy as np
import pandas as pd
import shapely
from pyproj import Geod
from shapely.geometry import box

from tatc.analysis import (
    collect_airborne_ground_track,
    collect_ground_track,
    compute_airborne_ground_track,
    find_coincident_ground_track,
)
from tatc.schemas import (
    Aircraft,
    Instrument,
    PointedInstrument,
    Satellite,
    TwoLineElements,
    Trajectory,
)

from ..schemas.air.common import START, straight_trajectory, surface_distance

GEOD = Geod(ellps="WGS84")
COLUMNS = ["time", "aircraft", "instrument", "valid_obs", "geometry"]
HALF_ANGLE = 25.0


def _scanner(name="scanner", **kwargs):
    """
    Gets a cross-track scanner with a half angle of `HALF_ANGLE` degrees.
    """
    return PointedInstrument(
        name=name,
        cross_track_field_of_view=2 * HALF_ANGLE,
        along_track_field_of_view=0.1,
        view_geometry="scan",
        is_swept=True,
        **kwargs,
    )


def _aircraft(data, instruments=None, name="aircraft"):
    """
    Gets an aircraft with a scanner flying a trajectory.
    """
    return Aircraft(
        name=name,
        trajectory=Trajectory.from_dataframe(data),
        instruments=[_scanner()] if instruments is None else instruments,
    )


def _offsets(results, data, row):
    """
    Gets the forward and rightward offsets (meters) of the exterior points of
    a footprint from the aircraft's position, relative to its track.
    """
    lon, lat = shapely.get_coordinates(results.geometry.iloc[row].exterior).T
    azimuth, _, distance = GEOD.inv(
        np.full(len(lon), data["longitude"][row]),
        np.full(len(lat), data["latitude"][row]),
        lon,
        lat,
    )
    angle = np.radians(np.asarray(azimuth) - data["track"][row])
    return distance * np.cos(angle), distance * np.sin(angle)


class TestAirborneTrack(unittest.TestCase):
    """
    Unit tests for airborne ground track analysis functions.
    """

    def test_columns_and_order(self):
        """
        Tests the columns of the results, sorted by time, for one or more
        aircraft and instruments.
        """
        data = straight_trajectory(n=10)
        aircraft = _aircraft(
            data, [_scanner("a"), Instrument(name="b", field_of_regard=30)]
        )
        results = collect_airborne_ground_track(aircraft)
        self.assertEqual(list(results.columns), COLUMNS)
        self.assertEqual(len(results), 10)
        self.assertTrue((results["instrument"] == "a").all())
        self.assertTrue(results["valid_obs"].all())
        self.assertTrue(results["time"].is_monotonic_increasing)
        pd.testing.assert_series_equal(results["time"], data["time"], check_names=False)
        combined = collect_airborne_ground_track(
            [aircraft, _aircraft(data, name="other")], instrument_index=None
        )
        self.assertEqual(len(combined), 30)
        self.assertEqual(
            combined.groupby(["aircraft", "instrument"]).size().to_dict(),
            {("aircraft", "a"): 10, ("aircraft", "b"): 10, ("other", "scanner"): 10},
        )
        self.assertTrue(combined["time"].is_monotonic_increasing)

    def test_empty(self):
        """
        Tests that times outside the trajectory are omitted.
        """
        aircraft = _aircraft(straight_trajectory(n=10))
        results = collect_airborne_ground_track(
            aircraft, times=[START - timedelta(minutes=1)]
        )
        self.assertTrue(results.empty)
        self.assertEqual(list(results.columns), COLUMNS)
        combined = collect_airborne_ground_track(
            [aircraft], times=[START + timedelta(hours=1)]
        )
        self.assertTrue(combined.empty)
        self.assertEqual(list(combined.columns), COLUMNS)

    def test_not_an_aircraft(self):
        """
        Tests that a data frame is not accepted in place of an aircraft.
        """
        with self.assertRaisesRegex(TypeError, "must be an Aircraft"):
            collect_airborne_ground_track(straight_trajectory(n=3))

    def test_level_swath(self):
        """
        Tests that a level scan is as wide as a scan over a flat surface and
        as long as the distance flown in a time step, centered on the aircraft.
        """
        data = straight_trajectory(n=20, altitude=8000.0, speed=150.0)
        results = collect_airborne_ground_track(_aircraft(data))
        forward, right = _offsets(results, data, 10)
        width = 2 * 8000.0 * np.tan(np.radians(HALF_ANGLE))
        self.assertAlmostEqual(right.max() - right.min(), width, delta=0.002 * width)
        self.assertAlmostEqual(right.max() + right.min(), 0, delta=1.0)
        # the time step of one second, centered on the time
        self.assertAlmostEqual(forward.max(), 75.0, delta=0.5)
        self.assertAlmostEqual(forward.min(), -75.0, delta=0.5)

    def test_time_step(self):
        """
        Tests that a time step sets the length of the swept footprints and
        that a gap between times does not stretch the footprints beside it.
        """
        data = straight_trajectory(n=20, speed=100.0)
        results = collect_airborne_ground_track(
            _aircraft(data), time_step=timedelta(seconds=3)
        )
        forward, _ = _offsets(results, data, 10)
        self.assertAlmostEqual(forward.max() - forward.min(), 300.0, delta=1.0)
        gapped = data.drop(index=range(5, 15)).reset_index(drop=True)
        results = collect_airborne_ground_track(_aircraft(gapped))
        for row in (4, 5):
            forward, _ = _offsets(results, gapped, row)
            self.assertAlmostEqual(forward.max() - forward.min(), 100.0, delta=1.0)

    def test_times(self):
        """
        Tests that times between the states of a trajectory are sampled.
        """
        data = straight_trajectory(n=11, speed=100.0, step=10.0)
        times = [START + timedelta(seconds=s) for s in np.arange(0, 100.5, 0.5)]
        results = collect_airborne_ground_track(_aircraft(data), times=times)
        self.assertEqual(len(results), len(times))
        swath = shapely.union_all(
            results.to_crs(results.estimate_utm_crs()).geometry.values
        )
        width = 2 * 10000.0 * np.tan(np.radians(HALF_ANGLE))
        # the swath from the first to the last state (trimmed to the
        # trajectory, not 0.25 seconds beyond it)
        self.assertAlmostEqual(
            swath.area, width * 100.0 * 100.0, delta=0.005 * swath.area
        )
        self.assertIsInstance(swath, shapely.Polygon)

    def test_sweeps_trimmed(self):
        """
        Tests that swept footprints are trimmed to the trajectory: they do not
        extend before the first state, after the last, or into a gap longer
        than the maximum gap.
        """
        data = straight_trajectory(n=20, speed=100.0)
        later = data.iloc[10:].copy()
        later["time"] += timedelta(minutes=5)
        gapped = pd.concat([data.iloc[:10], later], ignore_index=True)
        aircraft = _aircraft(gapped)
        results = collect_airborne_ground_track(
            aircraft, time_step=timedelta(seconds=10)
        )
        # the forward extents (meters) of the footprints of a 10-second step
        # at 100 meters per second, in parts from 0 to 9 and 10 to 19 seconds,
        # trimmed within 0.1 seconds (10 meters) inside the trajectory
        for row, extent in {
            0: (0, 500),
            5: (-500, 400),
            9: (-500, 0),
            10: (0, 500),
            19: (-500, 0),
        }.items():
            with self.subTest(row=row):
                forward, _ = _offsets(results, gapped, row)
                self.assertTrue(extent[0] - 1 <= forward.min() <= extent[0] + 10.5)
                self.assertTrue(extent[1] - 10.5 <= forward.max() <= extent[1] + 1)
        self.assertTrue(results["valid_obs"].all())
        # by default, half the time to the neighboring states
        forward, _ = _offsets(collect_airborne_ground_track(aircraft), gapped, 0)
        self.assertAlmostEqual(forward.min(), 0, delta=0.5)
        self.assertAlmostEqual(forward.max(), 50, delta=0.5)

    def test_times_order_and_repeats(self):
        """
        Tests that the results of unsorted times are sorted with their
        footprints, and that a single or repeated time has a footprint.
        """
        aircraft = _aircraft(straight_trajectory(n=10, speed=100.0))
        times = [START + timedelta(seconds=s) for s in (5, 3, 7, 5)]
        results = collect_airborne_ground_track(aircraft, times=times)
        ordered = collect_airborne_ground_track(aircraft, times=sorted(times))
        self.assertEqual(list(results["time"]), sorted(times))
        self.assertTrue(results.geometry.geom_equals_exact(ordered.geometry, 0).all())
        self.assertTrue(results["valid_obs"].all())
        for instrument in (_scanner(), Instrument(field_of_regard=10.0)):
            with self.subTest(instrument=instrument.name):
                single = collect_airborne_ground_track(
                    _aircraft(straight_trajectory(n=10), [instrument]), times=times[:1]
                )
                self.assertTrue(single["valid_obs"].all())
                self.assertFalse(single.geometry.is_empty.any())
        empty = collect_airborne_ground_track(aircraft, times=[])
        self.assertTrue(empty.empty)
        self.assertEqual(list(empty.columns), COLUMNS)

    def _assert_swath_area(self, data, track, latitude):
        """
        Asserts that the footprints of a straight, level flight at 10 km and
        200 meters per second tile a swath of the expected area (from their
        geodesic areas, which are not distorted by longitude and latitude).
        """
        results = collect_airborne_ground_track(_aircraft(data))
        self.assertTrue(results["valid_obs"].all())
        self.assertTrue(results.geometry.is_valid.all())
        area = sum(abs(GEOD.geometry_area_perimeter(g)[0]) for g in results.geometry)
        width = 2 * 10000.0 * np.tan(np.radians(HALF_ANGLE))
        length = surface_distance(200.0 * (len(data) - 1), latitude, track, 10000.0)
        self.assertAlmostEqual(area, width * length, delta=0.005 * area)
        return results

    def test_antimeridian(self):
        """
        Tests the footprints of a flight across the anti-meridian: split
        there, each part narrow in longitude, with the area of the swath.
        """
        data = straight_trajectory(n=20, track=90.0, longitude=179.99, latitude=10.0)
        self.assertTrue((data["longitude"] < 0).any() and (data["longitude"] > 0).any())
        results = self._assert_swath_area(data, 90.0, 10.0)
        for geometry in results.geometry:
            for part in getattr(geometry, "geoms", [geometry]):
                self.assertLess(part.bounds[2] - part.bounds[0], 0.01)
        coverage = compute_airborne_ground_track(_aircraft(data)).geometry.iloc[0]
        self.assertEqual(coverage.geom_type, "MultiPolygon")

    def test_over_pole(self):
        """
        Tests the footprints of a flight over the north pole.
        """
        data = straight_trajectory(n=20, track=0.0, latitude=89.98)
        # the track turns from north to south past the pole
        self.assertEqual(set(np.round(data["track"])), {0.0, 180.0})
        self._assert_swath_area(data, 0.0, 89.98)

    def test_times_in_gaps(self):
        """
        Tests that times within a gap between states longer than the maximum
        gap, or outside the trajectory, are omitted.
        """
        data = straight_trajectory(n=3, step=10.0)
        later = straight_trajectory(n=3, step=10.0, latitude=35.5)
        later["time"] += timedelta(hours=1)
        both = pd.concat([data, later], ignore_index=True)
        times = [START + timedelta(minutes=m) for m in (-1, 0, 0.25, 30, 60, 60.25, 61)]
        results = collect_airborne_ground_track(_aircraft(both), times=times)
        self.assertEqual(list(results["time"]), [times[i] for i in (1, 2, 4, 5)])
        self.assertTrue(results["valid_obs"].all())
        aircraft = Aircraft(
            name="aircraft",
            trajectory=Trajectory.from_dataframe(both, max_gap=None),
            instruments=[_scanner()],
        )
        self.assertEqual(len(collect_airborne_ground_track(aircraft, times=times)), 5)

    def test_altaz(self):
        """
        Tests the altitude/azimuth angles of the aircraft and the sun seen from
        the footprint centers.
        """
        data = straight_trajectory(n=5, roll=10.0, track=0.0)
        results = collect_airborne_ground_track(
            _aircraft(data, [Instrument(field_of_regard=10.0)]),
            aircraft_altaz=True,
            solar_altaz=True,
        )
        for column in (
            "aircraft_alt",
            "aircraft_az",
            "solar_alt",
            "solar_az",
            "glint_angle",
        ):
            self.assertIn(column, results.columns)
        # a generic instrument observes about nadir: the aircraft is overhead
        np.testing.assert_allclose(results["aircraft_alt"], 90.0, atol=1e-3)
        self.assertTrue(
            ((results["solar_alt"] > 0) & (results["solar_alt"] < 90)).all()
        )
        rolled = collect_airborne_ground_track(_aircraft(data), aircraft_altaz=True)
        # rolled right wing down, the scan centers to the left (west) of a
        # northward track, from which the aircraft is seen to the east
        np.testing.assert_allclose(rolled["aircraft_alt"], 80.0, atol=0.05)
        np.testing.assert_allclose(rolled["aircraft_az"], 90.0, atol=0.5)

    def test_compute(self):
        """
        Tests that the dissolved coverage is the union of the valid footprints.
        """
        data = straight_trajectory(n=20, altitude=5000.0)
        aircraft = _aircraft(data)
        track = collect_airborne_ground_track(aircraft)
        coverage = compute_airborne_ground_track(aircraft)
        self.assertEqual(len(coverage), 1)
        self.assertAlmostEqual(
            coverage.geometry.iloc[0].area,
            shapely.union_all(track.geometry.values).area,
            delta=1e-12,
        )
        self.assertTrue(compute_airborne_ground_track(aircraft, elevation=6000.0).empty)

    def test_underflight(self):
        """
        Tests that the footprints of an aircraft beneath a satellite's ground
        track are coincident with the satellite's.
        """
        satellite = Satellite(
            name="ISS",
            orbit=TwoLineElements(
                tle=[
                    "1 25544U 98067A   26074.50000000  .00016717  00000-0  10270-3 0  9994",
                    "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391 50000",
                ]
            ),
            instruments=[Instrument(name="imager", field_of_regard=20.0)],
        )
        times = [START + timedelta(seconds=s) for s in range(0, 6000, 10)]
        ground_track = collect_ground_track(satellite, times)
        # an aircraft flying beneath the satellite's nadir point at one time
        center = ground_track.geometry.iloc[0].centroid
        data = straight_trajectory(n=30, latitude=center.y, longitude=center.x)
        data["time"] = times[0] + pd.to_timedelta(np.arange(30) - 15, unit="s")
        airborne = collect_airborne_ground_track(_aircraft(data))
        pairs = find_coincident_ground_track(
            pd.concat([ground_track, airborne], ignore_index=True),
            max_time_difference=timedelta(minutes=1),
        )
        self.assertFalse(pairs.empty)
        self.assertTrue(
            (pairs["satellite_1"].fillna(pairs["satellite_2"]) == "ISS").all()
        )
        self.assertTrue(
            (pairs["aircraft_1"].fillna(pairs["aircraft_2"]) == "aircraft").all()
        )

    def test_roll_shifts_swath(self):
        """
        Tests that a positive roll angle (right wing down) shifts the swath to
        the left.
        """
        roll = 10.0
        data = straight_trajectory(n=20, roll=roll)
        _, right = _offsets(collect_airborne_ground_track(_aircraft(data)), data, 10)
        alpha, phi = np.radians([HALF_ANGLE, roll])
        self.assertAlmostEqual(
            right.min(), -10000.0 * np.tan(alpha + phi), delta=0.003 * 10000
        )
        self.assertAlmostEqual(
            right.max(), 10000.0 * np.tan(alpha - phi), delta=0.003 * 10000
        )

    def test_pitch_shifts_swath(self):
        """
        Tests that a positive pitch angle (nose up) shifts the swath forward.
        """
        data = straight_trajectory(n=20, speed=100.0, pitch=5.0)
        forward, _ = _offsets(collect_airborne_ground_track(_aircraft(data)), data, 10)
        ahead = 10000.0 * np.tan(np.radians(5.0))
        self.assertAlmostEqual(forward.max(), ahead + 50.0, delta=3.0)
        self.assertAlmostEqual(forward.min(), ahead - 50.0, delta=3.0)

    def test_heading_orients_scan(self):
        """
        Tests that the scan is across the heading, not the track, in a
        crosswind.
        """
        data = straight_trajectory(n=20, speed=1.0, track=30.0, heading=50.0)
        results = collect_airborne_ground_track(_aircraft(data))
        projected = results.to_crs(results.estimate_utm_crs()).geometry.iloc[10]
        rectangle = np.array(projected.minimum_rotated_rectangle.exterior.coords)
        sides = np.diff(rectangle, axis=0)
        longest = sides[np.argmax(np.hypot(*sides.T))]
        azimuth = np.mod(np.degrees(np.arctan2(*longest)), 180)
        # across the heading (50 degrees), less the grid convergence (about
        # 1 degree)
        self.assertAlmostEqual(azimuth, 140.0, delta=1.5)

    def test_elevation(self):
        """
        Tests that the swath at an elevation is as wide as a scan over a flat
        surface from the height above it, and that a level at or above the
        aircraft is not observed by a downward scan.
        """
        data = straight_trajectory(n=10, altitude=8000.0)
        aircraft = _aircraft(data)
        results = collect_airborne_ground_track(aircraft, elevation=3000.0)
        _, right = _offsets(results, data, 5)
        width = 2 * 5000.0 * np.tan(np.radians(HALF_ANGLE))
        self.assertAlmostEqual(right.max() - right.min(), width, delta=0.002 * width)
        # the footprints are at the elevation
        np.testing.assert_allclose(
            shapely.get_coordinates(results.geometry.iloc[5], include_z=True)[:, 2],
            3000.0,
        )
        above = collect_airborne_ground_track(
            aircraft, elevation=9000.0, aircraft_altaz=True, solar_altaz=True
        )
        self.assertFalse(above["valid_obs"].any())
        self.assertTrue(above.geometry.is_empty.all())
        # without a footprint, the angles at its center are not defined
        for column in (
            "aircraft_alt",
            "aircraft_az",
            "solar_alt",
            "solar_az",
            "glint_angle",
        ):
            self.assertTrue(above[column].isna().all(), column)

    def test_turn_has_no_gaps(self):
        """
        Tests that the swath of a banked turn has no gaps between consecutive
        footprints.
        """
        n, rate, speed = 60, 3.0, 130.0
        lon, lat, track = [-110.0], [35.0], [0.0]
        for _ in range(n - 1):
            x, y, _ = GEOD.fwd(lon[-1], lat[-1], track[-1] + rate / 2, speed)
            lon.append(x)
            lat.append(y)
            track.append(track[-1] + rate)
        data = pd.DataFrame(
            {
                "time": START + pd.to_timedelta(np.arange(n), unit="s"),
                "latitude": lat,
                "longitude": lon,
                "altitude": 6000.0,
                "track": np.mod(track, 360),
                "ground_speed": speed,
                "roll": np.degrees(np.arctan(speed * np.radians(rate) / 9.80665)),
            }
        )
        results = collect_airborne_ground_track(_aircraft(data))
        self.assertTrue(results["valid_obs"].all())
        swath = shapely.union_all(results.geometry.values)
        self.assertIsInstance(swath, shapely.Polygon)
        self.assertEqual(len(swath.interiors), 0)

    def test_mask(self):
        """
        Tests that the footprints are clipped to a mask.
        """
        data = straight_trajectory(n=30, track=0.0)
        mask = box(-111.0, data["latitude"][10], -109.0, data["latitude"][20])
        results = collect_airborne_ground_track(_aircraft(data), mask=mask)
        self.assertTrue(10 <= len(results) <= 12)
        self.assertTrue(results.geometry.within(mask.buffer(1e-9)).all())
        self.assertTrue(results["time"].is_monotonic_increasing)

    def test_validity(self):
        """
        Tests that observations are valid if they satisfy the instrument's
        requirements (here, a sunlit target).
        """
        data = straight_trajectory(n=3)
        for hour, sunlit in ((19, True), (7, False)):
            with self.subTest(hour=hour):
                data["time"] = pd.Timestamp(
                    f"2026-03-01T{hour:02d}:00Z"
                ) + pd.to_timedelta(np.arange(3), unit="s")
                results = collect_airborne_ground_track(
                    _aircraft(data, [_scanner(req_target_sunlit=True)]),
                    solar_altaz=True,
                )
                self.assertEqual(results["valid_obs"].tolist(), [sunlit] * 3)
                # an invalid observation with a footprint keeps its angles
                self.assertEqual((results["solar_alt"] > 0).tolist(), [sunlit] * 3)


if __name__ == "__main__":
    unittest.main()

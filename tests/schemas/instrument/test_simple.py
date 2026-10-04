"""
Unit tests for the Instrument schema.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

import numpy as np
from pydantic import ValidationError
from skyfield.api import EarthSatellite, wgs84
from skyfield.framelib import itrs

from tatc.constants import de421, timescale
from tatc.schemas import CircularOrbit, Instrument
from tatc.utils import geodesic_distance


class TestInstrument(unittest.TestCase):
    """
    Unit tests for the Instrument schema.
    """

    def setUp(self):
        noon_utc = datetime(2020, 3, 20, 12, tzinfo=timezone.utc)
        self.test_time = timescale.from_datetime(noon_utc)
        self.test_sat_1 = EarthSatellite.from_satrec(
            CircularOrbit(
                mean_altitude=400000,
                true_anomaly=0,
                epoch=noon_utc,
                inclination=0.0,
                right_ascension_ascending_node=0.0,
            )
            .to_gp_orbit()
            .elements[0]
            .to_satrec(),
            timescale,
        )
        self.test_sat_2 = EarthSatellite.from_satrec(
            CircularOrbit(
                mean_altitude=400000,
                true_anomaly=0,
                epoch=noon_utc,
                inclination=0.0,
                right_ascension_ascending_node=80.0,
            )
            .to_gp_orbit()
            .elements[0]
            .to_satrec(),
            timescale,
        )
        self.test_sat_3 = EarthSatellite.from_satrec(
            CircularOrbit(
                mean_altitude=400000,
                true_anomaly=0,
                epoch=noon_utc,
                inclination=0.0,
                right_ascension_ascending_node=100.0,
            )
            .to_gp_orbit()
            .elements[0]
            .to_satrec(),
            timescale,
        )
        self.test_sat_4 = EarthSatellite.from_satrec(
            CircularOrbit(
                mean_altitude=400000,
                true_anomaly=0,
                epoch=noon_utc,
                inclination=0.0,
                right_ascension_ascending_node=180.0,
            )
            .to_gp_orbit()
            .elements[0]
            .to_satrec(),
            timescale,
        )
        self.test_sat_5 = EarthSatellite.from_satrec(
            CircularOrbit(
                mean_altitude=400000,
                true_anomaly=0,
                epoch=noon_utc,
                inclination=45.0,
                right_ascension_ascending_node=0.0,
            )
            .to_gp_orbit()
            .elements[0]
            .to_satrec(),
            timescale,
        )

    def test_good_data(self):
        """
        Test that an Instrument can be created with valid data.
        """
        good_data = {
            "name": "Test Instrument",
            "field_of_regard": 20.0,
            "min_access_time": timedelta(seconds=10),
            "req_self_sunlit": None,
            "req_target_sunlit": None,
        }
        o = Instrument(**good_data)
        self.assertEqual(o.name, good_data.get("name"))
        self.assertEqual(o.field_of_regard, good_data.get("field_of_regard"))
        self.assertEqual(o.min_access_time, good_data.get("min_access_time"))
        self.assertEqual(o.req_self_sunlit, good_data.get("req_self_sunlit"))
        self.assertEqual(o.req_target_sunlit, good_data.get("req_target_sunlit"))

    def test_field_of_regard_bounds(self):
        """
        Test that field_of_regard must be in the interval (0, 360].
        """
        Instrument(name="Test Instrument", field_of_regard=360.0)
        with self.assertRaises(ValidationError):
            Instrument(name="Test Instrument", field_of_regard=0.0)
        with self.assertRaises(ValidationError):
            Instrument(name="Test Instrument", field_of_regard=360.1)
        with self.assertRaises(ValidationError):
            Instrument(name="Test Instrument", field_of_regard=-10.0)

    def test_access_time_fixed_default(self):
        """
        Test that access_time_fixed defaults to False.
        """
        o = Instrument(name="Test Instrument")
        self.assertFalse(o.access_time_fixed)
        o = Instrument(name="Test Instrument", access_time_fixed=True)
        self.assertTrue(o.access_time_fixed)

    def test_get_swath_width(self):
        """
        Test that the swath width can be computed from the field of regard.
        """
        o = Instrument(name="GMI", field_of_regard=15.0)
        self.assertAlmostEqual(o.get_swath_width(705000), 185815, delta=1.0)

    def test_get_min_elevation_angle(self):
        """
        Test that the minimum elevation angle can be computed from the field of regard.
        """
        o = Instrument(name="GMI", field_of_regard=15.0)
        self.assertAlmostEqual(o.get_min_elevation_angle(705000), 81.66446, delta=0.01)

    def test_valid_observation_no_constraints(self):
        """
        Test that an observation is valid when there are no constraints
        on sunlit conditions.
        """
        o = Instrument(name="Test Instrument")
        self.assertTrue(o.is_valid_observation(self.test_sat_1.at(self.test_time)).all())  # type: ignore
        self.assertTrue(o.is_valid_observation(self.test_sat_2.at(self.test_time)).all())  # type: ignore
        self.assertTrue(o.is_valid_observation(self.test_sat_3.at(self.test_time)).all())  # type: ignore
        self.assertTrue(o.is_valid_observation(self.test_sat_4.at(self.test_time)).all())  # type: ignore

    def test_valid_observation_self_sunlit(self):
        """
        Test that an observation is valid when the instrument requires
        self-sunlit conditions.
        """
        o = Instrument(name="Test Instrument", req_self_sunlit=True)
        self.assertTrue(o.is_valid_observation(self.test_sat_1.at(self.test_time)).all())  # type: ignore
        self.assertTrue(o.is_valid_observation(self.test_sat_2.at(self.test_time)).all())  # type: ignore
        self.assertTrue(o.is_valid_observation(self.test_sat_3.at(self.test_time)).all())  # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_4.at(self.test_time)).all())  # type: ignore

    def test_valid_observation_self_not_sunlit(self):
        """
        Test that an observation is valid when the instrument requires
        self-not-sunlit conditions.
        """
        o = Instrument(name="Test Instrument", req_self_sunlit=False)
        self.assertFalse(o.is_valid_observation(self.test_sat_1.at(self.test_time)).all())  # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_2.at(self.test_time)).all())  # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_3.at(self.test_time)).all())  # type: ignore
        self.assertTrue(o.is_valid_observation(self.test_sat_4.at(self.test_time)).all())  # type: ignore

    def test_valid_observation_target_sunlit(self):
        """
        Test that an observation is valid when the instrument requires
        target-sunlit conditions.
        """
        o = Instrument(name="Test Instrument", req_target_sunlit=True)
        self.assertTrue(o.is_valid_observation(self.test_sat_1.at(self.test_time)).all())  # type: ignore
        self.assertTrue(o.is_valid_observation(self.test_sat_2.at(self.test_time)).all())  # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_3.at(self.test_time)).all())  # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_4.at(self.test_time)).all())  # type: ignore

    def test_valid_observation_target_not_sunlit(self):
        """
        Test that an observation is valid when the instrument requires
        target-not-sunlit conditions.
        """
        o = Instrument(name="Test Instrument", req_target_sunlit=False)
        self.assertFalse(o.is_valid_observation(self.test_sat_1.at(self.test_time)).all())  # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_2.at(self.test_time)).all())  # type: ignore
        self.assertTrue(o.is_valid_observation(self.test_sat_3.at(self.test_time)).all())  # type: ignore
        self.assertTrue(o.is_valid_observation(self.test_sat_4.at(self.test_time)).all())  # type: ignore

    def test_valid_observation_self_sunlit_target_sunlit(self):
        """
        Test that an observation is valid when the instrument requires
        both self-sunlit and target-sunlit conditions."""
        o = Instrument(
            name="Test Instrument", req_self_sunlit=True, req_target_sunlit=True
        )
        self.assertTrue(o.is_valid_observation(self.test_sat_1.at(self.test_time)).all())  # type: ignore
        self.assertTrue(o.is_valid_observation(self.test_sat_2.at(self.test_time)).all())  # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_3.at(self.test_time)).all())  # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_4.at(self.test_time)).all())  # type: ignore

    def test_valid_observation_self_not_sunlit_target_sunlit(self):
        """
        Test that an observation is valid when the instrument requires
        self-not-sunlit and target-sunlit conditions.
        """
        o = Instrument(
            name="Test Instrument", req_self_sunlit=False, req_target_sunlit=True
        )
        self.assertFalse(o.is_valid_observation(self.test_sat_1.at(self.test_time)).all())  # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_2.at(self.test_time)).all())  # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_3.at(self.test_time)).all())  # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_4.at(self.test_time)).all())  # type: ignore

    def test_valid_observation_self_sunlit_target_not_sunlit(self):
        """
        Test that an observation is valid when the instrument requires
        self-sunlit and target-not-sunlit conditions.
        """
        o = Instrument(
            name="Test Instrument", req_self_sunlit=True, req_target_sunlit=False
        )
        self.assertFalse(o.is_valid_observation(self.test_sat_1.at(self.test_time)).all())  # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_2.at(self.test_time)).all())  # type: ignore
        self.assertTrue(o.is_valid_observation(self.test_sat_3.at(self.test_time)).all())  # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_4.at(self.test_time)).all())  # type: ignore

    def test_valid_observation_self_not_sunlit_target_not_sunlit(self):
        """
        Test that an observation is valid when the instrument requires
        both self-not-sunlit and target-not-sunlit conditions.
        """
        o = Instrument(
            name="Test Instrument", req_self_sunlit=False, req_target_sunlit=False
        )
        self.assertFalse(o.is_valid_observation(self.test_sat_1.at(self.test_time)).all())  # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_2.at(self.test_time)).all())  # type: ignore
        self.assertFalse(o.is_valid_observation(self.test_sat_3.at(self.test_time)).all())  # type: ignore
        self.assertTrue(o.is_valid_observation(self.test_sat_4.at(self.test_time)).all())  # type: ignore

    def test_valid_observation_self_sunlit_vector(self):
        """
        Test that an observation is valid when the instrument requires
        self-sunlit conditions for a vector of times.
        """
        o = Instrument(name="Test Instrument", req_self_sunlit=True)
        times = timescale.utc(2020, 3, 20, [11, 12, 13])  # type: ignore
        results = o.is_valid_observation(self.test_sat_1.at(times))  # type: ignore
        self.assertEqual(len(results), 3)
        self.assertFalse(results[0])
        self.assertTrue(results[1])
        self.assertFalse(results[2])

    def test_valid_observation_target_sunlit_vector(self):
        """
        Test that an observation is valid when the instrument requires
        target-sunlit conditions for a vector of times.
        """
        o = Instrument(name="Test Instrument", req_target_sunlit=True)
        times = timescale.utc(2020, 3, 20, [11, 12, 13])  # type: ignore
        results = o.is_valid_observation(self.test_sat_1.at(times))  # type: ignore
        self.assertEqual(len(results), 3)
        self.assertFalse(results[0])
        self.assertTrue(results[1])
        self.assertFalse(results[2])

    def test_valid_observation_self_sunlit_vector_inclined(self):
        """
        Test that an observation is valid when the instrument requires self-sunlit conditions for a vector of times with an inclined orbit.
        """
        o = Instrument(name="Test Instrument", req_self_sunlit=True)
        times = timescale.utc(2020, 3, 20, [11, 12, 13])  # type: ignore
        results = o.is_valid_observation(self.test_sat_5.at(times))  # type: ignore
        self.assertEqual(len(results), 3)
        self.assertFalse(results[0])
        self.assertTrue(results[1])
        self.assertFalse(results[2])

    def test_valid_observation_target_sunlit_vector_inclined(self):
        """
        Test that an observation is valid when the instrument requires target-sunlit conditions for a vector of times with an inclined orbit.
        """
        o = Instrument(name="Test Instrument", req_target_sunlit=True)
        times = timescale.utc(2020, 3, 20, [11, 12, 13])  # type: ignore
        results = o.is_valid_observation(self.test_sat_5.at(times))  # type: ignore
        self.assertEqual(len(results), 3)
        self.assertFalse(results[0])
        self.assertTrue(results[1])
        self.assertFalse(results[2])

    def test_target_solar_elevation_bounds(self):
        """
        Test that target solar elevation angles are bounded and that the
        minimum cannot exceed the maximum.
        """
        for kwargs in [
            {"min_target_solar_elevation": -91},
            {"max_target_solar_elevation": 91},
            {"min_target_solar_elevation": 10, "max_target_solar_elevation": 5},
        ]:
            with self.assertRaises(ValidationError):
                Instrument(name="Test Instrument", **kwargs)
        o = Instrument(
            name="Test Instrument",
            min_target_solar_elevation=5,
            max_target_solar_elevation=5,
        )
        self.assertEqual(o.min_target_solar_elevation, 5)

    def test_valid_observation_target_solar_elevation(self):
        """
        Test that an observation is valid when the target's solar elevation
        angle lies within the required range, for a vector of times over
        half a day with an inclined orbit, alone and combined with the
        target-sunlit requirement.
        """
        times = timescale.utc(2020, 3, 20, 6, range(0, 720, 5))  # type: ignore
        orbit_track = self.test_sat_5.at(times)  # type: ignore
        subpoint = wgs84.subpoint_of(orbit_track)
        solar_alt = (
            (de421["earth"] + subpoint)
            .at(times)
            .observe(de421["sun"])
            .apparent()
            .altaz()[0]
            .degrees
        )
        # the sample spans day and night, including low solar elevations
        self.assertTrue(np.any((solar_alt > 0) & (solar_alt < 10)))
        self.assertTrue(np.any(solar_alt < -12))
        for kwargs, expected in [
            ({"min_target_solar_elevation": 10}, solar_alt >= 10),
            ({"max_target_solar_elevation": -12}, solar_alt <= -12),
            (
                {"min_target_solar_elevation": 0, "max_target_solar_elevation": 30},
                (solar_alt >= 0) & (solar_alt <= 30),
            ),
            (
                {"req_target_sunlit": True, "max_target_solar_elevation": 30},
                (solar_alt > 0) & (solar_alt <= 30),
            ),
            (
                {"req_target_sunlit": False, "min_target_solar_elevation": 10},
                np.zeros_like(solar_alt, dtype=bool),
            ),
        ]:
            o = Instrument(name="Test Instrument", **kwargs)
            np.testing.assert_array_equal(
                o.is_valid_observation(orbit_track, subpoint), expected, str(kwargs)
            )

    def test_compute_footprint_center_geocentric_nadir(self):
        """
        Test that the footprint center of an instrument with a geocentric
        nadir reference lies on the line from the satellite to the Earth's
        center, away from the geodetic sub-satellite point for an inclined
        orbit off the equator.
        """
        o = Instrument(name="Test Instrument", nadir_reference="geocentric")
        orbit_track = self.test_sat_5.at(timescale.utc(2020, 3, 20, 12, 20))  # type: ignore
        center = o.compute_footprint_center(orbit_track)
        position = np.array(orbit_track.frame_xyz(itrs).m)
        ground = np.array(center.itrs_xyz.m)
        self.assertLess(
            np.linalg.norm(np.cross(position, ground))
            / np.linalg.norm(position)
            / np.linalg.norm(ground),
            1e-9,
        )
        subpoint = wgs84.subpoint_of(orbit_track)
        self.assertGreater(abs(subpoint.latitude.degrees), 20)
        self.assertGreater(
            geodesic_distance(
                center.longitude.degrees,
                center.latitude.degrees,
                subpoint.longitude.degrees,
                subpoint.latitude.degrees,
            ),
            500,
        )

    def test_compute_footprint_center_matches_subpoint(self):
        """
        Test that the footprint center of a nadir-pointing instrument
        (zero field of view, roll, and pitch) coincides with Skyfield's
        own WGS 84 sub-satellite point, for both equatorial and inclined
        orbits.
        """
        o = Instrument(name="Test Instrument")
        for sat in (self.test_sat_1, self.test_sat_2, self.test_sat_5):
            orbit_track = sat.at(self.test_time)  # type: ignore
            center = o.compute_footprint_center(orbit_track)
            subpoint = wgs84.subpoint_of(orbit_track)
            self.assertAlmostEqual(
                geodesic_distance(
                    center.longitude.degrees,
                    center.latitude.degrees,
                    subpoint.longitude.degrees,
                    subpoint.latitude.degrees,
                ),
                0,
                delta=1e-3,
            )
            self.assertAlmostEqual(center.elevation.m, subpoint.elevation.m, delta=1e-3)

    def test_compute_footprint_cross_track_extent_matches_swath_width(self):
        """
        Test that the cross-track extent of a computed footprint (the
        geodesic distance between the footprint edge points directly
        left and right of nadir) matches `get_swath_width`, evaluated at
        the satellite's actual height above the WGS 84 ellipsoid (which,
        due to Earth's oblateness, differs slightly from the orbit's
        nominal mean altitude away from the equator). This cross-checks
        the closed-form swath width formula (assumes a spherical Earth)
        against the WGS 84 ellipsoid footprint geometry.
        """
        o = Instrument(name="Test Instrument", field_of_regard=30.0)
        for sat in (self.test_sat_1, self.test_sat_5):
            orbit_track = sat.at(self.test_time)  # type: ignore
            height = wgs84.geographic_position_of(orbit_track).elevation.m
            expected_swath_width = o.get_swath_width(height)
            # request 5 points so the polygon samples exactly the
            # cross-track edges (angle=0 and angle=180 degrees)
            footprint = o.compute_footprint(orbit_track, number_points=5)
            coords = list(footprint[0].exterior.coords)
            # coords are sampled at angle = [0, 90, 180, 270, 360] degrees;
            # index 0 (angle=0) and index 2 (angle=180) are the cross-track edges
            cross_track_extent = geodesic_distance(
                coords[0][0], coords[0][1], coords[2][0], coords[2][1]
            )
            self.assertAlmostEqual(cross_track_extent, expected_swath_width, delta=50.0)

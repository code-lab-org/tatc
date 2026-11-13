"""
Tests of swept footprints (`PointedInstrument.is_swept`): the strips swept by
an instrument's cross-track line over the time steps of a ground track.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import json
import unittest
from datetime import datetime, timedelta, timezone

import numpy as np
import shapely
from pyproj import Geod, Transformer

from tatc.analysis import (
    collect_ground_track,
    collect_observations,
    collect_region_observations,
    compute_ground_track,
)
from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import (
    Attitude,
    GeneralPerturbationsOrbit,
    PointedInstrument,
    Satellite,
)
from tatc.utils import (
    VelocityFrame,
    ViewGeometry,
    compute_along_track_field_of_view,
    compute_swept_footprints,
)

# NOAA-20, in a sun-synchronous orbit at about 830 km
NOAA_20_TLE = [
    "1 43013U 17073A   26224.00000000  .00000058  00000-0  27779-4 0  2098",
    "2 43013 098.7781 162.9718 0002024 099.5641 067.3984 14.19520627452396",
]
START = datetime(2026, 8, 12, 3, tzinfo=timezone.utc)
DURATION = 600  # seconds
GEOD = Geod(ellps="WGS84")


def get_times(step: float, duration: float = DURATION, start: datetime = START):
    """Gets the times at the middle of each step of a duration."""
    return [
        start + timedelta(seconds=step * (k + 0.5))
        for k in range(int(round(duration / step)))
    ]


def get_satellite(step: float, attitude: Attitude | None = None, **parameters):
    """Gets a satellite with one pointed instrument, tailored to a time step."""
    orbit = GeneralPerturbationsOrbit.from_tle(NOAA_20_TLE)
    parameters.setdefault("cross_track_field_of_view", 20)
    parameters.setdefault("is_rectangular", True)
    instruments = parameters.pop("instruments", None) or [
        PointedInstrument(
            name="Test",
            field_of_regard=180,
            along_track_field_of_view=compute_along_track_field_of_view(
                830e3, step, 98.7
            ),
            **parameters,
        )
    ]
    return Satellite(
        name="Test", orbit=orbit, instruments=instruments, attitude=attitude
    )


def swept(satellite: Satellite) -> Satellite:
    """Gets a copy of a satellite whose instruments are swept."""
    return satellite.model_copy(
        update={
            "instruments": [
                instrument.model_copy(update={"is_swept": True})
                for instrument in satellite.instruments
            ]
        }
    )


def transform_geometry(transform, geometry):
    """
    Transforms the coordinates of a geometry by a function of their x and y
    arrays (such as a pyproj transform), as the deprecated
    `shapely.ops.transform` did.
    """
    return shapely.transform(
        geometry, lambda xy: np.column_stack(transform(xy[:, 0], xy[:, 1]))
    )


def iou(a, b) -> float:
    """Gets the intersection over union of two geometries, in a local projection."""
    center = a.centroid
    forward = Transformer.from_crs(
        "EPSG:4326",
        f"+proj=aeqd +lat_0={center.y} +lon_0={center.x} +units=m",
        always_xy=True,
    ).transform
    a, b = (transform_geometry(forward, g).buffer(0) for g in (a, b))
    return a.intersection(b).area / a.union(b).area


def swath(satellite: Satellite, step: float, **kwargs):
    """Gets the swath of a satellite over the duration, sampled at a time step."""
    track = compute_ground_track(satellite, get_times(step), **kwargs)
    return track.geometry.iloc[0]


class TestSweptFootprints(unittest.TestCase):
    """
    Tests swept footprints of pointed instruments.
    """

    def test_not_swept_by_default(self):
        """
        Test that instruments are not swept by default, so that their footprints
        are their views at each time.
        """
        satellite = get_satellite(10)
        self.assertFalse(satellite.instruments[0].is_swept)
        times = get_times(10)
        track = collect_ground_track(satellite, times)
        orbit_track = satellite.orbit.to_gp_orbit().get_orbit_track(times)
        views = satellite.get_instrument().compute_footprint(orbit_track)
        for footprint, view in zip(track.geometry, views):
            self.assertTrue(shapely.equals_exact(footprint, view, tolerance=1e-9))

    def test_round_trip(self):
        """
        Test that a swept instrument keeps its field through JSON.
        """
        satellite = swept(get_satellite(10))
        loaded = Satellite.model_validate(json.loads(satellite.model_dump_json()))
        self.assertTrue(loaded.instruments[0].is_swept)
        self.assertEqual(loaded, satellite)

    def test_independent_of_time_step(self):
        """
        Test that a swept swath does not depend on the time step, for views
        rolled, pitched and yawed by the body, with a roll profile, in scan
        geometry, and in the inertial frame.
        """
        cases = {
            "nadir": ({}, None),
            "rolled": ({"roll_angle": 25}, None),
            "roll profile": (
                {"roll_angle_profile": [(0, 20), (90, 30), (180, 20), (270, 30)]},
                None,
            ),
            "inertial, yawed and pitched": (
                {},
                Attitude(
                    velocity_frame=VelocityFrame.INERTIAL,
                    yaw_angle=4,
                    pitch_angle=-8,
                ),
            ),
            "scan": (
                {
                    "cross_track_field_of_view": 100,
                    "pitch_angle": 15,
                    "view_geometry": ViewGeometry.SCAN,
                },
                None,
            ),
        }
        for label, (parameters, attitude) in cases.items():
            with self.subTest(label):
                fine = swath(swept(get_satellite(5, attitude, **parameters)), 5)
                coarse = swath(swept(get_satellite(120, attitude, **parameters)), 120)
                self.assertGreater(iou(fine, coarse), 0.999)

    def test_matches_fine_views(self):
        """
        Test that a swept swath at a long time step matches the union of views
        tailored to a short one, which tile the ground track near nadir.
        """
        for parameters in ({}, {"roll_angle": 25}):
            with self.subTest(parameters):
                views = swath(get_satellite(1, **parameters), 1)
                strips = swath(swept(get_satellite(120, **parameters)), 120)
                self.assertGreater(iou(views, strips), 0.995)
                # the views tailored to the long time step reach beyond it
                long_views = swath(get_satellite(120, **parameters), 120)
                self.assertLess(iou(views, long_views), iou(views, strips))

    def test_union_of_instruments(self):
        """
        Test that the swath of several swept instruments is the union of their
        strips at any time step.
        """
        instruments = [
            PointedInstrument(
                name=f"Test {k}",
                field_of_regard=180,
                cross_track_field_of_view=10,
                along_track_field_of_view=1,
                roll_angle=roll,
                is_rectangular=True,
                is_swept=True,
            )
            for k, roll in enumerate((-15, 15))
        ]
        satellite = get_satellite(10, instruments=instruments)
        fine = swath(satellite, 5, instrument_index=None)
        coarse = swath(satellite, 120, instrument_index=None)
        self.assertGreater(iou(fine, coarse), 0.999)
        # the strips of the two instruments do not touch
        self.assertEqual(shapely.get_num_geometries(coarse), 2)

    def test_observed_points_are_in_the_swath(self):
        """
        Test that the points that `collect_observations` observes are those
        within the swept swath, away from its boundary.
        """
        satellite = swept(get_satellite(60, roll_angle=10))
        times = get_times(60, duration=300)
        strip = compute_ground_track(satellite, times).geometry.iloc[0]
        points = generate_points_fibonacci_lattice(50e3)
        points = points[points.intersects(strip.buffer(2))].reset_index(drop=True)
        end = START + timedelta(seconds=300)
        observations = collect_observations(points, satellite, START, end)
        # observations whose view sweeps over the point within the duration
        # (others' epochs are clipped to its bounds)
        observations = observations[
            (observations.epoch > START) & (observations.epoch < end)
        ]
        observed = set(
            np.flatnonzero(
                points.geometry.apply(shapely.to_wkb).isin(
                    observations.geometry.apply(shapely.to_wkb)
                )
            )
        )
        center = strip.centroid
        forward = Transformer.from_crs(
            "EPSG:4326",
            f"+proj=aeqd +lat_0={center.y} +lon_0={center.x} +units=m",
            always_xy=True,
        ).transform
        projected = transform_geometry(forward, strip)
        checked = 0
        for k, point in enumerate(points.geometry):
            p = transform_geometry(forward, point)
            if projected.exterior.distance(p) < 5e3:
                continue
            checked += 1
            self.assertEqual(projected.contains(p), k in observed, point)
        self.assertGreater(checked, 50)

    def test_crosses_the_pole_and_antimeridian(self):
        """
        Test that a wide swept swath over a pole is valid and does not depend on
        the time step.
        """
        attitude = Attitude(velocity_frame=VelocityFrame.INERTIAL)
        satellite = swept(get_satellite(120, attitude, cross_track_field_of_view=110))
        times = get_times(120, duration=3000)
        track = collect_ground_track(satellite, times)
        self.assertTrue(all(shapely.is_valid(track.geometry)))
        self.assertGreater(track.geometry.total_bounds[3], 89)
        areas = []
        for step in (10, 120):
            geometry = compute_ground_track(
                satellite, get_times(step, duration=3000)
            ).geometry.iloc[0]
            areas.append(abs(GEOD.geometry_area_perimeter(geometry)[0]))
        self.assertAlmostEqual(areas[0] / areas[1], 1, delta=2e-3)

    def test_trimmed_to_blind_intervals(self):
        """
        Test that swept footprints are trimmed to the times when the instrument
        observes, so that a swath does not depend on the time step, and a step
        whose time falls in a blind interval keeps its observing part.
        """
        # an interval within the duration (arguments of latitude from about
        # 85 to 121 deg)
        attitude = Attitude(blind_intervals=[(95, 97)])
        fine = swath(swept(get_satellite(2, attitude)), 2)
        coarse = swath(swept(get_satellite(120, attitude)), 120)
        self.assertGreater(iou(fine, coarse), 0.999)
        # the blind interval splits the swath
        self.assertEqual(shapely.get_num_geometries(coarse), 2)
        track = collect_ground_track(
            swept(get_satellite(120, attitude)), get_times(120)
        )
        blind = ~track.geometry.apply(lambda g: g.is_empty)
        self.assertTrue(all(track.valid_obs[blind]))

    def test_gaps_do_not_stretch_strips(self):
        """
        Test that a gap between times does not stretch the strips beside it.
        """
        satellite = swept(get_satellite(60))
        first = get_times(60, duration=300)
        second = get_times(60, duration=300, start=START + timedelta(hours=1))
        joint = compute_ground_track(
            satellite, first + second, dissolve_orbits=True
        ).geometry.iloc[0]
        apart = shapely.union_all(
            [
                compute_ground_track(satellite, times).geometry.iloc[0]
                for times in (first, second)
            ]
        )
        self.assertGreater(iou(joint, apart), 0.9999)

    def test_single_time_is_a_view(self):
        """
        Test that the footprint of a single time is the view at that time.
        """
        satellite = swept(get_satellite(10))
        track = collect_ground_track(satellite, [START])
        orbit_track = satellite.orbit.to_gp_orbit().get_orbit_track(START)
        view = satellite.get_instrument().compute_footprint(orbit_track)[0]
        self.assertTrue(
            shapely.equals_exact(track.geometry.iloc[0], view, tolerance=1e-9)
        )


class TestSweptTimeStep(unittest.TestCase):
    """
    Tests the time step that each time of a swept ground track represents.
    """

    def test_explicit_time_step(self):
        """
        Test that an explicit time step gives the strip of a single time, and
        matches the time step inferred from regular times.
        """
        satellite = swept(get_satellite(120))
        times = get_times(120)
        step = timedelta(seconds=120)
        inferred = compute_ground_track(satellite, times).geometry.iloc[0]
        explicit = compute_ground_track(satellite, times, time_step=step)
        self.assertGreater(iou(inferred, explicit.geometry.iloc[0]), 0.99999)
        # a single time sweeps over its time step
        single = compute_ground_track(satellite, [times[2]], time_step=step)
        middle = compute_ground_track(
            satellite, get_times(5, duration=120, start=times[2] - step / 2)
        ).geometry.iloc[0]
        self.assertGreater(iou(middle, single.geometry.iloc[0]), 0.999)


class TestSweptRegionObservations(unittest.TestCase):
    """
    Tests region observations by swept instruments.
    """

    def test_region_swath_matches_ground_track(self):
        """
        Test that the swath of a region observed by a swept instrument is the
        part of the region within its swept ground track, and that its period
        is that of its cross-track line, not of a view tailored to a long time
        step.
        """
        satellite = swept(get_satellite(120, roll_angle=10))
        strip = compute_ground_track(satellite, get_times(5)).geometry.iloc[0]
        # a region across the middle of the swath, wider than it
        center = strip.centroid
        region = shapely.box(center.x - 15, center.y - 1, center.x + 15, center.y + 1)
        end = START + timedelta(seconds=DURATION)
        observations = collect_region_observations(region, satellite, START, end)
        self.assertEqual(len(observations), 1)
        observation = observations.iloc[0]
        expected = region.intersection(strip)
        self.assertGreater(iou(expected, observation.geometry), 0.995)
        # the period is that in which the cross-track line meets the region
        line = satellite.get_instrument().model_copy(
            update={"along_track_field_of_view": 1e-4}
        )
        orbit = satellite.orbit.to_gp_orbit()
        for time, inside in (
            (observation.start - timedelta(seconds=1), False),
            (observation.start + timedelta(seconds=1), True),
            (observation.end - timedelta(seconds=1), True),
            (observation.end + timedelta(seconds=1), False),
        ):
            view = line.compute_footprint(orbit.get_orbit_track(time))[0]
            self.assertEqual(view.intersects(region), inside, time)
        # a view tailored to 120 s meets it for about a minute longer
        duration = (observation.end - observation.start).total_seconds()
        tailored = collect_region_observations(
            region, get_satellite(120, roll_angle=10), START, end
        ).iloc[0]
        self.assertGreater(
            (tailored.end - tailored.start).total_seconds(), duration + 60
        )


class TestComputeSweptFootprints(unittest.TestCase):
    """
    Tests the projection of swept footprints.
    """

    def setUp(self):
        orbit = GeneralPerturbationsOrbit.from_tle(NOAA_20_TLE)
        self.track = orbit.to_gp_orbit().get_orbit_track(get_times(5, duration=20))

    def test_counts(self):
        """
        Test that spans need at least two times, and counts that sum to the
        number of times.
        """
        with self.assertRaises(ValueError):
            compute_swept_footprints(self.track, [1, 3], 20)
        with self.assertRaises(ValueError):
            compute_swept_footprints(self.track, [2, 3], 20)
        self.assertEqual(compute_swept_footprints(self.track, [], 20), [])

    def test_strips_share_their_bounds(self):
        """
        Test that consecutive spans that share a time meet without a gap.
        """
        orbit = GeneralPerturbationsOrbit.from_tle(NOAA_20_TLE)
        times = get_times(5, duration=20)
        track = orbit.to_gp_orbit().get_orbit_track(times[:3] + times[2:])
        first, second = compute_swept_footprints(track, [3, 2], 20)
        whole = compute_swept_footprints(self.track, [4], 20)[0]
        self.assertGreater(iou(whole, shapely.union_all([first, second])), 0.9999)
        self.assertLess(iou(first, second), 1e-3)

"""
Unit tests for the tatc.analysis.tasking module.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

import geopandas as gpd
import numpy as np
from shapely.geometry import Point
from skyfield.api import wgs84

from tatc.analysis import collect_imaging_opportunities, schedule_tasks
from tatc.analysis.tasking import _get_direction, _Pointing
from tatc.constants import timescale
from tatc.schemas import (
    Agility,
    CircularOrbit,
    LimbSounder,
    PointedInstrument,
    Satellite,
    WalkerConstellation,
)

START = datetime(2026, 6, 1, tzinfo=timezone.utc)
END = START + timedelta(days=1)
IMAGER = PointedInstrument(
    name="imager",
    field_of_regard=60,
    cross_track_field_of_view=1,
    along_track_field_of_view=1,
)
SATELLITE = Satellite(
    name="agile",
    orbit=CircularOrbit(altitude=500e3, inclination=97.4, epoch=START),
    instruments=[IMAGER],
)


def _requests(count=30, seed=1, **columns):
    """Requests at random points over the western United States."""
    rng = np.random.default_rng(seed)
    points = [
        Point(x, y)
        for x, y in zip(rng.uniform(-120, -80, count), rng.uniform(30, 48, count))
    ]
    return gpd.GeoDataFrame(columns, geometry=points, crs="EPSG:4326")


def _off_nadir(roll, pitch):
    return np.degrees(np.arccos(np.cos(np.radians(roll)) * np.cos(np.radians(pitch))))


class TestOpportunities(unittest.TestCase):
    """
    Tests the opportunities to image requests.
    """

    def test_within_field_of_regard(self):
        """
        Tests that opportunities are the periods in which the target is within
        the field of regard: at its edge at their ends, and within it between.
        """
        opportunities = collect_imaging_opportunities(
            SATELLITE, _requests(), START, END
        )
        self.assertGreater(len(opportunities), 0)
        self.assertTrue(opportunities["start"].is_monotonic_increasing)
        for row in opportunities.itertuples():
            pointing = _Pointing(
                SATELLITE,
                IMAGER,
                row.geometry,
                row.start,
                row.end,
                timedelta(seconds=1),
            )
            off_nadir = _off_nadir(pointing.roll, pointing.pitch)
            self.assertLessEqual(off_nadir.max(), 30 + 0.01)
            self.assertAlmostEqual(off_nadir[0], 30, delta=0.01)
            self.assertAlmostEqual(off_nadir[-1], 30, delta=0.01)

    def test_clipped_to_requests(self):
        """
        Tests that opportunities lie between the onsets and expiry of their
        requests.
        """
        rng = np.random.default_rng(2)
        onset = [START + timedelta(hours=h) for h in rng.uniform(0, 12, 30)]
        expires = [time + timedelta(hours=6) for time in onset]
        requests = _requests(onset=onset, expires=expires)
        everything = collect_imaging_opportunities(SATELLITE, _requests(), START, END)
        clipped = collect_imaging_opportunities(SATELLITE, requests, START, END)
        self.assertLess(len(clipped), len(everything))
        for row in clipped.itertuples():
            self.assertGreaterEqual(row.start, requests.onset[row.request])
            self.assertLessEqual(row.end, requests.expires[row.request])

    def test_constellations(self):
        """
        Tests that constellations are their members, and that requests
        without opportunities give an empty frame.
        """
        walker = WalkerConstellation(
            name="walker",
            orbit=SATELLITE.orbit,
            instruments=[IMAGER],
            number_satellites=2,
            number_planes=2,
        )
        opportunities = collect_imaging_opportunities(walker, _requests(), START, END)
        self.assertEqual(
            set(opportunities["satellite"]),
            {member.name for member in walker.generate_members()},
        )
        empty = collect_imaging_opportunities(
            SATELLITE, _requests(), START, START + timedelta(minutes=1)
        )
        self.assertEqual(len(empty), 0)
        self.assertEqual(list(empty.columns), list(opportunities.columns))

    def test_sounders_are_not_pointed(self):
        """
        Tests that a sounder, which views no targets, is not pointed at them,
        and that the opportunities of an imager carried with one are its own.
        """
        sounder = LimbSounder(
            name="limb",
            scan_azimuth=0,
            scan_elevations=[10e3, 50e3],
            scan_duration=timedelta(seconds=20),
        )
        satellite = SATELLITE.model_copy(update={"instruments": [sounder, IMAGER]})
        end = START + timedelta(hours=6)
        with self.assertRaisesRegex(ValueError, "collect_limb_soundings"):
            collect_imaging_opportunities(satellite, _requests(), START, end)
        with self.assertRaisesRegex(ValueError, "collect_limb_soundings"):
            schedule_tasks(
                satellite,
                _requests(),
                START,
                end,
                Agility(max_rate=1),
                dwell=timedelta(seconds=10),
            )
        imaged = collect_imaging_opportunities(
            satellite, _requests(), START, end, instrument_index=1
        )
        expected = collect_imaging_opportunities(SATELLITE, _requests(), START, end)
        self.assertEqual(len(imaged), len(expected))
        self.assertGreater(len(imaged), 0)
        self.assertTrue(
            (imaged["start"].to_numpy() == expected["start"].to_numpy()).all()
        )


class TestPointing(unittest.TestCase):
    """
    Tests the pointing of views.
    """

    def test_directions(self):
        """
        Tests that the directions of views are unit vectors whose angle from
        nadir is the off-nadir angle of their roll and pitch.
        """
        roll, pitch = np.array([0, 20, -15, 10]), np.array([0, 0, 25, -30])
        directions = _get_direction(roll, pitch)
        np.testing.assert_allclose(np.linalg.norm(directions, axis=1), 1)
        np.testing.assert_allclose(
            np.degrees(np.arccos(directions[:, 2])), _off_nadir(roll, pitch)
        )


SLOW = Agility(max_rate=0.5, max_acceleration=0.1, settle_time=10)
DWELL = timedelta(seconds=5)
FAST = Agility(max_rate=3, max_acceleration=1, settle_time=5)


def _dense_requests(count=120, seed=3):
    """Requests crowded over a small region, which compete for time."""
    rng = np.random.default_rng(seed)
    points = [
        Point(x, y)
        for x, y in zip(rng.uniform(-115, -105, count), rng.uniform(33, 43, count))
    ]
    return gpd.GeoDataFrame(
        {"priority": rng.integers(1, 4, count)}, geometry=points, crs="EPSG:4326"
    )


class TestScheduleTasks(unittest.TestCase):
    """
    Tests the schedules of agile satellites.
    """

    def test_slews_fit(self):
        """
        Tests that each satellite has time to slew between its images, which
        lie within their opportunities and the field of regard, and that
        every request is either imaged or unscheduled.
        """
        requests = _dense_requests()
        dwell = timedelta(seconds=5)
        plan = schedule_tasks(SATELLITE, requests, START, END, SLOW, dwell=dwell)
        images = plan.observations
        self.assertGreater(len(images), 1)
        self.assertEqual(len(images) + len(plan.unscheduled), len(requests))
        self.assertFalse(set(images["request"]) & set(plan.unscheduled.index))
        self.assertTrue((images["end"] - images["start"] == dwell).all())
        self.assertLessEqual(images["off_nadir"].max(), 30 + 0.01)
        gaps = (images["start"].iloc[1:].array - images["end"].iloc[:-1].array) / (
            timedelta(seconds=1)
        )
        slews = images["slew_time"].iloc[1:] / timedelta(seconds=1)
        self.assertGreaterEqual(np.min(gaps - slews.to_numpy()), -1e-5)
        for column in ("target_hash", "satellite", "instrument", "epoch"):
            self.assertIn(column, images.columns)
        opportunities = collect_imaging_opportunities(SATELLITE, requests, START, END)
        for image in images.itertuples():
            windows = opportunities[opportunities["request"] == image.request]
            self.assertTrue(
                (
                    (windows["start"] <= image.start) & (image.end <= windows["end"])
                ).any()
            )

    def test_agility(self):
        """
        Tests that a more agile satellite images at least as many requests,
        and that agilities may be given by satellite.
        """
        requests = _dense_requests()
        slow = schedule_tasks(SATELLITE, requests, START, END, SLOW, dwell=DWELL)
        fast = schedule_tasks(SATELLITE, requests, START, END, FAST, dwell=DWELL)
        self.assertGreater(len(fast.observations), len(slow.observations))
        named = schedule_tasks(
            SATELLITE, requests, START, END, {"agile": FAST}, dwell=DWELL
        )
        self.assertEqual(
            list(named.observations["request"]), list(fast.observations["request"])
        )
        with self.assertRaises(ValueError):
            schedule_tasks(
                SATELLITE, requests, START, END, {"other": FAST}, dwell=DWELL
            )
        with self.assertRaises(ValueError):
            schedule_tasks(SATELLITE, requests, START, END, FAST, dwell=timedelta(0))

    def test_priorities(self):
        """
        Tests that the higher priority of two requests that compete for an
        opportunity is imaged, whatever their order.
        """
        requests = _requests(count=2, seed=5)
        opportunities = collect_imaging_opportunities(SATELLITE, requests, START, END)
        # images as long as the opportunities, so that only one fits
        dwell = (opportunities["end"] - opportunities["start"]).min() * 0.9
        for favored in (0, 1):
            prioritized = requests.assign(
                priority=[2 if i == favored else 1 for i in range(2)]
            )
            plan = schedule_tasks(SATELLITE, prioritized, START, END, SLOW, dwell=dwell)
            self.assertIn(favored, list(plan.observations["request"]))

    def test_earliest_opportunities(self):
        """
        Tests that, with nearly instant slews and images, each request is
        imaged at the start of its earliest opportunity.
        """
        requests = _requests()
        instant = Agility(max_rate=1e6)
        plan = schedule_tasks(
            SATELLITE, requests, START, END, instant, dwell=timedelta(microseconds=1)
        )
        opportunities = collect_imaging_opportunities(SATELLITE, requests, START, END)
        earliest = opportunities.groupby("request")["start"].min()
        self.assertEqual(len(plan.observations), len(earliest))
        for image in plan.observations.itertuples():
            self.assertLess(
                abs((image.start - earliest[image.request]).total_seconds()), 1e-3
            )

    def test_attitudes(self):
        """
        Tests that the attitude of a satellite points at each target during
        its image, and at nadir, without observing, between images.
        """
        plan = schedule_tasks(
            SATELLITE, _dense_requests(), START, END, FAST, dwell=timedelta(seconds=5)
        )
        schedule = plan.attitudes["agile"]
        self.assertEqual(len(schedule.segments), 1 + 2 * len(plan.observations))
        self.assertEqual(
            [segment.observing for segment in schedule.segments[:3]],
            [False, True, False],
        )
        (satellite,) = plan.satellites()
        self.assertEqual(satellite.attitude, schedule)
        first = plan.observations.iloc[0]
        target = schedule.get_attitude(first["epoch"].to_pydatetime()).targets[0]
        self.assertAlmostEqual(target.longitude, first.geometry.x)

    def test_attitudes_hold_the_view(self):
        """
        Tests that the attitudes of the plan hold each target in the narrow
        field of view, observing, at the epoch of its image, and do not
        observe while slewing to it.
        """
        plan = schedule_tasks(
            SATELLITE, _dense_requests(), START, END, FAST, dwell=DWELL
        )
        (tasked,) = plan.satellites()
        instrument = tasked.get_instrument(0)
        for image in plan.observations.itertuples():
            target = wgs84.latlon(image.geometry.y, image.geometry.x)
            times = [
                image.epoch.to_pydatetime(),
                (image.start - image.slew_time / 2).to_pydatetime(),
            ]
            track = tasked.get_track(timescale.from_datetimes(times))
            self.assertTrue(instrument.is_in_field_of_view(track, target)[0])
            self.assertEqual(list(instrument.is_observing(track)), [True, False])

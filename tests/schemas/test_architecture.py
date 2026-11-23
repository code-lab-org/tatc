"""
Unit tests for architecture schemas.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
import warnings
from datetime import datetime, timedelta, timezone

from pydantic import ValidationError

from tatc.schemas import (
    ActiveWindow,
    Architecture,
    CircularOrbit,
    GeosynchronousOrbit,
    GroundStation,
    Instrument,
    RadarStation,
    Relay,
    Satellite,
    SurfacePlatform,
    WalkerConstellation,
)

START = datetime(2026, 1, 1, tzinfo=timezone.utc)
DAY = timedelta(days=1)


def _satellite(name):
    return Satellite(
        name=name,
        orbit=CircularOrbit(altitude=500e3, inclination=60),
        instruments=[Instrument(name="i")],
    )


class TestActiveWindow(unittest.TestCase):
    """Tests the windows in which members are active."""

    def test_overlap(self):
        window = ActiveWindow(start=START + DAY, end=START + 3 * DAY)
        self.assertEqual(
            window.overlap(START, START + 2 * DAY), (START + DAY, START + 2 * DAY)
        )
        self.assertEqual(
            window.overlap(START + 2 * DAY, START + 5 * DAY),
            (START + 2 * DAY, START + 3 * DAY),
        )
        self.assertIsNone(window.overlap(START, START + DAY))
        self.assertEqual(
            ActiveWindow().overlap(START, START + DAY), (START, START + DAY)
        )

    def test_contains(self):
        window = ActiveWindow(start=START, end=START + DAY)
        self.assertTrue(window.contains(START))
        self.assertFalse(window.contains(START + DAY))
        # (a time without a time zone is UTC)
        self.assertTrue(window.contains(datetime(2026, 1, 1, 12)))
        self.assertTrue(ActiveWindow(end=START).contains(START - DAY))

    def test_order(self):
        with self.assertRaises(ValidationError):
            ActiveWindow(start=START, end=START)


class TestArchitecture(unittest.TestCase):
    """Tests architectures."""

    def setUp(self):
        self.station = GroundStation(name="g", latitude=70, longitude=20)
        self.radar = RadarStation(name="r", latitude=35, longitude=-97)
        self.tower = SurfacePlatform(
            name="t",
            latitude=40,
            longitude=-105,
            height=100,
            instruments=[Instrument(name="i")],
        )
        self.relay = Relay(
            name="relay",
            satellites=[Satellite(name="geo", orbit=GeosynchronousOrbit(longitude=0))],
        )
        self.architecture = Architecture(
            name="a",
            satellites=[_satellite("s1"), _satellite("s2")],
            surface_platforms=[self.tower],
            radar_stations=[self.radar],
            ground_stations=[self.station],
            relays=[self.relay],
            windows={"s2": ActiveWindow(start=START + DAY)},
        )

    def test_members(self):
        self.assertEqual(
            [m.name for m in self.architecture.platforms], ["s1", "s2", "t", "r"]
        )
        self.assertEqual(
            [m.name for m in self.architecture.members],
            ["s1", "s2", "t", "r", "g", "relay"],
        )
        self.assertIs(
            self.architecture.get_member("g"), self.architecture.ground_stations[0]
        )
        with self.assertRaises(KeyError):
            self.architecture.get_member("x")

    def test_windows(self):
        self.assertFalse(self.architecture.is_active("s2", START))
        self.assertTrue(self.architecture.is_active("s1", START))
        self.assertEqual(
            self.architecture.get_active_period("s2", START, START + 2 * DAY),
            (START + DAY, START + 2 * DAY),
        )
        self.assertIsNone(self.architecture.get_active_period("s2", START, START + DAY))

    def test_distinct_names(self):
        with self.assertRaisesRegex(ValidationError, "distinct names"):
            Architecture(satellites=[_satellite("g")], ground_stations=[self.station])
        with self.assertRaisesRegex(ValidationError, "not in the architecture"):
            Architecture(satellites=[_satellite("s")], windows={"x": ActiveWindow()})

    def test_constellations(self):
        walker = WalkerConstellation(
            name="w",
            orbit=CircularOrbit(altitude=500e3, inclination=60),
            configuration="delta",
            number_satellites=4,
            number_planes=2,
        )
        architecture = Architecture(satellites=[walker])
        self.assertEqual(
            [s.name for s in architecture.satellites], ["w 1", "w 2", "w 3", "w 4"]
        )

    def test_stations_deprecated(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            architecture = Architecture(name="a", stations=[self.station])
        self.assertEqual(architecture.ground_stations, [self.station])
        self.assertTrue(any(w.category is DeprecationWarning for w in caught))
        with self.assertRaises(ValidationError):
            Architecture(stations=[self.station], ground_stations=[self.station])

    def test_with_and_without(self):
        other = GroundStation(name="h", latitude=-70, longitude=20)
        added = self.architecture.with_members(
            other, [_satellite("s3")], windows={"h": ActiveWindow(end=START)}
        )
        self.assertEqual([m.name for m in added.ground_stations], ["g", "h"])
        self.assertEqual([m.name for m in added.satellites], ["s1", "s2", "s3"])
        self.assertEqual(set(added.windows), {"s2", "h"})
        removed = added.without_members("s2", ["h", "relay"])
        self.assertEqual([m.name for m in removed.members], ["s1", "s3", "t", "r", "g"])
        self.assertEqual(removed.windows, {})
        with self.assertRaises(KeyError):
            self.architecture.without_members("x")
        with self.assertRaises(TypeError):
            self.architecture.with_members("x")

    def test_merge(self):
        baseline = Architecture(
            satellites=[_satellite("b")], windows={"b": ActiveWindow(end=START)}
        )
        merged = self.architecture.merge(baseline)
        self.assertEqual(merged.name, "a")
        self.assertEqual([s.name for s in merged.satellites], ["s1", "s2", "b"])
        self.assertEqual(set(merged.windows), {"s2", "b"})

    def test_serialization(self):
        restored = Architecture.model_validate_json(self.architecture.model_dump_json())
        self.assertEqual(restored, self.architecture)


if __name__ == "__main__":
    unittest.main()

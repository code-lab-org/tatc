"""
Unit tests for building blocks of constellations.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, time, timezone

import numpy as np

from tatc.schemas import Instrument, Satellite
from tatc.tradespace import (
    Categorical,
    Choice,
    Derived,
    DesignSpace,
    Integer,
    Real,
    Uncertain,
    build_constellation,
    count_planes,
    count_satellites,
    hours_to_time,
    remove_failures,
    sso_planes_space,
    walker_space,
)

EPOCH = datetime(2026, 1, 1, tzinfo=timezone.utc)
IMAGER = Instrument(name="imager", field_of_regard=60)


def _space():
    return DesignSpace(
        altitude=Real(600e3, 700e3, step=100e3),
        layout=Choice(
            walker=walker_space(
                inclination=Real(45, 60, step=15),
                planes=Integer(1, 3),
                per_plane=Integer(1, 2),
            ),
            sso=sso_planes_space(
                crossing=Real(9, 15, step=3),
                satellites=Integer(1, 2),
                count=Integer(1, 2),
            ),
        ),
        constraints=[lambda d: count_satellites(d.layout) - 4],
    )


class TestFragments(unittest.TestCase):
    """
    Unit tests for fragments of design spaces.
    """

    def test_enumeration(self):
        """Fragments enumerate their layouts, with constraints."""
        designs = _space().enumerate()
        walker = [d for d in designs if d.layout.choice == "walker"]
        sso = [d for d in designs if d.layout.choice == "sso"]
        # 2 altitudes x 2 inclinations x (planes, per plane, spacing) with at most
        # 4 satellites: 1x1, 1x2 (1 each), 2x1, 2x2 (2 each), 3x1 (3)
        self.assertEqual(len(walker), 2 * 2 * (1 + 1 + 2 + 2 + 3))
        # 6 plane designs: one plane (6), or two with different crossings
        # (6 x 6 - same-crossing pairs) / 2 unordered = 12
        self.assertEqual(len(sso), 2 * (6 + 12))
        for design in sso:
            crossings = [plane.crossing for plane in design.layout.planes]
            self.assertEqual(len(set(crossings)), len(crossings))
        self.assertTrue(all(count_satellites(d.layout) <= 4 for d in designs))

    def test_counts(self):
        """Satellites and planes are counted for each kind of layout."""
        designs = _space().enumerate()
        walker = [d for d in designs if d.layout.choice == "walker"][-1]
        self.assertEqual(count_planes(walker.layout), walker.layout.planes)
        self.assertEqual(
            count_satellites(walker.layout),
            walker.layout.planes * walker.layout.per_plane,
        )
        sso = [d for d in designs if d.layout.choice == "sso"][-1]
        self.assertEqual(count_planes(sso.layout), sso.layout.planes.count)
        with self.assertRaises(ValueError):
            count_planes(sso)

    def test_defaults(self):
        """Fragments have default variables."""
        self.assertGreater(DesignSpace(layout=walker_space(altitude=500e3)).size(), 100)
        self.assertGreater(
            DesignSpace(layout=sso_planes_space(altitude=500e3)).size(), 100
        )


class TestBuild(unittest.TestCase):
    """
    Unit tests for building constellations.
    """

    def test_walker(self):
        """Walker layouts build the members of a Walker constellation."""
        space = DesignSpace(
            altitude=650e3,
            layout=walker_space(
                inclination=55, planes=3, per_plane=2, relative_spacing=1
            ),
        )
        design = space.enumerate()[0]
        satellites = build_constellation(
            design.layout, IMAGER, epoch=EPOCH, altitude=design.altitude, name="w"
        )
        self.assertEqual(len(satellites), 6)
        self.assertEqual(satellites[0].name, "w 1")
        self.assertEqual(satellites[0].instruments, [IMAGER])
        raans = {round(s.orbit.right_ascension_ascending_node, 6) for s in satellites}
        self.assertEqual(len(raans), 3)
        self.assertTrue(all(s.orbit.inclination == 55 for s in satellites))

    def test_sso(self):
        """Sun-synchronous layouts build evenly spaced satellites in their planes."""
        space = DesignSpace(
            layout=sso_planes_space(
                crossing=13.5,
                satellites=4,
                count=1,
                phase=Uncertain.uniform(0, 360),
                altitude=700e3,
            ),
        )
        design = space.enumerate()[0].realize(0)
        satellites = build_constellation(design.layout, [IMAGER], epoch=EPOCH)
        self.assertEqual(
            [s.name for s in satellites][:2], ["satellite 1-1", "satellite 1-2"]
        )
        anomalies = [s.orbit.true_anomaly for s in satellites]
        gaps = np.diff(np.unwrap(np.radians(anomalies)))
        np.testing.assert_allclose(np.degrees(gaps) % 360, 90)
        self.assertEqual(satellites[0].orbit.equator_crossing_time, time(13, 30))
        self.assertAlmostEqual(anomalies[0], design.layout.planes[0].phase)

    def test_derived_planes(self):
        """Plane values may be derived from their index and enclosing values."""
        space = DesignSpace(
            per_plane=Integer(1, 2),
            layout=sso_planes_space(
                crossing=Derived(lambda d: (13.5 + 12 * d.index / d.planes.count) % 24),
                satellites=Derived(lambda d: d.per_plane),
                count=Integer(1, 3),
                ordered=True,
            ),
        )
        design = [d for d in space.enumerate() if len(d.layout.planes) == 3][-1]
        satellites = build_constellation(
            design.layout, IMAGER, epoch=EPOCH, altitude=830e3
        )
        self.assertEqual(len(satellites), 3 * design.per_plane)
        times = sorted({s.orbit.equator_crossing_time for s in satellites})
        self.assertEqual(times, [time(13, 30), time(17, 30), time(21, 30)])

    def test_catalog_payloads(self):
        """Catalog entries are built with their validated attitude."""
        space = DesignSpace(
            sounder=Categorical(["atms", "cris"]),
            layout=sso_planes_space(crossing=13.5, satellites=1, count=1),
        )
        design = space.enumerate()[0]
        satellite = build_constellation(
            design.layout, design.sounder, epoch=EPOCH, altitude=830e3, time_step=30
        )[0]
        self.assertIsInstance(satellite, Satellite)
        self.assertIsNotNone(satellite.attitude)
        both = build_constellation(
            design.layout,
            ["atms:swath", "cris:swath"],
            epoch=EPOCH,
            altitude=830e3,
            time_step=30,
        )[0]
        self.assertEqual(len(both.instruments), 2)
        custom = build_constellation(
            design.layout,
            lambda name, orbit: Satellite(name=name, orbit=orbit, instruments=[IMAGER]),
            epoch=EPOCH,
            altitude=830e3,
        )[0]
        self.assertEqual(custom.instruments, [IMAGER])
        with self.assertRaises(ValueError):
            build_constellation(design.layout, IMAGER, epoch=EPOCH)

    def test_hours_to_time(self):
        """Hours convert to times of day, wrapping and rounding to seconds."""
        self.assertEqual(hours_to_time(13.5), time(13, 30))
        self.assertEqual(hours_to_time(15.9), time(15, 54))
        self.assertEqual(hours_to_time(24.25), time(0, 15))
        self.assertEqual(hours_to_time(10.9999999), time(11, 0))


class TestFailures(unittest.TestCase):
    """
    Unit tests for failures of satellites.
    """

    def test_remove_failures(self):
        """Failures depend on the seed and position, common to every design."""
        satellites = [
            Satellite(name=f"s{k}", orbit=_orbit(), instruments=[IMAGER])
            for k in range(200)
        ]
        survivors = remove_failures(satellites, 0.25, seed=7)
        self.assertAlmostEqual(len(survivors) / 200, 0.75, delta=0.08)
        # the first 10 satellites fail alike in a design of 10
        first = {s.name for s in remove_failures(satellites[:10], 0.25, seed=7)}
        self.assertEqual(first, {s.name for s in survivors if int(s.name[1:]) < 10})
        self.assertEqual(len(remove_failures(satellites, 0.25, seed=None)), 200)
        self.assertEqual(remove_failures(satellites, 1.0, seed=3), [])
        with self.assertRaises(ValueError):
            remove_failures(satellites, 1.5, seed=3)

    def test_seed(self):
        """Uncertain seeds are drawn for each replicate, common to designs."""
        space = DesignSpace(a=Integer(1, 2), seed=Uncertain.seed())
        a, b = space.enumerate()
        self.assertIsNone(a.seed)
        self.assertEqual(a.realize(3).seed, b.realize(3).seed)
        self.assertNotEqual(a.realize(3).seed, a.realize(4).seed)


def _orbit():
    from tatc.schemas import (
        SunSynchronousOrbit,
    )  # pylint: disable=import-outside-toplevel

    return SunSynchronousOrbit(
        altitude=700e3, equator_crossing_time="10:30", epoch=EPOCH
    )


def _walker(**kwargs):
    """The satellites of a Walker layout (on its own) of given values."""
    values = {
        "inclination": 55.0,
        "planes": 3,
        "per_plane": 4,
        "relative_spacing": 1,
        "configuration": "delta",
        **kwargs,
    }
    space = DesignSpace(layout=walker_space(**values))
    (design,) = space.enumerate()
    return design, build_constellation(
        design.layout, IMAGER, epoch=EPOCH, altitude=700e3
    )


def _orbits(satellites):
    """The (node, anomaly) of each satellite, as a set, ignoring names."""
    return {
        (
            round(s.orbit.right_ascension_ascending_node % 360, 6) % 360,
            round(s.orbit.true_anomaly % 360, 6) % 360,
        )
        for s in satellites
    }


class TestWalkerPhase(unittest.TestCase):
    """
    Unit tests for the phase and node of the lead satellite of Walker layouts.
    """

    def test_default(self):
        """The lead satellite is at the origin by default."""
        design, satellites = _walker()
        self.assertEqual(design.layout.true_anomaly, 0)
        self.assertEqual(design.layout.right_ascension, 0)
        self.assertEqual(satellites[0].orbit.true_anomaly, 0)
        self.assertEqual(satellites[0].orbit.right_ascension_ascending_node, 0)

    def test_fractions(self):
        """Fractions give angles of the intervals of distinct constellations."""
        design, satellites = _walker(phase_fraction=0.5, node_fraction=0.25)
        self.assertAlmostEqual(design.layout.true_anomaly, 0.5 * 360 / 4)
        self.assertAlmostEqual(design.layout.right_ascension, 0.25 * 360 / 3)
        self.assertAlmostEqual(satellites[0].orbit.true_anomaly, 45)
        self.assertAlmostEqual(satellites[0].orbit.right_ascension_ascending_node, 30)
        star, _ = _walker(configuration="star", node_fraction=0.25)
        self.assertAlmostEqual(star.layout.right_ascension, 90)

    def test_symmetries(self):
        """A whole interval of phase or node gives the same constellation."""
        _, base = _walker()
        _, phased = _walker(phase_fraction=1.0)
        self.assertEqual(_orbits(phased), _orbits(base))
        # a shift of node by 360 / planes is a shift of phase by
        # -360 * relative_spacing / satellites (a fraction of -1 / planes)
        _, noded = _walker(node_fraction=1.0)
        _, shifted = _walker(phase_fraction=-1 / 3)
        self.assertEqual(_orbits(noded), _orbits(shifted))
        self.assertNotEqual(_orbits(noded), _orbits(base))
        # fractions within the intervals give other constellations
        _, half = _walker(phase_fraction=0.5)
        self.assertNotEqual(_orbits(half), _orbits(base))

    def test_uncertain_and_designed(self):
        """Fractions are design variables or uncertain parameters."""
        space = DesignSpace(
            layout=walker_space(
                inclination=55.0,
                planes=Integer(1, 2),
                per_plane=2,
                relative_spacing=0,
                phase_fraction=Uncertain.uniform(0, 1),
                node_fraction=Real(0, 0.5, step=0.5),
            )
        )
        designs = space.enumerate()
        self.assertEqual(len(designs), 4)
        realized = designs[-1].realize(0)
        fraction = realized.layout.phase_fraction
        self.assertTrue(0 <= fraction < 1)
        self.assertAlmostEqual(realized.layout.true_anomaly, fraction * 180)
        self.assertAlmostEqual(realized.layout.right_ascension, 0.5 * 180)
        # the same draw in every design (common random numbers)
        self.assertEqual(designs[0].realize(0).layout.phase_fraction, fraction)

import unittest
import warnings
from datetime import datetime, timedelta, timezone

import numpy as np
import pandas as pd

from tatc.analysis import (
    collect_downlinks,
    collect_observations,
    collect_relay_contacts,
    compute_latencies,
)
from tatc.analysis.data_flow import (
    analyze_data_flow,
    compute_data_flow,
    compute_generation,
    is_data_flow_modeled,
)
from tatc.schemas import (
    CircularOrbit,
    DataFlow,
    DataGeneration,
    GeosynchronousOrbit,
    GroundStation,
    Instrument,
    InstrumentResources,
    Point,
    Relay,
    RelayCommunications,
    Satellite,
    SpacecraftCommunications,
    SpacecraftResources,
    StationCommunications,
    SunSynchronousOrbit,
    Transmitter,
)

START = datetime(2025, 1, 1, tzinfo=timezone.utc)
STEP = timedelta(seconds=10)


def _satellite(name, generation=1e6, downlink=10e6, storage=None, crosslink=None):
    transmitters = [Transmitter(band="X", max_data_rate=downlink)]
    if crosslink is not None:
        transmitters.append(
            Transmitter(band="Ka", link="crosslink", max_data_rate=crosslink)
        )
    return Satellite(
        name=name,
        orbit=CircularOrbit(altitude=500e3, inclination=50, epoch=START),
        resources=SpacecraftResources(storage=storage) if storage else None,
        communications=SpacecraftCommunications(
            transmitters=transmitters, generation=DataGeneration(rate=generation)
        ),
    )


def _station(name, antennas=None, min_access_time=timedelta(0)):
    return GroundStation(
        name=name,
        latitude=0,
        longitude=0,
        min_access_time=min_access_time,
        communications=(StationCommunications(antennas=antennas) if antennas else None),
    )


def _samples(contacts):
    """
    Builds samples (see `sample_link_geometry`) of contacts, each a tuple of
    the satellite, the ground station (or relay and relay satellite), and the
    first and last cells.
    """
    rows = []
    for index, (satellite, endpoint, first, last) in enumerate(contacts):
        station, relay, relay_satellite = (
            (endpoint, None, None) if isinstance(endpoint, str) else (None, *endpoint)
        )
        for cell in range(first, last + 1):
            rows.append(
                {
                    "contact": index,
                    "satellite": satellite,
                    "ground_station": station,
                    "relay": relay,
                    "relay_satellite": relay_satellite,
                    "cell": cell,
                    "time": START + STEP * (cell + 0.5),
                    "fraction": 1.0,
                    "distance": 1e6,
                    "elevation": 45.0 if relay is None else np.nan,
                }
            )
    return pd.DataFrame(rows)


def _analyze(contacts, satellites, stations=(), relays=(), cells=100, **model):
    return analyze_data_flow(
        pd.DataFrame(),
        satellites,
        stations,
        relays,
        start=START,
        end=START + STEP * cells,
        model=DataFlow(step=STEP, **model),
        samples=_samples(contacts),
    )


def _find(result, satellite, ends):
    observations = pd.DataFrame(
        {"satellite": satellite, "end": [START + STEP * end for end in ends]}
    )
    return result.find(observations)


class TestStorage(unittest.TestCase):
    def test_single_contact(self):
        """
        Tests a satellite that stores data until one contact delivers it, first
        in, first out.
        """
        # 1 Mbit/s for 100 cells of 10 s; 10 Mbit/s in cells 50 to 59
        result = _analyze([("a", "g", 50, 59)], [_satellite("a")], [_station("g")])
        summary = result.summary.iloc[0]
        self.assertAlmostEqual(summary["generated"], 1e9)
        # 500 Mbit stored at the contact's start, and 100 Mbit generated during
        # it, delivered at 10 Mbit/s for 60 s (600 Mbit)
        self.assertAlmostEqual(summary["delivered"], 6e8)
        self.assertAlmostEqual(summary["final_backlog"], 4e8)
        self.assertAlmostEqual(summary["peak_backlog"], 5e8 + 5e8 - 5e8, delta=1)
        # data observed at 100 s (cell 10) is the 100th Mbit: delivered 10 s
        # into the contact (at 510 s)
        found = _find(result, "a", [10, 49, 99])
        self.assertEqual(found["downlinked"][0], START + timedelta(seconds=510))
        self.assertEqual(found["ground_station"][0], "g")
        # the 490th Mbit, 49 s into the contact
        self.assertEqual(found["downlinked"][1], START + timedelta(seconds=549))
        # not delivered by the end
        self.assertTrue(pd.isna(found["downlinked"][2]))
        self.assertFalse(found["lost"].any())

    def test_streaming(self):
        """
        Tests that data generated during a contact with an empty storage is
        delivered as it is generated.
        """
        result = _analyze([("a", "g", 0, 99)], [_satellite("a")], [_station("g")])
        found = _find(result, "a", [10.5, 50])
        np.testing.assert_array_equal(
            found["downlinked"].to_numpy(),
            pd.to_datetime([START + STEP * 10.5, START + STEP * 50]).to_numpy(),
        )
        self.assertAlmostEqual(result.summary["final_backlog"][0], 0)

    def test_initial_backlog(self):
        """
        Tests that data stored at the start is delivered first.
        """
        result = _analyze(
            [("a", "g", 0, 99)],
            [_satellite("a")],
            [_station("g")],
            initial_backlog=90e6,
        )
        # 90 Mbit stored and 1 Mbit/s generated, delivered at 10 Mbit/s: the
        # data of 0 s is delivered at 9 s
        found = _find(result, "a", [0])
        self.assertEqual(found["downlinked"][0], START + timedelta(seconds=9))

    def test_overflow(self):
        """
        Tests that data generated while storage is full is lost.
        """
        result = _analyze(
            [("a", "g", 80, 99)],
            [_satellite("a", storage=300e6)],
            [_station("g")],
        )
        summary = result.summary.iloc[0]
        # storage fills at 300 s; 500 Mbit is lost before the contact
        self.assertAlmostEqual(summary["lost"], 5e8)
        self.assertAlmostEqual(summary["peak_backlog"], 3e8)
        self.assertAlmostEqual(
            summary["generated"],
            summary["delivered"] + summary["lost"] + summary["final_backlog"],
        )
        found = _find(result, "a", [10, 40, 70])
        np.testing.assert_array_equal(found["lost"].to_numpy(), [False, True, True])
        self.assertTrue(pd.notna(found["downlinked"][0]))

    def test_storage_series(self):
        """
        Tests the time series of storage.
        """
        result = _analyze([("a", "g", 50, 59)], [_satellite("a")], [_station("g")])
        storage = result.storage
        self.assertEqual(len(storage), 100)
        self.assertEqual(storage["time"].iloc[0], START + STEP)
        self.assertAlmostEqual(storage["backlog"].iloc[49], 5e8)
        self.assertAlmostEqual(storage["generated"].iloc[-1], 1e9)


class TestContention(unittest.TestCase):
    def setUp(self):
        self.satellites = [_satellite("a"), _satellite("b")]

    def test_one_antenna(self):
        """
        Tests that a ground station with one antenna serves one satellite at a
        time, the first to wait first, and the other once the first's storage
        empties.
        """
        result = _analyze(
            [("a", "g", 50, 99), ("b", "g", 50, 99)],
            self.satellites,
            [_station("g", antennas=1)],
        )
        schedule = result.schedule
        first = schedule.iloc[0]
        self.assertEqual(first["satellite"], "a")
        # a's 500 Mbit takes 500 / (10 - 1) s, so b starts in cell 56
        self.assertEqual(schedule.iloc[1]["satellite"], "b")
        self.assertEqual(schedule.iloc[1]["start"], START + STEP * 56)
        self.assertAlmostEqual(result.utilization["busy_fraction"][0], 0.5, delta=0.01)

    def test_hold_to_end(self):
        """
        Tests that with links released at the end of their contacts, the first
        satellite holds the antenna for its whole contact.
        """
        result = _analyze(
            [("a", "g", 50, 79), ("b", "g", 50, 99)],
            self.satellites,
            [_station("g", antennas=1)],
            release="end",
        )
        schedule = result.schedule
        self.assertEqual(list(schedule["satellite"]), ["a", "b"])
        self.assertEqual(schedule.iloc[0]["end"], START + STEP * 80)
        self.assertEqual(schedule.iloc[1]["start"], START + STEP * 80)

    def test_two_antennas(self):
        """
        Tests that a ground station with two antennas serves both at once.
        """
        result = _analyze(
            [("a", "g", 50, 99), ("b", "g", 50, 99)],
            self.satellites,
            [_station("g", antennas=2)],
        )
        np.testing.assert_array_equal(
            result.schedule["start"].to_numpy(),
            pd.to_datetime([START + STEP * 50] * 2).to_numpy(),
        )

    def test_unlimited(self):
        """
        Tests that a ground station without communications serves any number.
        """
        result = _analyze(
            [("a", "g", 50, 99), ("b", "g", 50, 99)], self.satellites, [_station("g")]
        )
        self.assertEqual(len(result.schedule), 2)
        self.assertTrue(np.isnan(result.utilization["busy_fraction"][0]))

    def test_one_link_per_transmitter(self):
        """
        Tests that a transmitter sends over one link at a time, handing over
        from one ground station to the next.
        """
        result = _analyze(
            [("a", "g", 50, 59), ("a", "h", 55, 69)],
            [_satellite("a", generation=10e6)],
            [_station("g"), _station("h")],
            release="end",
        )
        schedule = result.schedule
        self.assertEqual(list(schedule["ground_station"]), ["g", "h"])
        self.assertEqual(schedule.iloc[1]["start"], START + STEP * 60)

    def test_min_access_time(self):
        """
        Tests that a link starts only if its contact lasts the minimum access
        time.
        """
        station = _station("g", min_access_time=timedelta(seconds=150))
        result = _analyze(
            [("a", "g", 10, 19), ("a", "g", 50, 69)], [_satellite("a")], [station]
        )
        self.assertEqual(len(result.schedule), 1)
        self.assertEqual(result.schedule.iloc[0]["start"], START + STEP * 50)

    def test_min_backlog(self):
        """
        Tests that a satellite requests a link only once more than the minimum
        backlog is stored.
        """
        result = _analyze(
            [("a", "g", 0, 99)],
            [_satellite("a")],
            [_station("g")],
            min_backlog=100e6,
        )
        # 100 Mbit is stored after 100 s, then delivered in about 11 s
        self.assertEqual(result.schedule.iloc[0]["start"], START + STEP * 10)
        self.assertGreater(len(result.schedule), 1)

    def test_feeder_rate(self):
        """
        Tests that a relay satellite's feeder rate is shared by its links.
        """
        relay = Relay(
            name="r",
            satellites=[Satellite(name="r1", orbit=GeosynchronousOrbit(longitude=0))],
            communications=RelayCommunications(antennas=2, feeder_rate=12e6),
        )
        satellites = [
            _satellite("a", generation=10e6, crosslink=10e6),
            _satellite("b", generation=10e6, crosslink=10e6),
        ]
        result = _analyze(
            [("a", ("r", "r1"), 0, 99), ("b", ("r", "r1"), 0, 99)],
            satellites,
            relays=[relay],
        )
        delivered = result.summary["delivered"].sum()
        self.assertAlmostEqual(delivered, 12e6 * 1000, delta=1e6)
        self.assertEqual(set(result.utilization["relay_satellite"]), {"r1"})


class TestGeneration(unittest.TestCase):
    def test_activity(self):
        """
        Tests generation while always active, sunlit, or over daylight.
        """
        end = START + timedelta(days=1)
        satellite = Satellite(
            name="s",
            orbit=SunSynchronousOrbit(
                altitude=600e3, equator_crossing_time="10:30", epoch=START
            ),
            instruments=[
                Instrument(name="i", resources=InstrumentResources(data_rate=2e6))
            ],
        )
        always = compute_generation(satellite, START, end, timedelta(minutes=1))
        self.assertEqual(len(always), 1440)
        self.assertAlmostEqual(always.sum(), 2e6 * 86400)
        fractions = {}
        for activity in ["sunlit", "day"]:
            active = satellite.model_copy(
                update={
                    "communications": SpacecraftCommunications(
                        generation=DataGeneration(activity=activity, duty_cycle=0.5)
                    )
                }
            )
            volume = compute_generation(active, START, end, timedelta(minutes=1))
            fractions[activity] = volume.sum() / (0.5 * always.sum())
        # sunlit about 64% of a 600 km orbit at 10:30, over daylight about half
        self.assertGreater(fractions["sunlit"], 0.55)
        self.assertLess(fractions["sunlit"], 0.75)
        self.assertAlmostEqual(fractions["day"], 0.5, delta=0.05)

    def test_partial_last_cell(self):
        """
        Tests that the last cell ends at the end.
        """
        satellite = _satellite("a")
        volume = compute_generation(
            satellite, START, START + timedelta(seconds=25), STEP
        )
        np.testing.assert_allclose(volume, [1e7, 1e7, 5e6])


class TestComputeDataFlow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.end = START + timedelta(hours=12)
        orbit = SunSynchronousOrbit(
            altitude=600e3, equator_crossing_time="10:30", epoch=START
        )
        cls.satellite = Satellite(
            name="s",
            orbit=orbit,
            instruments=[
                Instrument(
                    name="i",
                    field_of_regard=60,
                    resources=InstrumentResources(data_rate=1e6),
                )
            ],
            resources=SpacecraftResources(downlink_rate=1e12),
        )
        cls.station = GroundStation(
            name="g", latitude=78.23, longitude=15.4, min_elevation_angle=5
        )
        cls.points = [
            Point(id=i, latitude=lat, longitude=lon)
            for i, (lat, lon) in enumerate([(0, 0), (40, -100), (-30, 140), (60, 30)])
        ]
        cls.observations = collect_observations(
            cls.points, cls.satellite, START, cls.end
        )
        cls.downlinks = collect_downlinks(
            [cls.station], [cls.satellite], START, cls.end
        )

    def test_modeled(self):
        """
        Tests which satellites' data flow can be modeled: from instrument data
        rates and a resources' downlink rate.
        """
        self.assertTrue(is_data_flow_modeled(self.satellite))
        self.assertFalse(
            is_data_flow_modeled(self.satellite.model_copy(update={"resources": None}))
        )

    def test_fast_links(self):
        """
        Tests that with fast links and no contention, data is delivered no
        later than by timing alone, delivering at the end of a contact.
        """
        latencies, result = compute_data_flow(
            self.observations,
            self.downlinks,
            [self.satellite],
            [self.station],
            start=START,
            end=self.end,
            model=DataFlow(step=timedelta(seconds=10)),
        )
        timing = compute_latencies(self.observations, self.downlinks)
        self.assertEqual(len(latencies), len(timing))
        self.assertIn("lost", latencies.columns)
        self.assertIn("relay", latencies.columns)
        merged = latencies.merge(
            timing[["observed", "latency"]], on="observed", suffixes=("", "_timing")
        )
        delivered = merged["latency"].notna()
        self.assertTrue(delivered.any())
        self.assertTrue(
            (
                merged.loc[delivered, "latency"]
                <= merged.loc[delivered, "latency_timing"] + timedelta(seconds=10)
            ).all()
        )
        # delivered within seconds of the start of its contact
        self.assertTrue((latencies["ground_station"].dropna() == "g").all())
        self.assertFalse(latencies["lost"].any())
        self.assertEqual(result.satellites, ["s"])

    def test_timing_fallback(self):
        """
        Tests that the latencies of a satellite whose data flow cannot be
        modeled are found by timing alone, or raise if strict.
        """
        satellite = self.satellite.model_copy(update={"resources": None})
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            latencies, result = compute_data_flow(
                self.observations,
                self.downlinks,
                [satellite],
                [self.station],
                start=START,
                end=self.end,
            )
        self.assertTrue(any("cannot be modeled" in str(w.message) for w in caught))
        timing = compute_latencies(self.observations, self.downlinks)
        pd.testing.assert_series_equal(
            latencies["latency"].reset_index(drop=True),
            timing["latency"].reset_index(drop=True),
        )
        self.assertEqual(result.satellites, [])
        with self.assertRaises(ValueError):
            compute_data_flow(
                self.observations,
                self.downlinks,
                [satellite],
                [self.station],
                start=START,
                end=self.end,
                model=DataFlow(strict=True),
            )

    def test_relay(self):
        """
        Tests delivery through a relay's crosslink, faster than waiting for a
        ground station.
        """
        relay = Relay(
            name="r",
            satellites=[
                Satellite(
                    name=f"r{i}",
                    orbit=GeosynchronousOrbit(longitude=longitude, epoch=START),
                )
                for i, longitude in enumerate([-41, 174, 62])
            ],
        )
        satellite = self.satellite.model_copy(
            update={
                "communications": SpacecraftCommunications(
                    transmitters=[
                        Transmitter(band="X", max_data_rate=100e6),
                        Transmitter(band="Ka", link="crosslink", max_data_rate=50e6),
                    ]
                )
            }
        )
        contacts = pd.concat(
            [
                self.downlinks.assign(relay=None),
                collect_relay_contacts([relay], [satellite], START, self.end),
            ],
            ignore_index=True,
        )
        latencies, _ = compute_data_flow(
            self.observations,
            contacts,
            [satellite],
            [self.station],
            [relay],
            start=START,
            end=self.end,
        )
        self.assertTrue((latencies["relay"] == "r").any())
        self.assertLess(latencies["latency"].max(), timedelta(minutes=10))

    def test_empty(self):
        """
        Tests no observations.
        """
        latencies, _ = compute_data_flow(
            self.observations.iloc[:0],
            self.downlinks,
            [self.satellite],
            [self.station],
            start=START,
            end=self.end,
        )
        self.assertTrue(latencies.empty)
        self.assertIn("lost", latencies.columns)

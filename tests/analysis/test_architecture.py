"""
Unit tests for the evaluation of architectures.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

import pandas as pd
from shapely.geometry import Point as ShapelyPoint
from shapely.geometry import box

from tatc.analysis import (
    collect_downlinks,
    collect_footprint_observations,
    collect_footprints,
    collect_observations,
    collect_platform_observations,
    collect_region_observations,
    collect_relay_contacts,
    collect_space_observations,
    compute_latencies,
    evaluate_architecture,
)
from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import (
    ActiveWindow,
    Aircraft,
    Architecture,
    DataFlow,
    DataGeneration,
    GeosynchronousOrbit,
    GroundStation,
    Instrument,
    RadarStation,
    Relay,
    RelayCommunications,
    Satellite,
    SpacecraftCommunications,
    SpacecraftResources,
    StationCommunications,
    SunSynchronousOrbit,
    SurfacePlatform,
    Trajectory,
    Transmitter,
)

START = datetime(2026, 1, 1, tzinfo=timezone.utc)
END = START + timedelta(hours=12)
HOUR = timedelta(hours=1)


def _satellite(name, crossing="10:30"):
    return Satellite(
        name=name,
        orbit=SunSynchronousOrbit(altitude=600e3, equator_crossing_time=crossing),
        instruments=[Instrument(name="imager", field_of_regard=90)],
    )


def _metrics(observations):
    return {
        "observations": len(observations),
        "targets": observations["target_hash"].nunique(),
        "latency": (
            observations["latency"].mean() if "latency" in observations else None
        ),
        "label": "x",
    }


def _sorted(frame):
    """A frame's rows in a canonical order, for comparisons."""
    columns = [c for c in ["start", "target_hash", "satellite"] if c in frame]
    return frame.sort_values(columns).reset_index(drop=True)


class TestPlatformObservations(unittest.TestCase):
    """Tests the observations of one platform of any kind."""

    def setUp(self):
        self.points = generate_points_fibonacci_lattice(count=100)
        self.satellite = _satellite("s")

    def test_satellite(self):
        observed = collect_platform_observations(
            self.points, self.satellite, START, END
        )
        expected = collect_observations(self.points, self.satellite, START, END)
        pd.testing.assert_frame_equal(
            _sorted(observed), _sorted(expected), check_like=True
        )

    def test_satellite_regions(self):
        region = box(-20, -20, 20, 20)
        targets = [*self.points.geometry[:5], region]
        observed = collect_platform_observations(targets, self.satellite, START, END)
        regions = collect_region_observations([region], self.satellite, START, END)
        points = collect_observations(
            list(self.points.geometry[:5]), self.satellite, START, END
        )
        self.assertEqual(len(observed), len(regions) + len(points))
        self.assertTrue(observed["start"].is_monotonic_increasing)

    def test_fixed_platform(self):
        tower = SurfacePlatform(
            name="t",
            latitude=35,
            longitude=-110,
            height=500,
            instruments=[Instrument(name="i", field_of_regard=170)],
        )
        point = ShapelyPoint(-110.001, 35.001)
        observed = collect_platform_observations([point], tower, START, END)
        # a fixed platform observes a target in view throughout the period
        self.assertEqual(len(observed), 1)
        self.assertEqual(observed["start"].iloc[0], pd.Timestamp(START))
        self.assertEqual(observed["end"].iloc[0], pd.Timestamp(END))
        stepped = collect_platform_observations(
            [point], tower, START, END, time_step=HOUR
        )
        times = list(pd.date_range(START, END, freq=HOUR).to_pydatetime())
        expected = collect_footprint_observations(
            [point], collect_footprints(tower, times, time_step=HOUR), platforms=tower
        )
        pd.testing.assert_frame_equal(stepped, expected)

    def test_radar_station(self):
        radar = RadarStation(name="r", latitude=35, longitude=-97)
        point = ShapelyPoint(-97.5, 35.5, 3000)
        with self.assertRaises(ValueError):
            collect_platform_observations([point], radar, START, END)
        observed = collect_platform_observations(
            [point], radar, START, END, elevation=3000
        )
        self.assertEqual(list(observed["radar_station"]), ["r"])

    def test_aircraft(self):
        trajectory = Trajectory.from_waypoints(
            [(-110, 35), (-109, 35)],
            START + HOUR,
            speed=150,
            altitude=8000,
            time_step=timedelta(seconds=10),
        )
        aircraft = Aircraft(
            name="a",
            trajectory=trajectory,
            instruments=[Instrument(name="i", field_of_regard=60)],
        )
        point = ShapelyPoint(-109.5, 35)
        observed = collect_platform_observations([point], aircraft, START, END)
        self.assertEqual(len(observed), 1)
        # (outside the period, none)
        self.assertEqual(
            len(
                collect_platform_observations(
                    [point], aircraft, START, START + HOUR / 2
                )
            ),
            0,
        )
        # sampled at time steps within the trajectory
        stepped = collect_platform_observations(
            [point], aircraft, START, END, time_step=timedelta(seconds=30)
        )
        self.assertEqual(len(stepped), 1)

    def test_invalid(self):
        with self.assertRaises(TypeError):
            collect_platform_observations(
                self.points,
                GroundStation(name="g", latitude=0, longitude=0),
                START,
                END,
            )


class TestRelayContacts(unittest.TestCase):
    """Tests the contacts of satellites with relays."""

    def setUp(self):
        self.satellite = _satellite("s")
        self.geo = [
            Satellite(name="east", orbit=GeosynchronousOrbit(longitude=-41)),
            Satellite(name="west", orbit=GeosynchronousOrbit(longitude=-174)),
        ]

    def test_line_of_sight(self):
        relay = Relay(name="r", satellites=self.geo, min_grazing_altitude=100e3)
        contacts = collect_relay_contacts(relay, self.satellite, START, END)
        sight = collect_space_observations(
            self.satellite,
            START,
            END,
            to_satellites=self.geo,
            min_grazing_altitude=100e3,
        )
        self.assertEqual(len(contacts), len(sight))
        self.assertTrue((contacts["relay"] == "r").all())
        self.assertTrue(contacts["ground_station"].isna().all())
        self.assertTrue(
            (
                contacts["epoch"]
                == contacts["start"] + (contacts["end"] - contacts["start"]) / 2
            ).all()
        )

    def test_ground_stations(self):
        # a station in view of one relay satellite only
        station = GroundStation(
            name="wsc", latitude=0, longitude=-10, min_elevation_angle=10
        )
        relay = Relay(
            name="r",
            satellites=self.geo,
            ground_stations=[station],
            min_grazing_altitude=100e3,
        )
        contacts = collect_relay_contacts(relay, self.satellite, START, END)
        self.assertGreater(len(contacts), 0)
        self.assertEqual(set(contacts["relay_satellite"]), {"east"})
        self.assertTrue((contacts["ground_station"] == "wsc").all())

    def test_min_access_time(self):
        relay = Relay(
            name="r", satellites=self.geo, min_access_time=timedelta(minutes=30)
        )
        contacts = collect_relay_contacts(relay, self.satellite, START, END)
        self.assertTrue(
            (contacts["end"] - contacts["start"] >= timedelta(minutes=30)).all()
        )

    def test_empty(self):
        relay = Relay(name="r", satellites=self.geo, max_range=1)
        contacts = collect_relay_contacts(relay, self.satellite, START, END)
        self.assertEqual(len(contacts), 0)
        self.assertIn("relay", contacts.columns)


class TestEvaluateArchitecture(unittest.TestCase):
    """Tests the evaluation of architectures."""

    def setUp(self):
        self.points = generate_points_fibonacci_lattice(count=100)
        self.satellites = [_satellite("s1"), _satellite("s2", "13:30")]
        self.stations = [
            GroundStation(name="svalbard", latitude=78.2, longitude=15.4),
            GroundStation(name="mcmurdo", latitude=-77.8, longitude=166.7),
        ]
        self.relay = Relay(
            name="relay",
            satellites=[
                Satellite(name="geo", orbit=GeosynchronousOrbit(longitude=-41))
            ],
            min_grazing_altitude=100e3,
        )
        self.architecture = Architecture(
            satellites=self.satellites,
            ground_stations=self.stations,
            relays=[self.relay],
            windows={
                "s2": ActiveWindow(start=START + 4 * HOUR),
                "mcmurdo": ActiveWindow(end=START + 6 * HOUR),
            },
        )

    def test_observations(self):
        evaluation = evaluate_architecture(
            Architecture(satellites=self.satellites), self.points, START, END
        )
        expected = collect_observations(self.points, self.satellites, START, END)
        pd.testing.assert_frame_equal(
            _sorted(evaluation.observations), _sorted(expected), check_like=True
        )
        self.assertIsNone(evaluation.metrics)
        for column in ["ground_station", "latency"]:
            self.assertNotIn(column, evaluation.observations.columns)

    def test_windows(self):
        evaluation = self.architecture.evaluate(self.points, START, END, _metrics)
        observations = evaluation.observations
        # a satellite observes over the part of the period in its window
        expected = collect_observations(
            self.points, self.satellites[1], START + 4 * HOUR, END
        )
        pd.testing.assert_frame_equal(
            _sorted(observations[observations["satellite"] == "s2"][expected.columns]),
            _sorted(expected),
            check_like=True,
        )
        # a ground station contacts satellites only in its window
        downlinks = evaluation.downlinks
        mcmurdo = downlinks[downlinks["ground_station"] == "mcmurdo"]
        self.assertGreater(len(mcmurdo), 0)
        self.assertTrue((mcmurdo["end"] <= pd.Timestamp(START + 6 * HOUR)).all())
        self.assertTrue(
            (
                downlinks[downlinks["satellite"] == "s2"]["start"]
                >= pd.Timestamp(START + 4 * HOUR)
            ).all()
        )

    def test_latencies(self):
        evaluation = self.architecture.evaluate(self.points, START, END, _metrics)
        observations = evaluation.observations
        latencies = compute_latencies(
            observations.drop(
                columns=["ground_station", "relay", "downlinked", "latency"]
            ),
            evaluation.downlinks,
        )
        expected = latencies.sort_values(
            ["observed", "target_hash", "satellite"]
        ).reset_index(drop=True)
        actual = observations.sort_values(
            ["epoch", "target_hash", "satellite"]
        ).reset_index(drop=True)
        pd.testing.assert_series_equal(
            actual["latency"], expected["latency"], check_names=False
        )
        pd.testing.assert_series_equal(
            actual["relay"].astype(object),
            expected["relay"].astype(object),
            check_names=False,
        )
        self.assertEqual(set(observations["relay"].dropna()), {"relay"})
        self.assertEqual(evaluation.metrics["latency"], observations["latency"].mean())
        # the downlinks are those of the ground stations and the relay
        direct = collect_downlinks(self.stations[:1], self.satellites[0], START, END)
        self.assertEqual(
            len(
                evaluation.downlinks.query(
                    "satellite == 's1' and ground_station == 'svalbard'"
                )
            ),
            len(direct),
        )

    def test_leave_one_out(self):
        evaluation = self.architecture.evaluate(self.points, START, END, _metrics)
        values = evaluation.leave_one_out()
        self.assertEqual(
            list(values.index), ["s1", "s2", "svalbard", "mcmurdo", "relay"]
        )
        self.assertEqual(
            list(values["kind"]),
            ["satellite", "satellite", "ground_station", "ground_station", "relay"],
        )
        for name in values.index:
            alone = evaluate_architecture(
                self.architecture.without_members(name),
                self.points,
                START,
                END,
                _metrics,
            )
            for metric in ["observations", "targets", "latency"]:
                self.assertEqual(
                    values.loc[name, metric], alone.metrics[metric], (name, metric)
                )
            self.assertEqual(
                values.loc[name, "observations contribution"],
                evaluation.metrics["observations"] - alone.metrics["observations"],
            )
        # the observations are those of the other members' analyses
        without = evaluation.without_members("s1")
        alone = evaluate_architecture(
            self.architecture.without_members("s1"), self.points, START, END
        )
        pd.testing.assert_frame_equal(without.observations, alone.observations)
        # (not numbers: no contribution)
        self.assertTrue(values["label contribution"].isna().all())
        # groups
        grouped = evaluation.leave_one_out(
            {"stations": ["svalbard", "mcmurdo"], "s1": "s1"}
        )
        self.assertEqual(list(grouped["kind"]), ["group", "satellite"])
        self.assertTrue(
            pd.isna(grouped.loc["stations", "latency"])
            or grouped.loc["stations", "latency"] >= timedelta(0)
        )

    def test_add_one(self):
        evaluation = self.architecture.evaluate(self.points, START, END, _metrics)
        station = GroundStation(name="fairbanks", latitude=64.8, longitude=-147.7)
        satellite = _satellite("s3", "06:00")
        values = evaluation.add_one([station, satellite])
        self.assertEqual(list(values.index), ["fairbanks", "s3"])
        for member in [station, satellite]:
            alone = evaluate_architecture(
                self.architecture.with_members(member),
                self.points,
                START,
                END,
                _metrics,
            )
            self.assertEqual(
                values.loc[member.name, "observations"], alone.metrics["observations"]
            )
            self.assertEqual(
                values.loc[member.name, "latency"], alone.metrics["latency"]
            )
        self.assertEqual(values.loc["fairbanks", "observations contribution"], 0)
        self.assertGreater(values.loc["s3", "observations contribution"], 0)
        # with a window
        windowed = evaluation.add_one(
            {"late": (satellite, ActiveWindow(start=END - HOUR))}
        )
        self.assertLess(
            windowed.loc["late", "observations"], values.loc["s3", "observations"]
        )

    def test_without_data_path(self):
        architecture = Architecture(
            satellites=self.satellites, ground_stations=self.stations
        )
        evaluation = architecture.evaluate(self.points, START, END, _metrics)
        # the columns of downlinks remain without ground stations
        without = evaluation.without_members("svalbard", "mcmurdo")
        self.assertTrue(without.observations["latency"].isna().all())
        self.assertTrue(pd.isna(without.metrics["latency"]))

    def test_invalid(self):
        with self.assertRaises(ValueError):
            self.architecture.evaluate(self.points, START, END).leave_one_out()
        with self.assertRaises(TypeError):
            evaluate_architecture(self.satellites, self.points, START, END)
        with self.assertRaises(ValueError):
            self.architecture.evaluate(self.points, START, END, during_contact="later")

    def test_empty(self):
        evaluation = Architecture().evaluate(self.points, START, END, _metrics)
        self.assertEqual(evaluation.metrics["observations"], 0)


def _flow_satellite(name, crossing="10:30", rate=20e6, crosslink=False):
    """A satellite that generates data at a rate and downlinks it at 150 Mbps."""
    transmitters = [Transmitter(band="X", max_data_rate=150e6)]
    if crosslink:
        transmitters.append(
            Transmitter(band="Ka", link="crosslink", max_data_rate=300e6)
        )
    return _satellite(name, crossing).model_copy(
        update={
            "resources": SpacecraftResources(storage=100e9),
            "communications": SpacecraftCommunications(
                transmitters=transmitters, generation=DataGeneration(rate=rate)
            ),
        }
    )


def _flow_metrics(observations):
    return {
        "latency": observations["latency"].mean(),
        "lost": int(observations["lost"].sum()) if "lost" in observations else 0,
    }


class TestDataFlowArchitecture(unittest.TestCase):
    """Tests the evaluation of architectures with a data flow model."""

    def setUp(self):
        self.points = generate_points_fibonacci_lattice(count=60)
        self.station = GroundStation(
            name="svalbard",
            latitude=78.2,
            longitude=15.4,
            min_elevation_angle=5,
            communications=StationCommunications(antennas=1),
        )
        self.satellites = [_flow_satellite("s1"), _flow_satellite("s2", "10:35")]
        self.architecture = Architecture(
            satellites=self.satellites, ground_stations=[self.station]
        )

    def test_default_unchanged(self):
        """
        Tests that communications do not change an evaluation by timing alone.
        """
        plain = Architecture(
            satellites=[_satellite("s1"), _satellite("s2", "10:35")],
            ground_stations=[self.station.model_copy(update={"communications": None})],
        )
        expected = evaluate_architecture(plain, self.points, START, END)
        evaluation = evaluate_architecture(self.architecture, self.points, START, END)
        pd.testing.assert_frame_equal(evaluation.observations, expected.observations)
        self.assertIsNone(evaluation.flow)
        self.assertNotIn("lost", evaluation.observations.columns)

    def test_data_flow(self):
        """
        Tests the downlinks of observations through the flow of data.
        """
        evaluation = evaluate_architecture(
            self.architecture,
            self.points,
            START,
            END,
            _flow_metrics,
            data_flow=DataFlow(),
        )
        observations = evaluation.observations
        self.assertIn("lost", observations.columns)
        self.assertIsNotNone(evaluation.flow)
        self.assertEqual(set(evaluation.flow.summary["satellite"]), {"s1", "s2"})
        delivered = observations["downlinked"].notna()
        self.assertTrue(delivered.any())
        self.assertTrue(
            (observations.loc[delivered, "ground_station"] == "svalbard").all()
        )
        # one antenna serves one satellite at a time
        schedule = evaluation.flow.schedule
        for _, run in schedule.iterrows():
            overlapping = schedule[
                (schedule["start"] < run["end"]) & (schedule["end"] > run["start"])
            ]
            self.assertEqual(len(overlapping), 1)

    def test_leave_one_out(self):
        """
        Tests that leaving out a satellite that competes for an antenna
        reduces the other's loss and latency.
        """
        evaluation = evaluate_architecture(
            self.architecture,
            self.points,
            START,
            END,
            _flow_metrics,
            data_flow=DataFlow(),
        )
        both = evaluation.flow.summary.set_index("satellite")
        alone = evaluation.without_members("s2")
        self.assertIsNotNone(alone.flow)
        one = alone.flow.summary.set_index("satellite")
        self.assertLessEqual(one.loc["s1", "lost"], both.loc["s1", "lost"])
        self.assertGreater(one.loc["s1", "delivered"], both.loc["s1", "delivered"])
        table = evaluation.leave_one_out(["s2"])
        self.assertIn("latency contribution", table.columns)

    def test_add_relay(self):
        """
        Tests that adding a relay with a crosslink reduces latency.
        """
        satellites = [
            _flow_satellite("s1", crosslink=True),
            _flow_satellite("s2", "10:35", crosslink=True),
        ]
        architecture = Architecture(
            satellites=satellites, ground_stations=[self.station]
        )
        relay = Relay(
            name="relay",
            satellites=[
                Satellite(name=f"geo{i}", orbit=GeosynchronousOrbit(longitude=lon))
                for i, lon in enumerate([-41, 174, 62])
            ],
            communications=RelayCommunications(antennas=1, feeder_rate=300e6),
        )
        evaluation = evaluate_architecture(
            architecture,
            self.points,
            START,
            END,
            _flow_metrics,
            data_flow=DataFlow(),
        )
        table = evaluation.add_one([relay])
        self.assertLess(table.loc["relay", "latency"], evaluation.metrics["latency"])
        self.assertLessEqual(table.loc["relay", "lost"], evaluation.metrics["lost"])

    def test_window(self):
        """
        Tests that a satellite generates data only in its window.
        """
        architecture = self.architecture.model_copy(
            update={"windows": {"s2": ActiveWindow(start=START + 6 * HOUR)}}
        )
        evaluation = evaluate_architecture(
            architecture, self.points, START, END, data_flow=DataFlow()
        )
        storage = evaluation.flow.storage
        s2 = storage[storage["satellite"] == "s2"]
        early = s2[s2["time"] <= pd.Timestamp(START + 6 * HOUR)]
        self.assertEqual(early["generated"].max(), 0)
        self.assertAlmostEqual(s2["generated"].iloc[-1], 20e6 * 6 * 3600)

    def test_during_contact_warning(self):
        """
        Tests a warning that during_contact is ignored with a data flow model.
        """
        with self.assertWarns(UserWarning):
            evaluate_architecture(
                self.architecture,
                self.points,
                START,
                START + HOUR,
                data_flow=DataFlow(),
                during_contact="next",
            )


if __name__ == "__main__":
    unittest.main()

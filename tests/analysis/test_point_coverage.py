"""
Unit tests for the point coverage analysis functions in tatc.analysis.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

from datetime import datetime, timedelta, timezone

import geopandas as gpd
import numpy as np
import pandas as pd
from skyfield.api import wgs84
from skyfield.framelib import itrs
from shapely.geometry import Point as ShapelyPoint
from shapely.geometry import LineString, box

from tatc.analysis import (
    aggregate_observations,
    collect_multi_observations,
    collect_observations,
    compute_access_periods,
    grid_observations,
    reduce_observations,
)
from tatc.analysis.observations import _refine_access_periods
from tatc.constants import timescale
from tatc.schemas import (
    ConicalInstrument,
    GeneralPerturbationsOrbit,
    Instrument,
    Point,
    PointedInstrument,
    Satellite,
)
from tatc.utils import compute_cone_and_azimuth, compute_view_tangents

from .common import IssConstellationTestCase


class TestCoverageAnalysis(IssConstellationTestCase):
    """
    Unit tests for the coverage analysis functions in tatc.analysis.
    """

    def setUp(self):
        super().setUp()
        self.point = Point(id=0, latitude=0, longitude=0)

    def test_collect_observations(self):
        """
        Test that observations can be collected for a single satellite and point.
        """
        collect_observations(
            self.point,
            self.satellite,
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 2, tzinfo=timezone.utc),
            instrument_index=0,
            omit_solar=True,
        )

    def test_collect_observations_with_solar(self):
        """
        Test that observations can be collected for a single satellite and point with solar constraints.
        """
        collect_observations(
            self.point,
            self.satellite,
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 2, tzinfo=timezone.utc),
            instrument_index=0,
            omit_solar=False,
        )

    def test_collect_observations_all_culminate(self):
        """
        Test that observations can be collected for a single satellite
        and point when all observations culminate.
        """
        start = datetime(2022, 6, 1, 0, 43, tzinfo=timezone.utc)
        end = datetime(2022, 6, 1, 0, 45, tzinfo=timezone.utc)
        results = collect_observations(
            self.point,
            self.satellite,
            start,
            end,
            instrument_index=0,
        )
        self.assertEqual(len(results), 1)
        self.assertEqual(results.iloc[0].start, start)
        self.assertEqual(results.iloc[0].end, end)

    def test_collect_observations_all_culminate_no_instrument_index(self):
        """
        Test that observations can be collected for a single satellite and p
        oint when all observations culminate and no instrument index is provided.
        """
        start = datetime(2022, 6, 1, 0, 43, tzinfo=timezone.utc)
        end = datetime(2022, 6, 1, 0, 45, tzinfo=timezone.utc)
        results = collect_observations(
            self.point,
            self.satellite,
            start,
            end,
        )
        self.assertEqual(len(results), 1)
        self.assertEqual(results.iloc[0].start, start)
        self.assertEqual(results.iloc[0].end, end)

    def test_collect_observations_miss_first_rise(self):
        """
        Test that observations can be collected for a single satellite and point
        when the first observation is missed due to the start time.
        """
        start = datetime(2022, 6, 1, 0, 43, tzinfo=timezone.utc)
        end = datetime(2022, 6, 1, 1, tzinfo=timezone.utc)
        results = collect_observations(
            self.point,
            self.satellite,
            start,
            end,
            instrument_index=0,
        )
        self.assertEqual(len(results), 1)
        self.assertEqual(results.iloc[0].start, start)

    def test_collect_observations_miss_last_set(self):
        """
        Test that observations can be collected for a single satellite and point
        when the last observation is missed due to the end time.
        """
        start = datetime(2022, 6, 1, tzinfo=timezone.utc)
        end = datetime(2022, 6, 1, 0, 45, tzinfo=timezone.utc)
        results = collect_observations(
            self.point,
            self.satellite,
            start,
            end,
            instrument_index=0,
        )
        self.assertEqual(len(results), 1)
        self.assertEqual(results.iloc[0].end, end)

    def test_collect_observations_narrow_window_fully_visible(self):
        """
        Regression test: a narrow window entirely contained within a longer
        visible pass, with no rise, set, or culmination event inside it
        (elevation angle stays continuously above the threshold and never
        reaches a local maximum in-window), must still be reported as one
        continuous observation spanning the whole window -- not as no
        observation at all.
        """
        start = datetime(2022, 6, 1, 0, 43, 0, tzinfo=timezone.utc)
        end = datetime(2022, 6, 1, 0, 43, 30, tzinfo=timezone.utc)
        results = collect_observations(
            self.point,
            self.satellite,
            start,
            end,
            instrument_index=0,
        )
        self.assertEqual(len(results), 1)
        self.assertEqual(results.iloc[0].start, start)
        self.assertEqual(results.iloc[0].end, end)

    def test_collect_observations_narrow_window_fully_invisible(self):
        """
        Test that a narrow window entirely outside any visible pass (also
        producing no rise/set/culmination events) correctly yields no
        observations, distinguishing this from the fully-visible case in
        `test_collect_observations_narrow_window_fully_visible` (both
        produce zero events, but only one is a true miss).
        """
        start = datetime(2022, 6, 1, 0, 30, 0, tzinfo=timezone.utc)
        end = datetime(2022, 6, 1, 0, 30, 30, tzinfo=timezone.utc)
        results = collect_observations(
            self.point,
            self.satellite,
            start,
            end,
            instrument_index=0,
        )
        self.assertTrue(results.empty)

    def _collect_pointed_observations(self, field_of_regard=100, **kwargs):
        """
        Collects one day of observations of a set of points by a pointed
        instrument with a wide, thin rectangular view (100 deg across track),
        together with those of a nadir instrument with a 100 deg field of
        regard.
        """
        instrument = PointedInstrument(
            name="Pointed",
            field_of_regard=field_of_regard,
            cross_track_field_of_view=100,
            along_track_field_of_view=1,
            is_rectangular=True,
            **kwargs,
        )
        satellite = self.satellite.model_copy(
            update={"instruments": [instrument, Instrument(field_of_regard=100)]}
        )
        start = datetime(2022, 6, 1, tzinfo=timezone.utc)
        points = [
            Point(id=i, latitude=latitude, longitude=longitude)
            for i, (latitude, longitude) in enumerate(
                [(0, 0), (20, 45), (-35, -100), (45, 170)]
            )
        ]
        return satellite, [
            pd.concat(
                [
                    collect_observations(
                        point, satellite, start, start + timedelta(days=1), index
                    )
                    for point in points
                ]
            ).reset_index(drop=True)
            for index in (0, 1)
        ]

    def test_collect_observations_pointed_epoch_at_view_crossing(self):
        """
        Test that a pointed instrument's observation epochs are the times
        its view sweeps over the point, for both velocity frames: the
        point's along-track view angle is zero and it lies in the field of
        view.
        """
        for velocity_frame in ("earth_fixed", "inertial"):
            satellite, (pointed, _) = self._collect_pointed_observations(
                velocity_frame=velocity_frame
            )
            self.assertGreater(len(pointed), 0)
            instrument = satellite.instruments[0]
            for _, observation in pointed.iterrows():
                orbit_track = satellite.orbit.to_gp_orbit().get_orbit_track(
                    observation.epoch
                )
                target = wgs84.latlon(observation.geometry.y, observation.geometry.x)
                along, _ = compute_view_tangents(orbit_track, target, velocity_frame)
                self.assertAlmostEqual(float(along), 0, delta=1e-5)
                self.assertTrue(
                    instrument.is_in_field_of_view(orbit_track, target).all()
                )
                self.assertTrue(
                    observation.start <= observation.epoch <= observation.end
                )

    def test_collect_observations_pointed_geocentric_nadir(self):
        """
        Test that, with a geocentric nadir reference, a pointed instrument's
        observation epochs are the times its geocentric view sweeps over the
        point, which differ from those of a geodetic view by a fraction of a
        second away from the equator.
        """
        satellite, (geocentric, _) = self._collect_pointed_observations(
            nadir_reference="geocentric"
        )
        _, (geodetic, _) = self._collect_pointed_observations()
        # the field of regard is also measured from the nadir reference, so
        # marginal passes can differ
        self.assertLessEqual(abs(len(geocentric) - len(geodetic)), 1)
        for _, observation in geocentric.iterrows():
            orbit_track = satellite.orbit.to_gp_orbit().get_orbit_track(
                observation.epoch
            )
            target = wgs84.latlon(observation.geometry.y, observation.geometry.x)
            along, _ = compute_view_tangents(
                orbit_track, target, nadir_reference="geocentric"
            )
            self.assertAlmostEqual(float(along), 0, delta=1e-5)
        matched = pd.merge_asof(
            geocentric.sort_values("epoch"),
            geodetic[["point_id", "epoch"]]
            .rename(columns={"epoch": "epoch_geodetic"})
            .sort_values("epoch_geodetic"),
            left_on="epoch",
            right_on="epoch_geodetic",
            by="point_id",
            direction="nearest",
            tolerance=pd.Timedelta(seconds=10),
        ).dropna(subset=["epoch_geodetic"])
        self.assertGreaterEqual(len(matched), min(len(geocentric), len(geodetic)) - 1)
        difference = np.abs(
            (matched.epoch - matched.epoch_geodetic).dt.total_seconds().to_numpy()
        )
        self.assertLess(difference.max(), 1)
        self.assertGreater(
            difference[matched.geometry.y.abs().to_numpy() > 30].max(), 0.05
        )

    def test_collect_observations_access_period_at_field_of_regard(self):
        """
        Test that the access period of a nadir instrument starts and ends
        when the point's angle from nadir equals half the field of regard,
        for points at low and high latitudes and both nadir references.
        """
        start = datetime(2022, 6, 1, tzinfo=timezone.utc)
        for reference in ("geodetic", "geocentric"):
            instrument = Instrument(
                name="Nadir", field_of_regard=100, nadir_reference=reference
            )
            satellite = self.satellite.model_copy(update={"instruments": [instrument]})
            for latitude in (0, 50):
                point = Point(id=0, latitude=latitude, longitude=30)
                observations = collect_observations(
                    point, satellite, start, start + timedelta(days=1)
                )
                self.assertGreater(len(observations), 0)
                target = wgs84.latlon(latitude, 30)
                for column in ("start", "end"):
                    times = observations[column].tolist()
                    angle, _ = compute_cone_and_azimuth(
                        self.satellite.orbit.to_gp_orbit().get_orbit_track(times),
                        target,
                        nadir_reference=reference,
                    )
                    np.testing.assert_allclose(angle, 50, atol=1e-3)

    def test_get_orbit_track_repeated(self):
        """
        Test that the orbit track has the Earth-fixed position and velocity
        of the orbit's element maintained on its repeat ground track at the
        times shifted by whole repeat cycles (here, two cycles after the
        epoch, and one cycle before it), expressed at the unshifted times;
        and that, without a repeat cycle, it is directly propagated.
        """
        orbit = GeneralPerturbationsOrbit.from_tle(
            [
                "1 39084U 13008A   26213.27824675  .00000294  00000+0  75333-4 0  9990",
                "2 39084  98.2277 282.8718 0001275  92.4910 267.6434 14.57104473704466",
            ],
            remove_drag=True,
            repeat_cycle="auto",
        )
        repeat_cycle = orbit.get_repeat_cycle()
        epoch = orbit.get_epoch()
        times = [epoch + timedelta(days=40), epoch - timedelta(days=20)]
        shifts = [2 * repeat_cycle, -repeat_cycle]
        repeated = orbit.get_orbit_track(times)
        maintained = orbit.get_repeat_element().to_skyfield()
        direct = maintained.at(
            timescale.from_datetimes([t - d for t, d in zip(times, shifts)])
        )
        self.assertEqual(repeated.t.utc_datetime().tolist(), times)
        repeated_position, repeated_velocity = repeated.frame_xyz_and_velocity(itrs)
        direct_position, direct_velocity = direct.frame_xyz_and_velocity(itrs)
        np.testing.assert_allclose(repeated_position.m, direct_position.m, atol=1e-3)
        np.testing.assert_allclose(
            repeated_velocity.m_per_s, direct_velocity.m_per_s, atol=1e-6
        )
        np.testing.assert_allclose(
            orbit.model_copy(update={"repeat_cycle": None})
            .get_orbit_track(times)
            .position.m,
            orbit.elements[0]
            .to_skyfield(remove_drag=True)
            .at(timescale.from_datetimes(times))
            .position.m,
        )

    def test_collect_observations_continuous_at_repeat_boundary(self):
        """
        Test that, when observation events are repeated, observations just
        before the end of a repeat cycle are at the times of the element
        maintained on its repeat ground track, like those just after it,
        rather than offset by the drift of the element's own mean motion
        over the cycle (about 10 s for this Landsat 8 element set).
        """
        orbit = GeneralPerturbationsOrbit.from_tle(
            [
                "1 39084U 13008A   26213.27824675  .00000294  00000+0  75333-4 0  9990",
                "2 39084  98.2277 282.8718 0001275  92.4910 267.6434 14.57104473704466",
            ],
            remove_drag=True,
            repeat_cycle="auto",
        )
        boundary = orbit.get_epoch() + orbit.get_repeat_cycle()
        satellite = Satellite(
            name="Landsat 8",
            orbit=orbit,
            instruments=[Instrument(name="nadir", field_of_regard=15)],
        )
        for minutes in (-3, 3):
            time = boundary + timedelta(minutes=minutes)
            sub = wgs84.subpoint_of(orbit.get_orbit_track(time))
            point = Point(
                id=0,
                latitude=float(sub.latitude.degrees),
                longitude=float(sub.longitude.degrees),
            )
            observations = collect_observations(
                point,
                satellite,
                boundary - timedelta(hours=3),
                boundary + timedelta(hours=3),
            )
            nearest = (observations.epoch - time).abs().min()
            self.assertLess(nearest, pd.Timedelta(seconds=1))

    def test_collect_observations_repeat_cycle_maintained_orbit(self):
        """
        Test that, when observation events are repeated with the orbit's
        repeat cycle, the observations of each repeated cycle are those of
        the first cycle shifted by whole repeat cycles (a maintained orbit),
        for a pointed instrument whose epochs are solved within the access
        periods.
        """
        orbit = GeneralPerturbationsOrbit.from_tle(
            [
                "1 39084U 13008A   26213.27824675  .00000294  00000+0  75333-4 0  9990",
                "2 39084  98.2277 282.8718 0001275  92.4910 267.6434 14.57104473704466",
            ],
            remove_drag=True,
            repeat_cycle="auto",
        )
        repeat_cycle = orbit.get_repeat_cycle()
        self.assertIsNotNone(repeat_cycle)
        satellite = Satellite(
            name="Landsat 8",
            orbit=orbit,
            instruments=[
                PointedInstrument(
                    name="Imager",
                    field_of_regard=30,
                    cross_track_field_of_view=15,
                    along_track_field_of_view=0.1,
                    is_rectangular=True,
                )
            ],
        )
        start = orbit.get_epoch()
        end = start + 4 * repeat_cycle
        point = Point(id=0, latitude=40, longitude=-105)
        observations = collect_observations(point, satellite, start, end)
        cycle = ((observations.epoch - start) // repeat_cycle).to_numpy()
        first = observations.epoch[cycle == 0].reset_index(drop=True)
        self.assertGreater(len(first), 0)
        for k in range(1, 4):
            repeated = observations.epoch[cycle == k].reset_index(drop=True)
            self.assertEqual(len(repeated), len(first))
            np.testing.assert_allclose(
                (repeated - first - k * repeat_cycle).dt.total_seconds(), 0, atol=1e-2
            )

    def test_collect_observations_repeat_cycle_anchored_at_epoch(self):
        """
        Test that observations ten repeat cycles after the orbit's epoch are
        those of the first two repeat cycles after the epoch, shifted by ten
        repeat cycles: the repeated cycle is anchored at the epoch, rather
        than propagated (with drag) to the start of the analysis period.
        """
        orbit = GeneralPerturbationsOrbit.from_tle(
            [
                "1 39084U 13008A   26213.27824675  .00000294  00000+0  75333-4 0  9990",
                "2 39084  98.2277 282.8718 0001275  92.4910 267.6434 14.57104473704466",
            ],
            remove_drag=True,
            repeat_cycle="auto",
        )
        repeat_cycle = orbit.get_repeat_cycle()
        self.assertIsNotNone(repeat_cycle)
        satellite = Satellite(
            name="Landsat 8",
            orbit=orbit,
            instruments=[
                PointedInstrument(
                    name="Imager",
                    field_of_regard=30,
                    cross_track_field_of_view=15,
                    along_track_field_of_view=0.1,
                    is_rectangular=True,
                )
            ],
        )
        epoch = orbit.get_epoch()
        point = Point(id=0, latitude=40, longitude=-105)
        near = collect_observations(point, satellite, epoch, epoch + 2 * repeat_cycle)
        far = collect_observations(
            point, satellite, epoch + 10 * repeat_cycle, epoch + 12 * repeat_cycle
        )
        self.assertGreater(len(near), 0)
        self.assertEqual(len(far), len(near))
        np.testing.assert_allclose(
            (far.epoch - near.epoch - 10 * repeat_cycle).dt.total_seconds(),
            0,
            atol=1e-2,
        )

    def test_collect_observations_pointed_matches_field_of_regard(self):
        """
        Test that a wide, thin pointed view observes the same passes as a
        nadir instrument with (about) the same field of regard, within
        seconds of the time of closest approach.
        """
        _, (pointed, nadir) = self._collect_pointed_observations()
        self.assertEqual(len(pointed), len(nadir))
        np.testing.assert_array_less(
            np.abs((pointed.epoch - nadir.epoch).dt.total_seconds()), 5
        )

    def test_collect_observations_pointed_pitch(self):
        """
        Test that a forward (aft) pitched view observes a point before
        (after) the time of closest approach, by about the time to travel
        420 km * tan(20 deg) = 153 km (22 s). The pitched view reaches
        slightly farther from nadir across track than the nadir instrument's
        field of regard, so it may observe additional passes.
        """
        for pitch_angle, sign in [(20, -1), (-20, 1)]:
            _, (pointed, nadir) = self._collect_pointed_observations(
                field_of_regard=140, pitch_angle=pitch_angle
            )
            matched = pd.merge_asof(
                nadir.sort_values("epoch"),
                pointed[["point_id", "epoch"]]
                .rename(columns={"epoch": "epoch_pointed"})
                .sort_values("epoch_pointed"),
                left_on="epoch",
                right_on="epoch_pointed",
                by="point_id",
                direction="nearest",
                tolerance=pd.Timedelta(minutes=1),
            )
            self.assertFalse(matched.epoch_pointed.isna().any())
            self.assertLessEqual(len(pointed) - len(nadir), 1)
            np.testing.assert_allclose(
                sign * (matched.epoch_pointed - matched.epoch).dt.total_seconds(),
                22,
                atol=2,
            )

    def test_collect_observations_pointed_pitch_scan_matches_frame(self):
        """
        Test that a pitched view in scan geometry, whose pitch angle tilts the
        scan plane, observes points at the same times as a pitched view in
        frame geometry: both views sweep the plane containing the cross-track
        axis and the tilted boresight.
        """
        _, (frame, _) = self._collect_pointed_observations(
            field_of_regard=140, pitch_angle=20
        )
        _, (scan, _) = self._collect_pointed_observations(
            field_of_regard=140, pitch_angle=20, view_geometry="scan"
        )
        self.assertEqual(len(frame), len(scan))
        np.testing.assert_allclose(
            (scan.epoch - frame.epoch).dt.total_seconds(), 0, atol=0.01
        )

    def test_collect_observations_pointed_pitch_profile(self):
        """
        Test that a view pitched forward over the northern half of the orbit
        and aft over the southern half (pitch_angle_profile) observes points
        in the northern hemisphere before the time of closest approach and
        points in the southern hemisphere after it, by about 22 s.
        """
        _, (pointed, nadir) = self._collect_pointed_observations(
            field_of_regard=140,
            pitch_angle_profile=[(5, 20), (175, 20), (185, -20), (355, -20)],
        )
        matched = pd.merge_asof(
            nadir.sort_values("epoch"),
            pointed[["point_id", "epoch"]]
            .rename(columns={"epoch": "epoch_pointed"})
            .sort_values("epoch_pointed"),
            left_on="epoch",
            right_on="epoch_pointed",
            by="point_id",
            direction="nearest",
            tolerance=pd.Timedelta(minutes=1),
        )
        # points away from the equator, where the pitch changes
        matched = matched[matched.geometry.y.abs() > 10]
        self.assertGreater(len(matched), 0)
        self.assertFalse(matched.epoch_pointed.isna().any())
        sign = np.where(matched.geometry.y > 0, -1, 1)
        np.testing.assert_allclose(
            sign * (matched.epoch_pointed - matched.epoch).dt.total_seconds(),
            22,
            atol=2,
        )

    def test_collect_observations_pointed_forward_wide_field_of_regard(self):
        """
        Test that a strongly pitched view with a field of regard reaching the
        horizon (so the target is behind the view at the ends of each access
        period) observes points at its view crossings, at least on every pass
        observed by a narrower nadir instrument.
        """
        satellite, (pointed, nadir) = self._collect_pointed_observations(
            field_of_regard=180, pitch_angle=45
        )
        self.assertGreaterEqual(len(pointed), len(nadir))
        instrument = satellite.instruments[0]
        for _, observation in pointed.iterrows():
            orbit_track = satellite.orbit.to_gp_orbit().get_orbit_track(
                observation.epoch
            )
            target = wgs84.latlon(observation.geometry.y, observation.geometry.x)
            along, _ = compute_view_tangents(
                orbit_track, target, instrument.velocity_frame, 0, 45
            )
            self.assertAlmostEqual(float(along), 0, delta=1e-5)
            self.assertTrue(instrument.is_in_field_of_view(orbit_track, target).all())

    def test_collect_observations_pointed_access_time_fixed(self):
        """
        Test that a fixed access time is centered on a pointed instrument's
        view crossing time.
        """
        _, (pointed, _) = self._collect_pointed_observations(
            min_access_time=timedelta(seconds=10), access_time_fixed=True
        )
        self.assertGreater(len(pointed), 0)
        self.assertTrue(
            (pointed.epoch - pointed.start == timedelta(seconds=5)).all()
            and (pointed.end - pointed.epoch == timedelta(seconds=5)).all()
        )

    def _collect_conical_observations(self, **kwargs):
        """
        Collects one day of observations of a set of points by a conical
        instrument with a 45 deg cone, together with those of a nadir
        instrument with a 90 deg field of regard.
        """
        instrument = ConicalInstrument(
            name="Conical", cone_angle=45, along_track_field_of_view=1, **kwargs
        )
        satellite = self.satellite.model_copy(
            update={"instruments": [instrument, Instrument(field_of_regard=90)]}
        )
        start = datetime(2022, 6, 1, tzinfo=timezone.utc)
        points = [
            Point(id=i, latitude=latitude, longitude=longitude)
            for i, (latitude, longitude) in enumerate(
                [(0, 0), (20, 45), (-35, -100), (45, 170)]
            )
        ]
        return satellite, [
            pd.concat(
                [
                    collect_observations(
                        point, satellite, start, start + timedelta(days=1), index
                    )
                    for point in points
                ]
            ).reset_index(drop=True)
            for index in (0, 1)
        ]

    def test_collect_observations_conical_epoch_on_cone(self):
        """
        Test that a conical instrument's observation epochs are times when
        the point lies on the cone within the scan sector.
        """
        for kwargs in [
            {},
            {"scan_half_width": 60},
            {"scan_center_azimuth": 180, "scan_half_width": 60},
        ]:
            satellite, (conical, _) = self._collect_conical_observations(**kwargs)
            self.assertGreater(len(conical), 0)
            instrument = satellite.instruments[0]
            for _, observation in conical.iterrows():
                orbit_track = satellite.orbit.to_gp_orbit().get_orbit_track(
                    observation.epoch
                )
                target = wgs84.latlon(observation.geometry.y, observation.geometry.x)
                cone, azimuth = compute_cone_and_azimuth(orbit_track, target)
                self.assertAlmostEqual(float(cone), 45, delta=1e-4)
                offset = (
                    float(azimuth) - instrument.scan_center_azimuth + 180
                ) % 360 - 180
                self.assertLessEqual(abs(offset), instrument.scan_half_width)

    def test_collect_observations_conical_fore_and_aft(self):
        """
        Test that a full rotation observes points twice per pass (entering
        and leaving the cone), a forward sector at the start of the nadir
        instrument's access period, and an aft sector at its end.
        """
        _, (full, nadir) = self._collect_conical_observations()
        self.assertEqual(len(full), 2 * len(nadir))
        _, (forward, nadir) = self._collect_conical_observations(scan_half_width=89)
        self.assertGreater(len(forward), 0)
        merged = pd.merge_asof(
            forward.sort_values("epoch"),
            nadir[["point_id", "start", "end"]].sort_values("start"),
            left_on="epoch",
            right_on="start",
            by="point_id",
            direction="nearest",
        )
        np.testing.assert_array_less(
            np.abs((merged.epoch - merged.start_y).dt.total_seconds()), 10
        )
        _, (aft, nadir) = self._collect_conical_observations(
            scan_center_azimuth=180, scan_half_width=89
        )
        merged = pd.merge_asof(
            aft.sort_values("epoch"),
            nadir[["point_id", "start", "end"]].sort_values("end"),
            left_on="epoch",
            right_on="end",
            by="point_id",
            direction="nearest",
        )
        np.testing.assert_array_less(
            np.abs((merged.epoch - merged.end_y).dt.total_seconds()), 10
        )

    def test_collect_observations_null(self):
        """
        Test that no observations are collected for a single satellite and point
        when the time window does not overlap with any observations.
        """
        start = datetime(2022, 6, 1, tzinfo=timezone.utc)
        end = datetime(2022, 6, 1, 0, 30, tzinfo=timezone.utc)
        results = collect_observations(
            self.point,
            self.satellite,
            start,
            end,
            instrument_index=0,
        )
        self.assertTrue(results.empty)

    def test_collect_observations_null_with_solar(self):
        """
        Test that no observations are collected for a single satellite and point
        when the time window does not overlap with any observations, even with solar constraints.
        """
        start = datetime(2022, 6, 1, tzinfo=timezone.utc)
        end = datetime(2022, 6, 1, 0, 30, tzinfo=timezone.utc)
        results = collect_observations(
            self.point, self.satellite, start, end, instrument_index=0, omit_solar=False
        )
        self.assertTrue(results.empty)

    def test_collect_multi_observations(self):
        """
        Test that observations can be collected for multiple satellites and a single point.
        """
        collect_multi_observations(
            self.point,
            self.constellation.generate_members(),
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 2, tzinfo=timezone.utc),
        )

    def test_collect_multi_observations_null(self):
        """
        Test that no observations are collected for multiple satellites and a single point
        when the time window does not overlap with any observations.
        """
        start = datetime(2022, 6, 1, 0, 10, tzinfo=timezone.utc)
        end = datetime(2022, 6, 1, 0, 30, tzinfo=timezone.utc)
        results = collect_multi_observations(
            self.point,
            self.constellation.generate_members(),
            start,
            end,
        )
        self.assertTrue(results.empty)

    def test_collect_multi_observations_no_satellites(self):
        """
        Regression test: an empty `satellites` list must return an empty
        DataFrame with the expected columns, not raise (there is nothing to
        concatenate when no satellite/instrument pair ever runs).
        """
        results = collect_multi_observations(
            self.point,
            [],
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 2, tzinfo=timezone.utc),
        )
        self.assertTrue(results.empty)
        self.assertIn("start", results.columns)

    @staticmethod
    def _make_observation(point_id, satellite, instrument, start, end):
        """
        Build a synthetic observation record matching the schema
        `collect_observations` produces, for direct, deterministic control
        over `aggregate_observations` inputs (bypassing orbital propagation).
        """
        return {
            "point_id": point_id,
            "geometry": ShapelyPoint(0, 0),
            "satellite": satellite,
            "instrument": instrument,
            "start": pd.Timestamp(start),
            "epoch": pd.Timestamp(start)
            + (pd.Timestamp(end) - pd.Timestamp(start)) / 2,
            "end": pd.Timestamp(end),
        }

    def test_aggregate_observations_merges_overlapping_and_nested_intervals(self):
        """
        Test the core interval-merging algorithm directly with synthetic,
        deliberately overlapping/nested/gapped observations: satellite B's
        window is fully nested inside A's, and C's window starts before A
        ends but after B ends (a case a naive "compare only to the previous
        row" check would wrongly split, since C.start > B.end even though
        C still overlaps A's still-ongoing window -- the running max via
        `.cummax()` is what keeps this correct). D is fully separate.
        """
        t0 = datetime(2022, 6, 1, tzinfo=timezone.utc)
        observations = gpd.GeoDataFrame(
            [
                self._make_observation(0, "A", "Test", t0, t0 + timedelta(minutes=10)),
                self._make_observation(
                    0, "B", "Test", t0 + timedelta(minutes=2), t0 + timedelta(minutes=4)
                ),
                self._make_observation(
                    0,
                    "C",
                    "Test",
                    t0 + timedelta(minutes=9),
                    t0 + timedelta(minutes=15),
                ),
                self._make_observation(
                    0,
                    "D",
                    "Test",
                    t0 + timedelta(minutes=20),
                    t0 + timedelta(minutes=25),
                ),
            ],
            crs="EPSG:4326",
        )
        results = aggregate_observations(observations)
        self.assertEqual(len(results.index), 2)
        self.assertEqual(results.iloc[0].satellite, "A, B, C")
        self.assertEqual(results.iloc[0].start, t0)
        self.assertEqual(results.iloc[0].end, t0 + timedelta(minutes=15))
        self.assertTrue(pd.isna(results.iloc[0].revisit))
        self.assertEqual(results.iloc[1].satellite, "D")
        self.assertEqual(results.iloc[1].start, t0 + timedelta(minutes=20))
        self.assertEqual(results.iloc[1].end, t0 + timedelta(minutes=25))
        self.assertEqual(results.iloc[1].revisit, timedelta(minutes=5))

    def test_aggregate_observations_isolates_point_ids(self):
        """
        Test that merging and revisit computation are scoped per point_id:
        a point_id=1 observation must not be merged with, or treated as a
        revisit predecessor for, a point_id=0 observation, even if their
        windows would otherwise overlap/abut.
        """
        t0 = datetime(2022, 6, 1, tzinfo=timezone.utc)
        observations = gpd.GeoDataFrame(
            [
                self._make_observation(0, "A", "Test", t0, t0 + timedelta(minutes=10)),
                self._make_observation(
                    1,
                    "B",
                    "Test",
                    t0 + timedelta(minutes=5),
                    t0 + timedelta(minutes=15),
                ),
            ],
            crs="EPSG:4326",
        )
        results = aggregate_observations(observations)
        self.assertEqual(len(results.index), 2)
        for i in range(len(results.index)):
            self.assertTrue(pd.isna(results.iloc[i].revisit))

    def test_aggregate_observations_epoch_is_merged_midpoint(self):
        """
        Test that the merged group's epoch is the midpoint of its (merged)
        start/end, not the mean of the constituent observations' own
        (pre-merge) epochs -- these differ here since B's epoch sits much
        earlier than the midpoint of the full A+B merged window.
        """
        t0 = datetime(2022, 6, 1, tzinfo=timezone.utc)
        observations = gpd.GeoDataFrame(
            [
                self._make_observation(0, "A", "Test", t0, t0 + timedelta(minutes=10)),
                self._make_observation(
                    0, "B", "Test", t0 + timedelta(minutes=1), t0 + timedelta(minutes=2)
                ),
            ],
            crs="EPSG:4326",
        )
        results = aggregate_observations(observations)
        self.assertEqual(len(results.index), 1)
        self.assertEqual(results.iloc[0].epoch, t0 + timedelta(minutes=5))

    def test_aggregate_observations_drops_per_observation_columns(self):
        """
        Test that per-observation columns that lose their meaning once
        merged across satellites (sat_alt, sat_az, sat_sunlit, solar_alt,
        solar_az, solar_time) are dropped, even if present on the input.
        """
        t0 = datetime(2022, 6, 1, tzinfo=timezone.utc)
        record = self._make_observation(0, "A", "Test", t0, t0 + timedelta(minutes=10))
        record["sat_alt"] = 45.0
        record["sat_az"] = 180.0
        record["solar_alt"] = 10.0
        observations = gpd.GeoDataFrame([record], crs="EPSG:4326")
        results = aggregate_observations(observations)
        for column in ("sat_alt", "sat_az", "solar_alt"):
            self.assertNotIn(column, results.columns)

    def test_aggregate_observations(self):
        """
        Test that observations can be aggregated for multiple satellites and a single point.
        """
        results = collect_multi_observations(
            self.point,
            self.constellation.generate_members(),
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 2, tzinfo=timezone.utc),
        )
        results = aggregate_observations(results)
        for i in range(len(results.index)):
            self.assertEqual(
                results.iloc[i].end - results.iloc[i].start, results.iloc[i].access
            )
        for i in range(1, len(results.index)):
            self.assertEqual(
                results.iloc[i].start - results.iloc[i - 1].end, results.iloc[i].revisit
            )

    def test_aggregate_observations_null(self):
        """
        Test that no observations are aggregated for multiple satellites and a single point
        when the time window does not overlap with any observations.
        """
        results = collect_multi_observations(
            self.point,
            self.constellation.generate_members(),
            datetime(2022, 6, 1, 0, 10, tzinfo=timezone.utc),
            datetime(2022, 6, 1, 0, 30, tzinfo=timezone.utc),
        )
        results = aggregate_observations(results)
        self.assertTrue(results.empty)

    @staticmethod
    def _make_aggregated_observation(point_id, access_minutes, revisit_minutes):
        """
        Build a synthetic aggregated-observation record (matching
        `aggregate_observations`'s output schema) for direct, deterministic
        control over `reduce_observations` inputs. `revisit_minutes=None`
        produces `pandas.NaT`, matching the first observation for a point.
        """
        return {
            "point_id": point_id,
            "geometry": ShapelyPoint(0, 0),
            "access": pd.Timedelta(minutes=access_minutes),
            "revisit": (
                pd.NaT
                if revisit_minutes is None
                else pd.Timedelta(minutes=revisit_minutes)
            ),
        }

    def test_reduce_observations_computes_mean_access_and_skips_first_revisit(self):
        """
        Test that access is averaged over every sample, but revisit is
        averaged only over the samples with a defined revisit -- the first
        sample's revisit is undefined (NaT, no prior observation), and must
        be skipped rather than treated as a zero-minute revisit (which
        would otherwise skew the mean down substantially).
        """
        observations = gpd.GeoDataFrame(
            [
                self._make_aggregated_observation(0, 5, None),
                self._make_aggregated_observation(0, 10, 55),
                self._make_aggregated_observation(0, 15, 50),
            ],
            crs="EPSG:4326",
        )
        results = reduce_observations(observations)
        self.assertEqual(len(results.index), 1)
        self.assertEqual(results.iloc[0].samples, 3)
        self.assertEqual(results.iloc[0].access, timedelta(minutes=10))
        self.assertEqual(results.iloc[0].revisit, timedelta(minutes=52.5))

    def test_reduce_observations_isolates_point_ids(self):
        """
        Test that statistics are computed independently per point_id, not
        pooled across points.
        """
        observations = gpd.GeoDataFrame(
            [
                self._make_aggregated_observation(0, 5, None),
                self._make_aggregated_observation(1, 5, None),
                self._make_aggregated_observation(1, 7, 30),
            ],
            crs="EPSG:4326",
        )
        results = reduce_observations(observations)
        self.assertEqual(len(results.index), 2)
        point_0 = results[results.point_id == 0].iloc[0]
        point_1 = results[results.point_id == 1].iloc[0]
        self.assertEqual(point_0.samples, 1)
        self.assertEqual(point_0.access, timedelta(minutes=5))
        self.assertTrue(pd.isna(point_0.revisit))
        self.assertEqual(point_1.samples, 2)
        self.assertEqual(point_1.access, timedelta(minutes=6))
        self.assertEqual(point_1.revisit, timedelta(minutes=30))

    def test_reduce_observations(self):
        """
        Test that observations can be reduced for multiple satellites and a single point.
        """
        results = collect_multi_observations(
            self.point,
            self.constellation.generate_members(),
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 10, tzinfo=timezone.utc),
        )
        aggregated_results = aggregate_observations(results)
        reduced_results = reduce_observations(aggregated_results)
        self.assertEqual(len(reduced_results.index), 1)
        self.assertAlmostEqual(
            reduced_results.iloc[0].access,
            aggregated_results[
                aggregated_results.point_id == reduced_results.iloc[0].point_id
            ].access.mean(),
            delta=timedelta(seconds=0.01),
        )
        self.assertAlmostEqual(
            reduced_results.iloc[0].revisit,
            aggregated_results[
                aggregated_results.point_id == reduced_results.iloc[0].point_id
            ].revisit.mean(),
            delta=timedelta(seconds=0.01),
        )

    def test_reduce_observations_null(self):
        """
        Test that no observations are reduced for multiple satellites and a single point
        when the time window does not overlap with any observations.
        """
        results = collect_multi_observations(
            self.point,
            self.constellation.generate_members(),
            datetime(2022, 6, 1, 0, 10, tzinfo=timezone.utc),
            datetime(2022, 6, 1, 0, 30, tzinfo=timezone.utc),
        )
        aggregated_results = aggregate_observations(results)
        reduced_results = reduce_observations(aggregated_results)
        self.assertTrue(reduced_results.empty)

    @staticmethod
    def _make_cell(cell_id, min_lon, min_lat, max_lon, max_lat):
        """
        Build a synthetic cell record matching the schema expected by
        `grid_observations` (a `cell_id` plus a polygon `geometry`).
        """
        return {"cell_id": cell_id, "geometry": box(min_lon, min_lat, max_lon, max_lat)}

    @staticmethod
    def _make_reduced_observation(
        point_id, lon, lat, access_seconds, revisit_seconds, samples
    ):
        """
        Build a synthetic reduced-observation record (matching
        `reduce_observations`'s output schema) for direct, deterministic
        control over `grid_observations` inputs.
        """
        return {
            "point_id": point_id,
            "geometry": ShapelyPoint(lon, lat),
            "access": pd.Timedelta(seconds=access_seconds),
            "revisit": pd.Timedelta(seconds=revisit_seconds),
            "samples": samples,
        }

    def test_grid_observations_empty_reduced_observations(self):
        """
        Test that an empty `reduced_observations` yields every cell with
        zero samples and undefined access/revisit, rather than an empty
        result.
        """
        cells = gpd.GeoDataFrame(
            [self._make_cell(0, 0, 0, 1, 1), self._make_cell(1, 2, 2, 3, 3)],
            crs="EPSG:4326",
        )
        reduced = gpd.GeoDataFrame(
            columns=["point_id", "geometry", "access", "revisit", "samples"],
            crs="EPSG:4326",
        )
        result = grid_observations(reduced, cells)
        self.assertEqual(len(result.index), 2)
        self.assertTrue((result.samples == 0).all())
        self.assertTrue(result.access.isna().all())
        self.assertTrue(result.revisit.isna().all())

    def test_grid_observations_single_point_passes_through_unchanged(self):
        """
        Test that a single point within a cell yields that point's own
        access/revisit/samples unchanged (the trivial case of a weighted
        mean over one value).
        """
        cells = gpd.GeoDataFrame([self._make_cell(0, 0, 0, 1, 1)], crs="EPSG:4326")
        reduced = gpd.GeoDataFrame(
            [self._make_reduced_observation(0, 0.5, 0.5, 5, 10, 100)],
            crs="EPSG:4326",
        )
        result = grid_observations(reduced, cells)
        self.assertEqual(len(result.index), 1)
        self.assertEqual(result.iloc[0].samples, 100)
        self.assertEqual(result.iloc[0].access, timedelta(seconds=5))
        self.assertEqual(result.iloc[0].revisit, timedelta(seconds=10))

    def test_grid_observations_uses_weighted_arithmetic_mean_for_access(self):
        """
        Test that access is combined across points in the same cell using
        a sample-weighted arithmetic mean.
        """
        cells = gpd.GeoDataFrame([self._make_cell(0, 0, 0, 1, 1)], crs="EPSG:4326")
        reduced = gpd.GeoDataFrame(
            [
                self._make_reduced_observation(0, 0.4, 0.4, 5, 10, 100),
                self._make_reduced_observation(1, 0.6, 0.6, 50, 100, 10),
            ],
            crs="EPSG:4326",
        )
        result = grid_observations(reduced, cells)
        self.assertEqual(len(result.index), 1)
        self.assertEqual(result.iloc[0].samples, 110)
        expected_access = (5 * 100 + 50 * 10) / 110
        self.assertAlmostEqual(
            result.iloc[0].access.total_seconds(), expected_access, places=6
        )

    def test_grid_observations_uses_weighted_harmonic_mean_for_revisit(self):
        """
        Test that revisit is combined across points in the same cell using
        a sample-weighted harmonic mean, not an arithmetic mean: revisit is
        a time-between-events (reciprocal-of-rate) quantity, so a naive
        arithmetic mean would under-weight the more frequently sampled
        point. The two means give clearly different results here (~10.9s
        harmonic vs. 55s arithmetic), so this distinguishes them concretely.
        """
        cells = gpd.GeoDataFrame([self._make_cell(0, 0, 0, 1, 1)], crs="EPSG:4326")
        reduced = gpd.GeoDataFrame(
            [
                self._make_reduced_observation(0, 0.4, 0.4, 5, 10, 100),
                self._make_reduced_observation(1, 0.6, 0.6, 50, 100, 10),
            ],
            crs="EPSG:4326",
        )
        result = grid_observations(reduced, cells)
        expected_harmonic_revisit = (100 + 10) / (100 / 10 + 10 / 100)
        naive_arithmetic_revisit = (10 * 100 + 100 * 10) / 110
        self.assertAlmostEqual(
            result.iloc[0].revisit.total_seconds(), expected_harmonic_revisit, places=6
        )
        self.assertNotAlmostEqual(
            result.iloc[0].revisit.total_seconds(), naive_arithmetic_revisit, places=1
        )

    def test_grid_observations_revisit_is_invariant_to_point_density(self):
        """
        Regression test: a cell's gridded revisit must not depend on how
        many (near-identical) points happen to fall inside it -- it should
        represent a typical point's revisit, not shrink just because the
        input point grid happened to be sampled more finely there. This
        specifically distinguishes the (correct) sample-weighted harmonic
        mean from summing raw per-point rates (1/revisit) unweighted by
        sample count, which would make revisit shrink roughly in
        proportion to the number of pooled points instead of staying
        constant.
        """
        cells = gpd.GeoDataFrame([self._make_cell(0, 0, 0, 1, 1)], crs="EPSG:4326")
        sparse = gpd.GeoDataFrame(
            [self._make_reduced_observation(0, 0.5, 0.5, 5, 100, 20)],
            crs="EPSG:4326",
        )
        dense = gpd.GeoDataFrame(
            [
                self._make_reduced_observation(i, 0.1 * i, 0.1 * i, 5, 100, 20)
                for i in range(1, 10)
            ],
            crs="EPSG:4326",
        )
        sparse_result = grid_observations(sparse, cells)
        dense_result = grid_observations(dense, cells)
        self.assertAlmostEqual(
            sparse_result.iloc[0].revisit.total_seconds(),
            dense_result.iloc[0].revisit.total_seconds(),
            places=6,
        )

    def test_grid_observations_cell_without_points_is_omitted(self):
        """
        Documents current behavior: unlike the fully-empty
        `reduced_observations` case (which zero-fills every cell), a cell
        with no points inside it is simply absent from a non-empty result,
        since the point-in-cell join is an inner join.
        """
        cells = gpd.GeoDataFrame(
            [self._make_cell(0, 0, 0, 1, 1), self._make_cell(1, 2, 2, 3, 3)],
            crs="EPSG:4326",
        )
        reduced = gpd.GeoDataFrame(
            [self._make_reduced_observation(0, 0.5, 0.5, 5, 10, 100)],
            crs="EPSG:4326",
        )
        result = grid_observations(reduced, cells)
        self.assertEqual(len(result.index), 1)
        self.assertEqual(result.iloc[0].cell_id, 0)


class TestCollectObservationsGeometry(IssConstellationTestCase):
    """
    Unit tests for collecting observations of shapely points.
    """

    def setUp(self):
        super().setUp()
        self.narrow = Instrument(name="Narrow", field_of_regard=60)
        self.narrow_satellite = Satellite(
            name="Narrow", orbit=self.orbit, instruments=[self.narrow]
        )
        self.start = datetime(2022, 6, 1, tzinfo=timezone.utc)
        self.end = self.start + timedelta(days=1)

    def test_shapely_point_matches_point(self):
        """
        Test that a shapely point is observed exactly as the equivalent
        TAT-C point (other than its identifier).
        """
        expected = collect_observations(
            Point(id=5, latitude=40.74, longitude=-74.03),
            self.satellite,
            self.start,
            self.end,
        )
        actual = collect_observations(
            ShapelyPoint(-74.03, 40.74),
            self.satellite,
            self.start,
            self.end,
        )
        columns = ["start", "end", "epoch", "sat_alt", "sat_az"]
        self.assertGreater(len(expected.index), 0)
        pd.testing.assert_frame_equal(expected[columns], actual[columns])

    def test_shapely_point_default_identifier(self):
        """
        Test that observations of a shapely point record a `point_id` of 0
        by default.
        """
        observations = collect_observations(
            ShapelyPoint(-74.03, 40.74), self.satellite, self.start, self.end
        )
        self.assertGreater(len(observations.index), 0)
        self.assertTrue((observations.point_id == 0).all())

    def test_invalid_geometry_type(self):
        """
        Test that a geometry other than a point (including a region, see
        collect_region_observations) raises a TypeError.
        """
        for geometry in [LineString([(0, 0), (1, 1)]), box(0, 0, 1, 1)]:
            with self.subTest(geometry=geometry.geom_type):
                with self.assertRaises(TypeError):
                    collect_observations(geometry, self.satellite, self.start, self.end)

    def test_refine_access_periods_splits_period(self):
        """
        Test that a period in which the residual changes sign several times
        is divided into the parts where it is not positive.
        """
        start = pd.Timestamp(self.start)

        def residual(orbit_track):
            seconds = (
                orbit_track.t.tt - timescale.from_datetime(self.start).tt
            ) * 86400
            return np.cos(2 * np.pi * np.asarray(seconds) / 600)

        periods = _refine_access_periods(
            residual,
            self.orbit,
            [pd.Interval(start, start + pd.Timedelta(seconds=1200))],
        )
        self.assertEqual(len(periods), 2)
        for period, (left, right) in zip(periods, [(150, 450), (750, 1050)]):
            self.assertAlmostEqual(
                (period.left - start).total_seconds(), left, delta=0.01
            )
            self.assertAlmostEqual(
                (period.right - start).total_seconds(), right, delta=0.01
            )

    def test_reductions_separate_points_sharing_identifier(self):
        """
        Test that aggregating and reducing observations of distinct shapely
        points that share the default identifier keeps the points apart,
        with the same results as TAT-C points with distinct identifiers.
        """
        points = [ShapelyPoint(-74.0, 40.7), ShapelyPoint(-118.2, 34.0)]
        shared = pd.concat(
            [
                collect_observations(point, self.satellite, self.start, self.end)
                for point in points
            ],
            ignore_index=True,
        )
        distinct = pd.concat(
            [
                collect_observations(
                    Point(id=i, latitude=point.y, longitude=point.x),
                    self.satellite,
                    self.start,
                    self.end,
                )
                for i, point in enumerate(points)
            ],
            ignore_index=True,
        )
        expected = reduce_observations(aggregate_observations(distinct))
        actual = reduce_observations(aggregate_observations(shared))
        self.assertTrue((actual.point_id == 0).all())
        self.assertEqual(len(actual.index), 2)
        self.assertTrue((actual.geometry.geom_type == "Point").all())
        actual = actual.set_index(actual.geometry.to_wkb())
        expected = expected.set_index(expected.geometry.to_wkb())
        for column in ["access", "revisit", "samples"]:
            pd.testing.assert_series_equal(
                actual[column], expected.loc[actual.index, column]
            )


class TestComputeAccessPeriods(IssConstellationTestCase):
    """
    Unit tests for `compute_access_periods`.
    """

    def setUp(self):
        super().setUp()
        self.start = datetime(2022, 6, 1, tzinfo=timezone.utc)
        self.end = self.start + timedelta(days=1)

    def test_elevation_angle_at_bounds(self):
        """
        Test that the satellite is at the minimum elevation angle at the
        start and end of each period that does not span the analysis
        period's bounds, and above it at its midpoint.
        """
        periods = compute_access_periods(
            Point(latitude=40.74, longitude=-74.03),
            self.satellite,
            self.start,
            self.end,
            10,
        )
        self.assertGreater(len(periods), 0)
        topos = wgs84.latlon(40.74, -74.03)
        for period in periods:
            times = [period.left, period.mid, period.right]
            track = self.orbit.get_orbit_track(times)
            altitude = (
                (track - topos.at(timescale.from_datetimes(times))).altaz()[0].degrees
            )
            np.testing.assert_allclose(altitude[[0, 2]], 10, atol=1e-3)
            self.assertGreater(altitude[1], 10)

    def test_shapely_point_matches_point(self):
        """
        Test that a shapely point has the same periods as the equivalent
        TAT-C point.
        """
        expected = compute_access_periods(
            Point(latitude=40.74, longitude=-74.03, elevation=100),
            self.satellite,
            self.start,
            self.end,
            10,
        )
        actual = compute_access_periods(
            ShapelyPoint(-74.03, 40.74, 100), self.satellite, self.start, self.end, 10
        )
        self.assertGreater(len(expected), 0)
        self.assertEqual(list(actual), list(expected))

    def test_rejects_constellation(self):
        """
        Test that a constellation raises a TypeError.
        """
        with self.assertRaises(TypeError):
            compute_access_periods(
                Point(latitude=0, longitude=0), self.constellation, self.start, self.end
            )

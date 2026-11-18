"""
Validates the ATMS catalog entry against the NOAA-20 ATMS geolocation product, using
the exact values of the entry through public TAT-C methods. The tests of pixels and
swath width follow `ValidateImagerATMS`, on whose day of data the entry's angles were
fitted, and those of sizes and spacing follow `ValidateSampleSizes`; the tests of
passes, epochs, and swath area are this script's own. Running on that day checks the
shared tests against the notebooks, while running on a day at least 30 days from it
tests the fitted angles on data held out from the fit.

Reads anonymous public data, and needs the `validation` extra:

    python validation/catalog/validate_catalog_atms.py <data directory> 2026-10-02
    python validation/catalog/validate_catalog_atms.py <data directory> 2026-08-12

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

NOAA-21's and Suomi NPP's ATMS have their own configurations, whose angles are fitted
with `--fit` (the notebooks' angles are NOAA-20's), and tested on a day at least 30
days from the fit:

    python validation/catalog/validate_catalog_atms.py <data directory> 2026-11-05 --platform noaa-21 --fit
    python validation/catalog/validate_catalog_atms.py <data directory> 2026-10-02 --platform noaa-21
    python validation/catalog/validate_catalog_atms.py <data directory> 2026-10-02 --platform suomi-npp --fit
    python validation/catalog/validate_catalog_atms.py <data directory> 2026-08-25 --platform suomi-npp

With `--transfer`, it tests NOAA-20's configurations on another platform instead, and
adds the results to the record as a transfer, which does not count toward the tier.
Without a day, it takes that of the entry's record (or, with `--fit`, of the fit),
and `--configuration` tests one configuration.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import concurrent.futures
import io
import pathlib
from datetime import datetime, timedelta

import h5py
import numpy as np
import pandas as pd
from skyfield.api import wgs84

from tatc.schemas import GeneralPerturbationsOrbit
import scanner
from validation_io import (
    ARCHIVES,
    EARTH_ROTATION,
    TO_GEODETIC,
    Archive,
    FIXTURES,
    ScanLattice,
    ecef,
    iet_to_utc,
    normalize,
    parse_run,
    up,
)

FIT_PLATFORM = "noaa-20"
# the configurations of each platform's fit and its fixture
CONFIGURATIONS = {
    "noaa-20": ("swath_noaa_20", "beams_noaa_20", FIXTURES / "atms_noaa_20.npz"),
    "noaa-21": ("swath_noaa_21", "beams_noaa_21", FIXTURES / "atms_noaa_21.npz"),
    "suomi-npp": (
        "swath_suomi_npp",
        "beams_suomi_npp",
        FIXTURES / "atms_suomi_npp.npz",
    ),
}
FIELDS = (
    "StartTime",
    "MidTime",
    "SCPosition",
    "SCVelocity",
    "Latitude",
    "Longitude",
    "Height",
    "QF1_ATMSSDRGEO",
)
N_BEAMS, SCAN_PERIOD = 96, 8 / 3
ORBIT_MARGIN = 128
# published footprints: instrument, altitude (km), inclination (deg), scan angle
# (deg), and sizes along and across track (km), from the TC PRIMED documentation
PUBLISHED_FOOTPRINTS = (
    ("ATMS channels 1-2", 824, 98.7, 0, 74.8, 74.8),
    ("ATMS channels 1-2", 824, 98.7, 52.725, 141.8, 323.1),
    ("ATMS channels 3-16", 824, 98.7, 0, 31.6, 31.6),
    ("ATMS channels 3-16", 824, 98.7, 52.725, 60.0, 136.7),
    ("ATMS channels 17-22", 824, 98.7, 0, 15.8, 15.8),
    ("ATMS channels 17-22", 824, 98.7, 52.725, 30.0, 68.4),
)


def get_setup(platform: str) -> scanner.Setup:
    """Gets what is validated with a platform's fit: its configurations and fixture."""
    swath, beams, fixture_path = CONFIGURATIONS[platform]
    return scanner.Setup(
        swath=swath,
        beams=beams,
        size_instrument="ATMS channels 3-16",
        published_footprints=PUBLISHED_FOOTPRINTS,
        fixture_path=fixture_path,
        reference="ATMS SDR geolocation",
        published_reference="Footprints of the TC PRIMED documentation; beam centers "
        "of the ATMS SDR geolocation",
        time_steps=(SCAN_PERIOD, 32, 128),
        pass_radius=60,
    )


def load_day(
    day: datetime, directory: pathlib.Path, archive: Archive
) -> dict[str, np.ndarray]:
    """Downloads and reduces the geolocation of a day, or loads it from a cache."""
    path = directory / f"atms_gatmo_{day:%Y%m%d}.npz"
    if not path.exists():

        def read(key):
            with h5py.File(io.BytesIO(archive.fetch(key))) as f:
                return {k: f["All_Data/ATMS-SDR-GEO_All/" + k][:] for k in FIELDS}

        keys = archive.list_keys(f"ATMS-SDR-GEO/{day:%Y/%m/%d}/")
        with concurrent.futures.ThreadPoolExecutor(8) as pool:
            granules = list(pool.map(read, keys))
        np.savez_compressed(
            path, **{k: np.concatenate([g[k] for g in granules]) for k in FIELDS}
        )
    return dict(np.load(path))


class Reference(ScanLattice):
    """
    A day of ATMS scans, with what is derived from them. The beams of a scan are
    sampled in sequence about its middle time.
    """

    def __init__(self, data: dict[str, np.ndarray]):
        good = (data["QF1_ATMSSDRGEO"] == 0) & np.all(
            np.abs(data["Latitude"]) <= 90, axis=1
        )
        self.mid_iet = data["MidTime"][good]
        self.time = iet_to_utc(self.mid_iet)
        self.start = iet_to_utc(data["StartTime"][good])
        self.position = data["SCPosition"][good].astype(float)
        self.velocity = data["SCVelocity"][good].astype(float)
        self.lat = data["Latitude"][good].astype(float)
        self.lon = data["Longitude"][good].astype(float)
        self.height = data["Height"][good].astype(float)
        self.unit = up(self.lon, self.lat)
        self.seconds = (self.time - self.time[0]).total_seconds().to_numpy()
        period = np.diff(self.seconds)
        self.is_next_contiguous = period < 1.5 * SCAN_PERIOD
        self.gaps = [
            (self.time[i], self.time[i + 1]) for i in np.nonzero(period > 60)[0]
        ]
        beam_period = (
            2 * np.median(self.mid_iet - data["StartTime"][good]) / 1e6 / N_BEAMS
        )
        self.beam_offsets = (np.arange(N_BEAMS) - (N_BEAMS - 1) / 2) * beam_period

    @property
    def sample_seconds(self) -> np.ndarray:
        """The time of each beam (s from the first scan), scans by beams."""
        return self.seconds[:, np.newaxis] + self.beam_offsets[np.newaxis, :]

    def pass_time(self, before: int, sample: int, fraction: float):
        """Gets the time of a pass between a beam of a scan and of the next."""
        seconds = (
            self.seconds[before]
            + self.beam_offsets[sample]
            + fraction * (self.seconds[before + 1] - self.seconds[before])
        )
        return self.time[0] + timedelta(seconds=float(seconds))

    @property
    def span(self) -> tuple[datetime, datetime]:
        """The start of the first scan and the end of the last."""
        return (
            self.start[0].to_pydatetime(),
            self.start[-1].to_pydatetime() + timedelta(seconds=SCAN_PERIOD),
        )

    def contiguous_scans(self) -> np.ndarray:
        """Gets the scans whose neighbors are one scan period before and after."""
        period = np.diff(self.seconds)
        regular = np.abs(period - SCAN_PERIOD) < 0.01
        return np.nonzero(regular[:-1] & regular[1:])[0] + 1

    def score_swath_width(self, instrument, orbit, setup, _attitude=None) -> float:
        """Gets the median absolute difference (km) from the width through the beams."""
        scans = np.arange(0, len(self.time), setup.width_stride)
        xyz = ecef(self.lon[scans], self.lat[scans], self.height[scans])
        observed = np.sum(np.linalg.norm(np.diff(xyz, axis=1), axis=-1), axis=1)
        altitude = wgs84.geographic_position_of(
            orbit.get_orbit_track(list(self.time[scans]))
        ).elevation.m
        predicted = np.vectorize(instrument.get_swath_width)(altitude)
        return float(np.median(np.abs(predicted - observed))) / 1e3

    def observe_area(self, orbit, grid, _setup):
        """Gets the cells compared over one orbit, and which of them ATMS observed."""
        return self.observe_orbit(orbit, grid, ORBIT_MARGIN)

    def measure_pixels(self, satellite, index: int, setup, _frame) -> np.ndarray:
        """
        Gets rows of seconds from the start of the day, pixel index, position
        difference (km) from the beam center, and the beam's latitude and longitude,
        for the pixels projected at the middle times of a sample of scans.
        """
        scans = np.arange(0, len(self.time), setup.pixel_stride)
        instrument = satellite.get_instrument(index)
        times = self.time[scans]
        track = satellite.orbit.to_gp_orbit().get_orbit_track(list(times))
        seconds = (times - times[0].floor("D")).total_seconds().to_numpy()
        rows = []
        for beam in range(N_BEAMS):
            # TAT-C numbers pixels from the right of the motion, and ATMS from the left
            pixel = N_BEAMS - 1 - beam
            position = instrument.compute_projected_pixel_position(track, pixel, 0)
            difference = np.linalg.norm(
                ecef(position.longitude.degrees, position.latitude.degrees)
                - ecef(self.lon[scans, beam], self.lat[scans, beam]),
                axis=-1,
            )
            rows.append(
                np.stack(
                    [
                        seconds,
                        np.full(len(scans), pixel),
                        difference / 1e3,
                        self.lat[scans, beam],
                        self.lon[scans, beam],
                    ],
                    axis=-1,
                )
            )
        return np.concatenate(rows)

    def measure_spacing(self, setup) -> tuple[pd.DatetimeIndex, np.ndarray]:
        """
        Measures the spacing of beams as `ValidateSampleSizes` does, from half the
        distance between the beams on either side, with each beam at the time it is
        sampled. Gets the times, and rows of latitude, longitude, along-track and
        cross-track spacing (km).
        """
        scans = self.contiguous_scans()[:: setup.spacing_stride]
        xyz = ecef(self.lon, self.lat)
        along = np.linalg.norm(xyz[scans + 1] - xyz[scans - 1], axis=-1) / 2e3
        cross = np.full(along.shape, np.nan)
        cross[:, 1:-1] = np.linalg.norm(xyz[scans, 2:] - xyz[scans, :-2], axis=-1) / 2e3
        times = iet_to_utc(
            (
                self.mid_iet[scans, np.newaxis] + 1e6 * self.beam_offsets[np.newaxis, :]
            ).ravel()
        )
        rows = np.stack(
            [
                self.lat[scans].ravel(),
                self.lon[scans].ravel(),
                along.ravel(),
                cross.ravel(),
            ],
            axis=-1,
        )
        keep = np.isfinite(rows[:, 3])
        return times[keep], rows[keep]


def fit(day: datetime, reference: Reference) -> None:
    """
    Prints the angles fitted to a day in the inertial frame: those of the outermost
    beams, and the spacing and roll of a line fitted to the beams' angles (with
    `scanner.fit`); and the slope of the scan line and the along-track angle
    at nadir, with the spacecraft placed at each beam's time, which removes the skew
    of beams sampled in sequence. The slope is that of `scanner.fit_yaw`, so
    that the attitude's yaw angle is its negative, and its pitch angle the
    along-track angle.
    """
    lon, lat, altitude = TO_GEODETIC.transform(*reference.position.T)
    scanner.fit(
        scanner.Reference(
            day,
            {
                "time": (reference.time - pd.Timestamp(day)).total_seconds().to_numpy(),
                "lat": reference.lat,
                "lon": reference.lon,
                "sc_lat": np.asarray(lat),
                "sc_lon": np.asarray(lon),
                "sc_alt": np.asarray(altitude),
            },
            SCAN_PERIOD,
        ),
        "inertial",
    )
    zenith = up(np.asarray(lon), np.asarray(lat))
    motion = reference.velocity + np.cross([0, 0, EARTH_ROTATION], reference.position)
    along_axis = normalize(motion - np.sum(motion * zenith, axis=-1)[:, None] * zenith)
    left_axis = np.cross(zenith, along_axis)
    # the spacecraft's position (Earth-fixed) at each beam's time
    spacecraft = (
        reference.position[:, np.newaxis]
        + reference.velocity[:, np.newaxis] * reference.beam_offsets[:, np.newaxis]
    )
    rays = ecef(reference.lon, reference.lat, reference.height) - spacecraft
    down = -np.sum(rays * zenith[:, np.newaxis], axis=-1)
    along = np.sum(rays * along_axis[:, np.newaxis], axis=-1) / down
    cross = np.sum(rays * left_axis[:, np.newaxis], axis=-1) / down
    centered = cross - cross.mean(axis=1, keepdims=True)
    slope = np.degrees(
        np.arctan(np.sum(centered * along, axis=1) / np.sum(centered**2, axis=1))
    )
    middle = np.degrees(
        np.arctan(along[:, N_BEAMS // 2 - 1 : N_BEAMS // 2 + 1].mean(axis=1))
    )
    print(
        f"inertial, with beams at their times: slope {np.median(slope):+.4f} deg "
        f"(p5 {np.percentile(slope, 5):+.4f}, p95 {np.percentile(slope, 95):+.4f}); "
        f"along-track angle at nadir {np.median(middle):+.4f} deg (p5 "
        f"{np.percentile(middle, 5):+.4f}, p95 {np.percentile(middle, 95):+.4f})"
    )


def main() -> None:
    """Validates the entry on a day of data, or fits its values."""
    run = parse_run("atms", __doc__)
    archive = ARCHIVES[run.platform]
    directory = (
        run.directory
        if run.platform == FIT_PLATFORM
        else run.directory / archive.prefix
    )
    directory.mkdir(parents=True, exist_ok=True)
    reference = Reference(load_day(run.day, directory, archive))
    if run.fit:
        fit(run.day, reference)
        return
    setup = get_setup(FIT_PLATFORM if run.transfer else run.platform)
    tle = run.element_set(directory=directory)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    print(
        f"{run.day:%Y-%m-%d}: {tle[0][:32]}; {len(reference.time)} scans, "
        f"{len(reference.gaps)} gaps over a minute; configurations "
        f"{', '.join(run.configurations)}",
        flush=True,
    )
    if setup.beams in run.configurations:
        scanner.validate_beams(setup, run, reference, orbit, tle)
    if setup.swath in run.configurations:
        scanner.validate_swath(setup, run, reference, orbit, tle)


if __name__ == "__main__":
    main()

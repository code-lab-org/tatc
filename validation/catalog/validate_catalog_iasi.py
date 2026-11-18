"""
Validates the IASI catalog entry against the pixel positions of Metop-B's sea surface
temperature product (OSI SAF, distributed by NASA PO.DAAC), which holds the latitude and
longitude of all four pixels of each of 30 fields of regard, over land as well as sea,
using the exact values of the entry through public TAT-C methods.

IASI steps across 30 fields of regard in about 6.3 s, every 8 s. The center of the four
pixels of a field of regard is compared with TAT-C's pixels, as the array of four turns
with the scan and no single rectangular view represents it. The time across the scan
tilts the scan line relative to its position at the middle of the scan, which the pixel
configuration represents with a constant yaw of the satellite, not a physical one. The
swath configuration has no yaw, as IASI samples a point when it is abeam.

The product's times are in whole seconds, and the phase of its scans changes between
days while staying fixed within a day, so they carry an offset that differs by day. They
are not calibrated against the orbit, which is what they test, so epochs and pixel
positions include this uncertainty.

* The swath is that between the outermost centers of fields of regard, whose width and
  area are compared as for ATMS, with the area over one orbit.
* Passes are those of lattice points inside the cells between centers, with the times
  those of the nearest field of regard, interpolated along track.

The entry's angles come from `--fit`, on a day of data. The product does not record
satellite positions, and no element set is near that day, so positions come from the
ground track of the central fields of regard and heights from a later element set, in
the same part of an orbit whose shape is maintained. Running on a day at least 30 days
from it tests the angles on data held out from the fit. That element set, with an
epoch of 2026-11-07, is used on both days (as registered in
`elements/element_sets.yml`).

With `--platform metop-c`, it tests the same values on Metop-C, whose pixel positions
come from the IASI carbon monoxide product (EO:EUM:DAT:METOP:IASIL2COX, in BUFR, of
the EUMETSAT Data Store), which holds every pixel of every scan, with times in whole
seconds; the results are recorded as a transfer, which does not count toward the tier.
It needs a login to the Data Store (`eumdac auth --interactive`).

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`:

    python validation/catalog/validate_catalog_iasi.py <data directory> 2026-10-03 --fit
    python validation/catalog/validate_catalog_iasi.py <data directory> 2026-11-04

Without a day, it takes that of the entry's record (or, with `--fit`, of the fit),
and `--configuration` tests one configuration.

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import concurrent.futures
import json
import pathlib
from datetime import datetime, timedelta

import earthaccess
import h5py
import numpy as np
import pandas as pd
import requests
from skyfield.framelib import itrs

from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import GeneralPerturbationsOrbit
import scanner
from validation_io import (
    EARTH_ROTATION,
    FIXTURES,
    TO_GEODETIC,
    ScanLattice,
    build_grid,
    download_eumetsat,
    ecef,
    login_earthdata,
    normalize,
    parse_run,
    pick_passes,
    recompute_passes,
    score_passes,
    score_swath_area,
    up,
    write_results,
)

PRODUCT = "IASI_SST_METOP_B-OSISAF-L2P-v1.0"
# the product of Metop-C's pixel positions, in the EUMETSAT Data Store
COX_COLLECTION = "EO:EUM:DAT:METOP:IASIL2COX"

N_FORS, N_PIXELS = 30, 4
SCAN_PERIOD = 8.0
FOR_PERIOD = 0.2181  # s, the slope with field of regard of the product's times
MIDDLE = (N_FORS - 1) / 2
ORBIT_MARGIN = 128
GATED = (
    "pass_recall",
    "pass_precision",
    "epoch_p95",
    "epoch_median",
    "swath_width_median",
    "footprint_iou_median",
)
# published footprints: instrument, altitude (km), inclination (deg), scan angle (deg),
# and sizes along and across track (km)
PUBLISHED_FOOTPRINTS = (
    ("IASI", 817, 98.7, 0, 12.0, 12.0),
    ("IASI", 817, 98.7, 48.3, 20.0, 39.0),
)
FIT_SCAN_STRIDE = 10
SETUP = scanner.Setup(
    swath="swath",
    beams="pixels",
    size_instrument="IASI",
    published_footprints=PUBLISHED_FOOTPRINTS,
    fixture_path=FIXTURES / "iasi.npz",
    reference="Pixel positions of the OSI SAF IASI sea surface temperature product",
    published_reference="Published footprints; pixel positions of the OSI SAF IASI "
    "sea surface temperature product",
    time_steps=(SCAN_PERIOD, 32, 128),
    pass_radius=130,
    pixel_stride=30,
    spacing_stride=40,
    width_stride=5,
)


def reduce_granule(url: str, session, directory: pathlib.Path) -> pathlib.Path:
    """Reduces a granule to its start time and pixel positions, or finds it in a cache."""
    path = directory / (pathlib.PurePosixPath(url).stem + ".npz")
    if path.exists():
        return path
    response = session.get(url, timeout=300)
    response.raise_for_status()
    temporary_file = path.with_suffix(".tmp.nc")
    temporary_file.write_bytes(response.content)
    with h5py.File(temporary_file, "r") as f:
        start = f.attrs["time_coverage_start"].decode()
        arrays = {
            "start": np.array(
                pd.Timestamp(start, tz="UTC").as_unit("ns").value, dtype="int64"
            ),
            "latitude": f["lat"][:].astype(np.float32),
            "longitude": f["lon"][:].astype(np.float32),
        }
    temporary_file.unlink()
    # written whole, so that an interrupted download leaves no partial granule
    temporary = path.with_suffix(".tmp.npz")
    np.savez_compressed(temporary, **arrays)
    temporary.replace(path)
    return path


def load_day(day: datetime, directory: pathlib.Path) -> dict:
    """Downloads and reduces the pixel positions of a day, or loads them from a cache."""
    path = directory / f"iasi_{day:%Y%m%d}.npz"
    if not path.exists():
        login_earthdata()
        session = earthaccess.get_requests_https_session()
        entries = requests.get(
            "https://cmr.earthdata.nasa.gov/search/granules.json",
            params={
                "short_name": PRODUCT,
                "temporal": f"{day:%Y-%m-%d}T00:00:00Z,{day:%Y-%m-%d}T23:59:59Z",
                "page_size": 2000,
            },
            timeout=120,
        ).json()["feed"]["entry"]
        urls = sorted(
            next(
                link["href"]
                for link in entry["links"]
                if link["href"].endswith(".nc") and link["href"].startswith("https")
            )
            for entry in entries
        )
        granules = directory / f"iasi_{day:%Y%m%d}"
        granules.mkdir(parents=True, exist_ok=True)
        with concurrent.futures.ThreadPoolExecutor(6) as pool:
            paths = list(
                pool.map(lambda url: reduce_granule(url, session, granules), urls)
            )
        starts, seconds, latitudes, longitudes = [], [], [], []
        for granule in sorted(paths):
            with np.load(granule) as data:
                n_scans = len(data["latitude"])
                starts.append(np.full(n_scans, data["start"]))
                seconds.append(np.arange(n_scans) * SCAN_PERIOD)
                latitudes.append(data["latitude"])
                longitudes.append(data["longitude"])
        np.savez_compressed(
            path,
            start=np.concatenate(starts),
            offset=np.concatenate(seconds),
            latitude=np.concatenate(latitudes),
            longitude=np.concatenate(longitudes),
        )
        print(f"{len(urls)} granules", flush=True)
    return dict(np.load(path))


def read_cox(path: pathlib.Path) -> list[tuple]:
    """
    Reads the scans of a file of the IASI carbon monoxide product (BUFR, a scan in
    each message): each scan's time and its pixels' latitudes and longitudes, in the
    order of the fields of view (the four pixels of each field of regard in turn).
    """
    # imported here, as only Metop-C's product needs it
    import eccodes  # pylint: disable=import-outside-toplevel

    scans = []
    with open(path, "rb") as f:
        while (message := eccodes.codes_bufr_new_from_file(f)) is not None:
            try:
                eccodes.codes_set(message, "unpack", 1)
                fov = eccodes.codes_get_array(message, "fieldOfViewNumber")
                if len(fov) != N_FORS * N_PIXELS:
                    continue
                order = np.argsort(fov)
                stamp = [
                    int(eccodes.codes_get_array(message, key)[0])
                    for key in ("year", "month", "day", "hour", "minute", "second")
                ]
                scans.append(
                    (
                        pd.Timestamp(datetime(*stamp), tz="UTC"),
                        eccodes.codes_get_array(message, "latitude")[order],
                        eccodes.codes_get_array(message, "longitude")[order],
                    )
                )
            finally:
                eccodes.codes_release(message)
    return scans


def load_day_cox(day: datetime, directory: pathlib.Path) -> dict:
    """
    Downloads and reduces Metop-C's pixel positions of a day, from the IASI carbon
    monoxide product, or loads them from a cache.
    """
    path = directory / f"iasi_metopc_{day:%Y%m%d}.npz"
    if not path.exists():
        files = download_eumetsat(
            COX_COLLECTION,
            day - timedelta(minutes=5),
            day + timedelta(days=1, minutes=5),
            directory / f"cox_{day:%Y%m%d}",
            match="METOPC+IASI",
            entry=".bin",
        )
        scans = {}
        for file in files:
            for time, latitude, longitude in read_cox(file):
                if day <= time.to_pydatetime() < day + timedelta(days=1):
                    scans[time] = (latitude, longitude)
        times = sorted(scans)
        np.savez_compressed(
            path,
            start=np.array([t.as_unit("ns").value for t in times], dtype="int64"),
            offset=np.zeros(len(times)),
            latitude=np.array([scans[t][0] for t in times], dtype=np.float32),
            longitude=np.array([scans[t][1] for t in times], dtype=np.float32),
        )
        print(f"{len(files)} files, {len(times)} scans", flush=True)
    return dict(np.load(path))


class Reference(ScanLattice):
    """
    A day of IASI scans, with what is derived from them. Each field of regard is
    at its time from the middle of its scan.
    """

    def __init__(self, data: dict[str, np.ndarray], day: datetime):
        self.day = day
        latitude = data["latitude"].astype(float)
        longitude = data["longitude"].astype(float)
        good = np.all(np.isfinite(latitude), axis=1) & np.all(
            np.isfinite(longitude), axis=1
        )
        self.n_dropped = int(np.sum(~good))
        start = pd.to_datetime(data["start"][good], unit="ns", utc=True)
        self.time = start + pd.to_timedelta(data["offset"][good], unit="s")
        pixels = ecef(longitude[good], latitude[good])
        n_scans = len(pixels)
        # the center of the pixels of each field of regard
        self.xyz = pixels.reshape(n_scans, N_FORS, N_PIXELS, 3).mean(axis=2)
        lon, lat, _ = TO_GEODETIC.transform(*self.xyz.reshape(-1, 3).T)
        self.lon = np.asarray(lon).reshape(n_scans, N_FORS)
        self.lat = np.asarray(lat).reshape(n_scans, N_FORS)
        self.unit = up(self.lon, self.lat)
        self.seconds = (self.time - self.time[0]).total_seconds().to_numpy()
        period = np.diff(self.seconds)
        self.is_next_contiguous = period < 1.5 * SCAN_PERIOD
        self.gaps = [
            (self.time[i], self.time[i + 1]) for i in np.nonzero(period > 60)[0]
        ]
        self.scan_period = float(np.median(period))
        # from the start of the scan, as the fitted constant of the product's times is
        # small; a scan's time here is that of its middle
        self.for_offsets = (np.arange(N_FORS) - MIDDLE) * FOR_PERIOD
        self.middle_time = self.time + pd.Timedelta(seconds=MIDDLE * FOR_PERIOD)
        self.middle_seconds = self.seconds + MIDDLE * FOR_PERIOD

    @property
    def sample_seconds(self) -> np.ndarray:
        """The time of each field of regard (s from the first scan)."""
        return self.middle_seconds[:, np.newaxis] + self.for_offsets[np.newaxis, :]

    def pass_time(self, before: int, sample: int, fraction: float):
        """Gets the time of a pass between a field of regard of a scan and the next."""
        seconds = (
            self.middle_seconds[before]
            + self.for_offsets[sample]
            + fraction * (self.middle_seconds[before + 1] - self.middle_seconds[before])
        )
        return self.time[0] + timedelta(seconds=float(seconds))

    @property
    def span(self) -> tuple[datetime, datetime]:
        """The first and last scans' times, widened by a scan period."""
        margin = timedelta(seconds=SCAN_PERIOD)
        return (
            self.time[0].floor("us").to_pydatetime() - margin,
            self.time[-1].floor("us").to_pydatetime() + margin,
        )

    def contiguous_scans(self) -> np.ndarray:
        """Gets the scans whose neighbors are one scan period before and after."""
        period = np.diff(self.seconds)
        regular = np.abs(period - self.scan_period) < 0.01
        return np.nonzero(regular[:-1] & regular[1:])[0] + 1

    def score_swath_width(self, instrument, orbit, setup, _attitude=None) -> float:
        """Gets the median absolute difference (km) from the width through the centers."""
        scans = np.arange(0, len(self.time), setup.width_stride)
        observed = np.sum(
            np.linalg.norm(np.diff(self.xyz[scans], axis=1), axis=-1), axis=1
        )
        track = orbit.to_gp_orbit().get_orbit_track(
            [t.floor("us").to_pydatetime() for t in self.middle_time[scans]]
        )
        altitude = np.asarray(TO_GEODETIC.transform(*track.frame_xyz(itrs).m)[2])
        predicted = np.vectorize(instrument.get_swath_width)(altitude)
        print(
            f"swath width: observed {np.median(observed) / 1e3:.2f} km, "
            f"predicted {np.median(predicted) / 1e3:.2f} km",
            flush=True,
        )
        return float(np.median(np.abs(predicted - observed))) / 1e3

    def observe_area(self, orbit, grid, _setup):
        """
        Gets the cells compared over one orbit from noon of the day of the middle
        scan, and which of them IASI observed, with the scans at their middle times.
        """
        return self.observe_orbit(
            orbit,
            grid,
            ORBIT_MARGIN,
            begin=self.time[len(self.time) // 2].floor("D") + pd.Timedelta(hours=12),
            times=self.middle_time,
            floor=True,
        )

    def measure_pixels(self, satellite, index: int, setup, _frame) -> np.ndarray:
        """
        Gets rows of seconds from the start of the day, pixel index, position
        difference (km) from the center of a field of regard, and that center's
        latitude and longitude, for pixels projected at the middle times of a sample
        of scans. Prints the median difference of some pixels.
        """
        scans = np.arange(0, len(self.time), setup.pixel_stride)
        instrument = satellite.get_instrument(index)
        times = self.middle_time[scans]
        track = satellite.orbit.to_gp_orbit().get_orbit_track(
            [t.floor("us").to_pydatetime() for t in times]
        )
        seconds = (times - pd.Timestamp(self.day)).total_seconds().to_numpy()
        rows = []
        for column in range(N_FORS):
            # TAT-C numbers pixels from the right of the motion, and IASI from the left
            pixel = N_FORS - 1 - column
            position = instrument.compute_projected_pixel_position(track, pixel, 0)
            difference = np.linalg.norm(
                ecef(position.longitude.degrees, position.latitude.degrees)
                - self.xyz[scans, column],
                axis=-1,
            )
            rows.append(
                np.stack(
                    [
                        seconds,
                        np.full(len(scans), pixel),
                        difference / 1e3,
                        self.lat[scans, column],
                        self.lon[scans, column],
                    ],
                    axis=-1,
                )
            )
        pixels = np.concatenate(rows)
        for pixel in (0, 7, 14, 22, 29):
            print(
                f"  pixel {pixel}: median "
                f"{np.median(pixels[pixels[:, 1] == pixel, 2]):.2f} km"
            )
        return pixels

    def measure_spacing(self, setup):
        """
        Measures the spacing of fields of regard, from half the distance between
        those on either side. Gets the times, and rows of latitude, longitude,
        along-track and cross-track spacing (km).
        """
        scans = self.contiguous_scans()[:: setup.spacing_stride]
        along = np.linalg.norm(self.xyz[scans + 1] - self.xyz[scans - 1], axis=-1) / 2e3
        cross = np.full(along.shape, np.nan)
        cross[:, 1:-1] = (
            np.linalg.norm(self.xyz[scans, 2:] - self.xyz[scans, :-2], axis=-1) / 2e3
        )
        times = self.time[0] + pd.to_timedelta(
            (
                self.middle_seconds[scans, np.newaxis] + self.for_offsets[np.newaxis, :]
            ).ravel(),
            unit="s",
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


def build_heights(orbit, day: datetime):
    """Gets geodetic heights of an orbit, by latitude, on its two sides."""
    period = 360 / orbit.elements[0].mean_motion
    seconds = np.arange(0, period, 10.0)
    track = orbit.to_gp_orbit().get_orbit_track(
        [day + timedelta(seconds=float(s)) for s in seconds]
    )
    _, lat, height = TO_GEODETIC.transform(*track.frame_xyz(itrs).m)
    lat, height = np.asarray(lat), np.asarray(height)
    ascending = np.gradient(lat) > 0
    tables = {}
    for side in (True, False):
        order = np.argsort(lat[ascending == side])
        tables[side] = (lat[ascending == side][order], height[ascending == side][order])
    return tables


def fit_geometry(reference: Reference, orbit, day: datetime) -> None:
    """Measures the angles of the scan from the ground track of its central positions."""
    tables = build_heights(orbit, day)
    scans = reference.contiguous_scans()[::FIT_SCAN_STRIDE]

    def nadir(rows):
        return (reference.xyz[rows, 14] + reference.xyz[rows, 15]) / 2

    lon, lat, _ = TO_GEODETIC.transform(*nadir(scans).T)
    lon, lat = np.asarray(lon), np.asarray(lat)
    _, after, _ = TO_GEODETIC.transform(*nadir(scans + 1).T)
    ascending = np.asarray(after) > lat
    height = np.where(
        ascending,
        np.interp(lat, *tables[True]),
        np.interp(lat, *tables[False]),
    )

    def position_of(rows):
        row_lon, row_lat, _ = TO_GEODETIC.transform(*nadir(rows).T)
        return ecef(np.asarray(row_lon), np.asarray(row_lat), height)

    position = position_of(scans)
    velocity = (position_of(scans + 1) - position_of(scans - 1)) / (
        reference.seconds[scans + 1] - reference.seconds[scans - 1]
    )[:, None]
    zenith = up(lon, lat)
    print(f"{len(scans)} scans, median height {np.median(height) / 1e3:.1f} km")
    for frame in ("earth_fixed", "inertial"):
        motion = velocity
        if frame == "inertial":
            motion = velocity + np.cross([0, 0, EARTH_ROTATION], position)
        along_axis = normalize(
            motion - np.sum(motion * zenith, axis=-1)[:, None] * zenith
        )
        left_axis = np.cross(zenith, along_axis)
        rays = reference.xyz[scans] - position[:, None]
        down = -np.sum(rays * zenith[:, None], axis=-1)
        along = np.sum(rays * along_axis[:, None], axis=-1) / down
        cross = np.sum(rays * left_axis[:, None], axis=-1) / down
        # the slope of along-track on cross-track tangents is the tangent of a yaw
        slope = np.sum(along * cross, axis=1) / np.sum(cross * cross, axis=1)
        for name, chosen in (("ascending", ascending), ("descending", ~ascending)):
            print(
                f"{frame}, {name}: first {np.median(np.degrees(np.arctan(cross[chosen, 0]))):+.3f}"
                f", last {np.median(np.degrees(np.arctan(cross[chosen, -1]))):+.3f} deg; "
                f"along at first {np.median(np.degrees(np.arctan(along[chosen, 0]))):+.3f}"
                f", middle {np.median(np.degrees(np.arctan(along[chosen, 14:16]))):+.3f}"
                f", last {np.median(np.degrees(np.arctan(along[chosen, -1]))):+.3f} deg; "
                f"slope {np.median(np.degrees(np.arctan(slope[chosen]))):+.4f} deg"
            )
        angles = np.degrees(np.arctan(cross))
        fit = np.polyfit(np.arange(N_FORS), np.median(angles, axis=0), 1)
        residual = np.median(angles, axis=0) - np.polyval(fit, np.arange(N_FORS))
        print(
            f"{frame}: spacing {-fit[0]:.4f} deg, roll "
            f"{np.polyval(fit, MIDDLE):+.4f} deg, largest residual "
            f"{np.abs(residual).max():.4f} deg; extent "
            f"{np.abs(np.median(angles[:, 0] - angles[:, -1])):.4f} deg; yaw slope "
            f"{np.median(np.degrees(np.arctan(slope))):+.4f} deg, pitch "
            f"{np.median(np.degrees(np.arctan(along[:, 14:16]))):+.4f} deg"
        )
    width = np.sum(
        np.linalg.norm(np.diff(reference.xyz[scans], axis=1), axis=-1), axis=1
    )
    print(
        f"swath width: median {np.median(width) / 1e3:.1f} km, "
        f"p5 {np.percentile(width, 5) / 1e3:.1f}, p95 {np.percentile(width, 95) / 1e3:.1f}"
    )


def validate_swath(setup: scanner.Setup, run, reference: Reference, orbit, tle):
    """
    Scores the swath configuration of a run, and writes its results if the run
    writes. As `scanner.validate_swath`, but its record also gives the median
    difference of epochs.
    """
    entry, day = run.entry, run.day

    def satellite_for(time_step: float):
        return entry.to_satellite(
            orbit, configuration=setup.swath, time_step=time_step, platform=run.platform
        )

    satellite = satellite_for(setup.time_steps[0])
    points = generate_points_fibonacci_lattice(setup.point_spacing)
    observed = reference.observe_passes(points, setup.pass_radius)
    start, end = reference.span
    results = score_passes(points, satellite, 0, observed, start, end, reference.gaps)
    results["swath_width_median"] = reference.score_swath_width(
        satellite.instruments[0], orbit, setup
    )
    print(json.dumps(results), flush=True)
    grid = build_grid()
    first, last, compared, covered = reference.observe_area(orbit, grid, setup)
    results["footprint_iou_median"] = score_swath_area(
        satellite_for, 0, first, last, setup.time_steps, grid, compared, covered
    )
    print(json.dumps(results), flush=True)
    if not run.write:
        return
    passes = pick_passes(points, observed, day, setup.fixture_passes)
    write_results(
        run,
        setup.swath,
        reference=setup.reference,
        metrics=[{"name": name, "value": float(results[name])} for name in GATED],
        fixture={
            "tle": np.array(tle),
            "day": np.array(day.isoformat()),
            "passes": passes,
        },
        fixture_path=setup.fixture_path,
        fixture_metrics=recompute_passes(satellite, 0, day, passes),
    )


def main() -> None:
    """Fits or validates the entry on a day of data."""
    run = parse_run("iasi", __doc__)
    tle = run.element_set()
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    data = (
        load_day_cox(run.day, run.directory)
        if run.platform == "metop-c"
        else load_day(run.day, run.directory)
    )
    reference = Reference(data, run.day)
    print(
        f"{run.day:%Y-%m-%d}: {len(reference.time)} scans, {reference.n_dropped} "
        f"dropped, {len(reference.gaps)} gaps over a minute, "
        f"median scan period {reference.scan_period:.4f} s",
        flush=True,
    )
    if run.fit:
        fit_geometry(reference, orbit, run.day)
        return
    if SETUP.beams in run.configurations:
        scanner.validate_beams(SETUP, run, reference, orbit, tle)
    if SETUP.swath in run.configurations:
        validate_swath(SETUP, run, reference, orbit, tle)


if __name__ == "__main__":
    main()

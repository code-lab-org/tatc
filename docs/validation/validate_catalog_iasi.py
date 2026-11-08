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
from it tests the angles on data held out from the fit.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`:

    python docs/validation/validate_catalog_iasi.py <data directory> 2026-10-03 --fit
    python docs/validation/validate_catalog_iasi.py <data directory> 2026-11-04

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import concurrent.futures
import json
import pathlib
from datetime import datetime, timedelta, timezone

import earthaccess
import h5py
import numpy as np
import pandas as pd
import requests
from scipy.spatial import cKDTree
from skyfield.framelib import itrs

import tatc
from tatc.catalog import get_entry
from tatc.catalog.metrics import summarize
from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import GeneralPerturbationsOrbit
from validation_io import (
    TO_GEODETIC,
    build_grid,
    compute_published_errors,
    compute_spacing_ratios,
    ecef,
    get_holdout_gap,
    in_lattice_cells,
    normalize,
    predict_passes,
    match_passes,
    recompute_passes,
    score_published,
    score_swath_area,
    up,
    update_fixture,
    update_record,
)

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
RECORD_PATH = REPOSITORY / "src" / "tatc" / "resources" / "validation" / "iasi.json"
FIXTURE_PATH = REPOSITORY / "tests" / "catalog" / "fixtures" / "iasi.npz"

PRODUCT = "IASI_SST_METOP_B-OSISAF-L2P-v1.0"
FIT_DAY = datetime(2026, 10, 3, tzinfo=timezone.utc)
PLATFORM = "metop-b"
ELEMENT_SET = [
    "1 38771U 12049A   26311.29297005  .00000108  00000+0  68954-4 0  9993",
    "2 38771  98.6325 359.8075 0001272 306.3444  53.7616 14.21470410733593",
]

SWATH, PIXELS = "swath", "pixels"
N_FORS, N_PIXELS = 30, 4
SCAN_PERIOD = 8.0
FOR_PERIOD = 0.2181  # s, the slope with field of regard of the product's times
MIDDLE = (N_FORS - 1) / 2
TIME_STEPS = (SCAN_PERIOD, 32, 128)
ORBIT_MARGIN = 128
PASS_RADIUS = 2 * np.sin(130 / 6371 / 2)  # chord, on the unit sphere
POINT_SPACING = 500e3
EARTH_ROTATION = np.array([0, 0, 7.2921150e-5])
# published footprints: instrument, altitude (km), inclination (deg), scan angle (deg),
# and sizes along and across track (km)
PUBLISHED_FOOTPRINTS = (
    ("IASI", 817, 98.7, 0, 12.0, 12.0),
    ("IASI", 817, 98.7, 48.3, 20.0, 39.0),
)
FIT_SCAN_STRIDE, PIXEL_SCAN_STRIDE, SPACING_SCAN_STRIDE, WIDTH_SCAN_STRIDE = (
    10,
    30,
    40,
    5,
)
FIXTURE_PASSES, FIXTURE_PIXELS, FIXTURE_SPACING = 60, 200, 200


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
        if earthaccess.login(strategy="environment") is None:
            raise RuntimeError("Earthdata login failed: is EARTHDATA_TOKEN set?")
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


class Reference:
    """A day of IASI scans, with what is derived from them."""

    def __init__(self, data: dict[str, np.ndarray]):
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

    def in_cells(self, points: np.ndarray, first: int, last: int) -> np.ndarray:
        """Gets whether unit vectors lie in cells on a range of scans."""
        return in_lattice_cells(
            self.unit[first : last + 1], self.is_next_contiguous[first:last], points
        )

    def observe_passes(self, points) -> pd.DataFrame:
        """Finds the passes on which IASI observed points, and their times."""
        point_unit = up(points.geometry.x.to_numpy(), points.geometry.y.to_numpy())
        flat = self.unit.reshape(-1, 3)
        seconds_of = (
            self.middle_seconds[:, np.newaxis] + self.for_offsets[np.newaxis, :]
        ).ravel()
        tree = cKDTree(flat)
        rows = []
        for p, hits in enumerate(
            tree.query_ball_point(point_unit, r=PASS_RADIUS, workers=-1)
        ):
            if not hits:
                continue
            hits = np.array(hits)
            hits = hits[np.argsort(seconds_of[hits])]
            for group in np.split(
                hits, np.nonzero(np.diff(seconds_of[hits]) > 1200)[0] + 1
            ):
                k, j = np.divmod(
                    group[
                        np.argmin(np.linalg.norm(flat[group] - point_unit[p], axis=-1))
                    ],
                    N_FORS,
                )
                first, last = max(k - 1, 0), min(k + 1, len(self.time) - 1)
                if not self.in_cells(point_unit[[p]], first, last)[0]:
                    continue
                before = (
                    k - 1
                    if np.sum(
                        (point_unit[p] - self.unit[k, j])
                        * (self.unit[k, j] - self.unit[k - 1, j])
                    )
                    < 0
                    else k
                )
                before = int(np.clip(before, 0, len(self.time) - 2))
                a, b = self.unit[before, j], self.unit[before + 1, j]
                fraction = np.clip(
                    np.sum((point_unit[p] - a) * (b - a)) / np.sum((b - a) ** 2), 0, 1
                )
                seconds = (
                    self.middle_seconds[before]
                    + self.for_offsets[j]
                    + fraction
                    * (self.middle_seconds[before + 1] - self.middle_seconds[before])
                )
                rows.append(
                    {
                        "point_id": p,
                        "time": self.time[0] + timedelta(seconds=float(seconds)),
                    }
                )
        return pd.DataFrame(rows).astype({"point_id": "int64"})

    def contiguous_scans(self) -> np.ndarray:
        """Gets the scans whose neighbors are one scan period before and after."""
        period = np.diff(self.seconds)
        regular = np.abs(period - self.scan_period) < 0.01
        return np.nonzero(regular[:-1] & regular[1:])[0] + 1


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
            motion = velocity + np.cross(EARTH_ROTATION, position)
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


def score_swath_width(reference: Reference, instrument, orbit) -> float:
    """Gets the median absolute difference (km) from the width through the centers."""
    scans = np.arange(0, len(reference.time), WIDTH_SCAN_STRIDE)
    observed = np.sum(
        np.linalg.norm(np.diff(reference.xyz[scans], axis=1), axis=-1), axis=1
    )
    track = orbit.to_gp_orbit().get_orbit_track(
        [t.floor("us").to_pydatetime() for t in reference.middle_time[scans]]
    )
    altitude = np.asarray(TO_GEODETIC.transform(*track.frame_xyz(itrs).m)[2])
    predicted = np.vectorize(instrument.get_swath_width)(altitude)
    print(
        f"swath width: observed {np.median(observed) / 1e3:.2f} km, "
        f"predicted {np.median(predicted) / 1e3:.2f} km",
        flush=True,
    )
    return float(np.median(np.abs(predicted - observed))) / 1e3


def observe_orbit(reference: Reference, orbit, grid):
    """Gets the cells compared over one orbit, and which of them IASI observed."""
    lon, lat, _ = grid
    unit = up(lon, lat)
    begin = reference.time[len(reference.time) // 2].floor("D") + pd.Timedelta(hours=12)
    period = pd.Timedelta(days=1 / orbit.elements[0].mean_motion * 360 / 86400)
    scans = np.nonzero((reference.time >= begin) & (reference.time < begin + period))[0]
    first, last = reference.middle_time[scans[0]], reference.middle_time[scans[-1]]
    distance, nearest = cKDTree(reference.unit[scans].reshape(-1, 3)).query(
        unit, workers=-1
    )
    elapsed = (reference.middle_time[scans] - first).total_seconds().to_numpy()
    nearest_time = elapsed[nearest // N_FORS]
    compared = (
        (distance < 2 * np.sin(300 / 6371 / 2))
        & (nearest_time > ORBIT_MARGIN)
        & (nearest_time < (last - first).total_seconds() - ORBIT_MARGIN)
    )
    observed = np.zeros(len(unit), dtype=bool)
    observed[compared] = reference.in_cells(unit[compared], scans[0], scans[-1])
    return (
        first.floor("us").to_pydatetime(),
        last.floor("us").to_pydatetime(),
        compared,
        observed,
    )


def measure_pixels(reference: Reference, satellite, day: datetime) -> np.ndarray:
    """
    Gets rows of seconds from the start of the day, pixel index, position difference
    (km) from the center of a field of regard, and that center's latitude and longitude,
    for pixels projected at the middle times of a sample of scans.
    """
    scans = np.arange(0, len(reference.time), PIXEL_SCAN_STRIDE)
    instrument = satellite.get_instrument(0)
    times = reference.middle_time[scans]
    track = satellite.orbit.to_gp_orbit().get_orbit_track(
        [t.floor("us").to_pydatetime() for t in times]
    )
    seconds = (times - pd.Timestamp(day)).total_seconds().to_numpy()
    rows = []
    for column in range(N_FORS):
        # TAT-C numbers pixels from the right of the motion, and IASI from the left
        pixel = N_FORS - 1 - column
        position = instrument.compute_projected_pixel_position(track, pixel, 0)
        difference = np.linalg.norm(
            ecef(position.longitude.degrees, position.latitude.degrees)
            - reference.xyz[scans, column],
            axis=-1,
        )
        rows.append(
            np.stack(
                [
                    seconds,
                    np.full(len(scans), pixel),
                    difference / 1e3,
                    reference.lat[scans, column],
                    reference.lon[scans, column],
                ],
                axis=-1,
            )
        )
    return np.concatenate(rows)


def measure_spacing(reference: Reference):
    """
    Measures the spacing of fields of regard, from half the distance between those on
    either side. Gets the times, and rows of latitude, longitude, along-track and
    cross-track spacing (km).
    """
    scans = reference.contiguous_scans()[::SPACING_SCAN_STRIDE]
    along = (
        np.linalg.norm(reference.xyz[scans + 1] - reference.xyz[scans - 1], axis=-1)
        / 2e3
    )
    cross = np.full(along.shape, np.nan)
    cross[:, 1:-1] = (
        np.linalg.norm(reference.xyz[scans, 2:] - reference.xyz[scans, :-2], axis=-1)
        / 2e3
    )
    times = reference.time[0] + pd.to_timedelta(
        (
            reference.middle_seconds[scans, np.newaxis]
            + reference.for_offsets[np.newaxis, :]
        ).ravel(),
        unit="s",
    )
    rows = np.stack(
        [
            reference.lat[scans].ravel(),
            reference.lon[scans].ravel(),
            along.ravel(),
            cross.ravel(),
        ],
        axis=-1,
    )
    keep = np.isfinite(rows[:, 3])
    return times[keep], rows[keep]


def validate_pixels(day, tle, entry, orbit, reference, write) -> None:
    """Scores the sample size configuration, and writes its results if requested."""
    within, largest = score_published(
        compute_published_errors(entry, PIXELS, PUBLISHED_FOOTPRINTS)
    )
    print(
        f"published dimensions: {within:.1f}% within 3%, largest error {largest:.2f}%",
        flush=True,
    )
    satellite = entry.to_satellite(orbit, PIXELS, platform=PLATFORM)
    pixels = measure_pixels(reference, satellite, day)
    median = summarize(pixels[:, 2])["median"]
    p95 = summarize(pixels[:, 2])["p95"]
    print(f"pixel position: median {median:.3f} km, 95th percentile {p95:.3f} km")
    for pixel in (0, 7, 14, 22, 29):
        print(
            f"  pixel {pixel}: median "
            f"{np.median(pixels[pixels[:, 1] == pixel, 2]):.2f} km"
        )
    times, rows = measure_spacing(reference)
    names = ["IASI"] * len(rows)
    along, cross = compute_spacing_ratios(satellite, names, times, rows)
    print(f"spacing, predicted / measured: {along:.4f} along, {cross:.4f} across")
    if not write:
        return
    rng = np.random.default_rng(3)
    picks = rng.choice(len(rows), FIXTURE_SPACING, replace=False)
    fixture_along, fixture_cross = compute_spacing_ratios(
        satellite, ["IASI"] * len(picks), times[picks], rows[picks]
    )
    reduced = pixels[rng.choice(len(pixels), FIXTURE_PIXELS, replace=False)]
    update_fixture(
        FIXTURE_PATH,
        tle=np.array(tle),
        day=np.array(day.isoformat()),
        pixels=reduced,
        spacing=rows[picks],
        spacing_seconds=(times[picks] - pd.Timestamp(day)).total_seconds().to_numpy(),
    )
    update_record(
        RECORD_PATH,
        "iasi",
        PIXELS,
        {
            "digest": entry.configurations[PIXELS].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": "docs/validation/validate_catalog_iasi.py",
            "reference": "Published footprints; pixel positions of the OSI SAF IASI "
            "sea surface temperature product",
            "period": f"{day:%Y-%m-%d}",
            "platform": PLATFORM,
            "holdout_gap_days": get_holdout_gap(day, FIT_DAY),
            "metrics": [
                {"name": "sample_size_within_3pct", "value": within},
                {"name": "sample_size_largest_error", "value": largest},
                {"name": "pixel_position_median", "value": median},
                {"name": "pixel_position_p95", "value": p95},
                {"name": "spacing_along_ratio_median", "value": along},
                {"name": "spacing_cross_ratio_median", "value": cross},
            ],
            "fixture_metrics": {
                "sample_size_within_3pct": within,
                "pixel_position_median": summarize(reduced[:, 2])["median"],
                "spacing_along_ratio_median": fixture_along,
                "spacing_cross_ratio_median": fixture_cross,
            },
        },
    )


def validate_swath(day, tle, entry, orbit, reference, write) -> None:
    """Scores the swath configuration, and writes its results if requested."""

    def satellite_for(time_step: float):
        return entry.to_satellite(orbit, SWATH, time_step=time_step, platform=PLATFORM)

    satellite = satellite_for(SCAN_PERIOD)
    points = generate_points_fibonacci_lattice(POINT_SPACING)
    observed = reference.observe_passes(points)
    margin = timedelta(seconds=SCAN_PERIOD)
    start = reference.time[0].floor("us").to_pydatetime() - margin
    end = reference.time[-1].floor("us").to_pydatetime() + margin
    predicted = predict_passes(points, satellite, 0, start, end, reference.gaps)
    results = match_passes(predicted, observed)
    results["swath_width_median"] = score_swath_width(
        reference, satellite.instruments[0], orbit
    )
    print(json.dumps(results), flush=True)
    grid = build_grid()
    first, last, compared, covered = observe_orbit(reference, orbit, grid)
    results["footprint_iou_median"] = score_swath_area(
        satellite_for, 0, first, last, TIME_STEPS, grid, compared, covered
    )
    print(json.dumps(results), flush=True)
    if not write:
        return
    rng = np.random.default_rng(1)
    picks = observed.iloc[rng.choice(len(observed), FIXTURE_PASSES, replace=False)]
    passes = np.stack(
        [
            points.geometry.y.to_numpy()[picks.point_id],
            points.geometry.x.to_numpy()[picks.point_id],
            (picks.time - pd.Timestamp(day)).dt.total_seconds().to_numpy(),
        ],
        axis=-1,
    )
    update_fixture(
        FIXTURE_PATH, tle=np.array(tle), day=np.array(day.isoformat()), passes=passes
    )
    gated = (
        "pass_recall",
        "pass_precision",
        "epoch_p95",
        "epoch_median",
        "swath_width_median",
        "footprint_iou_median",
    )
    update_record(
        RECORD_PATH,
        "iasi",
        SWATH,
        {
            "digest": entry.configurations[SWATH].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": "docs/validation/validate_catalog_iasi.py",
            "reference": "Pixel positions of the OSI SAF IASI sea surface temperature "
            "product",
            "period": f"{day:%Y-%m-%d}",
            "platform": PLATFORM,
            "holdout_gap_days": get_holdout_gap(day, FIT_DAY),
            "metrics": [
                {"name": name, "value": float(results[name])} for name in gated
            ],
            "fixture_metrics": recompute_passes(satellite, 0, day, passes),
        },
    )


def main() -> None:
    """Fits or validates the entry on a day of data."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("day")
    parser.add_argument("--fit", action="store_true")
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)

    day = datetime.fromisoformat(arguments.day).replace(tzinfo=timezone.utc)
    orbit = GeneralPerturbationsOrbit.from_tle(ELEMENT_SET)
    reference = Reference(load_day(day, arguments.directory))
    print(
        f"{arguments.day}: {len(reference.time)} scans, {reference.n_dropped} dropped, "
        f"{len(reference.gaps)} gaps over a minute, "
        f"median scan period {reference.scan_period:.4f} s",
        flush=True,
    )
    if arguments.fit:
        fit_geometry(reference, orbit, day)
        return
    entry = get_entry("iasi")
    validate_pixels(day, ELEMENT_SET, entry, orbit, reference, arguments.write)
    validate_swath(day, ELEMENT_SET, entry, orbit, reference, arguments.write)
    if arguments.write:
        print(f"wrote {RECORD_PATH}")


if __name__ == "__main__":
    main()

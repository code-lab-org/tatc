"""
Validates the MODIS catalog entry against the 1 km geolocation products of Terra and
Aqua, using the exact values of the entry through public TAT-C methods. The tests of
passes and pixels follow `ValidateDriftMODIS`, whose day of Terra data the entry's
angles were fitted to, and those of sizes and spacing follow `ValidateSampleSizes`:
running on that day checks this script against the notebooks, while running on a day
at least 30 days from it tests the fitted angles on data held out from the fit.

The notebooks do not test the swath's width or area, which the standard requires. The
observed width is the length through the frames of a scan, midway between its first
and last detectors, and the observed area is the union of the cells between those
lines over one orbit. Unlike the notebook, the reference time of a pass includes the
time within a scan of its frame, since MODIS sweeps the Earth in about 0.45 s.

The entry is fitted on Terra. On Aqua, results are reported as a transfer to a sibling
platform, under metric names that do not count toward the entry's validation.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`:

    python docs/validation/validate_catalog_modis.py <data directory> 2026-10-04
    python docs/validation/validate_catalog_modis.py <data directory> 2026-11-05
    python docs/validation/validate_catalog_modis.py <data directory> 2026-11-05 \
        --platform aqua

With `--write`, it writes the entry's validation record and, on Terra, the reduced
fixture from which unit tests recompute the results.

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

import tatc
from tatc.catalog import get_entry
from tatc.catalog.metrics import summarize
from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import GeneralPerturbationsOrbit
from validation_io import (
    TO_GEODETIC,
    RemoteFile,
    add_transfer_metrics,
    build_grid,
    compute_published_errors,
    compute_spacing_ratios,
    ecef,
    get_holdout_gap,
    in_lattice_cells,
    match_passes,
    normalize,
    predict_passes,
    recompute_passes,
    score_published,
    score_swath_area,
    up,
    update_fixture,
    update_record,
)

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
RECORD_PATH = REPOSITORY / "src" / "tatc" / "resources" / "validation" / "modis.json"
FIXTURE_PATH = REPOSITORY / "tests" / "catalog" / "fixtures" / "modis.npz"

FIT_DAY = datetime(2026, 10, 4, tzinfo=timezone.utc)
FIT_PLATFORM = "terra"
PRODUCTS = {"terra": "MOD03", "aqua": "MYD03"}
# element sets with epochs about two days after the day of data
ELEMENT_SETS = {
    ("terra", "2026-10-04"): [
        "1 25994U 99068A   26279.29978334  .00000220  00000+0  53719-4 0  9994",
        "2 25994  97.9333 324.3740 0001348 204.7527 222.4948 14.61170513425833",
    ],
    ("terra", "2026-11-05"): [
        "1 25994U 99068A   26311.28017572  .00000398  00000+0  89808-4 0  9992",
        "2 25994  97.9278 355.0236 0003362 122.2605 309.2823 14.61191265430504",
    ],
    ("aqua", "2026-11-05"): [
        "1 27424U 02022A   26311.28865478  .00000831  00000+0  17465-3 0  9993",
        "2 27424  98.4481 283.5818 0000983 104.9545   2.1913 14.62267061304185",
    ],
}
GEO = "HDFEOS/SWATHS/MODIS_Swath_Type_GEO/"
FRAMES = np.r_[np.arange(0, 1354, 8), 1353]
N_DETECTORS = 10
FRAME_PERIOD = 0.45 / 1354
CENTER_FRAME = 676.5
SCAN_PERIOD = 1.477
# seconds of International Atomic Time since 1993-01-01, when it led UTC by 27 s, not
# the current 37 s
TAI93_EPOCH = pd.Timestamp("1993-01-01", tz="UTC") - pd.Timedelta(seconds=10)
TIME_STEPS = (SCAN_PERIOD, 30, 120)
ORBIT_MARGIN = 120
PASS_RADIUS = 30e3
SWATH, PIXELS = "swath", "pixels"
# published footprints: instrument, altitude (km), inclination (deg), scan angle
# (deg), and sizes along and across track (km), from NASA LAADS
PUBLISHED_FOOTPRINTS = (
    ("MODIS", 705, 98.2, 0, 1.0, 1.0),
    ("MODIS", 705, 98.2, 55.0, 2.0, 4.8),
)
PIXEL_SCAN_STRIDE, SPACING_SCAN_STRIDE, WIDTH_SCAN_STRIDE = 200, 240, 10
FIXTURE_PASSES, FIXTURE_PIXELS, FIXTURE_SPACING = 60, 200, 200


def reduce_granule(url: str, session, directory: pathlib.Path) -> pathlib.Path:
    """
    Reduces a granule to its scan times, navigation, and first and last detectors of
    every scan at every eighth frame, or finds it in a cache.
    """
    path = directory / (pathlib.PurePosixPath(url).stem + ".npz")
    if path.exists():
        return path
    with h5py.File(RemoteFile(url, session, block_size=2**20), "r") as f:
        n_scans = len(f["EV center time"])
        rows = np.sort(
            np.r_[
                np.arange(n_scans) * N_DETECTORS,
                np.arange(n_scans) * N_DETECTORS + N_DETECTORS - 1,
            ]
        )

        def grid(name: str) -> np.ndarray:
            values = f[GEO + name][:][rows][:, FRAMES].astype(np.float32)
            values[values <= -999] = np.nan
            return values.reshape(n_scans, 2, len(FRAMES))

        arrays = {
            "time": f["EV center time"][:],
            "position": f["orb_pos"][:],
            "latitude": grid("Geolocation Fields/Latitude"),
            "longitude": grid("Geolocation Fields/Longitude"),
            "height": grid("Data Fields/Height"),
        }
    # written whole, so that an interrupted download leaves no partial granule
    temporary = path.with_suffix(".tmp.npz")
    np.savez_compressed(temporary, **arrays)
    temporary.replace(path)
    return path


def load_day(platform: str, day: datetime, directory: pathlib.Path) -> dict:
    """Downloads and reduces the geolocation of a day, or loads it from a cache."""
    product = PRODUCTS[platform]
    path = directory / f"{product}_{day:%Y%m%d}.npz"
    if not path.exists():
        if earthaccess.login(strategy="environment") is None:
            raise RuntimeError("Earthdata login failed: is EARTHDATA_TOKEN set?")
        session = earthaccess.get_requests_https_session()
        items = requests.get(
            "https://cmr.earthdata.nasa.gov/search/granules.umm_json",
            params={
                "short_name": product,
                "version": "7",
                "temporal": f"{day:%Y-%m-%d}T00:00:00Z,{day:%Y-%m-%d}T23:59:59Z",
                "page_size": 300,
            },
            timeout=120,
        ).json()["items"]
        urls = sorted(
            next(
                link["URL"]
                for link in item["umm"]["RelatedUrls"]
                if link["URL"].startswith("https://data.laads")
            )
            for item in items
        )
        granules = directory / f"{product}_{day:%Y%m%d}"
        granules.mkdir(parents=True, exist_ok=True)
        with concurrent.futures.ThreadPoolExecutor(6) as pool:
            paths = list(
                pool.map(lambda url: reduce_granule(url, session, granules), urls)
            )
        arrays: dict[str, list] = {}
        for granule in sorted(paths):
            with np.load(granule) as data:
                for key in data.files:
                    arrays.setdefault(key, []).append(data[key])
        np.savez_compressed(path, **{k: np.concatenate(v) for k, v in arrays.items()})
        print(f"{product}: {len(urls)} granules", flush=True)
    return dict(np.load(path))


class Reference:
    """A day of MODIS scans, with what is derived from them."""

    def __init__(self, data: dict[str, np.ndarray]):
        position = data["position"].astype(float)
        good = (
            (np.linalg.norm(position, axis=1) > 6e6)
            & (data["time"] > 0)
            & np.all(np.isfinite(data["latitude"]), axis=(1, 2))
            & np.all(np.isfinite(data["longitude"]), axis=(1, 2))
        )
        self.n_dropped = int(np.sum(~good))
        self.position = position[good]
        self.time = TAI93_EPOCH + pd.to_timedelta(data["time"][good], unit="s")
        self.lat = data["latitude"][good].astype(float)
        self.lon = data["longitude"][good].astype(float)
        self.height = np.nan_to_num(data["height"][good].astype(float))
        # midway between the first and last detectors, to which scan times refer
        self.center = (ecef(self.lon, self.lat).sum(axis=1)) / 2
        self.unit = normalize(up(self.lon, self.lat).sum(axis=1))
        self.seconds = (self.time - self.time[0]).total_seconds().to_numpy()
        period = np.diff(self.seconds)
        self.is_next_contiguous = period < 1.5 * SCAN_PERIOD
        self.gaps = [
            (self.time[i], self.time[i + 1]) for i in np.nonzero(period > 60)[0]
        ]
        self.scan_period = float(np.median(period))

    def in_cells(self, points: np.ndarray, first: int, last: int) -> np.ndarray:
        """Gets whether unit vectors lie in cells on a range of scans."""
        return in_lattice_cells(
            self.unit[first : last + 1], self.is_next_contiguous[first:last], points
        )

    def observe_passes(self, points) -> pd.DataFrame:
        """
        Finds the passes on which MODIS observed points, as `ValidateDriftMODIS` does:
        the nearest sampled pixel of the first detectors, excluding the outermost
        frames, within a distance. Gets their times, and those at the scan's center.
        """
        n_scans, n_frames = self.lat.shape[0], len(FRAMES)
        xyz = ecef(self.lon[:, 0], self.lat[:, 0], self.height[:, 0]).reshape(-1, 3)
        point_xyz = ecef(points.geometry.x.to_numpy(), points.geometry.y.to_numpy())
        rows = []
        tree = cKDTree(xyz)
        hits_of_points = tree.query_ball_point(point_xyz, PASS_RADIUS, workers=-1)
        for p, hits in enumerate(hits_of_points):
            if not hits:
                continue
            hits = np.array(hits)
            hits = hits[np.argsort(self.seconds[hits // n_frames])]
            seconds = self.seconds[hits // n_frames]
            for group in np.split(hits, np.nonzero(np.diff(seconds) > 600)[0] + 1):
                nearest = group[
                    np.argmin(np.linalg.norm(xyz[group] - point_xyz[p], axis=-1))
                ]
                scan, column = np.divmod(nearest, n_frames)
                if column in (0, n_frames - 1):
                    continue
                before = int(min(scan, n_scans - 2))
                a, b = self.center[before, column], self.center[before + 1, column]
                fraction = np.dot(point_xyz[p] - a, b - a) / np.dot(b - a, b - a)
                center = self.seconds[before] + fraction * (
                    self.seconds[before + 1] - self.seconds[before]
                )
                offset = (FRAMES[column] - CENTER_FRAME) * FRAME_PERIOD
                rows.append(
                    {
                        "point_id": p,
                        "time": self.time[0]
                        + timedelta(seconds=float(center + offset)),
                        "center_time": self.time[0] + timedelta(seconds=float(center)),
                    }
                )
        return pd.DataFrame(rows).astype({"point_id": "int64"})

    def contiguous_scans(self) -> np.ndarray:
        """Gets the scans whose neighbors are one scan period before and after."""
        period = np.diff(self.seconds)
        regular = np.abs(period - self.scan_period) < 0.01
        return np.nonzero(regular[:-1] & regular[1:])[0] + 1


def score_swath_width(reference: Reference, instrument) -> float:
    """Gets the median absolute difference (km) from the width through the frames."""
    scans = np.arange(0, len(reference.time), WIDTH_SCAN_STRIDE)
    center = reference.center[scans]
    observed = np.sum(np.linalg.norm(np.diff(center, axis=1), axis=-1), axis=1)
    altitude = TO_GEODETIC.transform(*reference.position[scans].T)[2]
    predicted = np.vectorize(instrument.get_swath_width)(altitude)
    print(
        f"swath width: observed {np.median(observed) / 1e3:.2f} km, "
        f"predicted {np.median(predicted) / 1e3:.2f} km",
        flush=True,
    )
    return float(np.median(np.abs(predicted - observed))) / 1e3


def observe_orbit(reference: Reference, orbit, grid):
    """Gets the cells compared over one orbit, and which of them MODIS observed."""
    lon, lat, _ = grid
    unit = up(lon, lat)
    # a day's granules can begin on the previous day
    begin = reference.time[len(reference.time) // 2].floor("D") + pd.Timedelta(hours=12)
    period = pd.Timedelta(days=1 / orbit.elements[0].mean_motion * 360 / 86400)
    scans = np.nonzero((reference.time >= begin) & (reference.time < begin + period))[0]
    first, last = reference.time[scans[0]], reference.time[scans[-1]]
    distance, nearest = cKDTree(reference.unit[scans].reshape(-1, 3)).query(
        unit, workers=-1
    )
    elapsed = (reference.time[scans] - first).total_seconds().to_numpy()
    nearest_time = elapsed[nearest // len(FRAMES)]
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
    Gets rows of seconds from the start of the day, frame, along-track pixel, terrain
    height (m), position difference (km), and the pixel's latitude and longitude, for
    the first and last detectors of a sample of scans, at the center times of scans.
    """
    scans = np.arange(0, len(reference.time), PIXEL_SCAN_STRIDE)
    instrument = satellite.get_instrument(0)
    times = reference.time[scans]
    orbit = satellite.orbit.to_gp_orbit()
    track = orbit.get_orbit_track(list(times))
    tracks = [orbit.get_orbit_track(time) for time in times]
    seconds = (times - pd.Timestamp(day)).total_seconds().to_numpy()
    rows = []
    for j, frame in enumerate(FRAMES):
        # TAT-C numbers pixels from the right of the motion, as MODIS numbers frames,
        # and along track from the most forward, which is the last detector
        for row, pixel in ((0, N_DETECTORS - 1), (1, 0)):
            height = reference.height[scans, row, j]
            lat, lon = reference.lat[scans, row, j], reference.lon[scans, row, j]
            position = instrument.compute_projected_pixel_position(
                track, int(frame), pixel
            )
            predicted_lon = position.longitude.degrees
            predicted_lat = position.latitude.degrees
            # elevation is a single value, so pixels on terrain are projected singly
            for k in np.nonzero(height)[0]:
                position = instrument.compute_projected_pixel_position(
                    tracks[k], int(frame), pixel, elevation=float(height[k])
                )
                predicted_lon[k] = position.longitude.degrees
                predicted_lat[k] = position.latitude.degrees
            difference = np.linalg.norm(
                ecef(predicted_lon, predicted_lat) - ecef(lon, lat),
                axis=-1,
            )
            rows.append(
                np.stack(
                    [
                        seconds,
                        np.full(len(scans), frame),
                        np.full(len(scans), pixel),
                        height,
                        difference / 1e3,
                        lat,
                        lon,
                    ],
                    axis=-1,
                )
            )
    return np.concatenate(rows)


def measure_spacing(reference: Reference) -> tuple[pd.DatetimeIndex, np.ndarray]:
    """
    Measures the spacing of the first detectors of scans as `ValidateSampleSizes` does:
    across the scan, from half the distance between the frames on either side, and
    along track, the detector pitch or, if larger, the advance per detector. Gets the
    times, and rows of latitude, longitude, along-track and cross-track spacing (km).
    """
    scans = reference.contiguous_scans()[::SPACING_SCAN_STRIDE]
    xyz = ecef(reference.lon, reference.lat)
    first, last, following = xyz[scans, 0], xyz[scans, 1], xyz[scans + 1, 0]
    pitch = np.linalg.norm(last - first, axis=-1) / (N_DETECTORS - 1) / 1e3
    advance = np.linalg.norm(following - first, axis=-1) / N_DETECTORS / 1e3
    cross = np.full(pitch.shape, np.nan)
    cross[:, 1:-1] = (
        np.linalg.norm(first[:, 2:] - first[:, :-2], axis=-1)
        / (FRAMES[2:] - FRAMES[:-2])
        / 1e3
    )
    rows = np.stack(
        [
            reference.lat[scans, 0].ravel(),
            reference.lon[scans, 0].ravel(),
            np.fmax(pitch, advance).ravel(),
            cross.ravel(),
        ],
        axis=-1,
    )
    times = pd.DatetimeIndex(np.repeat(reference.time[scans], len(FRAMES)))
    keep = np.isfinite(rows[:, 3])
    return times[keep], rows[keep]


def validate_pixels(day, tle, entry, orbit, reference, platform, write) -> None:
    """Scores the sample size configuration, and writes its results if requested."""
    within, largest = score_published(
        compute_published_errors(entry, PIXELS, PUBLISHED_FOOTPRINTS)
    )
    print(
        f"published dimensions: {within:.1f}% within 3%, largest error {largest:.2f}%"
    )
    satellite = entry.to_satellite(orbit, PIXELS, platform=platform)
    pixels = measure_pixels(reference, satellite, day)
    median = summarize(pixels[:, 4])["median"]
    p95 = summarize(pixels[:, 4])["p95"]
    print(f"pixel position: median {median:.3f} km, 95th percentile {p95:.3f} km")
    for frame in (0, 680, 1353):
        chosen = pixels[:, 1] == min(FRAMES, key=lambda f: abs(f - frame))
        print(f"  frame {frame}: median {np.median(pixels[chosen, 4]):.2f} km")
    times, rows = measure_spacing(reference)
    names = ["MODIS"] * len(rows)
    along, cross = compute_spacing_ratios(satellite, names, times, rows)
    print(f"spacing, predicted / measured: {along:.4f} along, {cross:.4f} across")
    if not write:
        return
    if platform != FIT_PLATFORM:
        add_transfer_metrics(
            RECORD_PATH,
            PIXELS,
            f"{platform}, {day:%Y-%m-%d}",
            {
                "pixel_position_median": median,
                "pixel_position_p95": p95,
                "spacing_along_ratio_median": along,
                "spacing_cross_ratio_median": cross,
            },
        )
        return
    rng = np.random.default_rng(3)
    picks = rng.choice(len(rows), FIXTURE_SPACING, replace=False)
    fixture_along, fixture_cross = compute_spacing_ratios(
        satellite, ["MODIS"] * len(picks), times[picks], rows[picks]
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
        "modis",
        PIXELS,
        {
            "digest": entry.configurations[PIXELS].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": "docs/validation/validate_catalog_modis.py",
            "reference": "Footprints published by NASA LAADS; MOD03 geolocation",
            "period": f"{day:%Y-%m-%d}",
            "platform": platform,
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
                "pixel_position_median": summarize(reduced[:, 4])["median"],
                "spacing_along_ratio_median": fixture_along,
                "spacing_cross_ratio_median": fixture_cross,
            },
        },
    )


def validate_swath(day, tle, entry, orbit, reference, platform, write) -> None:
    """Scores the swath configuration, and writes its results if requested."""

    def satellite_for(time_step: float):
        return entry.to_satellite(orbit, SWATH, time_step=time_step, platform=platform)

    satellite = satellite_for(SCAN_PERIOD)
    points = generate_points_fibonacci_lattice(500e3)
    observed = reference.observe_passes(points)
    start = reference.time[0].floor("us").to_pydatetime() - timedelta(
        seconds=SCAN_PERIOD
    )
    end = reference.time[-1].floor("us").to_pydatetime() + timedelta(
        seconds=SCAN_PERIOD
    )
    predicted = predict_passes(points, satellite, 0, start, end, reference.gaps)
    results = match_passes(predicted, observed)
    at_center = match_passes(
        predicted, observed.drop(columns="time").rename(columns={"center_time": "time"})
    )
    print(
        "at scan centers, as in the notebook: "
        f"median {at_center['epoch_median']:.3f} s, p95 {at_center['epoch_p95']:.3f} s",
        flush=True,
    )
    results["swath_width_median"] = score_swath_width(
        reference, satellite.instruments[0]
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
    gated = (
        "pass_recall",
        "pass_precision",
        "epoch_p95",
        "swath_width_median",
        "footprint_iou_median",
    )
    if platform != FIT_PLATFORM:
        add_transfer_metrics(
            RECORD_PATH,
            SWATH,
            f"{platform}, {day:%Y-%m-%d}",
            {name: results[name] for name in gated},
        )
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
    update_record(
        RECORD_PATH,
        "modis",
        SWATH,
        {
            "digest": entry.configurations[SWATH].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": "docs/validation/validate_catalog_modis.py",
            "reference": "MOD03 geolocation",
            "period": f"{day:%Y-%m-%d}",
            "platform": platform,
            "holdout_gap_days": get_holdout_gap(day, FIT_DAY),
            "metrics": [
                {"name": name, "value": float(results[name])} for name in gated
            ],
            "fixture_metrics": recompute_passes(satellite, 0, day, passes),
        },
    )


def main() -> None:
    """Validates the entry on a day of data."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("day")
    parser.add_argument("--platform", choices=sorted(PRODUCTS), default=FIT_PLATFORM)
    parser.add_argument("--sample-sizes-only", action="store_true")
    parser.add_argument("--download-only", action="store_true")
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)

    day = datetime.fromisoformat(arguments.day).replace(tzinfo=timezone.utc)
    tle = ELEMENT_SETS[(arguments.platform, arguments.day)]
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    entry = get_entry("modis")
    data = load_day(arguments.platform, day, arguments.directory)
    if arguments.download_only:
        return
    reference = Reference(data)
    print(
        f"{arguments.platform}, {arguments.day}: {tle[0][:32]}; "
        f"{len(reference.time)} scans, {reference.n_dropped} dropped, "
        f"{len(reference.gaps)} gaps over a minute, "
        f"median scan period {reference.scan_period:.4f} s",
        flush=True,
    )
    validate_pixels(
        day, tle, entry, orbit, reference, arguments.platform, arguments.write
    )
    if not arguments.sample_sizes_only:
        validate_swath(
            day, tle, entry, orbit, reference, arguments.platform, arguments.write
        )
    if arguments.write:
        print(f"wrote {RECORD_PATH}")


if __name__ == "__main__":
    main()

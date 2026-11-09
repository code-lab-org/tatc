"""
Helpers shared by the scripts that validate cross-track scanners against the swath
products of the GPM mission at NASA GES DISC (the Level 2 DPR product and the Level
1C MHS products), whose scans give the position of each beam's center on the
ellipsoid: reading the products, fitting the scan's angles, and scoring.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import concurrent.futures
import json
import pathlib
from dataclasses import dataclass
from datetime import datetime, timedelta

import earthaccess
import h5py
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from skyfield.api import wgs84

import tatc
from tatc.catalog.metrics import summarize
from tatc.generation import generate_points_fibonacci_lattice
from validation_io import (
    TO_GEODETIC,
    RemoteFile,
    build_grid,
    compute_published_errors,
    compute_spacing_ratios,
    ecef,
    get_holdout_gap,
    get_with_retries,
    in_lattice_cells,
    normalize,
    recompute_passes,
    score_passes,
    score_published,
    score_swath_area,
    up,
    update_fixture,
    update_record,
)

EARTH_ROTATION = 7.2921150e-5  # rad/s
ORBIT_MARGIN = 120  # s


@dataclass(frozen=True)
class Product:
    """Where a product stores the fields of a scanner's swath."""

    short_name: str
    version: str
    group: str
    spacecraft: tuple[str, str, str]
    quality: str | None = None


def load_day(product: Product, day: datetime, directory: pathlib.Path) -> dict:
    """
    Reads the scans of a day from the granules overlapping it, with range requests
    for only the needed fields, or loads them from a cache. Keeps each scan once, in
    time order, within the day.
    """
    cache = directory / f"{product.short_name.lower()}_{day:%Y%m%d}.npz"
    if cache.exists():
        return dict(np.load(cache))
    if earthaccess.login(strategy="environment") is None:
        raise RuntimeError("Earthdata login failed: is EARTHDATA_TOKEN set?")
    session = earthaccess.get_requests_https_session()
    items = get_with_retries(
        "https://cmr.earthdata.nasa.gov/search/granules.umm_json",
        params={
            "short_name": product.short_name,
            "version": product.version,
            "temporal": f"{day - timedelta(hours=2):%Y-%m-%dT%H:%M:%SZ},"
            f"{day + timedelta(days=1, hours=2):%Y-%m-%dT%H:%M:%SZ}",
            "page_size": 100,
        },
    ).json()["items"]
    urls = [
        next(
            u["URL"]
            for u in item["umm"]["RelatedUrls"]
            if u["URL"].startswith("https://data.gesdisc")
            and u["URL"].endswith((".nc", ".HDF5", ".h5"))
        )
        for item in items
    ]
    group = product.group

    def read(url):
        with h5py.File(RemoteFile(url, session)) as f:
            g = f[group]
            date = pd.to_datetime(
                pd.DataFrame(
                    {
                        "year": g["ScanTime/Year"][:],
                        "month": g["ScanTime/Month"][:],
                        "day": g["ScanTime/DayOfMonth"][:],
                    }
                ),
                utc=True,
            )
            seconds = (date - pd.Timestamp(day)).dt.total_seconds().to_numpy() + g[
                "ScanTime/SecondOfDay"
            ][:].astype(float)
            part = {
                "time": seconds,
                "lat": g["Latitude"][:],
                "lon": g["Longitude"][:],
                "sc_lat": g[product.spacecraft[0]][:],
                "sc_lon": g[product.spacecraft[1]][:],
                "sc_alt": g[product.spacecraft[2]][:],
            }
            if product.quality is not None:
                quality = g[product.quality][:]
                part["quality"] = quality if quality.ndim == 1 else quality[:, 0]
            return part

    with concurrent.futures.ThreadPoolExecutor(4) as pool:
        parts = list(pool.map(read, urls))
    reduced = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    _, first = np.unique(np.round(reduced["time"], 3), return_index=True)
    keep = first[(reduced["time"][first] >= 0) & (reduced["time"][first] < 86400)]
    np.savez_compressed(
        cache, granules=len(urls), **{k: v[keep] for k, v in reduced.items()}
    )
    return dict(np.load(cache))


class Reference:
    """A day of scans of a cross-track scanner, with what is derived from them."""

    def __init__(self, day: datetime, data: dict, scan_period: float):
        self.day = day
        valid = np.all(np.abs(data["lat"]) <= 90, axis=1) & np.all(
            np.abs(data["lon"]) <= 180, axis=1
        )
        if "quality" in data:
            valid &= data["quality"] == 0
        altitude = data["sc_alt"][valid].astype(float)
        # altitudes are in km in some products and m in others
        self.altitude = altitude * 1e3 if np.median(altitude) < 1e4 else altitude
        self.seconds = data["time"][valid].astype(float)
        self.time = pd.DatetimeIndex(
            [pd.Timestamp(day) + pd.Timedelta(seconds=s) for s in self.seconds]
        )
        self.lat = data["lat"][valid].astype(float)
        self.lon = data["lon"][valid].astype(float)
        self.position = ecef(
            data["sc_lon"][valid].astype(float),
            data["sc_lat"][valid].astype(float),
            self.altitude,
        )
        self.scan_period = scan_period
        self.beams = self.lat.shape[1]
        self.unit = up(self.lon, self.lat)
        period = np.diff(self.seconds)
        self.is_next_contiguous = period < 1.5 * scan_period
        self.gaps = [
            (self.time[i].to_pydatetime(), self.time[i + 1].to_pydatetime())
            for i in np.nonzero(period > 60)[0]
        ]

    def velocity(self, frame: str) -> np.ndarray:
        """Gets the spacecraft's velocity (m/s), Earth-fixed or inertial, by scan."""
        velocity = np.gradient(self.position, self.seconds, axis=0)
        if frame == "inertial":
            velocity = velocity + np.cross(
                [0, 0, EARTH_ROTATION], self.position, axis=-1
            )
        return velocity

    def cross_track_angles(self, frame: str) -> np.ndarray:
        """
        Gets each beam's angle (deg) from the geodetic nadir, across track (positive
        to the left of the motion in the frame), by scan and beam.
        """
        lon, lat, _ = TO_GEODETIC.transform(*self.position.T)
        nadir = -up(lon, lat)
        left = normalize(np.cross(self.velocity(frame), nadir))
        sight = normalize(ecef(self.lon, self.lat) - self.position[:, np.newaxis, :])
        return np.degrees(
            np.arctan2(
                np.sum(sight * left[:, np.newaxis, :], axis=-1),
                np.sum(sight * nadir[:, np.newaxis, :], axis=-1),
            )
        )

    def in_cells(self, points: np.ndarray, first: int, last: int) -> np.ndarray:
        """Gets whether unit vectors lie in cells of beams on a range of scans."""
        return in_lattice_cells(
            self.unit[first : last + 1], self.is_next_contiguous[first:last], points
        )

    def observe_passes(self, points, radius: float) -> pd.DataFrame:
        """
        Finds the passes on which the scanner observed points, and their times, with
        each beam at the time of its scan.
        """
        point_unit = up(points.geometry.x.to_numpy(), points.geometry.y.to_numpy())
        flat = self.unit.reshape(-1, 3)
        beam_seconds = np.repeat(self.seconds, self.beams)
        tree = cKDTree(flat)
        rows = []
        for p, hits in enumerate(
            tree.query_ball_point(
                point_unit, r=2 * np.sin(radius / 6371 / 2), workers=-1
            )
        ):
            if not hits:
                continue
            hits = np.array(hits)
            hits = hits[np.argsort(beam_seconds[hits])]
            for group in np.split(
                hits, np.nonzero(np.diff(beam_seconds[hits]) > 1200)[0] + 1
            ):
                k, j = np.divmod(
                    group[
                        np.argmin(np.linalg.norm(flat[group] - point_unit[p], axis=-1))
                    ],
                    self.beams,
                )
                first, last = max(k - 1, 0), min(k + 1, len(self.time) - 1)
                if not self.in_cells(point_unit[[p]], first, last)[0]:
                    continue
                # interpolate between the scans before and after the point
                before = (
                    k - 1
                    if np.sum(
                        (point_unit[p] - self.unit[k, j])
                        * (self.unit[k, j] - self.unit[max(k - 1, 0), j])
                    )
                    < 0
                    else k
                )
                before = int(np.clip(before, 0, len(self.time) - 2))
                a, b = self.unit[before, j], self.unit[before + 1, j]
                fraction = np.clip(
                    np.sum((point_unit[p] - a) * (b - a)) / np.sum((b - a) ** 2), 0, 1
                )
                seconds = self.seconds[before] + fraction * (
                    self.seconds[before + 1] - self.seconds[before]
                )
                rows.append(
                    {
                        "point_id": p,
                        "time": pd.Timestamp(self.day)
                        + pd.Timedelta(seconds=float(seconds)),
                    }
                )
        return pd.DataFrame(rows).astype({"point_id": "int64"})

    def contiguous_scans(self) -> np.ndarray:
        """Gets the scans whose neighbors are one scan period before and after."""
        regular = np.abs(np.diff(self.seconds) - self.scan_period) < 0.05
        return np.nonzero(regular[:-1] & regular[1:])[0] + 1


def fit(reference: Reference, frame: str) -> None:
    """
    Prints the fitted angles of the scan: the outermost beams' angles on either side,
    the field of regard between them (centered on nadir), the spacing and offset of a
    line fitted to the beams' angles by index, and the scan period.
    """
    angles = reference.cross_track_angles(frame)
    # orient each scan so that its angles increase with the beam index (from right
    # to left), as the side of the first beam changes if the spacecraft flies backward
    flipped = angles[:, 0] > angles[:, -1]
    angles[flipped] = angles[flipped, ::-1]
    median = np.median(angles, axis=0)
    index = np.arange(reference.beams) - (reference.beams - 1) / 2
    spacing, offset = np.polyfit(index, median, 1)
    right, left = -median[0], median[-1]
    print(
        f"{frame}: outermost beams {left:.4f} deg left, {right:.4f} deg right; "
        f"field of regard {left + right:.4f} deg; "
        f"spacing {spacing:.5f} deg, roll {offset:.4f} deg; "
        f"residual p95 {np.percentile(np.abs(median - (spacing * index + offset)), 95):.4f} deg; "
        f"first beam on the left in {np.mean(flipped):.1%} of scans"
    )
    period = np.diff(reference.seconds)
    print(f"scan period: median {np.median(period):.6f} s")


def score_swath_width(reference: Reference, instrument, orbit, stride: int) -> float:
    """Gets the median absolute difference (km) from the width through the beams."""
    scans = np.arange(0, len(reference.time), stride)
    xyz = ecef(reference.lon[scans], reference.lat[scans])
    observed = np.sum(np.linalg.norm(np.diff(xyz, axis=1), axis=-1), axis=1)
    altitude = wgs84.geographic_position_of(
        orbit.get_orbit_track(list(reference.time[scans].to_pydatetime()))
    ).elevation.m
    predicted = np.vectorize(instrument.get_swath_width)(altitude)
    return float(np.median(np.abs(predicted - observed))) / 1e3


def observe_orbit(reference: Reference, orbit, grid, radius: float):
    """Gets the cells compared over one orbit, and which of them were observed."""
    lon, lat, _ = grid
    unit = up(lon, lat)
    begin = pd.Timestamp(reference.day) + pd.Timedelta(hours=12)
    period = pd.Timedelta(days=1 / orbit.elements[0].mean_motion)
    # the first orbit from noon without a gap of more than a minute
    while True:
        scans = np.nonzero(
            (reference.time >= begin) & (reference.time < begin + period)
        )[0]
        if len(scans) and np.all(np.diff(reference.seconds[scans]) < 60):
            break
        begin += pd.Timedelta(minutes=10)
        if begin > pd.Timestamp(reference.day) + pd.Timedelta(hours=22):
            raise RuntimeError("no orbit without a gap")
    first, last = reference.time[scans[0]], reference.time[scans[-1]]
    distance, nearest = cKDTree(reference.unit[scans].reshape(-1, 3)).query(
        unit, workers=-1
    )
    elapsed = reference.seconds[scans] - reference.seconds[scans[0]]
    nearest_time = elapsed[nearest // reference.beams]
    compared = (
        (distance < 2 * np.sin(radius / 6371 / 2))
        & (nearest_time > ORBIT_MARGIN)
        & (nearest_time < (last - first).total_seconds() - ORBIT_MARGIN)
    )
    observed = np.zeros(len(unit), dtype=bool)
    observed[compared] = reference.in_cells(unit[compared], scans[0], scans[-1])
    return first.to_pydatetime(), last.to_pydatetime(), compared, observed


def measure_pixels(reference: Reference, satellite, index: int, stride: int, frame):
    """
    Gets rows of seconds from the start of the day, pixel index, position difference
    (km) from the beam center, and the beam's latitude and longitude, for the pixels
    projected at the times of a sample of scans. TAT-C numbers pixels from the right
    of the motion, so each scan's beams are ordered by their side, which changes when
    the spacecraft flies backward.
    """
    scans = np.arange(0, len(reference.time), stride)
    angles = reference.cross_track_angles(frame)[scans]
    # beam 0 is on the right when the angles increase with the index
    increasing = angles[:, 0] < angles[:, -1]
    instrument = satellite.get_instrument(index)
    orbit = satellite.orbit.to_gp_orbit()
    rows = []
    for beam in range(reference.beams):
        pixel = np.where(increasing, beam, reference.beams - 1 - beam)
        for value in np.unique(pixel):
            subset = np.nonzero(pixel == value)[0]
            track = orbit.get_orbit_track(
                list(reference.time[scans[subset]].to_pydatetime())
            )
            position = instrument.compute_projected_pixel_position(track, int(value), 0)
            difference = np.linalg.norm(
                np.atleast_2d(
                    ecef(position.longitude.degrees, position.latitude.degrees)
                )
                - ecef(
                    reference.lon[scans[subset], beam],
                    reference.lat[scans[subset], beam],
                ),
                axis=-1,
            )
            rows.append(
                np.stack(
                    [
                        reference.seconds[scans[subset]],
                        np.full(len(subset), value),
                        difference / 1e3,
                        reference.lat[scans[subset], beam],
                        reference.lon[scans[subset], beam],
                    ],
                    axis=-1,
                )
            )
    return np.concatenate(rows)


def measure_spacing(reference: Reference, stride: int):
    """
    Measures the spacing of beams as `ValidateSampleSizes` does, from half the distance
    between the beams on either side. Gets the times, and rows of latitude, longitude,
    along-track and cross-track spacing (km).
    """
    scans = reference.contiguous_scans()[::stride]
    xyz = ecef(reference.lon, reference.lat)
    along = np.linalg.norm(xyz[scans + 1] - xyz[scans - 1], axis=-1) / 2e3
    cross = np.full(along.shape, np.nan)
    cross[:, 1:-1] = np.linalg.norm(xyz[scans, 2:] - xyz[scans, :-2], axis=-1) / 2e3
    times = pd.DatetimeIndex(np.repeat(reference.time[scans], reference.beams))
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


@dataclass(frozen=True)
class Setup:
    """What a script validates, and where it writes its results."""

    entry_id: str
    platform: str
    swath: str
    beams: str
    size_instrument: str
    published_footprints: tuple
    record_path: pathlib.Path
    fixture_path: pathlib.Path
    source: str
    reference: str
    published_reference: str
    fit_day: datetime
    time_steps: tuple[float, ...]
    pass_radius: float
    point_spacing: float = 500e3
    pixel_stride: int = 100
    spacing_stride: int = 150
    width_stride: int = 5
    fixture_passes: int = 60
    fixture_pixels: int = 200
    fixture_spacing: int = 200


def validate_beams(setup: Setup, day, tle, entry, orbit, reference, write: bool):
    """Scores the sample size configuration, and writes its results if requested."""
    within, largest = score_published(
        compute_published_errors(entry, setup.beams, setup.published_footprints)
    )
    print(
        f"published dimensions: {within:.1f}% within 3%, largest error {largest:.2f}%"
    )
    satellite = entry.to_satellite(orbit, setup.beams, platform=setup.platform)
    names = [instrument.name for instrument in satellite.instruments]
    index = names.index(setup.size_instrument)
    pixels = measure_pixels(
        reference,
        satellite,
        index,
        setup.pixel_stride,
        entry.configurations[setup.beams].velocity_frame.value,
    )
    median = summarize(pixels[:, 2])["median"]
    p95 = summarize(pixels[:, 2])["p95"]
    print(f"pixel position: median {median:.3f} km, 95th percentile {p95:.3f} km")
    times, rows = measure_spacing(reference, setup.spacing_stride)
    along, cross = compute_spacing_ratios(
        satellite, [setup.size_instrument] * len(rows), times, rows
    )
    print(f"spacing, predicted / measured: {along:.4f} along, {cross:.4f} across")
    if not write:
        return
    rng = np.random.default_rng(3)
    picks = rng.choice(len(rows), setup.fixture_spacing, replace=False)
    fixture_along, fixture_cross = compute_spacing_ratios(
        satellite, [setup.size_instrument] * len(picks), times[picks], rows[picks]
    )
    reduced = pixels[rng.choice(len(pixels), setup.fixture_pixels, replace=False)]
    update_fixture(
        setup.fixture_path,
        tle=np.array(tle),
        day=np.array(day.isoformat()),
        pixels=reduced,
        spacing=rows[picks],
        spacing_seconds=(times[picks] - pd.Timestamp(day)).total_seconds().to_numpy(),
    )
    update_record(
        setup.record_path,
        setup.entry_id,
        setup.beams,
        {
            "digest": entry.configurations[setup.beams].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": setup.source,
            "reference": setup.published_reference,
            "period": f"{day:%Y-%m-%d}",
            "platform": setup.platform,
            "holdout_gap_days": get_holdout_gap(day, setup.fit_day),
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


def validate_swath(setup: Setup, day, tle, entry, orbit, reference, write: bool):
    """Scores the swath configuration, and writes its results if requested."""

    def satellite_for(time_step: float):
        return entry.to_satellite(
            orbit, setup.swath, time_step=time_step, platform=setup.platform
        )

    satellite = satellite_for(setup.time_steps[0])
    points = generate_points_fibonacci_lattice(setup.point_spacing)
    observed = reference.observe_passes(points, setup.pass_radius)
    start = reference.time[0].to_pydatetime()
    end = reference.time[-1].to_pydatetime() + timedelta(seconds=reference.scan_period)
    results = score_passes(points, satellite, 0, observed, start, end, reference.gaps)
    results["swath_width_median"] = score_swath_width(
        reference, satellite.instruments[0], orbit, setup.width_stride
    )
    print(json.dumps(results), flush=True)
    grid = build_grid()
    first, last, compared, covered = observe_orbit(
        reference, orbit, grid, 2 * setup.pass_radius
    )
    results["footprint_iou_median"] = score_swath_area(
        satellite_for, 0, first, last, setup.time_steps, grid, compared, covered
    )
    print(json.dumps(results), flush=True)
    if not write:
        return
    rng = np.random.default_rng(1)
    picks = observed.iloc[
        rng.choice(len(observed), setup.fixture_passes, replace=False)
    ]
    passes = np.stack(
        [
            points.geometry.y.to_numpy()[picks.point_id],
            points.geometry.x.to_numpy()[picks.point_id],
            (picks.time - pd.Timestamp(day)).dt.total_seconds().to_numpy(),
        ],
        axis=-1,
    )
    update_fixture(
        setup.fixture_path,
        tle=np.array(tle),
        day=np.array(day.isoformat()),
        passes=passes,
    )
    update_record(
        setup.record_path,
        setup.entry_id,
        setup.swath,
        {
            "digest": entry.configurations[setup.swath].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": setup.source,
            "reference": setup.reference,
            "period": f"{day:%Y-%m-%d}",
            "platform": setup.platform,
            "holdout_gap_days": get_holdout_gap(day, setup.fit_day),
            "metrics": [
                {"name": name, "value": float(results[name])}
                for name in (
                    "pass_recall",
                    "pass_precision",
                    "epoch_p95",
                    "swath_width_median",
                    "footprint_iou_median",
                )
            ],
            "fixture_metrics": recompute_passes(satellite, 0, day, passes),
        },
    )

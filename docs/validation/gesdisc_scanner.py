"""
Helpers shared by the scripts that validate cross-track scanners against swath
products at NASA GES DISC (the GPM mission's Level 2 DPR product and Level 1C MHS
products, and the Level 1B AIRS and AMSU-A products of Aqua), whose scans give the
position of each beam's center on the ellipsoid: reading the products, fitting the
scan's angles, and scoring.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import concurrent.futures
import json
import pathlib
from collections.abc import Callable
from dataclasses import dataclass
from datetime import datetime, timedelta

import earthaccess
import h5py
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from pyproj import Geod
from skyfield.api import wgs84

import tatc
from tatc.catalog.metrics import summarize
from tatc.generation import generate_points_fibonacci_lattice
from tatc.utils import compute_projected_ray_position
from validation_io import (
    add_transfer_metrics,
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
GEOD = Geod(ellps="WGS84")
GATED = (
    "pass_recall",
    "pass_precision",
    "epoch_p95",
    "swath_width_median",
    "footprint_iou_median",
)
TAI93_EPOCH = pd.Timestamp("1993-01-01", tz="UTC")
# TAI - UTC was 27 s at the start of 1993, and has been 37 s since 2017
TAI93_LEAP_SECONDS = 10


@dataclass(frozen=True)
class Product:
    """
    Where a product stores the fields of a scanner's swath. Products of the GPM
    mission are read from their group's fields; others give a function that reads a
    granule, from its URL, an Earthdata session, the day, and a directory.
    """

    short_name: str
    version: str
    group: str = ""
    spacecraft: tuple[str, str, str] = ("", "", "")
    quality: str | None = None
    read: Callable[..., dict] | None = None
    host: str = "https://data.gesdisc"


def read_gpm(product: Product, url: str, session, day: datetime, _directory) -> dict:
    """
    Reads the fields of a granule of a GPM product, with range requests: each scan's
    time (seconds from the start of the day), the beams' positions, the spacecraft's
    position, and any quality flag.
    """
    with h5py.File(RemoteFile(url, session)) as f:
        g = f[product.group]
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


def tai93_to_seconds(tai93, day: datetime) -> np.ndarray:
    """Converts TAI seconds since 1993 to UTC seconds from the start of a day."""
    return (
        np.asarray(tai93, dtype=float)
        - TAI93_LEAP_SECONDS
        - (pd.Timestamp(day) - TAI93_EPOCH).total_seconds()
    )


def find_granules(product: Product, day: datetime) -> list[str]:
    """Gets the URLs of a product's granules from two hours before to after a day."""
    items = get_with_retries(
        "https://cmr.earthdata.nasa.gov/search/granules.umm_json",
        params={
            "short_name": product.short_name,
            "version": product.version,
            "temporal": f"{day - timedelta(hours=2):%Y-%m-%dT%H:%M:%SZ},"
            f"{day + timedelta(days=1, hours=2):%Y-%m-%dT%H:%M:%SZ}",
            "page_size": 2000,
        },
    ).json()["items"]
    return [
        next(
            u["URL"]
            for u in item["umm"]["RelatedUrls"]
            if u["URL"].startswith(product.host)
            and u["URL"].endswith((".nc", ".HDF5", ".h5", ".hdf", ".zip"))
        )
        for item in items
    ]


def load_day(product: Product, day: datetime, directory: pathlib.Path) -> dict:
    """
    Reads the scans of a day from the granules overlapping it, with range requests
    for only the needed fields where the format allows, or loads them from a cache.
    Keeps each scan once, in time order, within the day.
    """
    cache = directory / f"{product.short_name.lower()}_{day:%Y%m%d}.npz"
    if cache.exists():
        return dict(np.load(cache))
    if earthaccess.login(strategy="environment") is None:
        raise RuntimeError("Earthdata login failed: is EARTHDATA_TOKEN set?")
    session = earthaccess.get_requests_https_session()
    urls = find_granules(product, day)
    read = product.read or read_gpm

    def read_granule(url):
        return read(product, url, session, day, directory)

    with concurrent.futures.ThreadPoolExecutor(4) as pool:
        parts = list(pool.map(read_granule, urls))
    reduced = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    _, first = np.unique(np.round(reduced["time"], 3), return_index=True)
    keep = first[(reduced["time"][first] >= 0) & (reduced["time"][first] < 86400)]
    np.savez_compressed(
        cache, granules=len(urls), **{k: v[keep] for k, v in reduced.items()}
    )
    return dict(np.load(cache))


class Reference:
    """
    A day of scans of a cross-track scanner, with what is derived from them. Each scan
    has a time; products that also give each beam's time (`beam_time`, in seconds from
    the start of the day) time the passes by beam. A day reduced to some of the
    product's lines and samples gives the samples' indices in the product (`pixel`),
    their number (`pixels`), and the lines from one kept scan to the next
    (`line_stride`).
    """

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
        self.beam_seconds = (
            data["beam_time"][valid].astype(float)
            if "beam_time" in data
            else np.repeat(self.seconds[:, np.newaxis], data["lat"].shape[1], axis=1)
        )
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
        self.pixel = np.asarray(data.get("pixel", np.arange(self.beams)))
        self.pixels = int(data.get("pixels", self.beams))
        self.line_stride = int(data.get("line_stride", 1))
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
        each beam at its time (or the time of its scan, if the product gives no time
        by beam).
        """
        point_unit = up(points.geometry.x.to_numpy(), points.geometry.y.to_numpy())
        flat = self.unit.reshape(-1, 3)
        beam_seconds = self.beam_seconds.ravel()
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
                seconds = self.beam_seconds[before, j] + fraction * (
                    self.beam_seconds[before + 1, j] - self.beam_seconds[before, j]
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
    index = reference.pixel - (reference.pixels - 1) / 2
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


def fit_yaw(reference: Reference, frame: str) -> None:
    """
    Prints the yaw of the scan line: the angle whose tangent is the slope of the beams'
    along-track tangents on their cross-track tangents (positive if the left of the
    scan lies ahead of the right, as a yaw of the body to the right, a negative
    `yaw_angle`, gives), with the along-track angle at nadir, overall and by latitude
    and direction.
    """
    lon, lat, _ = TO_GEODETIC.transform(*reference.position.T)
    zenith = up(lon, lat)
    motion = reference.velocity(frame)
    along_axis = normalize(motion - np.sum(motion * zenith, axis=-1)[:, None] * zenith)
    left_axis = np.cross(zenith, along_axis)
    rays = ecef(reference.lon, reference.lat) - reference.position[:, np.newaxis]
    down = -np.sum(rays * zenith[:, np.newaxis], axis=-1)
    along = np.sum(rays * along_axis[:, np.newaxis], axis=-1) / down
    cross = np.sum(rays * left_axis[:, np.newaxis], axis=-1) / down
    centered = cross - cross.mean(axis=1, keepdims=True)
    slope = np.sum(centered * along, axis=1) / np.sum(centered**2, axis=1)
    middle = along[:, np.argmin(np.abs(np.median(cross, axis=0)))]
    yaw = np.degrees(np.arctan(slope))
    pitch = np.degrees(np.arctan(middle))
    ascending = np.gradient(np.asarray(lat)) > 0
    print(
        f"{frame}: yaw {np.median(yaw):+.4f} deg (p5 {np.percentile(yaw, 5):+.4f}, "
        f"p95 {np.percentile(yaw, 95):+.4f}); along-track angle at nadir "
        f"{np.median(pitch):+.4f} deg"
    )
    for name, chosen in (("ascending", ascending), ("descending", ~ascending)):
        bands = []
        for low in range(-90, 90, 30):
            band = chosen & (np.asarray(lat) >= low) & (np.asarray(lat) < low + 30)
            if np.any(band):
                bands.append(f"{low:+d}..{low + 30:+d}: {np.median(yaw[band]):+.3f}")
        print(f"  {name} yaw by latitude: " + ", ".join(bands))


def score_swath_width(reference: Reference, instrument, orbit, stride: int) -> float:
    """
    Gets the median absolute difference (km) from the width through the beams. The
    swath width of a pointed instrument follows its field of regard, centered on
    nadir, so that of a rolled or pitched view is found instead from the projections
    of its edges, and compared with the distance along the surface between the
    outermost beams.
    """
    scans = np.arange(0, len(reference.time), stride)
    track = orbit.get_orbit_track(list(reference.time[scans].to_pydatetime()))
    roll = getattr(instrument, "roll_angle", 0)
    if getattr(instrument, "pitch_angle_profile", None) or getattr(
        instrument, "pitch_angle", 0
    ):
        # the left and right ends of the view at each scan's pitch
        pitch = np.broadcast_to(instrument.get_pitch_angle(track), (len(scans),))
        ends = [
            [
                compute_projected_ray_position(
                    track[k],
                    instrument.cross_track_field_of_view,
                    1e-3,
                    roll_angle=roll,
                    pitch_angle=float(pitch[k]),
                    is_rectangular=True,
                    angle=angle,
                    velocity_frame=instrument.velocity_frame,
                )
                for angle in (0, 180)
            ]
            for k in range(len(scans))
        ]
        predicted = GEOD.inv(
            [float(e[0].longitude.degrees) for e in ends],
            [float(e[0].latitude.degrees) for e in ends],
            [float(e[1].longitude.degrees) for e in ends],
            [float(e[1].latitude.degrees) for e in ends],
        )[2]
        observed = GEOD.inv(
            reference.lon[scans, 0],
            reference.lat[scans, 0],
            reference.lon[scans, -1],
            reference.lat[scans, -1],
        )[2]
        return float(np.median(np.abs(predicted - observed))) / 1e3
    if roll:
        half = instrument.cross_track_field_of_view / 2
        edges = [
            compute_projected_ray_position(
                track,
                0,
                0,
                roll_angle=roll + side * half,
                velocity_frame=instrument.velocity_frame,
            )
            for side in (-1, 1)
        ]
        predicted = GEOD.inv(
            edges[0].longitude.degrees,
            edges[0].latitude.degrees,
            edges[1].longitude.degrees,
            edges[1].latitude.degrees,
        )[2]
        observed = GEOD.inv(
            reference.lon[scans, 0],
            reference.lat[scans, 0],
            reference.lon[scans, -1],
            reference.lat[scans, -1],
        )[2]
        return float(np.median(np.abs(predicted - observed))) / 1e3
    xyz = ecef(reference.lon[scans], reference.lat[scans])
    observed = np.sum(np.linalg.norm(np.diff(xyz, axis=1), axis=-1), axis=1)
    altitude = wgs84.geographic_position_of(track).elevation.m
    predicted = np.vectorize(instrument.get_swath_width)(altitude)
    return float(np.median(np.abs(predicted - observed))) / 1e3


def longest_span(reference: Reference, period: float) -> np.ndarray:
    """
    Gets the scans of the longest span of data, of at most an orbit's period (s),
    without a gap of more than a minute, for products that do not cover whole orbits.
    """
    breaks = np.nonzero(np.diff(reference.seconds) >= 60)[0] + 1
    spans = np.split(np.arange(len(reference.seconds)), breaks)
    longest = max(spans, key=lambda s: np.ptp(reference.seconds[s]) if len(s) else 0)
    return longest[reference.seconds[longest] < reference.seconds[longest[0]] + period]


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
            scans = longest_span(reference, period.total_seconds())
            print(
                "no orbit without a gap; the longest span without one lasts "
                f"{np.ptp(reference.seconds[scans]) / 60:.1f} min",
                flush=True,
            )
            break
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
        pixel = np.where(
            increasing,
            reference.pixel[beam],
            reference.pixels - 1 - reference.pixel[beam],
        )
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
    between the beams on either side (divided by the lines or samples between them on
    a reduced day). Gets the times, and rows of latitude, longitude, along-track and
    cross-track spacing (km).
    """
    scans = reference.contiguous_scans()[::stride]
    xyz = ecef(reference.lon, reference.lat)
    along = np.linalg.norm(xyz[scans + 1] - xyz[scans - 1], axis=-1) / (
        2e3 * reference.line_stride
    )
    cross = np.full(along.shape, np.nan)
    cross[:, 1:-1] = np.linalg.norm(xyz[scans, 2:] - xyz[scans, :-2], axis=-1) / (
        1e3 * (reference.pixel[2:] - reference.pixel[:-2])
    )
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


def validate_beams(
    setup: Setup, day, tle, entry, orbit, reference, write: bool, platform=None
):
    """
    Scores the sample size configuration, and writes its results if requested: on a
    platform other than the configuration's (`platform`), as transfer metrics.
    """
    platform = platform or setup.platform
    within, largest = score_published(
        compute_published_errors(entry, setup.beams, setup.published_footprints)
    )
    print(
        f"published dimensions: {within:.1f}% within 3%, largest error {largest:.2f}%"
    )
    satellite = entry.to_satellite(orbit, setup.beams, platform=platform)
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
    if platform != setup.platform:
        add_transfer_metrics(
            setup.record_path,
            setup.beams,
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


def validate_swath(
    setup: Setup, day, tle, entry, orbit, reference, write: bool, platform=None
):
    """
    Scores the swath configuration, and writes its results if requested: on a
    platform other than the configuration's (`platform`), as transfer metrics.
    """
    platform = platform or setup.platform

    def satellite_for(time_step: float):
        return entry.to_satellite(
            orbit, setup.swath, time_step=time_step, platform=platform
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
    if platform != setup.platform:
        add_transfer_metrics(
            setup.record_path,
            setup.swath,
            f"{platform}, {day:%Y-%m-%d}",
            {name: float(results[name]) for name in GATED},
        )
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

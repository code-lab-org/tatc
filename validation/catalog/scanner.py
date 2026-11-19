"""
Helpers shared by the scripts that validate cross-track scanners against swath
products whose scans give the position of each sample's center on the ellipsoid
(of NASA's DAACs, NOAA, and EUMETSAT): reading products of NASA's DAACs, fitting a
scan's angles, and scoring a sample size configuration (`validate_beams`) and a
swath configuration (`validate_swath`).

A product's scans are a `Reference`, a lattice of samples (`ScanLattice`). Products
read otherwise (such as those of the JPSS satellites) give their own reference,
with the methods that `validate_beams` and `validate_swath` call.

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
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric

from tatc.catalog.metrics import summarize
from tatc.generation import generate_points_fibonacci_lattice
from tatc.constants import timescale
from tatc.schemas import Attitude
from tatc.utils import (
    VelocityFrame,
    compute_projected_ray_position,
    compute_view_tangents,
)
from validation_io import (
    EARTH_ROTATION,
    TO_GEODETIC,
    RemoteFile,
    ScanLattice,
    build_grid,
    compute_published_errors,
    compute_spacing_ratios,
    ecef,
    get_with_retries,
    login_earthdata,
    normalize,
    pick_passes,
    recompute_passes,
    score_passes,
    score_published,
    score_swath_area,
    up,
    write_results,
)

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
    login_earthdata()
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


class Reference(ScanLattice):
    """
    A day of scans of a cross-track scanner, with what is derived from them. Each scan
    has a time; products that also give each beam's time (`beam_time`, in seconds from
    the start of the day) time the passes by beam. A day reduced to some of the
    product's lines and samples gives the samples' indices in the product (`pixel`),
    their number (`pixels`), and the lines from one kept scan to the next
    (`line_stride`).
    """

    def __init__(
        self, day: datetime, data: dict, scan_period: float, gap_threshold: float = 60
    ):
        # breaks longer than the threshold (s) are periods without data
        self.day = day
        self.gap_threshold = gap_threshold
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
            for i in np.nonzero(period > gap_threshold)[0]
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

    @property
    def sample_seconds(self) -> np.ndarray:
        """The time of each beam (s from the start of the day), scans by beams."""
        return self.beam_seconds

    def pass_time(self, before: int, sample: int, fraction: float) -> pd.Timestamp:
        """
        Gets the time of a pass between a beam of a scan and of the next, with each
        beam at its time (or the time of its scan, if the product gives no time by
        beam).
        """
        seconds = self.beam_seconds[before, sample] + fraction * (
            self.beam_seconds[before + 1, sample] - self.beam_seconds[before, sample]
        )
        return pd.Timestamp(self.day) + pd.Timedelta(seconds=float(seconds))

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


def fit_view(
    reference: Reference,
    frame: str,
    yaw_angle: float,
    pitch_angle: float,
    lines: int = 3000,
    iterations: int = 8,
    instrument_pitch: bool = False,
) -> None:
    """
    Prints the view that TAT-C centers on the outermost beams: given the body's yaw
    and pitch, the instrument's roll that centers the first and last beams across
    track, the body's pitch that centers the middle beam along track, and the field
    of view between the first and last beams, from TAT-C's view tangents, with the
    spacecraft's position and velocity of the product. A roll measured from the
    geodetic nadir (`fit`) does not center a view that is also pitched and yawed.
    With `instrument_pitch`, the pitch is the instrument's (applied after its roll),
    starting from `pitch_angle`, and the body is only yawed.
    """
    contiguous = (
        np.r_[False, reference.is_next_contiguous]
        & np.r_[reference.is_next_contiguous, False]
    )
    candidates = np.flatnonzero(contiguous)
    pick = candidates[np.linspace(0, len(candidates) - 1, lines).astype(int)]
    t = timescale.from_datetimes(list(reference.time[pick].to_pydatetime()))
    rotation = itrs.rotation_at(t)  # takes inertial (GCRS) coordinates to ITRS
    au_m = 149597870700.0
    position = np.einsum("jin,nj->in", rotation, reference.position[pick]) / au_m
    velocity = (
        np.einsum("jin,nj->in", rotation, reference.velocity("inertial")[pick])
        / au_m
        * 86400
    )
    track = Geocentric(position, velocity, t, center=399)
    velocity_frame = VelocityFrame(frame)
    middle = int(np.argmin(np.abs(np.median(reference.cross_track_angles(frame), 0))))
    targets = [
        wgs84.latlon(reference.lat[pick, k], reference.lon[pick, k])
        for k in (0, middle, reference.beams - 1)
    ]
    roll, pitch = 0.0, pitch_angle

    def tangents(roll, pitch):
        body = 0.0 if instrument_pitch else pitch
        attitude = Attitude(
            velocity_frame=velocity_frame, yaw_angle=yaw_angle, pitch_angle=body
        )
        return [
            compute_view_tangents(
                track,
                target,
                velocity_frame,
                roll,
                pitch if instrument_pitch else 0,
                attitude=attitude,
            )
            for target in targets
        ]

    for _ in range(iterations):
        (a0, c0), (am, _), (a1, c1) = tangents(roll, pitch)
        roll += float(np.degrees(np.arctan(np.median((c0 + c1) / 2))))
        pitch += float(np.degrees(np.arctan(np.median(am))))
    (a0, c0), (am, _), (a1, c1) = tangents(roll, pitch)
    width = 2 * np.degrees(np.arctan(np.median(np.abs(c0 - c1)) / 2))
    slope = np.degrees(np.arctan(np.median((a0 - a1) / (c0 - c1))))
    print(
        f"{frame} view, with yaw {yaw_angle:.4f} deg: roll {roll:.4f} deg, "
        f"{'instrument' if instrument_pitch else 'body'} pitch "
        f"{pitch:.4f} deg, cross-track field of view {width:.4f} deg; remaining "
        f"across tangents of the ends p5..p95 {np.percentile(c0, 5):+.5f}..."
        f"{np.percentile(c0, 95):+.5f} and {np.percentile(c1, 5):+.5f}..."
        f"{np.percentile(c1, 95):+.5f}; middle along p5..p95 "
        f"{np.percentile(am, 5):+.5f}..{np.percentile(am, 95):+.5f}; slope of the "
        f"line along on across {slope:+.4f} deg; {len(pick)} lines"
    )


def score_swath_width(
    reference: Reference, instrument, orbit, stride: int, attitude=None
) -> float:
    """
    Gets the median absolute difference (km) from the width through the beams. The
    swath width of a pointed instrument follows its field of regard, centered on
    nadir, so that of a rolled or pitched view is found instead from the projections
    of its edges, and compared with the distance along the surface between the
    outermost beams. A view is pitched by the instrument or by the body (the
    satellite's attitude, whose yaw and roll are taken as zero here).
    """
    scans = np.arange(0, len(reference.time), stride)
    track = orbit.get_orbit_track(list(reference.time[scans].to_pydatetime()))
    roll = getattr(instrument, "roll_angle", 0)
    pitch = np.zeros(len(scans))
    if hasattr(instrument, "get_pitch_angle"):
        pitch = pitch + instrument.get_pitch_angle(track)
    if attitude is not None:
        pitch = pitch + attitude.get_body_angles(track)[2]
    if np.any(pitch != 0):
        # the left and right ends of the view at each scan's pitch
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
    Gets the scans of the longest span of data, of at most a period (s), without a gap
    of more than a minute.
    """
    breaks = np.nonzero(np.diff(reference.seconds) >= 60)[0] + 1
    spans = np.split(np.arange(len(reference.seconds)), breaks)
    longest = max(spans, key=lambda s: np.ptp(reference.seconds[s]) if len(s) else 0)
    return longest[reference.seconds[longest] < reference.seconds[longest[0]] + period]


def observe_orbit(reference: Reference, orbit, grid, radius: float):
    """
    Gets the cells compared over a span of the day's data, and which of them were
    observed. The span is not one orbit but the rest of the data from noon (or from
    the first later start, in steps of 10 minutes until 22:00, after which no gap of a
    minute or more remains), or else the longest span without such a gap.
    """
    lon, lat, _ = grid
    unit = up(lon, lat)
    begin = pd.Timestamp(reference.day) + pd.Timedelta(hours=12)
    # the span is bounded by 1 / mean motion days; TAT-C's mean motion is in deg/s,
    # so this bound (about 15 days) exceeds the data, and every record scored with
    # this function compares the span described above, rather than one orbit
    period = pd.Timedelta(days=1 / orbit.elements[0].mean_motion)
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
                "no span from noon without a gap; the longest span without one "
                f"lasts {np.ptp(reference.seconds[scans]) / 60:.1f} min",
                flush=True,
            )
            break
    first, last = reference.time[scans[0]], reference.time[scans[-1]]
    chord = 2 * np.sin(radius / 6371 / 2)
    # bounding the search keeps it fast for cells far from the span's thin swath,
    # which are not compared (their distance is infinite, their index out of range)
    vertices = reference.unit[scans].reshape(-1, 3)
    distance, nearest = cKDTree(vertices).query(
        unit, distance_upper_bound=chord, workers=-1
    )
    elapsed = reference.seconds[scans] - reference.seconds[scans[0]]
    nearest_scan = np.minimum(nearest, len(vertices) - 1) // reference.beams
    nearest_time = elapsed[nearest_scan]
    compared = (
        (distance < chord)
        & (nearest_time > ORBIT_MARGIN)
        & (nearest_time < (last - first).total_seconds() - ORBIT_MARGIN)
    )
    if reference.gap_threshold < 60:
        # cells nearest the scans on either side of a shorter period without data
        # (as the strip it leaves unimaged) are not compared
        breaks = np.nonzero(np.diff(elapsed) > reference.gap_threshold)[0]
        bordering = np.zeros(len(scans), dtype=bool)
        bordering[breaks] = bordering[breaks + 1] = True
        compared &= ~bordering[nearest_scan]
        print(
            f"{len(breaks)} periods without data over {reference.gap_threshold} s in "
            "the span compared",
            flush=True,
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
    """
    What a script validates on a platform: the configurations of its swath and
    beams (or samples), the instrument and published footprints of its sample
    sizes, its references, its fixture, and the parameters of its tests.
    """

    swath: str
    beams: str | None
    size_instrument: str | None
    published_footprints: tuple
    fixture_path: pathlib.Path
    reference: str
    published_reference: str
    time_steps: tuple[float, ...]
    pass_radius: float
    point_spacing: float = 500e3
    pixel_stride: int = 100
    spacing_stride: int = 150
    width_stride: int = 5
    fixture_passes: int = 60
    fixture_pixels: int = 200
    fixture_spacing: int = 200


def validate_beams(setup: Setup, run, reference, orbit, tle) -> None:
    """
    Scores the sample size configuration of a run (see `validation_io.parse_run`),
    and writes its results if the run writes: as a transfer, on a platform on which
    the configuration was not tested. A reference read otherwise than as a
    `Reference` gives its own `measure_pixels(satellite, index, setup, frame)` and
    `measure_spacing(setup)`.
    """
    entry = run.entry
    within, largest = score_published(
        compute_published_errors(entry, setup.beams, setup.published_footprints)
    )
    print(
        f"published dimensions: {within:.1f}% within 3%, largest error {largest:.2f}%"
    )
    satellite = entry.build(
        configuration=setup.beams, orbit=orbit, platform=run.platform
    ).to_satellite(orbit, name=entry.name)
    names = [instrument.name for instrument in satellite.instruments]
    index = names.index(setup.size_instrument)
    frame = entry.configurations[setup.beams].velocity_frame.value
    if hasattr(reference, "measure_pixels"):
        pixels = reference.measure_pixels(satellite, index, setup, frame)
    else:
        pixels = measure_pixels(reference, satellite, index, setup.pixel_stride, frame)
    median = summarize(pixels[:, 2])["median"]
    p95 = summarize(pixels[:, 2])["p95"]
    print(f"pixel position: median {median:.3f} km, 95th percentile {p95:.3f} km")
    if hasattr(reference, "measure_spacing"):
        times, rows = reference.measure_spacing(setup)
    else:
        times, rows = measure_spacing(reference, setup.spacing_stride)
    along, cross = compute_spacing_ratios(
        satellite, [setup.size_instrument] * len(rows), times, rows
    )
    print(f"spacing, predicted / measured: {along:.4f} along, {cross:.4f} across")
    if not run.write:
        return
    fixture = fixture_metrics = None
    if not run.transfer:
        rng = np.random.default_rng(3)
        picks = rng.choice(len(rows), setup.fixture_spacing, replace=False)
        fixture_along, fixture_cross = compute_spacing_ratios(
            satellite, [setup.size_instrument] * len(picks), times[picks], rows[picks]
        )
        reduced = pixels[rng.choice(len(pixels), setup.fixture_pixels, replace=False)]
        fixture = {
            "tle": np.array(tle),
            "day": np.array(run.day.isoformat()),
            "pixels": reduced,
            "spacing": rows[picks],
            "spacing_seconds": (times[picks] - pd.Timestamp(run.day))
            .total_seconds()
            .to_numpy(),
        }
        fixture_metrics = {
            "sample_size_within_3pct": within,
            "pixel_position_median": summarize(reduced[:, 2])["median"],
            "spacing_along_ratio_median": fixture_along,
            "spacing_cross_ratio_median": fixture_cross,
        }
    write_results(
        run,
        setup.beams,
        reference=setup.published_reference,
        metrics=[
            {"name": "sample_size_within_3pct", "value": within},
            {"name": "sample_size_largest_error", "value": largest},
            {"name": "pixel_position_median", "value": median},
            {"name": "pixel_position_p95", "value": p95},
            {"name": "spacing_along_ratio_median", "value": along},
            {"name": "spacing_cross_ratio_median", "value": cross},
        ],
        transfer_metrics={
            "pixel_position_median": median,
            "pixel_position_p95": p95,
            "spacing_along_ratio_median": along,
            "spacing_cross_ratio_median": cross,
        },
        fixture=fixture,
        fixture_path=setup.fixture_path,
        fixture_metrics=fixture_metrics,
    )


def validate_swath(
    setup: Setup, run, reference, orbit, tle, *, fixture_extras: dict | None = None
) -> None:
    """
    Scores the swath configuration of a run (see `validation_io.parse_run`), and
    writes its results if the run writes: as a transfer, on a platform on which the
    configuration was not tested. A reference read otherwise than as a `Reference`
    gives its own `span` (the start and end of its data),
    `score_swath_width(instrument, orbit, setup, attitude)`, and
    `observe_area(orbit, grid, setup)`. Other arrays that the tests of the fixture
    need are given as `fixture_extras`.
    """
    entry, day = run.entry, run.day

    def satellite_for(time_step: float):
        # the period of the day, for a configuration that schedules its attitude
        return entry.to_satellite(
            orbit,
            configuration=setup.swath,
            time_step=time_step,
            platform=run.platform,
            start=day,
            end=day + timedelta(days=1),
        )

    satellite = satellite_for(setup.time_steps[0])
    points = generate_points_fibonacci_lattice(setup.point_spacing)
    observed = reference.observe_passes(points, setup.pass_radius)
    if hasattr(reference, "span"):
        start, end = reference.span
    else:
        start = reference.time[0].to_pydatetime()
        end = reference.time[-1].to_pydatetime() + timedelta(
            seconds=reference.scan_period
        )
    results = score_passes(points, satellite, 0, observed, start, end, reference.gaps)
    if hasattr(reference, "score_swath_width"):
        results["swath_width_median"] = reference.score_swath_width(
            satellite.instruments[0], orbit, setup, satellite.attitude
        )
    else:
        results["swath_width_median"] = score_swath_width(
            reference,
            satellite.instruments[0],
            orbit,
            setup.width_stride,
            satellite.attitude,
        )
    print(json.dumps(results), flush=True)
    grid = build_grid()
    if hasattr(reference, "observe_area"):
        first, last, compared, covered = reference.observe_area(orbit, grid, setup)
    else:
        first, last, compared, covered = observe_orbit(
            reference, orbit, grid, 2 * setup.pass_radius
        )
    results["footprint_iou_median"] = score_swath_area(
        satellite_for, 0, first, last, setup.time_steps, grid, compared, covered
    )
    print(json.dumps(results), flush=True)
    if not run.write:
        return
    fixture = fixture_metrics = None
    if not run.transfer:
        passes = pick_passes(points, observed, day, setup.fixture_passes)
        fixture = {
            "tle": np.array(tle),
            "day": np.array(day.isoformat()),
            "passes": passes,
            **(fixture_extras or {}),
        }
        fixture_metrics = recompute_passes(satellite, 0, day, passes)
    write_results(
        run,
        setup.swath,
        reference=setup.reference,
        metrics=[{"name": name, "value": float(results[name])} for name in GATED],
        transfer_metrics={name: float(results[name]) for name in GATED},
        fixture=fixture,
        fixture_path=setup.fixture_path,
        fixture_metrics=fixture_metrics,
    )

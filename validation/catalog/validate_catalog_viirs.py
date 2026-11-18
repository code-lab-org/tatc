"""
Validates the VIIRS catalog entry against the NOAA-20 VIIRS moderate-resolution
geolocation product, using the exact values of the entry through public TAT-C
methods. The tests of pixels and swath width follow `ValidateImagerVIIRS`, on whose
day of data the entry's angles were fitted, and those of sizes and spacing follow
`ValidateSampleSizes`; the tests of passes, epochs, and swath area are this script's
own. Running on that day checks the shared tests against the notebooks, while running
on a day at least 30 days from it tests the fitted angles on data held out from the
fit.

The notebooks do not test the swath's area, which the standard requires. The observed
area over one orbit is the union of the outlines (G-rings) of the granules and of the
strips between consecutive granules, on whose outermost pixel centers the outlines lie.

Reads anonymous public data, and needs the `validation` extra:

    python validation/catalog/validate_catalog_viirs.py <data directory> 2026-10-02
    python validation/catalog/validate_catalog_viirs.py <data directory> 2026-08-12

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

NOAA-21's and Suomi NPP's VIIRS have their own configurations, whose values are
fitted with `--fit` (the notebooks' are NOAA-20's): the swath's extent, and the zones
as `ValidateSampleSizes` measured them, which the fit reproduces on NOAA-20's day.
They are tested on a day at least 30 days from the fit:

    python validation/catalog/validate_catalog_viirs.py <data directory> 2026-11-05 --platform noaa-21 --fit
    python validation/catalog/validate_catalog_viirs.py <data directory> 2026-10-02 --platform noaa-21
    python validation/catalog/validate_catalog_viirs.py <data directory> 2026-10-02 --platform suomi-npp --fit
    python validation/catalog/validate_catalog_viirs.py <data directory> 2026-08-25 --platform suomi-npp

With `--transfer`, it tests NOAA-20's configurations on another platform instead, and
adds the results to the record as a transfer, which does not count toward the tier.
Without a day, it takes that of the entry's record (or, with `--fit`, of the fit),
and `--configuration` tests one configuration (such as `zones_noaa_20`, only the
zones).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import concurrent.futures
import pathlib
from datetime import datetime

import h5py
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from tatc.catalog.metrics import summarize
from tatc.schemas import GeneralPerturbationsOrbit
import scanner
from validation_io import (
    ARCHIVES,
    FIXTURES,
    TO_GEODETIC,
    Archive,
    RemoteFile,
    compute_published_errors,
    compute_spacing_ratios,
    ecef,
    iet_to_utc,
    inertial_velocity,
    normalize,
    parse_run,
    score_published,
    up,
    write_results,
)

FIT_PLATFORM = "noaa-20"
# the configurations of each platform's fit and its fixture
CONFIGURATIONS = {
    "noaa-20": ("swath_noaa_20", "zones_noaa_20", FIXTURES / "viirs_noaa_20.npz"),
    "noaa-21": ("swath_noaa_21", "zones_noaa_21", FIXTURES / "viirs_noaa_21.npz"),
    "suomi-npp": (
        "swath_suomi_npp",
        "zones_suomi_npp",
        FIXTURES / "viirs_suomi_npp.npz",
    ),
}
GEO = "All_Data/VIIRS-MOD-GEO-TC_All/"
GRANULE = "Data_Products/VIIRS-MOD-GEO-TC/VIIRS-MOD-GEO-TC_Gran_0"
N_SCANS, N_DETECTORS, N_SAMPLES = 48, 16, 3200
PIXEL_GRANULES = 10
COLUMNS = np.unique(np.r_[np.arange(0, N_SAMPLES, 8), N_SAMPLES - 1])
ROWS = np.array([0, 7, 8, 15])
# the first sample of each zone, from right to left, and the end of the last
ZONE_BOUNDS = np.array([0, 640, 1007, 2191, 2559, 3200])
ORBIT_MARGIN = 120
# published footprints: instrument, altitude (km), inclination (deg), scan angle
# (deg), and sizes along and across track (km), from NOAA/CIMSS for band M16
PUBLISHED_FOOTPRINTS = (
    ("VIIRS zone 2", 824, 98.7, 0, 0.742, 0.776),
    ("VIIRS zone 4", 824, 98.7, 56.0, 1.60, 1.58),
)


def get_setup(platform: str, scan_period: float) -> scanner.Setup:
    """
    Gets what is validated with a platform's fit: its configurations and fixture,
    with the swath's area compared at time steps from the day's scan period.
    """
    swath, zones, fixture_path = CONFIGURATIONS[platform]
    return scanner.Setup(
        swath=swath,
        beams=zones,
        # each zone is an instrument of its own
        size_instrument=None,
        published_footprints=PUBLISHED_FOOTPRINTS,
        fixture_path=fixture_path,
        reference="VIIRS moderate-resolution geolocation",
        published_reference="Footprints published by NOAA/CIMSS; samples of the "
        "VIIRS moderate-resolution geolocation",
        time_steps=(scan_period, 30, 120),
        # passes are found in the granules' outlines, not within a radius
        pass_radius=0,
    )


def open_granule(url: str) -> h5py.File:
    """Opens a remote granule, first getting the metadata at its start and end."""
    remote = RemoteFile(url)
    remote.fetch(0, remote.block_size)
    remote.fetch(
        (remote.size - 2**19) // remote.block_size * remote.block_size, remote.size
    )
    return h5py.File(remote)


def read_scans(url: str) -> dict:
    """Reads the scans and outline of a granule."""
    with open_granule(url) as f:
        granule = f[GRANULE]
        return {
            "gring_lat": granule.attrs["G-Ring_Latitude"].ravel(),
            "gring_lon": granule.attrs["G-Ring_Longitude"].ravel(),
            "n_scans": int(granule.attrs["N_Number_Of_Scans"].ravel()[0]),
            **{
                field: f[GEO + field][:]
                for field in ("StartTime", "MidTime", "SCPosition", "SCVelocity")
            },
        }


def read_pixels(url: str) -> dict:
    """Reads a sample of the pixels of a granule."""
    with open_granule(url) as f:
        lat, lon, height = (
            f[GEO + field][:].reshape(N_SCANS, N_DETECTORS, N_SAMPLES)
            for field in ("Latitude", "Longitude", "Height")
        )
    center = N_SCANS // 2
    return {
        "lat": lat[:, ROWS][..., COLUMNS],
        "lon": lon[:, ROWS][..., COLUMNS],
        "height": height[:, ROWS][..., COLUMNS],
        "lat_full": lat[center, 7:9],
        "lon_full": lon[center, 7:9],
        "height_full": height[center, 7:9],
    }


def load_day(
    day: datetime, directory: pathlib.Path, archive: Archive
) -> tuple[dict, dict]:
    """Downloads and reduces the geolocation of a day, or loads it from a cache."""
    scans_path = directory / f"viirs_gmtco_scans_{day:%Y%m%d}.npz"
    pixels_path = directory / f"viirs_gmtco_pixels_{day:%Y%m%d}.npz"
    if not (scans_path.exists() and pixels_path.exists()):
        keys = archive.list_keys(f"VIIRS-MOD-GEO-TC/{day:%Y/%m/%d}/")
        urls = [archive.url + key for key in keys]
        with concurrent.futures.ThreadPoolExecutor(8) as pool:
            scans = list(pool.map(read_scans, urls))
            pixels = list(pool.map(read_pixels, urls[::PIXEL_GRANULES]))
        np.savez_compressed(
            scans_path,
            **{field: np.stack([s[field] for s in scans]) for field in scans[0]},
        )
        np.savez_compressed(
            pixels_path,
            granule=np.arange(len(keys))[::PIXEL_GRANULES],
            **{field: np.stack([p[field] for p in pixels]) for field in pixels[0]},
        )
    return dict(np.load(scans_path)), dict(np.load(pixels_path))


class Reference:
    """
    A day of VIIRS scans and pixels, with what is derived from them. The observed
    swath is that of the granules' outlines (G-rings).
    """

    def __init__(self, scans: dict, pixels: dict):
        self.scans, self.pixels = scans, pixels
        valid = np.arange(N_SCANS)[np.newaxis, :] < scans["n_scans"][:, np.newaxis]
        self.time = iet_to_utc(scans["MidTime"][valid])
        self.position = scans["SCPosition"][valid].astype(float)
        self.velocity = scans["SCVelocity"][valid].astype(float)
        self.scan_period = float(np.median(np.diff(scans["StartTime"][valid]))) / 1e6
        self.granules = pixels["granule"]
        self.ok = np.all(np.abs(pixels["lat"]) <= 90, axis=(2, 3)) & (
            np.arange(N_SCANS)[np.newaxis, :]
            < scans["n_scans"][self.granules, np.newaxis]
        )
        self.ring = up(scans["gring_lon"], scans["gring_lat"])
        self.first = iet_to_utc(scans["MidTime"][:, 0])
        self.last = iet_to_utc(
            scans["MidTime"][np.arange(len(self.ring)), scans["n_scans"] - 1]
        )
        seconds = (self.time - self.time[0]).total_seconds().to_numpy()
        self.gaps = [
            (self.time[i], self.time[i + 1])
            for i in np.nonzero(np.diff(seconds) > 180)[0]
        ]

    def observe_passes(self, points, _radius=None) -> pd.DataFrame:
        """Finds the granules in whose outlines points lie, by pass, and their times."""
        point_xyz = up(points.geometry.x.to_numpy(), points.geometry.y.to_numpy())
        first_nadir = normalize(self.scans["SCPosition"][:, 0].astype(float))
        rows = []
        for g, vertices in enumerate(self.ring):
            edges = np.cross(vertices, np.roll(vertices, -1, axis=0))
            orientation = np.sign(edges @ vertices.mean(axis=0))
            inside = np.all((point_xyz @ edges.T) * orientation >= 0, axis=1)
            if not inside.any():
                continue
            lines = [
                normalize(np.cross(vertices[2], vertices[4])),
                normalize(np.cross(vertices[6], vertices[0])),
            ]
            if np.abs(lines[1] @ first_nadir[g]) < np.abs(lines[0] @ first_nadir[g]):
                lines = lines[::-1]
            to_first, to_last = (np.abs(point_xyz[inside] @ line) for line in lines)
            times = self.first[g] + (self.last[g] - self.first[g]) * (
                to_first / (to_first + to_last)
            )
            rows += [
                {"point_id": p, "time": t}
                for p, t in zip(np.flatnonzero(inside), times)
            ]
        observed = (
            pd.DataFrame(rows)
            .astype({"point_id": "int64"})
            .sort_values(["point_id", "time"])
        )
        new_pass = observed.groupby("point_id").time.diff().fillna(
            pd.Timedelta(days=1)
        ) > pd.Timedelta(minutes=20)
        return observed[new_pass.to_numpy()].reset_index(drop=True)

    def bridge(self, g: int) -> list[np.ndarray]:
        """
        Gets the outlines of the strip between a granule and the next, if contiguous.
        The strip is about as wide as a scan, so a bend in the scan line would put the
        center of a single outline on the wrong side of an edge: it has two parts.
        """
        gap = (self.first[g + 1] - self.last[g]).total_seconds()
        if gap > 1.5 * self.scan_period:
            return []
        lines = ((2, 3, 4), (6, 7, 0))
        after, before = self.ring[g + 1].mean(axis=0), self.ring[g].mean(axis=0)
        near = max(lines, key=lambda line: np.sum(self.ring[g][list(line)] @ after))
        other = max(
            lines, key=lambda line: np.sum(self.ring[g + 1][list(line)] @ before)
        )
        a, b = self.ring[g][list(near)], self.ring[g + 1][list(other)]
        if a[0] @ b[0] < a[0] @ b[2]:
            b = b[::-1]
        return [
            np.array([a[0], a[1], b[1], b[0]]),
            np.array([a[1], a[2], b[2], b[1]]),
        ]

    @property
    def span(self) -> tuple[datetime, datetime]:
        """The first and last scans' middle times."""
        return self.time[0].to_pydatetime(), self.time[-1].to_pydatetime()

    def score_swath_width(self, instrument, _orbit, _setup, _attitude=None) -> float:
        """Gets the median absolute difference (km) from the width through the samples."""
        ok = self.ok
        lat, lon, height = (self.pixels[k][ok] for k in ("lat", "lon", "height"))
        center = (
            ecef(lon[:, 1], lat[:, 1], height[:, 1])
            + ecef(lon[:, 2], lat[:, 2], height[:, 2])
        ) / 2
        observed = np.sum(np.linalg.norm(np.diff(center, axis=1), axis=-1), axis=1)
        position = self.scans["SCPosition"][self.granules].astype(float)[ok]
        altitude = TO_GEODETIC.transform(*position.T)[2]
        predicted = np.vectorize(instrument.get_swath_width)(altitude)
        return float(np.median(np.abs(predicted - observed))) / 1e3

    def observe_area(self, orbit, grid, _setup):
        """Gets the cells compared over one orbit, and which of them VIIRS observed."""
        return observe_orbit(self, orbit, grid)


def in_outline(vertices: np.ndarray, points: np.ndarray) -> np.ndarray:
    """Gets whether unit vectors lie in an outline connected by great circles."""
    edges = np.cross(vertices, np.roll(vertices, -1, axis=0))
    orientation = np.sign(edges @ vertices.mean(axis=0))
    return np.all((points @ edges.T) * orientation >= 0, axis=1)


def observe_orbit(reference: Reference, orbit, grid):
    """Gets the cells compared over one orbit, and which of them VIIRS observed."""
    lon, lat, _ = grid
    unit = up(lon, lat)
    begin = reference.time[0].floor("D") + pd.Timedelta(hours=12)
    period = pd.Timedelta(days=1 / orbit.elements[0].mean_motion * 360 / 86400)
    granules = np.nonzero(
        (reference.first >= begin) & (reference.last < begin + period)
    )[0]
    first, last = reference.first[granules[0]], reference.last[granules[-1]]
    # points along the granules' edges, to measure distance from the swath
    fractions = np.linspace(0, 1, 20)[:, np.newaxis, np.newaxis, np.newaxis]
    rings = reference.ring[granules]
    dense = normalize(
        (1 - fractions) * rings[np.newaxis] + fractions * np.roll(rings, -1, axis=1)
    )
    dense = dense.transpose(1, 0, 2, 3).reshape(len(granules), -1, 3)
    distance, nearest = cKDTree(dense.reshape(-1, 3)).query(unit, workers=-1)
    elapsed = (reference.first[granules] - first).total_seconds().to_numpy()
    nearest_time = elapsed[nearest // dense.shape[1]]
    compared = (
        (distance < 2 * np.sin(300 / 6371 / 2))
        & (nearest_time > ORBIT_MARGIN)
        & (nearest_time < (last - first).total_seconds() - ORBIT_MARGIN)
    )
    candidates = unit[compared]
    inside = np.zeros(len(candidates), dtype=bool)
    for g in granules:
        inside |= in_outline(reference.ring[g], candidates)
        if g != granules[-1]:
            for outline in reference.bridge(g):
                inside |= in_outline(outline, candidates)
    observed = np.zeros(len(unit), dtype=bool)
    observed[compared] = inside
    return first.to_pydatetime(), last.to_pydatetime(), compared, observed


def measure_pixels(reference: Reference, satellite) -> np.ndarray:
    """
    Gets rows of seconds from the start of the day, zone, pixel index in the zone,
    position difference (km) from the sample, and the sample's latitude and longitude,
    for the center detectors of the middle scans of a sample of granules, every eighth
    sample across the scan.
    """
    center = N_SCANS // 2
    pixels = reference.pixels
    ok = np.all(np.abs(pixels["lat_full"]) <= 90, axis=(1, 2)) & (
        center < reference.scans["n_scans"][reference.granules]
    )
    times = iet_to_utc(reference.scans["MidTime"][reference.granules[ok], center])
    observed = (
        ecef(pixels["lon_full"][ok, 0], pixels["lat_full"][ok, 0])
        + ecef(pixels["lon_full"][ok, 1], pixels["lat_full"][ok, 1])
    ) / 2
    lat = pixels["lat_full"][ok].astype(float).mean(axis=1)
    lon = pixels["lon_full"][ok].astype(float).mean(axis=1)
    track = satellite.orbit.to_gp_orbit().get_orbit_track(list(times))
    seconds = (times - times[0].floor("D")).total_seconds().to_numpy()
    rows = []
    for sample in COLUMNS:
        zone = int(np.searchsorted(ZONE_BOUNDS, sample, side="right") - 1)
        pixel = int(sample - ZONE_BOUNDS[zone])
        instrument = satellite.get_instrument(zone)
        fore, aft = (
            instrument.compute_projected_pixel_position(track, pixel, detector)
            for detector in (7, 8)
        )
        predicted = (
            ecef(fore.longitude.degrees, fore.latitude.degrees)
            + ecef(aft.longitude.degrees, aft.latitude.degrees)
        ) / 2
        difference = np.linalg.norm(predicted - observed[:, sample], axis=-1) / 1e3
        rows.append(
            np.stack(
                [
                    seconds,
                    np.full(len(seconds), zone),
                    np.full(len(seconds), pixel),
                    difference,
                    lat[:, sample],
                    lon[:, sample],
                ],
                axis=-1,
            )
        )
    return np.concatenate(rows)


def measure_spacing(reference: Reference) -> tuple[pd.DatetimeIndex, np.ndarray]:
    """
    Measures the spacing of samples as `ValidateSampleSizes` does, within zones. Gets
    the times, and rows of latitude, longitude, along-track and cross-track spacing
    (km), and zone: along track, the larger of the detector pitch and the advance
    between scans per detector.
    """
    center = N_SCANS // 2
    pixels = reference.pixels
    ok = np.all(np.abs(pixels["lat_full"]) <= 90, axis=(1, 2)) & (
        center + 1 < reference.scans["n_scans"][reference.granules]
    )
    ok &= np.all(np.abs(pixels["lat"][:, center : center + 2]) <= 90, axis=(1, 2, 3))
    lat = pixels["lat_full"][ok].astype(float).mean(axis=1)
    lon = pixels["lon_full"][ok].astype(float).mean(axis=1)
    xyz = ecef(lon, lat)
    cross = np.full(lat.shape, np.nan)
    cross[:, 1:-1] = np.linalg.norm(xyz[:, 2:] - xyz[:, :-2], axis=-1) / 2e3
    first, last, following = (
        ecef(pixels["lon"][ok, scan, row], pixels["lat"][ok, scan, row])
        for scan, row in ((center, 0), (center, 3), (center + 1, 0))
    )
    pitch = np.linalg.norm(last - first, axis=-1) / 15e3
    advance = np.linalg.norm(following - first, axis=-1) / 16e3
    along = np.fmax(pitch, advance)
    boundary = np.r_[ZONE_BOUNDS[1:-1] - 1, ZONE_BOUNDS[1:-1]]
    keep = ~np.isin(COLUMNS, boundary) & (COLUMNS > 0) & (COLUMNS < N_SAMPLES - 1)
    columns = COLUMNS[keep]
    zones = np.searchsorted(ZONE_BOUNDS, columns, side="right") - 1
    rows = np.stack(
        [
            lat[:, columns].ravel(),
            lon[:, columns].ravel(),
            along[:, keep].ravel(),
            cross[:, columns].ravel(),
            np.tile(zones, len(lat)),
        ],
        axis=-1,
    )
    times = np.repeat(
        iet_to_utc(reference.scans["MidTime"][reference.granules[ok], center]),
        len(columns),
    )
    return pd.DatetimeIndex(times), rows


def fit(day: datetime, reference: Reference) -> None:
    """
    Prints the angles fitted to a day in the inertial frame (with
    `scanner.fit` and `scanner.fit_yaw`), from the middle of each
    sampled scan, between its two middle detectors: those of the outermost samples,
    and the yaw and along-track angle of the scan line.
    """
    pixels, ok = reference.pixels, reference.ok
    lat = (pixels["lat"][:, :, 1] + pixels["lat"][:, :, 2]) / 2
    # the mean of two longitudes, across the antimeridian
    lon = np.degrees(
        np.angle(
            np.exp(1j * np.radians(pixels["lon"][:, :, 1]))
            + np.exp(1j * np.radians(pixels["lon"][:, :, 2]))
        )
    )
    position = reference.scans["SCPosition"][reference.granules][ok].astype(float)
    sc_lon, sc_lat, sc_alt = TO_GEODETIC.transform(*position.T)
    seconds = (
        (
            iet_to_utc(reference.scans["MidTime"][reference.granules][ok])
            - pd.Timestamp(day)
        )
        .total_seconds()
        .to_numpy()
    )
    order = np.argsort(seconds)
    sampled = scanner.Reference(
        day,
        {
            "time": seconds[order],
            "lat": lat[ok][order],
            "lon": lon[ok][order],
            "sc_lat": np.asarray(sc_lat)[order],
            "sc_lon": np.asarray(sc_lon)[order],
            "sc_alt": np.asarray(sc_alt)[order],
            "pixel": COLUMNS,
            "pixels": N_SAMPLES,
        },
        reference.scan_period,
    )
    scanner.fit(sampled, "inertial")
    scanner.fit_yaw(sampled, "inertial")
    fit_zones(reference)


def fit_zones(reference: Reference) -> None:
    """
    Prints the values of the zones configuration measured as `ValidateSampleSizes`
    measured them, from every sample of the middle scans of a sample of granules
    (between their two middle detectors), with the median scan angles (positive left
    of the inertial velocity, from the geodetic nadir): each zone of `ZONE_BOUNDS`,
    by which the tests assign samples to zones, has its pixels, a field of view of
    its pixels times its median spacing, and a roll angle at the middle of its
    outermost samples. The bounds of runs of constant spacing are printed for
    comparison, as a step across two zones can place one a sample to either side.
    """
    center = N_SCANS // 2
    pixels = reference.pixels
    ok = np.all(np.abs(pixels["lat_full"]) <= 90, axis=(1, 2)) & (
        center < reference.scans["n_scans"][reference.granules]
    )
    observed = (
        ecef(pixels["lon_full"][ok, 0], pixels["lat_full"][ok, 0])
        + ecef(pixels["lon_full"][ok, 1], pixels["lat_full"][ok, 1])
    ) / 2
    position = reference.scans["SCPosition"][reference.granules[ok], center].astype(
        float
    )
    velocity = inertial_velocity(
        position,
        reference.scans["SCVelocity"][reference.granules[ok], center].astype(float),
    )
    lon, lat, _ = TO_GEODETIC.transform(*position.T)
    nadir = -up(lon, lat)
    left = normalize(np.cross(velocity, nadir))
    sight = normalize(observed - position[:, np.newaxis])
    angle = np.median(
        np.degrees(
            np.arctan2(
                np.sum(sight * left[:, np.newaxis], axis=-1),
                np.sum(sight * nadir[:, np.newaxis], axis=-1),
            )
        ),
        axis=0,
    )
    step = np.abs(np.diff(angle))
    aggregation = np.round(step / np.median(step[:100])).astype(int)
    runs = np.r_[0, np.nonzero(np.diff(aggregation))[0] + 1, N_SAMPLES]
    # the detector pitch: the angle between the first and last detectors (rows 0 and
    # 3 of those kept) of the middle scans, over 15 pitches, within 10 deg of nadir
    near = np.abs(angle[COLUMNS]) < 10
    first, last = (
        normalize(
            ecef(pixels["lon"][ok, center, row], pixels["lat"][ok, center, row])
            - position[:, np.newaxis]
        )[:, near]
        for row in (0, 3)
    )
    pitch = np.degrees(np.arccos(np.clip(np.sum(first * last, axis=-1), -1, 1))) / 15
    print(
        f"zones from {ok.sum()} scans; bounds {ZONE_BOUNDS.tolist()} (runs of "
        f"constant spacing: {runs.tolist()}); detector pitch {np.median(pitch):.5f} "
        f"deg, along_track_field_of_view {16 * np.median(pitch):.5f}; outermost "
        f"zone edges {angle[-1] + step[-1] / 2:.4f} deg left, "
        f"{-angle[0] + step[0] / 2:.4f} deg right"
    )
    for k, (a, b) in enumerate(zip(ZONE_BOUNDS[:-1], ZONE_BOUNDS[1:])):
        spacing = float(np.median(step[a : min(b, N_SAMPLES - 1)]))
        print(
            f"  zone {k}: {b - a} pixels, cross_track_field_of_view "
            f"{(b - a) * spacing:.6f} (spacing {spacing:.7f} deg), roll_angle "
            f"{(angle[a] + angle[b - 1]) / 2:.3f}; scan angles {angle[a]:.4f} to "
            f"{angle[b - 1]:.4f}"
        )


def validate_zones(setup: scanner.Setup, run, reference: Reference, orbit, tle):
    """
    Scores the sample size configuration of a run, and writes its results if the run
    writes: as a transfer, on a platform on which the configuration was not tested.
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
    names = np.array([instrument.name for instrument in satellite.instruments])
    pixels = measure_pixels(reference, satellite)
    median, p95 = (summarize(pixels[:, 3])[key] for key in ("median", "p95"))
    print(f"pixel position: median {median:.3f} km, 95th percentile {p95:.3f} km")
    times, rows = measure_spacing(reference)
    along, cross = compute_spacing_ratios(
        satellite, names[rows[:, 4].astype(int)], times, rows
    )
    print(f"spacing, predicted / measured: {along:.4f} along, {cross:.4f} across")
    if not run.write:
        return
    fixture = fixture_metrics = None
    if not run.transfer:
        rng = np.random.default_rng(3)
        picks = rng.choice(len(rows), setup.fixture_spacing, replace=False)
        fixture_along, fixture_cross = compute_spacing_ratios(
            satellite, names[rows[picks, 4].astype(int)], times[picks], rows[picks]
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
            "pixel_position_median": summarize(reduced[:, 3])["median"],
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


def main() -> None:
    """Validates the entry on a day of data, or fits its values."""
    run = parse_run("viirs", __doc__)
    archive = ARCHIVES[run.platform]
    directory = (
        run.directory
        if run.platform == FIT_PLATFORM
        else run.directory / archive.prefix
    )
    directory.mkdir(parents=True, exist_ok=True)
    reference = Reference(*load_day(run.day, directory, archive))
    if run.fit:
        fit(run.day, reference)
        return
    setup = get_setup(
        FIT_PLATFORM if run.transfer else run.platform, reference.scan_period
    )
    tle = run.element_set(directory=directory)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    print(
        f"{run.day:%Y-%m-%d}: {tle[0][:32]}; {len(reference.time)} scans, "
        f"{len(reference.gaps)} gaps over three minutes; configurations "
        f"{', '.join(run.configurations)}",
        flush=True,
    )
    if setup.beams in run.configurations:
        validate_zones(setup, run, reference, orbit, tle)
    if setup.swath in run.configurations:
        # the tests of the swath's passes also need the day's scan period
        scanner.validate_swath(
            setup,
            run,
            reference,
            orbit,
            tle,
            fixture_extras={"scan_period": np.array(reference.scan_period)},
        )


if __name__ == "__main__":
    main()

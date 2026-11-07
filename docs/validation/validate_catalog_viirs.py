"""
Validates the VIIRS catalog entry against the NOAA-20 VIIRS moderate-resolution
geolocation product, using the exact values of the entry through public TAT-C
methods. The tests follow `ValidateImagerVIIRS`, whose day of data the entry's angles
were fitted to, and `ValidateSampleSizes`: running on that day checks this script
against the notebooks, while running on a day at least 30 days from it tests the
fitted angles on data held out from the fit.

The notebooks do not test the swath's area, which the standard requires. The observed
area over one orbit is the union of the outlines (G-rings) of the granules and of the
strips between consecutive granules, on whose outermost pixel centers the outlines lie.

Reads anonymous public data, and needs the `validation` extra:

    python docs/validation/validate_catalog_viirs.py <data directory> 2026-10-02
    python docs/validation/validate_catalog_viirs.py <data directory> 2026-08-12

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import concurrent.futures
import json
import pathlib
from datetime import datetime, timezone

import h5py
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import tatc
from tatc.catalog import get_entry
from tatc.catalog.metrics import summarize
from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import GeneralPerturbationsOrbit
from validation_io import (
    BUCKET,
    TO_GEODETIC,
    RemoteFile,
    build_grid,
    compute_published_errors,
    compute_spacing_ratios,
    ecef,
    get_holdout_gap,
    iet_to_utc,
    list_keys,
    load_element_set,
    normalize,
    recompute_passes,
    score_passes,
    score_published,
    score_swath_area,
    up,
    update_fixture,
    update_record,
)

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
RECORD_PATH = REPOSITORY / "src" / "tatc" / "resources" / "validation" / "viirs.json"
FIXTURE_PATH = REPOSITORY / "tests" / "catalog" / "fixtures" / "viirs.npz"

FIT_DAY = datetime(2026, 10, 2, tzinfo=timezone.utc)
GEO = "All_Data/VIIRS-MOD-GEO-TC_All/"
GRANULE = "Data_Products/VIIRS-MOD-GEO-TC/VIIRS-MOD-GEO-TC_Gran_0"
N_SCANS, N_DETECTORS, N_SAMPLES = 48, 16, 3200
PIXEL_GRANULES = 10
COLUMNS = np.unique(np.r_[np.arange(0, N_SAMPLES, 8), N_SAMPLES - 1])
ROWS = np.array([0, 7, 8, 15])
SWATH, ZONES = "swath", "zones"
# the first sample of each zone, from right to left, and the end of the last
ZONE_BOUNDS = np.array([0, 640, 1007, 2191, 2559, 3200])
ORBIT_MARGIN = 120
# published footprints: instrument, altitude (km), inclination (deg), scan angle
# (deg), and sizes along and across track (km), from NOAA/CIMSS for band M16
PUBLISHED_FOOTPRINTS = (
    ("VIIRS zone 2", 824, 98.7, 0, 0.742, 0.776),
    ("VIIRS zone 4", 824, 98.7, 56.0, 1.60, 1.58),
)
FIXTURE_PASSES, FIXTURE_PIXELS, FIXTURE_SPACING = 60, 200, 200


def open_granule(key: str) -> h5py.File:
    """Opens a remote granule, first getting the metadata at its start and end."""
    remote = RemoteFile(BUCKET + key)
    remote.fetch(0, remote.block_size)
    remote.fetch(
        (remote.size - 2**19) // remote.block_size * remote.block_size, remote.size
    )
    return h5py.File(remote)


def read_scans(key: str) -> dict:
    """Reads the scans and outline of a granule."""
    with open_granule(key) as f:
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


def read_pixels(key: str) -> dict:
    """Reads a sample of the pixels of a granule."""
    with open_granule(key) as f:
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


def load_day(day: datetime, directory: pathlib.Path) -> tuple[dict, dict]:
    """Downloads and reduces the geolocation of a day, or loads it from a cache."""
    scans_path = directory / f"viirs_gmtco_scans_{day:%Y%m%d}.npz"
    pixels_path = directory / f"viirs_gmtco_pixels_{day:%Y%m%d}.npz"
    if not (scans_path.exists() and pixels_path.exists()):
        keys = list_keys(f"VIIRS-MOD-GEO-TC/{day:%Y/%m/%d}/")
        with concurrent.futures.ThreadPoolExecutor(8) as pool:
            scans = list(pool.map(read_scans, keys))
            pixels = list(pool.map(read_pixels, keys[::PIXEL_GRANULES]))
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
    """A day of VIIRS scans and pixels, with what is derived from them."""

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

    def observe_passes(self, points) -> pd.DataFrame:
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


def score_swath_width(reference: Reference, instrument) -> float:
    """Gets the median absolute difference (km) from the width through the samples."""
    ok = reference.ok
    lat, lon, height = (reference.pixels[k][ok] for k in ("lat", "lon", "height"))
    center = (
        ecef(lon[:, 1], lat[:, 1], height[:, 1])
        + ecef(lon[:, 2], lat[:, 2], height[:, 2])
    ) / 2
    observed = np.sum(np.linalg.norm(np.diff(center, axis=1), axis=-1), axis=1)
    position = reference.scans["SCPosition"][reference.granules].astype(float)[ok]
    altitude = TO_GEODETIC.transform(*position.T)[2]
    predicted = np.vectorize(instrument.get_swath_width)(altitude)
    return float(np.median(np.abs(predicted - observed))) / 1e3


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


def validate_zones(day, tle, entry, orbit, reference, write: bool) -> None:
    """Scores the sample size configuration, and writes its results if requested."""
    within, largest = score_published(
        compute_published_errors(entry, ZONES, PUBLISHED_FOOTPRINTS)
    )
    print(
        f"published dimensions: {within:.1f}% within 3%, largest error {largest:.2f}%"
    )
    satellite = entry.to_satellite(orbit, ZONES, platform="noaa-20")
    names = np.array([instrument.name for instrument in satellite.instruments])
    pixels = measure_pixels(reference, satellite)
    median, p95 = (summarize(pixels[:, 3])[key] for key in ("median", "p95"))
    print(f"pixel position: median {median:.3f} km, 95th percentile {p95:.3f} km")
    times, rows = measure_spacing(reference)
    along, cross = compute_spacing_ratios(
        satellite, names[rows[:, 4].astype(int)], times, rows
    )
    print(f"spacing, predicted / measured: {along:.4f} along, {cross:.4f} across")
    if not write:
        return
    rng = np.random.default_rng(3)
    picks = rng.choice(len(rows), FIXTURE_SPACING, replace=False)
    fixture_along, fixture_cross = compute_spacing_ratios(
        satellite, names[rows[picks, 4].astype(int)], times[picks], rows[picks]
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
        "viirs",
        ZONES,
        {
            "digest": entry.configurations[ZONES].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": "docs/validation/validate_catalog_viirs.py",
            "reference": "Footprints published by NOAA/CIMSS; samples of the VIIRS "
            "moderate-resolution geolocation",
            "period": f"{day:%Y-%m-%d}",
            "platform": "noaa-20",
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
                "pixel_position_median": summarize(reduced[:, 3])["median"],
                "spacing_along_ratio_median": fixture_along,
                "spacing_cross_ratio_median": fixture_cross,
            },
        },
    )


def validate_swath(day, tle, entry, orbit, reference, write: bool) -> None:
    """Scores the swath configuration, and writes its results if requested."""

    def satellite_for(time_step: float):
        return entry.to_satellite(orbit, SWATH, time_step=time_step, platform="noaa-20")

    satellite = satellite_for(reference.scan_period)
    points = generate_points_fibonacci_lattice(500e3)
    observed = reference.observe_passes(points)
    start, end = reference.time[0].to_pydatetime(), reference.time[-1].to_pydatetime()
    results = score_passes(points, satellite, 0, observed, start, end, reference.gaps)
    results["swath_width_median"] = score_swath_width(
        reference, satellite.instruments[0]
    )
    print(json.dumps(results), flush=True)
    grid = build_grid()
    first, last, compared, covered = observe_orbit(reference, orbit, grid)
    results["footprint_iou_median"] = score_swath_area(
        satellite_for,
        0,
        first,
        last,
        (reference.scan_period, 30, 120),
        grid,
        compared,
        covered,
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
        FIXTURE_PATH,
        tle=np.array(tle),
        day=np.array(day.isoformat()),
        passes=passes,
        scan_period=np.array(reference.scan_period),
    )
    update_record(
        RECORD_PATH,
        "viirs",
        SWATH,
        {
            "digest": entry.configurations[SWATH].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": "docs/validation/validate_catalog_viirs.py",
            "reference": "VIIRS moderate-resolution geolocation",
            "period": f"{day:%Y-%m-%d}",
            "platform": "noaa-20",
            "holdout_gap_days": get_holdout_gap(day, FIT_DAY),
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


def main() -> None:
    """Validates the entry on a day of data."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("day")
    parser.add_argument("--sample-sizes-only", action="store_true")
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)

    day = datetime.fromisoformat(arguments.day).replace(tzinfo=timezone.utc)
    tle = load_element_set(day, arguments.directory)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    entry = get_entry("viirs")
    reference = Reference(*load_day(day, arguments.directory))
    print(
        f"{arguments.day}: {tle[0][:32]}; {len(reference.time)} scans, "
        f"{len(reference.gaps)} gaps over three minutes",
        flush=True,
    )
    validate_zones(day, tle, entry, orbit, reference, arguments.write)
    if not arguments.sample_sizes_only:
        validate_swath(day, tle, entry, orbit, reference, arguments.write)
    if arguments.write:
        print(f"wrote {RECORD_PATH} and {FIXTURE_PATH}")


if __name__ == "__main__":
    main()

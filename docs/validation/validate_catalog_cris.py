"""
Validates the CrIS catalog entry against the NOAA-20 CrIS geolocation product, using
the exact values of the entry through public TAT-C methods. The tests follow
`ValidateImagerCrIS`, whose day of data the entry's angles were fitted to, and
`ValidateSampleSizes`: running on that day checks this script against the notebooks,
while running on a day at least 30 days from it tests the fitted angles on data held
out from the fit.

Reads anonymous public data, and needs the `validation` extra:

    python docs/validation/validate_catalog_cris.py <data directory> 2026-10-02
    python docs/validation/validate_catalog_cris.py <data directory> 2026-08-12

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import concurrent.futures
import io
import json
import pathlib
from datetime import datetime, timedelta, timezone

import h5py
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from skyfield.api import wgs84

import tatc
from tatc.catalog import get_entry
from tatc.catalog.metrics import summarize
from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import GeneralPerturbationsOrbit
from validation_io import (
    build_grid,
    compute_published_errors,
    ecef,
    fetch,
    get_holdout_gap,
    iet_to_utc,
    in_lattice_cells,
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
RECORD_PATH = REPOSITORY / "src" / "tatc" / "resources" / "validation" / "cris.json"
FIXTURE_PATH = REPOSITORY / "tests" / "catalog" / "fixtures" / "cris.npz"

FIT_DAY = datetime(2026, 10, 2, tzinfo=timezone.utc)
FIELDS = (
    "StartTime",
    "MidTime",
    "FORTime",
    "Latitude",
    "Longitude",
    "Height",
    "QF1_CRISSDRGEO",
)
N_FOR, N_FOV, SCAN_PERIOD = 30, 9, 8
TIME_STEPS = (SCAN_PERIOD, 32, 128)
ORBIT_MARGIN = 128
SWATH, PIXELS = "swath", "pixels"
# published footprints: instrument, altitude (km), inclination (deg), scan angle
# (deg), and sizes along and across track (km), from the NWP SAF, as used by
# ValidateSampleSizes (the edge is the outermost field of view, at 48.33 + 1.1 deg)
PUBLISHED_FOOTPRINTS = (
    ("CrIS", 824, 98.7, 0, 14.0, 14.0),
    ("CrIS", 824, 98.7, 49.43, 24.0, 48.0),
)
PIXEL_SCAN_STRIDE, WIDTH_SCAN_STRIDE = 100, 5
FIXTURE_PASSES, FIXTURE_PIXELS = 60, 200
NADIR_FORS = (14, 15)


def load_day(day: datetime, directory: pathlib.Path) -> dict[str, np.ndarray]:
    """Downloads and reduces the geolocation of a day, or loads it from a cache."""
    path = directory / f"cris_gcrso_{day:%Y%m%d}.npz"
    if not path.exists():

        def read(key):
            with h5py.File(io.BytesIO(fetch(key))) as f:
                return {k: f["All_Data/CrIS-SDR-GEO_All/" + k][:] for k in FIELDS}

        keys = list_keys(f"CrIS-SDR-GEO/{day:%Y/%m/%d}/")
        with concurrent.futures.ThreadPoolExecutor(16) as pool:
            granules = list(pool.map(read, keys))
        np.savez_compressed(
            path, **{k: np.concatenate([g[k] for g in granules]) for k in FIELDS}
        )
    return dict(np.load(path))


class Reference:
    """A day of CrIS scans, with what is derived from them."""

    def __init__(self, data: dict[str, np.ndarray]):
        good = (data["QF1_CRISSDRGEO"] == 0) & np.all(
            np.abs(data["Latitude"]) <= 90, axis=(1, 2)
        )
        self.time = iet_to_utc(data["MidTime"][good])
        self.start = iet_to_utc(data["StartTime"][good])
        self.fov_lat = data["Latitude"][good].astype(float)
        self.fov_lon = data["Longitude"][good].astype(float)
        self.fov_height = data["Height"][good].astype(float)
        # the center of each field of regard, the mean of its fields of view
        self.unit = normalize(up(self.fov_lon, self.fov_lat).mean(axis=2))
        self.lat = np.degrees(np.arcsin(self.unit[..., 2]))
        self.lon = np.degrees(np.arctan2(self.unit[..., 1], self.unit[..., 0]))
        self.height = self.fov_height.mean(axis=2)
        day = self.time[0].floor("D")
        self.seconds = (self.time - day).total_seconds().to_numpy()
        self.for_seconds = (
            (iet_to_utc(data["FORTime"][good].ravel()) - day)
            .total_seconds()
            .to_numpy()
            .reshape(-1, N_FOR)
        )
        period = np.diff(self.seconds)
        self.is_next_contiguous = period < 1.5 * SCAN_PERIOD
        self.gaps = [
            (self.time[i], self.time[i + 1]) for i in np.nonzero(period > 60)[0]
        ]

    def in_cells(self, points: np.ndarray, first: int, last: int) -> np.ndarray:
        """Gets whether unit vectors lie in cells of fields of regard on scans."""
        return in_lattice_cells(
            self.unit[first : last + 1], self.is_next_contiguous[first:last], points
        )

    def observe_passes(self, points) -> pd.DataFrame:
        """
        Finds the passes on which CrIS observed points, between the centers of its
        fields of regard, and the times at which it sampled them.
        """
        point_unit = up(points.geometry.x.to_numpy(), points.geometry.y.to_numpy())
        flat = self.unit.reshape(-1, 3)
        flat_seconds = self.for_seconds.ravel()
        tree = cKDTree(flat)
        day = self.time[0].floor("D")
        rows = []
        for p, hits in enumerate(
            tree.query_ball_point(point_unit, r=2 * np.sin(150 / 6371 / 2), workers=-1)
        ):
            if not hits:
                continue
            hits = np.array(hits)
            hits = hits[np.argsort(flat_seconds[hits])]
            for group in np.split(
                hits, np.nonzero(np.diff(flat_seconds[hits]) > 1200)[0] + 1
            ):
                k, j = np.divmod(
                    group[
                        np.argmin(np.linalg.norm(flat[group] - point_unit[p], axis=-1))
                    ],
                    N_FOR,
                )
                first, last = max(k - 1, 0), min(k + 1, len(self.time) - 1)
                if not self.in_cells(point_unit[[p]], first, last)[0]:
                    continue
                # interpolate between the scans before and after the point
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
                seconds = self.for_seconds[before, j] + fraction * (
                    self.for_seconds[before + 1, j] - self.for_seconds[before, j]
                )
                rows.append(
                    {"point_id": p, "time": day + timedelta(seconds=float(seconds))}
                )
        return pd.DataFrame(rows).astype({"point_id": "int64"})


def score_swath_width(reference: Reference, instrument, orbit) -> float:
    """
    Gets the median absolute difference (km) from the width through the centers of
    the fields of regard.
    """
    scans = np.arange(0, len(reference.time), WIDTH_SCAN_STRIDE)
    xyz = ecef(reference.lon[scans], reference.lat[scans], reference.height[scans])
    observed = np.sum(np.linalg.norm(np.diff(xyz, axis=1), axis=-1), axis=1)
    altitude = wgs84.geographic_position_of(
        orbit.get_orbit_track(list(reference.time[scans]))
    ).elevation.m
    predicted = np.vectorize(instrument.get_swath_width)(altitude)
    return float(np.median(np.abs(predicted - observed))) / 1e3


def observe_orbit(reference: Reference, orbit, grid):
    """Gets the cells compared over one orbit, and which of them CrIS observed."""
    lon, lat, _ = grid
    unit = up(lon, lat)
    day = reference.time[0].floor("D")
    begin = day + pd.Timedelta(hours=12)
    period = pd.Timedelta(days=1 / orbit.elements[0].mean_motion * 360 / 86400)
    scans = np.nonzero((reference.time >= begin) & (reference.time < begin + period))[0]
    first, last = reference.time[scans[0]], reference.time[scans[-1]]
    distance, nearest = cKDTree(reference.unit[scans].reshape(-1, 3)).query(
        unit, workers=-1
    )
    elapsed = (reference.time[scans] - first).total_seconds().to_numpy()
    nearest_time = elapsed[nearest // N_FOR]
    compared = (
        (distance < 2 * np.sin(300 / 6371 / 2))
        & (nearest_time > ORBIT_MARGIN)
        & (nearest_time < (last - first).total_seconds() - ORBIT_MARGIN)
    )
    observed = np.zeros(len(unit), dtype=bool)
    observed[compared] = reference.in_cells(unit[compared], scans[0], scans[-1])
    return first.to_pydatetime(), last.to_pydatetime(), compared, observed


def match_pixels(reference: Reference, instrument, orbit) -> tuple[np.ndarray, ...]:
    """
    Gets the TAT-C pixel (cross-track and along-track index) of each field of view
    of each field of regard. The order of the fields of regard is that in which the
    scan crosses the track, and the pixels of the array are those closest to the
    fields of view of the fields of regard nearest nadir, where the array is not
    turned, so that matching does not hide the errors elsewhere in the scan.
    """
    scans = np.arange(0, len(reference.time), PIXEL_SCAN_STRIDE)
    track = orbit.get_orbit_track(list(reference.time[scans]))

    def position(cross: int, along: int) -> np.ndarray:
        p = instrument.compute_projected_pixel_position(track, cross, along)
        return ecef(p.longitude.degrees, p.latitude.degrees)

    def distance(cross: int, along: int, j: int, m: int) -> float:
        measured = ecef(reference.fov_lon[scans, j, m], reference.fov_lat[scans, j, m])
        return float(
            np.median(np.linalg.norm(position(cross, along) - measured, axis=-1))
        )

    # TAT-C numbers pixels from the right of the motion
    first_on_left = distance(3 * (N_FOR - 1) + 1, 1, 0, 4) < distance(1, 1, 0, 4)
    for_position = N_FOR - 1 - np.arange(N_FOR) if first_on_left else np.arange(N_FOR)
    cross_offset, along_index = np.zeros(N_FOV, int), np.zeros(N_FOV, int)
    for m in range(N_FOV):
        candidates = [(c, a) for c in range(3) for a in range(3)]
        scores = [
            sum(distance(3 * for_position[j] + c, a, j, m) for j in NADIR_FORS)
            for c, a in candidates
        ]
        cross_offset[m], along_index[m] = candidates[int(np.argmin(scores))]
    if len(set(zip(cross_offset, along_index))) != N_FOV:
        raise RuntimeError("fields of view do not match distinct pixels")
    cross_index = 3 * for_position[:, np.newaxis] + cross_offset[np.newaxis, :]
    return cross_index, np.broadcast_to(along_index, cross_index.shape)


def measure_pixels(reference: Reference, satellite) -> np.ndarray:
    """
    Gets rows of seconds from the start of the day, cross-track and along-track
    pixel index, position difference (km) from the field of view, and the field of
    view's latitude and longitude, for the pixels projected at the middle times of a
    sample of scans.
    """
    instrument = satellite.get_instrument(0)
    orbit = satellite.orbit.to_gp_orbit()
    cross_index, along_index = match_pixels(reference, instrument, orbit)
    scans = np.arange(0, len(reference.time), PIXEL_SCAN_STRIDE)
    track = orbit.get_orbit_track(list(reference.time[scans]))
    rows = []
    for j in range(N_FOR):
        for m in range(N_FOV):
            position = instrument.compute_projected_pixel_position(
                track, int(cross_index[j, m]), int(along_index[j, m])
            )
            lat, lon = reference.fov_lat[scans, j, m], reference.fov_lon[scans, j, m]
            difference = np.linalg.norm(
                ecef(position.longitude.degrees, position.latitude.degrees)
                - ecef(lon, lat),
                axis=-1,
            )
            rows.append(
                np.stack(
                    [
                        reference.seconds[scans],
                        np.full(len(scans), cross_index[j, m]),
                        np.full(len(scans), along_index[j, m]),
                        difference / 1e3,
                        lat,
                        lon,
                    ],
                    axis=-1,
                )
            )
    return np.concatenate(rows)


def validate_pixels(day, tle, entry, orbit, reference, write: bool) -> None:
    """Scores the sample size configuration, and writes its results if requested."""
    within, largest = score_published(
        compute_published_errors(entry, PIXELS, PUBLISHED_FOOTPRINTS)
    )
    print(
        f"published dimensions: {within:.1f}% within 3%, largest error {largest:.2f}%"
    )
    satellite = entry.to_satellite(orbit, PIXELS, platform="noaa-20")
    pixels = measure_pixels(reference, satellite)
    median = summarize(pixels[:, 3])["median"]
    p95 = summarize(pixels[:, 3])["p95"]
    central = pixels[(pixels[:, 1] % 3 == 1) & (pixels[:, 2] == 1)]
    print(
        f"pixel position: median {median:.3f} km, 95th percentile {p95:.3f} km; "
        f"central fields of view: median {np.median(central[:, 3]):.3f} km, "
        f"95th percentile {summarize(central[:, 3])['p95']:.3f} km"
    )
    if not write:
        return
    rng = np.random.default_rng(3)
    reduced = pixels[rng.choice(len(pixels), FIXTURE_PIXELS, replace=False)]
    update_fixture(
        FIXTURE_PATH, tle=np.array(tle), day=np.array(day.isoformat()), pixels=reduced
    )
    update_record(
        RECORD_PATH,
        "cris",
        PIXELS,
        {
            "digest": entry.configurations[PIXELS].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": "docs/validation/validate_catalog_cris.py",
            "reference": "Footprints published by the NWP SAF; fields of view of the "
            "CrIS SDR geolocation",
            "period": f"{day:%Y-%m-%d}",
            "platform": "noaa-20",
            "holdout_gap_days": get_holdout_gap(day, FIT_DAY),
            "metrics": [
                {"name": "sample_size_within_3pct", "value": within},
                {"name": "sample_size_largest_error", "value": largest},
                {"name": "pixel_position_median", "value": median},
                {"name": "pixel_position_p95", "value": p95},
                {
                    "name": "pixel_position_central_median",
                    "value": float(np.median(central[:, 3])),
                },
            ],
            "fixture_metrics": {
                "sample_size_within_3pct": within,
                "pixel_position_median": summarize(reduced[:, 3])["median"],
            },
        },
    )


def validate_swath(day, tle, entry, orbit, reference, write: bool) -> None:
    """Scores the swath configuration, and writes its results if requested."""

    def satellite_for(time_step: float):
        return entry.to_satellite(orbit, SWATH, time_step=time_step, platform="noaa-20")

    satellite = satellite_for(SCAN_PERIOD)
    points = generate_points_fibonacci_lattice(500e3)
    observed = reference.observe_passes(points)
    start = reference.start[0].to_pydatetime()
    end = reference.start[-1].to_pydatetime() + timedelta(seconds=SCAN_PERIOD)
    results = score_passes(points, satellite, 0, observed, start, end, reference.gaps)
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
    update_record(
        RECORD_PATH,
        "cris",
        SWATH,
        {
            "digest": entry.configurations[SWATH].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": "docs/validation/validate_catalog_cris.py",
            "reference": "CrIS SDR geolocation",
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
    entry = get_entry("cris")
    reference = Reference(load_day(day, arguments.directory))
    print(
        f"{arguments.day}: {tle[0][:32]}; {len(reference.time)} scans, "
        f"{len(reference.gaps)} gaps over a minute",
        flush=True,
    )
    validate_pixels(day, tle, entry, orbit, reference, arguments.write)
    if not arguments.sample_sizes_only:
        validate_swath(day, tle, entry, orbit, reference, arguments.write)
    if arguments.write:
        print(f"wrote {RECORD_PATH} and {FIXTURE_PATH}")


if __name__ == "__main__":
    main()

"""
Validates the ATMS catalog entry against the NOAA-20 ATMS geolocation product, using
the exact values of the entry through public TAT-C methods. The tests follow
`ValidateImagerATMS`, whose day of data the entry's angles were fitted to, and
`ValidateSampleSizes`: running on that day checks this script against the notebooks,
while running on a day at least 30 days from it tests the fitted angles on data held
out from the fit.

Reads anonymous public data, and needs the `validation` extra:

    python docs/validation/validate_catalog_atms.py <data directory> 2026-10-02
    python docs/validation/validate_catalog_atms.py <data directory> 2026-08-12

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
    compute_spacing_ratios,
    ecef,
    fetch,
    get_holdout_gap,
    iet_to_utc,
    list_keys,
    load_element_set,
    recompute_passes,
    score_passes,
    score_published,
    score_swath_area,
    up,
    update_fixture,
    update_record,
)

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
RECORD_PATH = REPOSITORY / "src" / "tatc" / "resources" / "validation" / "atms.json"
FIXTURE_PATH = REPOSITORY / "tests" / "catalog" / "fixtures" / "atms.npz"

FIT_DAY = datetime(2026, 10, 2, tzinfo=timezone.utc)
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
TIME_STEPS = (SCAN_PERIOD, 32, 128)
ORBIT_MARGIN = 128
SWATH, BEAMS = "swath", "beams"
SIZE_INSTRUMENT = "ATMS channels 3-16"
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
PIXEL_SCAN_STRIDE, SPACING_SCAN_STRIDE, WIDTH_SCAN_STRIDE = 100, 150, 5
FIXTURE_PASSES, FIXTURE_PIXELS, FIXTURE_SPACING = 60, 200, 200


def load_day(day: datetime, directory: pathlib.Path) -> dict[str, np.ndarray]:
    """Downloads and reduces the geolocation of a day, or loads it from a cache."""
    path = directory / f"atms_gatmo_{day:%Y%m%d}.npz"
    if not path.exists():

        def read(key):
            with h5py.File(io.BytesIO(fetch(key))) as f:
                return {k: f["All_Data/ATMS-SDR-GEO_All/" + k][:] for k in FIELDS}

        keys = list_keys(f"ATMS-SDR-GEO/{day:%Y/%m/%d}/")
        with concurrent.futures.ThreadPoolExecutor(8) as pool:
            granules = list(pool.map(read, keys))
        np.savez_compressed(
            path, **{k: np.concatenate([g[k] for g in granules]) for k in FIELDS}
        )
    return dict(np.load(path))


class Reference:
    """A day of ATMS scans, with what is derived from them."""

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

    def in_cells(self, points: np.ndarray, first: int, last: int) -> np.ndarray:
        """Gets whether unit vectors lie in cells of beams on a range of scans."""
        cells = self.unit[first : last + 1]
        _, nearest = cKDTree(cells.reshape(-1, 3)).query(points, k=2, workers=-1)
        inside = np.zeros(len(points), dtype=bool)
        for column in range(2):
            k, j = np.divmod(nearest[:, column], N_BEAMS)
            for dk in (-1, 0):
                for dj in (-1, 0):
                    kk, jj = k + dk, j + dj
                    valid = (
                        (kk >= 0)
                        & (kk < len(cells) - 1)
                        & (jj >= 0)
                        & (jj < N_BEAMS - 1)
                    )
                    kk = np.clip(kk, 0, len(cells) - 2)
                    jj = np.clip(jj, 0, N_BEAMS - 2)
                    valid &= self.is_next_contiguous[first + kk]
                    corners = [
                        cells[kk, jj],
                        cells[kk, jj + 1],
                        cells[kk + 1, jj + 1],
                        cells[kk + 1, jj],
                    ]
                    signs = np.stack(
                        [
                            np.sum(np.cross(a, b) * points, axis=-1)
                            for a, b in zip(corners, corners[1:] + corners[:1])
                        ]
                    )
                    inside |= valid & (
                        np.all(signs >= 0, axis=0) | np.all(signs <= 0, axis=0)
                    )
        return inside

    def observe_passes(self, points) -> pd.DataFrame:
        """Finds the passes on which ATMS observed points, and their times."""
        point_unit = up(points.geometry.x.to_numpy(), points.geometry.y.to_numpy())
        flat = self.unit.reshape(-1, 3)
        beam_seconds = (
            self.seconds[:, np.newaxis] + self.beam_offsets[np.newaxis, :]
        ).ravel()
        tree = cKDTree(flat)
        rows = []
        for p, hits in enumerate(
            tree.query_ball_point(point_unit, r=2 * np.sin(60 / 6371 / 2), workers=-1)
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
                    N_BEAMS,
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
                seconds = (
                    self.seconds[before]
                    + self.beam_offsets[j]
                    + fraction * (self.seconds[before + 1] - self.seconds[before])
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
        regular = np.abs(period - SCAN_PERIOD) < 0.01
        return np.nonzero(regular[:-1] & regular[1:])[0] + 1


def score_swath_width(reference: Reference, instrument, orbit) -> float:
    """Gets the median absolute difference (km) from the width through the beams."""
    scans = np.arange(0, len(reference.time), WIDTH_SCAN_STRIDE)
    xyz = ecef(reference.lon[scans], reference.lat[scans], reference.height[scans])
    observed = np.sum(np.linalg.norm(np.diff(xyz, axis=1), axis=-1), axis=1)
    altitude = wgs84.geographic_position_of(
        orbit.get_orbit_track(list(reference.time[scans]))
    ).elevation.m
    predicted = np.vectorize(instrument.get_swath_width)(altitude)
    return float(np.median(np.abs(predicted - observed))) / 1e3


def observe_orbit(reference: Reference, orbit, grid):
    """Gets the cells compared over one orbit, and which of them ATMS observed."""
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
    nearest_time = elapsed[nearest // N_BEAMS]
    compared = (
        (distance < 2 * np.sin(300 / 6371 / 2))
        & (nearest_time > ORBIT_MARGIN)
        & (nearest_time < (last - first).total_seconds() - ORBIT_MARGIN)
    )
    observed = np.zeros(len(unit), dtype=bool)
    observed[compared] = reference.in_cells(unit[compared], scans[0], scans[-1])
    return first.to_pydatetime(), last.to_pydatetime(), compared, observed


def measure_pixels(reference: Reference, satellite, index: int) -> np.ndarray:
    """
    Gets rows of seconds from the start of the day, pixel index, position difference
    (km) from the beam center, and the beam's latitude and longitude, for the pixels
    projected at the middle times of a sample of scans.
    """
    scans = np.arange(0, len(reference.time), PIXEL_SCAN_STRIDE)
    instrument = satellite.get_instrument(index)
    times = reference.time[scans]
    track = satellite.orbit.to_gp_orbit().get_orbit_track(list(times))
    seconds = (times - times[0].floor("D")).total_seconds().to_numpy()
    rows = []
    for beam in range(N_BEAMS):
        # TAT-C numbers pixels from the right of the motion, and ATMS from the left
        pixel = N_BEAMS - 1 - beam
        position = instrument.compute_projected_pixel_position(track, pixel, 0)
        difference = np.linalg.norm(
            ecef(position.longitude.degrees, position.latitude.degrees)
            - ecef(reference.lon[scans, beam], reference.lat[scans, beam]),
            axis=-1,
        )
        rows.append(
            np.stack(
                [
                    seconds,
                    np.full(len(scans), pixel),
                    difference / 1e3,
                    reference.lat[scans, beam],
                    reference.lon[scans, beam],
                ],
                axis=-1,
            )
        )
    return np.concatenate(rows)


def measure_spacing(reference: Reference) -> tuple[pd.DatetimeIndex, np.ndarray]:
    """
    Measures the spacing of beams as `ValidateSampleSizes` does, from half the distance
    between the beams on either side, with each beam at the time it is sampled. Gets
    the times, and rows of latitude, longitude, along-track and cross-track spacing (km).
    """
    scans = reference.contiguous_scans()[::SPACING_SCAN_STRIDE]
    xyz = ecef(reference.lon, reference.lat)
    along = np.linalg.norm(xyz[scans + 1] - xyz[scans - 1], axis=-1) / 2e3
    cross = np.full(along.shape, np.nan)
    cross[:, 1:-1] = np.linalg.norm(xyz[scans, 2:] - xyz[scans, :-2], axis=-1) / 2e3
    times = iet_to_utc(
        (
            reference.mid_iet[scans, np.newaxis]
            + 1e6 * reference.beam_offsets[np.newaxis, :]
        ).ravel()
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


def validate_beams(day, tle, entry, orbit, reference, write: bool) -> None:
    """Scores the sample size configuration, and writes its results if requested."""
    within, largest = score_published(
        compute_published_errors(entry, BEAMS, PUBLISHED_FOOTPRINTS)
    )
    print(
        f"published dimensions: {within:.1f}% within 3%, largest error {largest:.2f}%"
    )
    satellite = entry.to_satellite(orbit, BEAMS, platform="noaa-20")
    names = [instrument.name for instrument in satellite.instruments]
    index = names.index(SIZE_INSTRUMENT)
    pixels = measure_pixels(reference, satellite, index)
    median = summarize(pixels[:, 2])["median"]
    p95 = summarize(pixels[:, 2])["p95"]
    print(f"pixel position: median {median:.3f} km, 95th percentile {p95:.3f} km")
    times, rows = measure_spacing(reference)
    along, cross = compute_spacing_ratios(
        satellite, [SIZE_INSTRUMENT] * len(rows), times, rows
    )
    print(f"spacing, predicted / measured: {along:.4f} along, {cross:.4f} across")
    if not write:
        return
    rng = np.random.default_rng(3)
    picks = rng.choice(len(rows), FIXTURE_SPACING, replace=False)
    fixture_along, fixture_cross = compute_spacing_ratios(
        satellite, [SIZE_INSTRUMENT] * len(picks), times[picks], rows[picks]
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
        "atms",
        BEAMS,
        {
            "digest": entry.configurations[BEAMS].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": "docs/validation/validate_catalog_atms.py",
            "reference": "Footprints of the TC PRIMED documentation; beam centers of "
            "the ATMS SDR geolocation",
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
                "pixel_position_median": summarize(reduced[:, 2])["median"],
                "spacing_along_ratio_median": fixture_along,
                "spacing_cross_ratio_median": fixture_cross,
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
        "atms",
        SWATH,
        {
            "digest": entry.configurations[SWATH].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": "docs/validation/validate_catalog_atms.py",
            "reference": "ATMS SDR geolocation",
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
    entry = get_entry("atms")
    reference = Reference(load_day(day, arguments.directory))
    print(
        f"{arguments.day}: {tle[0][:32]}; {len(reference.time)} scans, "
        f"{len(reference.gaps)} gaps over a minute",
        flush=True,
    )
    validate_beams(day, tle, entry, orbit, reference, arguments.write)
    if not arguments.sample_sizes_only:
        validate_swath(day, tle, entry, orbit, reference, arguments.write)
    if arguments.write:
        print(f"wrote {RECORD_PATH} and {FIXTURE_PATH}")


if __name__ == "__main__":
    main()

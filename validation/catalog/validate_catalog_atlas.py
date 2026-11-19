"""
Validates the ICESat-2 ATLAS catalog entry against the ATLAS/ICESat-2 Land and
Vegetation Height Quick Look product (ATL08QL), using the exact values of the entry
through public TAT-C methods. It follows `ValidateLidarICESat2`, over Antarctica
(south of 60 deg S), where ICESat-2 points its beams along their reference ground
tracks rather than off them.

The entry models each of the three pairs of beams as a narrow view. Their centers
and widths are fitted to the cross-track positions of the pairs' segments relative
to TAT-C's ground track on one day (run with `--fit`): each view is centered on the
median offset of its pair, and is as wide as twice the 99.5th percentile of the
offsets' distance from it. The widths include the error of the element set, so they
hold for element sets of a similar age (a few days). Running on a day at least 30
days from the fit tests the fitted values on data held out from the fit.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`:

    python validation/catalog/validate_catalog_atlas.py <data directory> 2026-10-01 --fit
    python validation/catalog/validate_catalog_atlas.py <data directory> 2026-10-01
    python validation/catalog/validate_catalog_atlas.py <data directory> 2026-11-04

Without a day, it takes that of the entry's record (or, with `--fit`, of the fit).
The element set of each day is registered in `elements/element_sets.yml`: public
element sets from CelesTrak, with epochs a few days after each day of data (that of
2026-10-01 is the one of `ValidateLidarICESat2`).

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results. A record of the day of the fit has no
holdout, which the entry's validation reports as a gap.

With `--same-day-element-set`, it tests instead with the element set nearest the middle
of the day from elements/icesat2_2026-11-04.3le, for comparison with the views'
intended use with an element set a few days old; it does not write a record:

    python validation/catalog/validate_catalog_atlas.py <data directory> 2026-11-04 --same-day-element-set

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import concurrent.futures
import json
import pathlib
from datetime import datetime, timedelta, timezone

import earthaccess
import h5py
import numpy as np
import pandas as pd
from pyproj import Geod
from skyfield.api import wgs84

from tatc.analysis import collect_observations
from tatc.catalog.metrics import summarize
from tatc.constants import EARTH_MEAN_RADIUS, timescale
from tatc.schemas import GeneralPerturbationsOrbit, Point
from validation_io import (
    ELEMENT_SET_FOLDER,
    RemoteFile,
    get_element_set,
    login_earthdata,
    parse_run,
    write_results,
)

# the file of element sets of `--same-day-element-set`
SAME_DAY_ELEMENT_SET_FILE = ELEMENT_SET_FOLDER / "icesat2_2026-11-04.3le"
CONFIGURATION = "beam_pairs"
PAIRS = ("gt1", "gt2", "gt3")
BEAMS = ("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r")
ALTITUDE = 481e3  # m, nominal, to convert between angles and ground distances
LATITUDE_LIMIT = -60
GPS_EPOCH = datetime(1980, 1, 6, tzinfo=timezone.utc)
GPS_UTC = 18  # s
SEGMENT_STRIDE = 10
SAMPLE_PER_PAIR = 600
TIME_STEPS = (1, 10, 30)
SPACING_BIN = 10  # s
FIT_PERCENTILE = 99.5
FIXTURE_PASSES = 60
GEOD = Geod(ellps="WGS84")


def load_day(day: datetime, directory: pathlib.Path) -> pd.DataFrame:
    """
    Reads every 10th land segment of each beam south of 60 deg S from the day's
    files, with range requests, or loads them from a cache.
    """
    cache = directory / f"atl08ql_{day:%Y%m%d}.csv"
    if not cache.exists():
        login_earthdata()
        session = earthaccess.get_requests_https_session()
        results = earthaccess.search_data(
            short_name="ATL08QL",
            version="007",
            temporal=(
                f"{day:%Y-%m-%d}T00:00:00Z",
                f"{day:%Y-%m-%d}T23:59:59Z",
            ),
        )
        urls = [
            next(
                link["URL"]
                for link in r["umm"]["RelatedUrls"]
                if link["URL"].startswith("https") and link["URL"].endswith(".h5")
            )
            for r in results
        ]

        def read(url):
            frames = []
            with h5py.File(RemoteFile(url, session), "r") as f:
                epoch = float(f["ancillary_data/atlas_sdp_gps_epoch"][0])
                for beam in BEAMS:
                    if f"{beam}/land_segments" not in f:
                        continue
                    group = f[f"{beam}/land_segments"]
                    latitude = group["latitude"][::SEGMENT_STRIDE]
                    if not np.any(latitude < LATITUDE_LIMIT):
                        continue
                    frames.append(
                        pd.DataFrame(
                            {
                                "file": url.rsplit("/", 1)[-1],
                                "beam": beam,
                                "latitude": latitude,
                                "longitude": group["longitude"][::SEGMENT_STRIDE],
                                "gps_seconds": epoch
                                + group["delta_time"][::SEGMENT_STRIDE],
                            }
                        )
                    )
            return frames

        with concurrent.futures.ThreadPoolExecutor(8) as pool:
            parts = [frame for frames in pool.map(read, urls) for frame in frames]
        segments = pd.concat(parts, ignore_index=True)
        segments = segments[segments.latitude < LATITUDE_LIMIT]
        segments.to_csv(cache, index=False)
    segments = pd.read_csv(cache)
    segments["time"] = pd.to_datetime(GPS_EPOCH) + pd.to_timedelta(
        segments.gps_seconds - GPS_UTC, unit="s"
    )
    segments = segments[
        (segments.time >= pd.Timestamp(day))
        & (segments.time < pd.Timestamp(day) + pd.Timedelta(days=1))
    ]
    return segments.assign(pair=segments.beam.str[:3]).reset_index(drop=True)


def angle(distance):
    """Off-nadir angle (deg) of a signed ground distance (m) from nadir."""
    theta = np.asarray(distance) / EARTH_MEAN_RADIUS
    return np.degrees(
        np.arctan2(
            EARTH_MEAN_RADIUS * np.sin(theta),
            EARTH_MEAN_RADIUS + ALTITUDE - EARTH_MEAN_RADIUS * np.cos(theta),
        )
    )


def distance(off_nadir):
    """Signed ground distance (m) from nadir of an off-nadir angle (deg)."""
    a = np.radians(np.asarray(off_nadir))
    theta = np.arcsin((EARTH_MEAN_RADIUS + ALTITUDE) / EARTH_MEAN_RADIUS * np.sin(a))
    return EARTH_MEAN_RADIUS * (theta - a)


def track_bearings(samples: pd.DataFrame, orbit) -> tuple[np.ndarray, np.ndarray]:
    """
    Gets the bearings (deg, clockwise from the heading of the ground track) and
    distances (m) of samples from the ground track's point at their times.
    """
    satellite = orbit.elements[0].to_skyfield()
    times = timescale.from_datetimes(list(samples.time.dt.floor("us")))
    now = wgs84.subpoint_of(satellite.at(times))
    later = wgs84.subpoint_of(satellite.at(timescale.tt_jd(times.tt + 1 / 86400)))
    heading = GEOD.inv(
        now.longitude.degrees,
        now.latitude.degrees,
        later.longitude.degrees,
        later.latitude.degrees,
    )[0]
    azimuth, _, length = GEOD.inv(
        now.longitude.degrees,
        now.latitude.degrees,
        samples.longitude.to_numpy(),
        samples.latitude.to_numpy(),
    )
    return azimuth - heading, length


def add_offsets(segments: pd.DataFrame, orbit) -> pd.DataFrame:
    """Adds the cross-track offsets (km, positive to the left) from the ground track."""
    bearing, length = track_bearings(segments, orbit)
    return segments.assign(cross_km=-length * np.sin(np.radians(bearing)) / 1e3)


def fit(segments: pd.DataFrame) -> None:
    """Prints the fitted roll angle and field of view of each pair's view."""
    for pair in PAIRS:
        offsets = segments.cross_km[segments.pair == pair].to_numpy()
        center = np.median(offsets)
        half = np.percentile(np.abs(offsets - center), FIT_PERCENTILE)
        near, far = angle(1e3 * (center - half)), angle(1e3 * (center + half))
        print(
            f"{pair}: {len(offsets)} segments, center {center:.3f} km, half-width "
            f"{half:.3f} km; roll {(near + far) / 2:.4f} deg, cross-track field of "
            f"view {far - near:.4f} deg"
        )


def score_views(
    run, configuration: str, orbit, samples: pd.DataFrame, views, per_view: int
) -> tuple[dict, np.ndarray]:
    """
    Scores a sample of the segments (or shots) of each view (by `pair`): whether the
    view observes them within two minutes of their times, and the epochs'
    differences, at each time step. Gets the results by time step, and the sample as
    rows of latitude, longitude, seconds from the start of the day, and the index of
    the view.
    """
    rng = np.random.default_rng(1)
    sample = pd.concat(
        [
            group.iloc[rng.choice(len(group), min(per_view, len(group)), replace=False)]
            for _, group in samples.groupby("pair")
        ],
        ignore_index=True,
    )
    results = {}
    for time_step in TIME_STEPS:
        satellite = run.entry.to_satellite(
            orbit,
            configuration=configuration,
            time_step=time_step,
            platform=run.platform,
        )
        differences = find_epochs(satellite, sample, views)
        found = np.isfinite(differences)
        results[time_step] = {
            "segment_recall": 100 * float(found.mean()),
            "epoch_p95": float(summarize(differences[found])["p95"]),
            "epoch_median": float(summarize(differences[found])["median"]),
        }
        print(
            f"  time step {time_step} s: {json.dumps(results[time_step])}", flush=True
        )
    passes = np.stack(
        [
            sample.latitude.to_numpy(),
            sample.longitude.to_numpy(),
            (sample.time - pd.Timestamp(run.day)).dt.total_seconds().to_numpy(),
            sample.pair.map(views.index).to_numpy(),
        ],
        axis=-1,
    )
    return results, passes


def find_epochs(satellite, sample: pd.DataFrame, views) -> np.ndarray:
    """
    Gets the difference (s) of each sample's nearest predicted epoch by the view
    (`pair`, of `views`), within two minutes of its time, or NaN.
    """

    def find(row):
        time = row.time.floor("us").to_pydatetime()
        found = collect_observations(
            Point(latitude=row.latitude, longitude=row.longitude),
            satellite,
            time - timedelta(minutes=2),
            time + timedelta(minutes=2),
            views.index(row.pair),
        )
        if found.empty:
            return np.nan
        return min(
            ((e - pd.Timestamp(time)).total_seconds() for e in found.epoch), key=abs
        )

    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        return np.array(list(pool.map(find, sample.itertuples())), dtype=float)


def write_views(
    run,
    configuration: str,
    orbit,
    tle,
    scored: tuple[dict, np.ndarray],
    spacing: float,
    views,
    reference: str,
) -> None:
    """
    Writes the results of the views (see `score_views`) and the spacing error, with
    the reduced fixture of a sample of the passes, rescored at the first time step.
    """
    results, passes = scored
    satellite = run.entry.to_satellite(
        orbit,
        configuration=configuration,
        time_step=TIME_STEPS[0],
        platform=run.platform,
    )
    rng = np.random.default_rng(2)
    reduced = passes[rng.choice(len(passes), FIXTURE_PASSES, replace=False)]
    fixture_differences = find_epochs(
        satellite,
        pd.DataFrame(
            {
                "latitude": reduced[:, 0],
                "longitude": reduced[:, 1],
                "time": pd.Timestamp(run.day)
                + pd.to_timedelta(reduced[:, 2], unit="s"),
                "pair": [views[int(i)] for i in reduced[:, 3]],
            }
        ),
        views,
    )
    found = np.isfinite(fixture_differences)
    write_results(
        run,
        configuration,
        reference=reference,
        metrics=[
            {
                "name": name,
                "value": values[name],
                "scope": f"time step {time_step} s",
            }
            for time_step, values in results.items()
            for name in ("segment_recall", "epoch_p95", "epoch_median")
        ]
        + [{"name": "beam_spacing_error_median", "value": spacing}],
        fixture={
            "tle": np.array(tle),
            "day": np.array(run.day.isoformat()),
            "passes": reduced,
        },
        fixture_metrics={
            "passes_found": 100 * float(found.mean()),
            "epoch_p95": float(summarize(fixture_differences[found])["p95"]),
        },
    )


def score_spacing(entry, segments: pd.DataFrame) -> float:
    """
    Gets the median error (km) of the spacing across track between adjacent pairs,
    measured from the pairs' median offsets in bins of 10 s of each file, against
    the spacing of the views' centers.
    """
    rolls = [
        spec.parameters["roll_angle"]
        for spec in entry.configurations[CONFIGURATION].instruments
    ]
    centers = distance(rolls) / 1e3
    binned = (
        segments.assign(bin=(segments.gps_seconds // SPACING_BIN).astype("int64"))
        .groupby(["file", "bin", "pair"])
        .cross_km.median()
        .unstack("pair")
    )
    errors = np.concatenate(
        [
            np.abs(
                (binned[PAIRS[i]] - binned[PAIRS[i + 1]]).dropna().to_numpy()
                - (centers[i] - centers[i + 1])
            )
            for i in range(len(PAIRS) - 1)
        ]
    )
    print(f"  spacing: {len(errors)} bins, median error {np.median(errors):.4f} km")
    return float(np.median(errors))


def main() -> None:
    """Validates the entry on a day of data, or fits its views."""
    run = parse_run(
        "atlas",
        __doc__,
        add_arguments=lambda parser: parser.add_argument(
            "--same-day-element-set",
            action="store_true",
            help="test with the element set nearest the middle of the day from "
            f"elements/{SAME_DAY_ELEMENT_SET_FILE.name}, without writing",
        ),
    )
    if run.same_day_element_set and (run.write or run.fit):
        raise SystemExit("--same-day-element-set neither fits nor writes a record")
    tle = (
        get_element_set(run.platform, run.day, file=SAME_DAY_ELEMENT_SET_FILE)
        if run.same_day_element_set
        else run.element_set()
    )
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    segments = add_offsets(load_day(run.day, run.directory), orbit)
    print(
        f"{run.day:%Y-%m-%d}: {tle[0][:32]}; {len(segments)} segments in "
        f"{segments.file.nunique()} files",
        flush=True,
    )
    if run.fit:
        fit(segments)
        return
    scored = score_views(run, CONFIGURATION, orbit, segments, PAIRS, SAMPLE_PER_PAIR)
    spacing = score_spacing(run.entry, segments)
    if run.write:
        write_views(
            run,
            CONFIGURATION,
            orbit,
            tle,
            scored,
            spacing,
            PAIRS,
            "ATL08QL version 007 (land segments south of 60 deg S)",
        )


if __name__ == "__main__":
    main()

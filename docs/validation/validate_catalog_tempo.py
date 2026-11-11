"""
Validates the TEMPO catalog entry against the TEMPO Level 1 radiance product, using
the exact values of the entry through public TAT-C methods. TEMPO is a geostationary
imaging spectrometer whose north-south slit (2,048 pixels) a mirror steps from east
to west across a fixed region of North America, about hourly in daylight. A full scan
(nine granules) images the region once; the area it images is the outline of its
outermost pixels' corners (the first and last pixels of each mirror step, and the
first and last mirror steps), in the ultraviolet band.

The view's angles are fitted to one full scan (run with `--fit`), from the median
tangents of each edge of the area seen from the satellite (whose Earth-fixed position
the product gives), iterated until the view centers them. Running on a scan at least
30 days from it tests the fitted values on data held out from the fit, with a public
element set: the field of regard is compared with the area imaged on a grid of 0.05
deg, by the intersection over union weighted by area, and a sample of grid points is
checked with `collect_observations`.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`. Only the
geometry is read, with range requests:

    python docs/validation/validate_catalog_tempo.py <data directory> 2026-10-06 --fit
    python docs/validation/validate_catalog_tempo.py <data directory> 2026-11-08

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
from shapely import contains_xy, prepare
from shapely.geometry import Polygon
from skyfield.api import wgs84
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric

import tatc
from tatc.analysis import collect_observations
from tatc.catalog import STANDARD_VERSION, get_entry
from tatc.catalog.metrics import compute_mask_iou
from tatc.constants import timescale
from tatc.schemas import GeneralPerturbationsOrbit, Point
from tatc.utils import VelocityFrame, compute_view_tangents
from validation_io import RemoteFile, ecef, get_holdout_gap, get_with_retries, up

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
RECORD_PATH = REPOSITORY / "src" / "tatc" / "resources" / "validation" / "tempo.json"
FIXTURE_PATH = REPOSITORY / "tests" / "catalog" / "fixtures" / "tempo.npz"
SHORT_NAME, VERSION = "TEMPO_RAD_L1", "V04"
FIT_DAY = datetime(2026, 10, 6, tzinfo=timezone.utc)
PLATFORM = "intelsat-40e"
CONFIGURATION = "field_of_regard"
# public element sets from CelesTrak, with epochs after each day of data
ELEMENTS = {
    "2026-11-08": (
        "1 56174U 23052A   26313.52440553 -.00000165  00000+0  00000+0 0  9999",
        "2 56174   0.0146 267.4057 0001991 300.1798 298.9210  1.00271930 13283",
    ),
}
GPS_EPOCH = datetime(1980, 1, 6, tzinfo=timezone.utc)
BAND = "band_290_490_nm"
GRANULES = 9  # of a full scan
SCAN_HOUR = 18  # UTC, the full scan nearest which is compared
GRID_STEP, GRID_MARGIN = 0.05, 3.0  # deg
CHECK_SAMPLES = 400
FIXTURE_SAMPLES = 400
FIT_ITERATIONS = 8
AU_KM = 149597870.7
FILL = 1e30  # bounds at or beyond are fill values


def find_scan(day: datetime) -> list[str]:
    """Gets the URLs of the granules of the full scan nearest SCAN_HOUR of a day."""
    items = get_with_retries(
        "https://cmr.earthdata.nasa.gov/search/granules.umm_json",
        params={
            "short_name": SHORT_NAME,
            "version": VERSION,
            "temporal": f"{day:%Y-%m-%d}T00:00:00Z,{day:%Y-%m-%d}T23:59:59Z",
            "page_size": 500,
        },
    ).json()["items"]
    scans = {}
    for item in items:
        url = next(
            u["URL"]
            for u in item["umm"]["RelatedUrls"]
            if u["URL"].startswith("https://data.asdc") and u["URL"].endswith(".nc")
        )
        scan = url.split("_")[-1][:4]
        scans.setdefault(scan, []).append(url)
    full = {k: sorted(v) for k, v in scans.items() if len(v) == GRANULES}
    noon = day + timedelta(hours=SCAN_HOUR)

    def start(urls):
        stamp = urls[0].split("_")[-2]
        return datetime.strptime(stamp, "%Y%m%dT%H%M%SZ").replace(tzinfo=timezone.utc)

    return min(
        full.values(), key=lambda urls: abs((start(urls) - noon).total_seconds())
    )


def read_granule(url: str, session, first: bool, last: bool) -> dict:
    """
    Reads a granule's outermost pixels' corners (north and south ends of each mirror
    step, and all of the first or last step of the scan), times, and the satellite's
    Earth-fixed positions (km).
    """
    with h5py.File(RemoteFile(url, session, block_size=2**16)) as f:
        g = f[BAND]
        result = {
            "north": np.stack(
                [g["latitude_bounds"][:, 0, :], g["longitude_bounds"][:, 0, :]], -1
            ),
            "south": np.stack(
                [g["latitude_bounds"][:, -1, :], g["longitude_bounds"][:, -1, :]], -1
            ),
            "time": f["time"][:],
            "position": np.stack(
                [f[f"geometry/satellite_{c}"][:] for c in "XYZ"], axis=-1
            ),
        }
        if first:
            result["east"] = np.stack(
                [g["latitude_bounds"][0], g["longitude_bounds"][0]], -1
            )
        if last:
            result["west"] = np.stack(
                [g["latitude_bounds"][-1], g["longitude_bounds"][-1]], -1
            )
    return result


def load_scan(day: datetime, directory: pathlib.Path) -> dict:
    """Reads the outline of the full scan nearest SCAN_HOUR, or loads it from a cache."""
    cache = directory / f"tempo_{day:%Y%m%d}.npz"
    if not cache.exists():
        urls = find_scan(day)
        if earthaccess.login(strategy="environment") is None:
            raise RuntimeError("Earthdata login failed: is EARTHDATA_TOKEN set?")
        session = earthaccess.get_requests_https_session()
        with concurrent.futures.ThreadPoolExecutor(GRANULES) as pool:
            parts = list(
                pool.map(
                    lambda k: read_granule(
                        urls[k], session, k == 0, k == len(urls) - 1
                    ),
                    range(len(urls)),
                )
            )
        np.savez_compressed(
            cache,
            north=np.concatenate([p["north"] for p in parts]),
            south=np.concatenate([p["south"] for p in parts]),
            east=parts[0]["east"],
            west=parts[-1]["west"],
            time=np.concatenate([p["time"] for p in parts]),
            position=np.concatenate([p["position"] for p in parts]),
            granules=np.array([u.split("/")[-1] for u in urls]),
        )
    with np.load(cache, allow_pickle=False) as data:
        return {key: data[key] for key in data.files}


def outline(scan: dict) -> Polygon:
    """
    Gets the outline of the area imaged by a scan (longitude, latitude), from the
    outer corners of its outermost pixels: east to west along the north edge, north
    to south along the west edge, west to east along the south edge, and south to
    north along the east edge.
    """

    # corners of a pixel's bounds: the two with the most northern (or southern,
    # eastern, western) of its four corners, ordered along the edge
    def side(corners, key, reverse):
        points = []
        for pixel in corners:
            order = np.argsort(key(pixel))
            pair = pixel[order[-2:]]
            points.extend(pair[np.argsort(pair[:, 1])[:: -1 if reverse else 1]])
        return points

    def valid(corners):
        return corners[np.all(np.abs(corners) < FILL, axis=(1, 2))]

    north = side(valid(scan["north"]), lambda c: c[:, 0], True)
    west = side(valid(scan["west"]), lambda c: -c[:, 1], False)
    south = side(valid(scan["south"][::-1]), lambda c: -c[:, 0], False)
    east = side(valid(scan["east"][::-1]), lambda c: c[:, 1], False)
    ring = [(lon, lat) for lat, lon in (*north, *west, *south, *east)]
    return Polygon(ring).buffer(0)


def scan_time(scan: dict) -> datetime:
    """Gets the middle time of a scan."""
    return GPS_EPOCH + timedelta(seconds=float(np.median(scan["time"])))


def product_track(scan: dict) -> Geocentric:
    """
    Gets the satellite's position and inertial velocity at the middle of a scan, from
    its Earth-fixed position in the product.
    """
    time = scan_time(scan)
    fixed = np.median(scan["position"], axis=0)
    t = timescale.from_datetimes(
        [time - timedelta(seconds=1), time, time + timedelta(seconds=1)]
    )
    # the rotation takes inertial (GCRS) coordinates to Earth-fixed ones
    inertial = np.einsum("jin,j->in", itrs.rotation_at(t), fixed) / AU_KM
    velocity = (inertial[:, 2] - inertial[:, 0]) / 2 * 86400
    return Geocentric(
        inertial[:, 1], velocity, timescale.from_datetime(time), center=399
    )


def repeat_track(track: Geocentric, _count: int) -> Geocentric:
    """Shapes a position at one time to broadcast against many targets."""
    return Geocentric(
        track.position.au.reshape(3, 1),
        track.velocity.au_per_d.reshape(3, 1),
        timescale.tt_jd(np.atleast_1d(track.t.tt)[:1]),
        center=399,
    )


def edge_points(scan: dict) -> dict[str, np.ndarray]:
    """Gets the outer corners of each edge's pixels (latitude, longitude)."""

    def corners(values, extreme, coordinate):
        values = values[np.all(np.abs(values) < FILL, axis=(1, 2))]
        index = np.argmax(extreme * values[..., coordinate], axis=-1)
        return values[np.arange(len(values)), index]

    return {
        "north": corners(scan["north"], 1, 0),
        "south": corners(scan["south"], -1, 0),
        "east": corners(scan["east"], 1, 1),
        "west": corners(scan["west"], -1, 1),
    }


def tangents(track, points, roll, pitch) -> tuple[np.ndarray, np.ndarray]:
    """Gets the along- and cross-track tangents of points from a view's center."""
    target = wgs84.latlon(points[:, 0], points[:, 1])
    along, cross = compute_view_tangents(
        repeat_track(track, len(points)), target, VelocityFrame.INERTIAL, roll, pitch
    )
    return np.atleast_1d(along), np.atleast_1d(cross)


def fit(scan: dict) -> None:
    """
    Prints the view that centers the median tangents of the scan's edges: its roll
    and pitch angles, and its fields of view across track (north-south) and along
    track (east-west).
    """
    track = product_track(scan)
    edges = edge_points(scan)
    roll = pitch = 0.0
    for _ in range(FIT_ITERATIONS):
        medians = {}
        for name, points in edges.items():
            along, cross = tangents(track, points, roll, pitch)
            medians[name] = np.median(cross if name in ("north", "south") else along)
        cross_mid = (medians["north"] + medians["south"]) / 2
        along_mid = (medians["east"] + medians["west"]) / 2
        roll += np.degrees(np.arctan(cross_mid))
        pitch += np.degrees(np.arctan(along_mid))
    cross_fov = np.degrees(np.arctan(abs(medians["north"] - medians["south"]) / 2)) * 2
    along_fov = np.degrees(np.arctan(abs(medians["east"] - medians["west"]) / 2)) * 2
    print(f"scan {scan['granules'][0]} at {scan_time(scan):%Y-%m-%dT%H:%M:%SZ}")
    print(f"edge tangents after centering: {medians}")
    print(
        f"roll {roll:.4f} deg, pitch {pitch:.4f} deg, cross-track field of view "
        f"{cross_fov:.4f} deg, along-track field of view {along_fov:.4f} deg"
    )
    for name, points in edges.items():
        along, cross = tangents(track, points, roll, pitch)
        values = cross if name in ("north", "south") else along
        print(
            f"  {name}: {len(points)} corners, tangent p5 {np.percentile(values, 5):.5f}, "
            f"median {np.median(values):.5f}, p95 {np.percentile(values, 95):.5f}"
        )


def build_grid(area: Polygon) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Gets the points of a grid around an area, and their weights by area."""
    west, south, east, north = area.bounds
    lon, lat = np.meshgrid(
        np.arange(west - GRID_MARGIN, east + GRID_MARGIN, GRID_STEP),
        np.arange(south - GRID_MARGIN, north + GRID_MARGIN, GRID_STEP),
    )
    lon, lat = lon.ravel(), lat.ravel()
    return lon, lat, np.cos(np.radians(lat))


def predict_grid(satellite, track, lon, lat) -> np.ndarray:
    """Gets whether the instrument observes grid points, in view and above the horizon."""
    instrument = satellite.get_instrument(0)
    target = wgs84.latlon(lat, lon)
    in_view = instrument.is_in_field_of_view(repeat_track(track, len(lon)), target)
    position = np.asarray(track.frame_xyz(itrs).m).reshape(3)
    points = ecef(lon, lat)
    line = position - points
    return in_view & (np.sum(line * up(lon, lat), axis=-1) > 0)


def check_points(satellite, time, lon, lat, rng) -> tuple[np.ndarray, np.ndarray]:
    """
    Gets a sample of grid points near the edges and whether `collect_observations`
    observes each within ten minutes of the scan's middle.
    """

    def observe(k):
        found = collect_observations(
            Point(latitude=float(lat[k]), longitude=float(lon[k])),
            satellite,
            time - timedelta(minutes=10),
            time + timedelta(minutes=10),
        )
        return not found.empty

    picks = rng.choice(len(lon), CHECK_SAMPLES, replace=False)
    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        return picks, np.array(list(pool.map(observe, picks)), dtype=bool)


def validate(day, tle, scan, write) -> None:
    """Scores the configuration on a held-out scan, and writes its results if requested."""
    entry = get_entry("tempo")
    orbit = GeneralPerturbationsOrbit.from_tle(list(tle))
    satellite = entry.to_satellite(orbit, platform=PLATFORM)
    time = scan_time(scan)
    track = orbit.to_gp_orbit().get_orbit_track([time])[0]
    area = outline(scan)
    prepare(area)
    lon, lat, weights = build_grid(area)
    observed = contains_xy(area, lon, lat)
    predicted = predict_grid(satellite, track, lon, lat)
    iou = compute_mask_iou(observed, predicted, weights)
    results = {
        "field_of_regard_iou": iou,
        "observed_predicted": 100
        * float(np.sum(weights[observed & predicted]) / np.sum(weights[observed])),
        "predicted_observed": 100
        * float(np.sum(weights[observed & predicted]) / np.sum(weights[predicted])),
    }
    print(f"scan {scan['granules'][0]} at {time:%Y-%m-%dT%H:%M:%SZ}")
    print(json.dumps(results), flush=True)
    # a sample near the boundary, where the methods could disagree
    rng = np.random.default_rng(1)
    near = np.flatnonzero(
        np.abs(
            np.asarray(contains_xy(area.buffer(1.0), lon, lat), dtype=int)
            - np.asarray(contains_xy(area.buffer(-1.0), lon, lat), dtype=int)
        )
    )
    picks, collected = check_points(satellite, time, lon[near], lat[near], rng)
    agreement = 100 * float(np.mean(collected == predicted[near][picks]))
    print(
        f"collect_observations agrees with the grid on {agreement:.1f}% of {len(picks)} points near the edges"
    )
    if not write:
        return
    keep = rng.choice(len(lon), FIXTURE_SAMPLES, replace=False)
    sample = np.concatenate([near[picks], keep])
    fixture = np.stack([lat[sample], lon[sample], observed[sample]], axis=-1)
    np.savez_compressed(
        FIXTURE_PATH,
        tle=np.array(tle),
        time=np.array(time.isoformat()),
        points=fixture,
    )
    fixture_predicted = predict_grid(satellite, track, lon[sample], lat[sample])
    record = {
        "entry_id": "tempo",
        "standard_version": STANDARD_VERSION,
        "configurations": {
            CONFIGURATION: {
                "digest": entry.configurations[CONFIGURATION].digest,
                "evidence": "rerun",
                "tatc_version": tatc.__version__,
                "source": "docs/validation/validate_catalog_tempo.py",
                "reference": f"{SHORT_NAME} {VERSION} (outer corners of the outermost "
                "pixels of a full scan, ultraviolet band)",
                "period": f"{time:%Y-%m-%dT%H:%MZ}",
                "platform": PLATFORM,
                "holdout_gap_days": get_holdout_gap(day, FIT_DAY),
                "metrics": [
                    {"name": "field_of_regard_iou", "value": float(iou)},
                    {
                        "name": "imaged_area_predicted",
                        "value": results["observed_predicted"],
                    },
                    {
                        "name": "predicted_area_imaged",
                        "value": results["predicted_observed"],
                    },
                    {"name": "collect_observations_agreement", "value": agreement},
                ],
                "fixture_metrics": {
                    "points_agreeing": 100
                    * float(np.mean(fixture_predicted == fixture[:, 2].astype(bool)))
                },
            }
        },
    }
    RECORD_PATH.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {RECORD_PATH} and {FIXTURE_PATH}")


def main() -> None:
    """Fits or validates the entry on a scan."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("day")
    parser.add_argument("--fit", action="store_true")
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)
    day = datetime.fromisoformat(arguments.day).replace(tzinfo=timezone.utc)
    scan = load_scan(day, arguments.directory)
    if arguments.fit:
        fit(scan)
        return
    validate(day, ELEMENTS[arguments.day], scan, arguments.write)


if __name__ == "__main__":
    main()

"""
Validates the C-SAR catalog entry against the Level-1 Ground Range Detected (GRD)
products of Sentinel-1C, using the exact values of the entry through public TAT-C
methods. The annotation of each product holds a geolocation grid, with the zero-Doppler
time of each of its lines, and the restituted orbit:

* The swath edges are compared with the near-range and far-range points of the grid
  lines, against the rays at the edges of the view at the times of the lines.
* The swath area is compared with the union of consecutive products of an acquisition,
  from the middle of the first to the middle of the last, so that the arbitrary ends of
  products are left out, at several time steps.
* Passes are those of lattice points inside the grids.
* Epochs are compared with the zero-Doppler times of interior grid points.

Each mode is a configuration (`iw` and `ew`), whose field of view and roll profile (the
roll steering, as a function of the argument of latitude) come from `--fit`, on a day
of products. Running on a day at least 30 days from it tests them on data held out from
the fit.

Needs the `validation` extra:

    python validation/catalog/validate_catalog_csar.py <data directory> 2026-10-02 --fit
    python validation/catalog/validate_catalog_csar.py <data directory> 2026-08-12 \
        --configuration iw
    python validation/catalog/validate_catalog_csar.py <data directory> 2026-11-10 \
        --configuration ew

Without a day, it takes that of the record of the configuration given with
`--configuration` (or, with `--fit`, that of the fit). With `--write`, it writes the
entry's validation record and the reduced fixture from which unit tests recompute the
results.

With `--platform sentinel-1d`, it tests the same values on a sibling satellite, and
adds the results to the record as a transfer, which does not count toward the tier.

It also validates the entry's resolution configuration (`resolution_iw`, run with
`--configuration resolution_iw`): the resolutions of the IW sub-swaths, derived from
the look bandwidths and windows of the ground range detected product, are compared
with the resolutions of the Sentinel-1 Product Definition (Table 7-6), which the
gate tests.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import concurrent.futures
import json
import pathlib
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta

import numpy as np
import pandas as pd
from pyproj import Transformer
from shapely.geometry import Polygon
from shapely.ops import unary_union

from tatc.analysis import compute_ground_track
from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import GeneralPerturbationsOrbit, PointedInstrument
import sar
from sar import (
    CHUNK_PRODUCTS,
    MAX_PRODUCT_GAP,
    MIN_CHUNK_PRODUCTS,
    TIME_STEPS,
    to_utc,
)
from validation_io import (
    compute_published_errors,
    ecef,
    get_with_retries,
    parse_run,
    post_with_retries,
    recompute_passes,
    score_published,
    write_results,
)

# published resolutions of the sub-swaths of the IW mode's ground range detected high
# resolution product: instrument, altitude (km), inclination (deg), roll (deg, the
# look angle of the sub-swath's incidence angle), and resolutions in azimuth and
# ground range (km), from Table 7-6 of the Sentinel-1 Product Definition
# (S1-RS-MDA-52-7440 v2.13)
RESOLUTIONS = {
    "resolution_iw": (
        ("C-SAR IW1", 693, 98.18, -29.3332, 0.0225, 0.0204),
        ("C-SAR IW2", 693, 98.18, -33.9851, 0.0226, 0.0203),
        ("C-SAR IW3", 693, 98.18, -38.0423, 0.0226, 0.0205),
    ),
}
STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1/search"
TOKEN_URL = "https://planetarycomputer.microsoft.com/api/sas/v1/token/sentinel-1-grd"

REFERENCE = (
    "Geolocation grids and restituted orbits of the annotations "
    + "of Sentinel-1 Level-1 GRD products"
)
MODES = {"iw": "IW", "ew": "EW"}
SPEED_OF_LIGHT = 299792458


def search_items(day: datetime, platform: str) -> list[dict]:
    """Gets the catalog items of a platform's IW and EW products of a day."""
    body = {
        "collections": ["sentinel-1-grd"],
        "datetime": f"{day:%Y-%m-%d}T00:00:00Z/{day:%Y-%m-%d}T23:59:59Z",
        "limit": 1000,
        "query": {"platform": {"eq": platform}},
    }
    items = []
    while True:
        page = post_with_retries(STAC_URL, body).json()
        items += page["features"]
        following = [
            link for link in page.get("links", []) if link.get("rel") == "next"
        ]
        if not following or not page["features"]:
            break
        body = following[0].get("body", body)
    return [
        item
        for item in items
        if item["properties"]["sar:instrument_mode"] in MODES.values()
    ]


def read_annotation(item: dict, token: str, day: datetime) -> dict:
    """Reads the geolocation grid and orbit of a product from its annotation."""
    base = item["assets"]["safe-manifest"]["href"].rsplit("/", 1)[0]
    mode = item["properties"]["sar:instrument_mode"]
    polarization = item["properties"]["sar:polarizations"][0].lower()
    response = get_with_retries(
        f"{base}/annotation/{mode.lower()}-{polarization}.xml?{token}"
    )
    root = ET.fromstring(response.content)

    def seconds(text):
        return (pd.Timestamp(text, tz="UTC") - pd.Timestamp(day)).total_seconds()

    points = root.findall(".//geolocationGridPoint")
    lines = sorted({int(p.find("line").text) for p in points})
    samples = sorted({int(p.find("pixel").text) for p in points})

    def grid(tag, convert=float):
        return np.array([convert(p.find(tag).text) for p in points]).reshape(
            len(lines), len(samples)
        )

    orbits = root.findall(".//orbitList/orbit")
    return {
        "id": item["id"],
        "mode": mode,
        "datatake": str(item["properties"]["s1:datatake_id"]),
        "time": grid("azimuthTime", seconds),
        "slant_range": grid("slantRangeTime")[0] * SPEED_OF_LIGHT / 2,
        "lon": grid("longitude"),
        "lat": grid("latitude"),
        "orbit_time": np.array([seconds(o.find("time").text) for o in orbits]),
        "orbit_position": np.array(
            [[float(o.find(f"position/{c}").text) for c in "xyz"] for o in orbits]
        ),
        "orbit_velocity": np.array(
            [[float(o.find(f"velocity/{c}").text) for c in "xyz"] for o in orbits]
        ),
    }


def load_products(day: datetime, directory: pathlib.Path, platform: str) -> list[dict]:
    """
    Gets the reduced annotations of a platform's IW and EW products of a day, in time
    order.
    """
    path = directory / f"csar_{platform}_{day:%Y%m%d}.npz"
    if not path.exists():
        items = search_items(day, platform)
        print(f"reading the annotations of {len(items)} products", flush=True)
        token = get_with_retries(TOKEN_URL).json()["token"]
        with concurrent.futures.ThreadPoolExecutor(16) as pool:
            records = list(pool.map(lambda i: read_annotation(i, token, day), items))
        np.savez_compressed(
            path,
            id=[r["id"] for r in records],
            mode=[r["mode"] for r in records],
            datatake=[r["datatake"] for r in records],
            n_line=[r["time"].shape[0] for r in records],
            n_sample=[r["time"].shape[1] for r in records],
            n_orbit=[len(r["orbit_time"]) for r in records],
            **{
                k: np.concatenate([r[k].ravel() for r in records])
                for k in ("time", "lon", "lat")
            },
            **{
                k: np.concatenate([r[k] for r in records])
                for k in ("orbit_time", "orbit_position", "orbit_velocity")
            },
        )
    with np.load(path, allow_pickle=False) as data:
        cache = {key: data[key] for key in data.files}
    size = cache["n_line"] * cache["n_sample"]

    def cut(key, counts):
        return np.split(cache[key], np.cumsum(counts)[:-1])

    products = []
    for k, (t, lon, lat, ot, op, ov) in enumerate(
        zip(
            cut("time", size),
            cut("lon", size),
            cut("lat", size),
            cut("orbit_time", cache["n_orbit"]),
            cut("orbit_position", cache["n_orbit"]),
            cut("orbit_velocity", cache["n_orbit"]),
            strict=True,
        )
    ):
        shape = (int(cache["n_line"][k]), int(cache["n_sample"][k]))
        products.append(
            {
                "id": str(cache["id"][k]),
                "mode": str(cache["mode"][k]),
                "datatake": str(cache["datatake"][k]),
                # the samples of a line differ by microseconds
                "time": t.reshape(shape).mean(axis=1),
                "lon": lon.reshape(shape),
                "lat": lat.reshape(shape),
                "orbit_time": ot,
                "orbit_position": op,
                "orbit_velocity": ov,
            }
        )
    products.sort(key=lambda p: p["time"][0])
    print(
        f"{len(products)} products: "
        + ", ".join(
            f"{sum(p['mode'] == mode for p in products)} {mode}"
            for mode in MODES.values()
        ),
        flush=True,
    )
    return products


def fit(products, orbit, day) -> None:
    """Measures the field of view and roll profile of each mode."""
    lines = sar.measure_lines(products, orbit, day)
    for key, mode in MODES.items():
        chosen = lines[lines["mode"] == mode]
        field_of_view = float(chosen.field_of_view.median())
        roll = float(chosen.roll.median())
        print(
            f"{key}: {chosen['product'].nunique()} products, {len(chosen)} lines; near "
            f"{chosen.near.median():.4f} deg, far {chosen.far.median():.4f} deg, roll "
            f"{roll:.4f} deg (p5 {chosen.roll.quantile(0.05):.4f}, p95 "
            f"{chosen.roll.quantile(0.95):.4f}), field of view {field_of_view:.4f} deg "
            f"(p5 {chosen.field_of_view.quantile(0.05):.4f}, p95 "
            f"{chosen.field_of_view.quantile(0.95):.4f})"
        )
        profile = sar.fit_profile(chosen)
        print(f"  roll profile: {json.dumps(profile)}")
        print(f"  far look angle, largest: {-chosen.far.min():.4f} deg")
        # a check within the fit day: profiles of even products, on odd products
        even = chosen[chosen["product"] % 2 == 0]
        odd = sorted(set(chosen["product"][chosen["product"] % 2 == 1]))
        for label, parameters in (
            ("fixed roll", {"roll_angle": roll}),
            ("roll profile", {"roll_angle_profile": sar.fit_profile(even)}),
        ):
            instrument = PointedInstrument(
                name="check",
                field_of_regard=120,
                cross_track_field_of_view=field_of_view,
                along_track_field_of_view=1,
                is_rectangular=True,
                **parameters,
            )
            edges = sar.score_edges(
                [products[k] for k in odd], orbit, day, instrument, "right"
            )
            print(
                f"  {label}, odd products: edges p95 near "
                f"{np.percentile(edges['near'], 95):.2f} km, far "
                f"{np.percentile(edges['far'], 95):.2f} km"
            )


def outline(product: dict, first: int, last: int) -> np.ndarray:
    """Gets the Earth-fixed outline of grid lines from a first to a last."""
    lon, lat = product["lon"][first : last + 1], product["lat"][first : last + 1]
    ring_lon = np.concatenate([lon[0], lon[1:, -1], lon[-1, -2::-1], lon[-2:0:-1, 0]])
    ring_lat = np.concatenate([lat[0], lat[1:, -1], lat[-1, -2::-1], lat[-2:0:-1, 0]])
    return ecef(ring_lon, ring_lat)


def get_chunks(products: list[dict]) -> list[list[int]]:
    """Divides the consecutive products of each acquisition into chunks."""
    sequences, current = [], []
    for k, product in enumerate(products):
        if current:
            previous = products[current[-1]]
            if (
                product["datatake"] != previous["datatake"]
                or product["time"][0] - previous["time"][-1] > MAX_PRODUCT_GAP
            ):
                sequences.append(current)
                current = []
        current.append(k)
    sequences.append(current)
    chunks = []
    for sequence in sequences:
        for k in range(0, len(sequence), CHUNK_PRODUCTS):
            chunk = sequence[k : k + CHUNK_PRODUCTS]
            if len(chunk) >= MIN_CHUNK_PRODUCTS:
                chunks.append(chunk)
    return chunks


def score_chunk(chunk, products, day, satellite_for) -> dict[str, float]:
    """Gets the intersection over union of a chunk, by time step."""
    first, last = products[chunk[0]], products[chunk[-1]]
    middle_first, middle_last = len(first["time"]) // 2, len(last["time"]) // 2
    center = products[chunk[len(chunk) // 2]]
    lon0, lat0 = (
        float(np.mean(center["lon"])),
        float(np.mean(center["lat"])),
    )
    project = Transformer.from_crs(
        "EPSG:4978",
        f"+proj=aeqd +lat_0={lat0} +lon_0={lon0} +units=m",
        always_xy=True,
    )
    to_plane = Transformer.from_crs(
        "EPSG:4326",
        f"+proj=aeqd +lat_0={lat0} +lon_0={lon0} +units=m",
        always_xy=True,
    )

    def plane(xyz):
        x, y, _ = project.transform(*xyz.T)
        return Polygon(np.column_stack([x, y])).buffer(0)

    parts = []
    for position, k in enumerate(chunk):
        product = products[k]
        n = len(product["time"])
        if position == 0:
            parts.append(plane(outline(product, middle_first, n - 1)))
        elif position == len(chunk) - 1:
            parts.append(plane(outline(product, 0, middle_last)))
        else:
            parts.append(plane(outline(product, 0, n - 1)))
    observed = unary_union(parts)
    start = first["time"][middle_first]
    duration = last["time"][middle_last] - start
    results = {"products": len(chunk), "latitude": lat0}
    for time_step in TIME_STEPS:
        n = max(1, int(round(duration / time_step)))
        step = duration / n
        times = to_utc(day, start + step * (np.arange(n) + 0.5))
        track = compute_ground_track(
            satellite_for(step), times, time_step=timedelta(seconds=float(step))
        )
        geometry = track.geometry.iloc[0]
        predicted = unary_union(
            [
                Polygon(
                    np.column_stack(
                        to_plane.transform(*np.array(part.exterior.coords).T[:2])
                    )
                ).buffer(0)
                for part in getattr(geometry, "geoms", [geometry])
            ]
        )
        results[f"iou_{time_step}"] = (
            predicted.intersection(observed).area / predicted.union(observed).area
        )
        if time_step == TIME_STEPS[0]:
            results["observed_only"] = (
                observed.difference(predicted).area / observed.area
            )
            results["predicted_only"] = (
                predicted.difference(observed).area / predicted.area
            )
    return results


def score_area(products, day, satellite_for) -> pd.DataFrame:
    """Gets the intersection over union of each chunk of products, by time step."""
    chunks = get_chunks(products)
    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        table = pd.DataFrame(
            list(
                pool.map(
                    lambda chunk: score_chunk(chunk, products, day, satellite_for),
                    chunks,
                )
            )
        )
    print(f"swath area over {len(table)} chunks of 3 to {CHUNK_PRODUCTS} products:")
    for time_step in TIME_STEPS:
        print(
            f"  IoU at {time_step} s: median {table[f'iou_{time_step}'].median():.5f}, "
            f"p5 {table[f'iou_{time_step}'].quantile(0.05):.5f}"
        )
    print(
        f"  observed only {100 * table.observed_only.median():.2f}%, "
        f"predicted only {100 * table.predicted_only.median():.2f}%",
        flush=True,
    )
    return table


def validate_mode(run, key, tle, orbit, products) -> None:
    """Scores the configuration of a mode, and writes its results if requested."""
    day = run.day
    chosen = [p for p in products if p["mode"] == MODES[key]]
    print(f"--- {key}: {len(chosen)} products", flush=True)

    def satellite_for(time_step: float):
        return run.entry.to_satellite(
            orbit, configuration=key, time_step=time_step, platform=run.platform
        )

    satellite = satellite_for(TIME_STEPS[0])
    # right-looking views have negative roll angles, with near range on the left
    edges = sar.score_edges(chosen, orbit, day, satellite.instruments[0], "right")
    sar.print_edges(edges)
    points = sar.epoch_points(chosen)
    differences = sar.score_epochs(points, satellite, day)
    epochs = sar.summarize_epochs(differences)
    lattice = generate_points_fibonacci_lattice(sar.POINT_SPACING)
    passes = sar.score_passes(lattice, satellite, chosen, day)
    print(json.dumps(passes), flush=True)
    area = score_area(chosen, day, satellite_for)
    metrics = sar.get_metrics(passes, edges, area, epochs)
    print(json.dumps(metrics), flush=True)
    if not run.write:
        return
    if run.transfer:
        write_results(run, key, reference=REFERENCE, metrics=metrics)
        return
    fixture = sar.pick_points(points, differences)
    write_results(
        run,
        key,
        reference=REFERENCE,
        metrics=metrics,
        # each mode keeps its own element set and day, since modes may be tested on
        # different days
        fixture={
            f"tle_{key}": np.array(tle),
            f"day_{key}": np.array(day.isoformat()),
            f"passes_{key}": fixture,
        },
        fixture_metrics=recompute_passes(satellite, 0, day, fixture),
    )


def validate_resolution(run, key: str) -> None:
    """
    Scores a resolution configuration against the published resolutions, and writes
    its results if the run writes.
    """
    within, largest = score_published(
        compute_published_errors(run.entry, key, RESOLUTIONS[key])
    )
    print(
        f"published resolutions: {within:.1f}% within 3%, largest error {largest:.2f}%"
    )
    if run.write:
        write_results(
            run,
            key,
            reference="Resolutions of the Sentinel-1 Product Definition "
            "(S1-RS-MDA-52-7440 v2.13), Table 7-6 (IW GRDH)",
            metrics=[
                {"name": "sample_size_within_3pct", "value": within},
                {"name": "sample_size_largest_error", "value": largest},
            ],
        )


def main() -> None:
    """Fits or validates the entry on a day of data."""
    run = parse_run(
        "csar",
        __doc__,
        add_arguments=lambda parser: parser.add_argument(
            "--download-only",
            action="store_true",
            help="read and cache the reference data, without scoring",
        ),
    )
    products = load_products(run.day, run.directory, run.platform)
    if run.download_only:
        return
    tle = run.element_set()
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    print(f"{run.platform}, {run.day:%Y-%m-%d}: {tle[0][:32]}", flush=True)
    if run.fit:
        fit(products, orbit, run.day)
        return
    for key in run.configurations:
        if key in RESOLUTIONS:
            # compared with published tables, so not a transfer to a platform
            if not run.transfer:
                validate_resolution(run, key)
        else:
            validate_mode(run, key, tle, orbit, products)


if __name__ == "__main__":
    main()

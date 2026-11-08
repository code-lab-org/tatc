"""
Validates the MSI catalog entry against the Level-2A products of Sentinel-2A, 2B, and
2C, using the exact values of the entry through public TAT-C methods. Products are
tiles of a fixed grid, cut from acquisitions (datatakes) whose ends are ragged, so:

* The swath width is compared with the width of the tiles' valid data across the track,
  at times inside acquisitions.
* The swath area is compared with the union of the tiles of an acquisition, between
  lines across the track at the start and end of chunks inside it.
* Passes are those of lattice points inside tiles.
* Epochs are compared with the times at which the red band imaged points of the viewing
  angle grids of a sample of tiles, recovered with the satellites' GPS ephemerides.

The entry's angles come from `--fit`, on a day of products of the three satellites.
Running on a day at least 30 days from it tests them on data held out from the fit.

Needs the `validation` extra:

    python docs/validation/validate_catalog_msi.py <data directory> 2026-10-03 --fit
    python docs/validation/validate_catalog_msi.py <data directory> 2026-08-12

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import concurrent.futures
import json
import pathlib
import re
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta, timezone

import numpy as np
import pandas as pd
from pyproj import Geod, Transformer
from scipy.spatial import cKDTree
from shapely import STRtree, points as make_points
from shapely.geometry import LineString, Polygon, shape
from shapely.ops import unary_union
from skyfield.framelib import itrs

import tatc
from tatc.analysis import collect_observations, compute_ground_track
from tatc.catalog import get_entry
from tatc.catalog.metrics import summarize
from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import GeneralPerturbationsOrbit, Point
from tatc.utils import compute_projected_ray_position
from validation_io import (
    TO_GEODETIC,
    ecef,
    get_holdout_gap,
    get_with_retries,
    match_passes,
    normalize,
    post_with_retries,
    predict_passes,
    read_element_set_file,
    recompute_passes,
    up,
    update_fixture,
    update_record,
)

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
RECORD_PATH = REPOSITORY / "src" / "tatc" / "resources" / "validation" / "msi.json"
FIXTURE_PATH = REPOSITORY / "tests" / "catalog" / "fixtures" / "msi.npz"
ELEMENT_SET_FILE = REPOSITORY / "docs" / "validation" / "sentinel_2026-08-12.3le"

STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1/search"
TOKEN_URL = "https://planetarycomputer.microsoft.com/api/sas/v1/token/sentinel-2-l2a"

FIT_DAY = datetime(2026, 10, 3, tzinfo=timezone.utc)
PLATFORMS = {
    "sentinel-2a": ("Sentinel-2A", 40697),
    "sentinel-2b": ("Sentinel-2B", 42063),
    "sentinel-2c": ("Sentinel-2C", 60989),
}
# the element sets of the validation notebook, for the fit day
ELEMENT_SETS = {
    "2026-10-03": {
        "sentinel-2a": [
            "1 40697U 15028A   26277.51535578  .00000205  00000+0  94862-4 0  9994",
            "2 40697  98.5646 350.7732 0001083  89.0049 271.1258 14.30820153589369",
        ],
        "sentinel-2b": [
            "1 42063U 17013A   26277.50820593  .00000133  00000+0  67483-4 0  9991",
            "2 42063  98.5691 350.7118 0001149  89.0927 271.0388 14.30818478500277",
        ],
        "sentinel-2c": [
            "1 60989U 24157A   26277.47326514  .00000038  00000+0  31260-4 0  9993",
            "2 60989  98.5703 350.6835 0001156  81.8213 278.3101 14.30817151108611",
        ],
    },
}
KEY = "swath"
FIXTURE_PLATFORM = "sentinel-2b"
BAND = "3"  # the identifier of the red band (B04) in the tile metadata
GPS_UTC = 18  # s, by which GPS time is ahead of UTC
GRID_TILE_STRIDE = 5
EPOCH_POINTS = 1000

END_MARGIN = 10  # s, left out at each end of an acquisition
SECTION_STEP = 20  # s, between the cross sections that measure the swath width
SECTION_HALF_LENGTH = 250e3  # m
SECTION_SPACING = 1e3  # m
MIN_SECTION_COMPLETENESS = 0.98
CHUNK_SECONDS = 150
MIN_CHUNK_SECONDS = 60
BAND_HALF_WIDTH = 400e3  # m, of the band between the ends of a chunk
TIME_STEPS = (10, 30, 120)
POINT_SPACING = 100e3
SEARCH_WINDOW = timedelta(minutes=10)
FIXTURE_PASSES = 300

GEOD = Geod(ellps="WGS84")
EARTH_RADIUS = 6371.0088e3


def get_element_set(day: datetime, platform: str) -> list[str]:
    """Gets the element set of a platform on a day."""
    key = f"{day:%Y-%m-%d}"
    if key in ELEMENT_SETS:
        return ELEMENT_SETS[key][platform]
    return read_element_set_file(ELEMENT_SET_FILE, PLATFORMS[platform][1], day)


def seconds_of(day: datetime, times) -> np.ndarray:
    """Gets seconds from the start of a day."""
    delta = pd.to_datetime(pd.Series(times), utc=True) - pd.Timestamp(day)
    return delta.dt.total_seconds().to_numpy()


def to_utc(day: datetime, seconds) -> list[datetime]:
    """Converts seconds from the start of a day to times."""
    return [day + timedelta(seconds=float(s)) for s in np.atleast_1d(seconds)]


def strip_namespaces(text: str) -> str:
    """Removes the namespace prefixes of the tags of an XML document."""
    return re.sub(r"<(/?)[A-Za-z0-9_]+:", r"<\1", text)


def load_tiles(day: datetime, directory: pathlib.Path) -> pd.DataFrame:
    """Gets the catalog entries of a day's tiles."""
    path = directory / f"msi_items_{day:%Y%m%d}.json"
    if not path.exists():
        body = {
            "collections": ["sentinel-2-l2a"],
            "datetime": f"{day:%Y-%m-%d}T00:00:00Z/{day:%Y-%m-%d}T23:59:59Z",
            "limit": 1000,
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
        path.write_text(
            json.dumps(
                [
                    {
                        "id": i["id"],
                        "geometry": i["geometry"],
                        "platform": i["properties"]["platform"],
                        "datatake": i["properties"]["s2:datatake_id"],
                        "datastrip": i["properties"]["s2:datastrip_id"],
                        "granule_metadata": i["assets"]["granule-metadata"]["href"],
                        "datastrip_metadata": i["assets"]["datastrip-metadata"]["href"],
                    }
                    for i in items
                ]
            ),
            encoding="utf-8",
        )
    tiles = pd.DataFrame(json.loads(path.read_text(encoding="utf-8")))
    names = {name: key for key, (name, _) in PLATFORMS.items()}
    tiles["platform"] = tiles.platform.map(names)
    tiles["geometry"] = [shape(g) for g in tiles.geometry]
    print(
        f"{len(tiles)} tiles: "
        + ", ".join(
            f"{(tiles.platform == key).sum()} {key}" for key in sorted(PLATFORMS)
        ),
        flush=True,
    )
    return tiles


def load_datastrips(day, directory, tiles) -> pd.DataFrame:
    """Gets the sensing start and stop times of the datastrips of the tiles."""
    path = directory / f"msi_datastrips_{day:%Y%m%d}.csv"
    if not path.exists():
        token = get_with_retries(TOKEN_URL).json()["token"]

        def read_times(row):
            head = get_with_retries(
                f"{row.datastrip_metadata}?{token}",
                headers={"Range": "bytes=0-16383"},
            ).text
            return {
                "datastrip": row.datastrip,
                **{
                    tag: re.search(
                        f"<DATASTRIP_SENSING_{tag.upper()}>(.*?)</", head
                    ).group(1)
                    for tag in ("start", "stop")
                },
            }

        first = tiles.drop_duplicates("datastrip")
        with concurrent.futures.ThreadPoolExecutor(16) as pool:
            pd.DataFrame(list(pool.map(read_times, first.itertuples()))).to_csv(
                path, index=False
            )
    return pd.read_csv(path)


def get_datatakes(day, tiles, datastrips) -> pd.DataFrame:
    """Gets the platform, start, and stop (s from the start of the day) of datatakes."""
    strips = datastrips.merge(
        tiles.drop_duplicates("datastrip")[["datastrip", "datatake", "platform"]]
    )
    strips["start"] = seconds_of(day, strips.start)
    strips["stop"] = seconds_of(day, strips.stop)
    datatakes = (
        strips.groupby("datatake")
        .agg(
            platform=("platform", "first"), start=("start", "min"), stop=("stop", "max")
        )
        .join(tiles.groupby("datatake").size().rename("tiles"))
        .reset_index()
        .sort_values("start")
        .reset_index(drop=True)
    )
    return datatakes


def load_reference_times(day, directory, tiles) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Gets the GPS ephemerides and the red band's viewing angle grids of every fifth tile
    of the datatake with the most tiles of each platform.
    """
    ephemeris_path = directory / f"msi_ephemeris_{day:%Y%m%d}.csv"
    grids_path = directory / f"msi_grids_{day:%Y%m%d}.csv"
    if not (ephemeris_path.exists() and grids_path.exists()):
        token = get_with_retries(TOKEN_URL).json()["token"]
        largest = tiles.groupby(["platform", "datatake"]).size().reset_index(name="n")
        selected = largest.sort_values("n").groupby("platform").tail(1).datatake
        chosen = tiles[tiles.datatake.isin(selected)]

        def read_ephemeris(row):
            root = ET.fromstring(
                strip_namespaces(
                    get_with_retries(f"{row.datastrip_metadata}?{token}").text
                )
            )
            points = root.findall(".//GPS_Point")
            frame = pd.DataFrame(
                {
                    "platform": row.platform,
                    "time": [
                        pd.Timestamp(p.find("GPS_TIME").text, tz="UTC")
                        - pd.Timedelta(seconds=GPS_UTC)
                        for p in points
                    ],
                }
            )
            for kind in ("position", "velocity"):
                values = np.array(
                    [
                        [
                            float(v)
                            for v in p.find(f"{kind.upper()}_VALUES").text.split()
                        ]
                        for p in points
                    ]
                )
                for k, axis in enumerate("xyz"):
                    frame[f"{kind}_{axis}"] = values[:, k] / 1e3
            return frame

        def read_grids(row):
            root = ET.fromstring(
                strip_namespaces(
                    get_with_retries(f"{row.granule_metadata}?{token}").text
                )
            )
            epsg = int(root.find(".//HORIZONTAL_CS_CODE").text.split(":")[1])
            geoposition = next(
                g for g in root.findall(".//Geoposition") if g.get("resolution") == "10"
            )
            ulx = float(geoposition.find("ULX").text)
            uly = float(geoposition.find("ULY").text)
            to_lonlat = Transformer.from_crs(epsg, 4326, always_xy=True)
            frames = []
            for grid in root.findall(".//Viewing_Incidence_Angles_Grids"):
                if grid.get("bandId") != BAND:
                    continue
                zenith, azimuth = (
                    np.array(
                        [
                            [float(v) for v in r.text.split()]
                            for r in grid.find(f"{name}/Values_List")
                        ]
                    )
                    for name in ("Zenith", "Azimuth")
                )
                i, j = np.nonzero(np.isfinite(zenith))
                lon, lat = to_lonlat.transform(ulx + j * 5000.0, uly - i * 5000.0)
                frames.append(
                    pd.DataFrame(
                        {
                            "platform": row.platform,
                            "detector": int(grid.get("detectorId")),
                            "lon": lon,
                            "lat": lat,
                            "zenith": zenith[i, j],
                            "azimuth": azimuth[i, j],
                        }
                    )
                )
            return pd.concat(frames) if frames else None

        with concurrent.futures.ThreadPoolExecutor(16) as pool:
            ephemeris = pd.concat(
                list(
                    pool.map(
                        read_ephemeris, chosen.drop_duplicates("datastrip").itertuples()
                    )
                )
            )
            sampled = pd.concat(
                [
                    group.sort_values("id").iloc[::GRID_TILE_STRIDE]
                    for _, group in chosen.groupby("datatake")
                ]
            )
            grids = pd.concat(
                [g for g in pool.map(read_grids, sampled.itertuples()) if g is not None]
            )
        ephemeris.to_csv(ephemeris_path, index=False)
        grids.to_csv(grids_path, index=False)
    ephemeris = pd.read_csv(ephemeris_path)
    ephemeris["t"] = seconds_of(day, ephemeris.time)
    return ephemeris, pd.read_csv(grids_path)


def recover_imaging_times(grids: pd.DataFrame, ephemeris: pd.DataFrame) -> pd.DataFrame:
    """
    Gets the times (s from the start of the day) at which grid points were imaged: each
    point's direction to the satellite, followed to the satellite's geocentric radius,
    locates the satellite, whose time is that of its nearest approach in the ephemeris.
    """
    frames = []
    for platform, group in grids.groupby("platform"):
        e = (
            ephemeris[ephemeris.platform == platform]
            .sort_values("t")
            .drop_duplicates("t")
        )
        t_e = e.t.to_numpy()
        position_e = e[[f"position_{a}" for a in "xyz"]].to_numpy()
        velocity_e = e[[f"velocity_{a}" for a in "xyz"]].to_numpy()

        def interpolate(t, values, t_e=t_e):
            return np.stack(
                [np.interp(t, t_e, values[:, k]) for k in range(3)], axis=-1
            )

        lon, lat = group.lon.to_numpy(), group.lat.to_numpy()
        point = ecef(lon, lat)
        vertical = up(lon, lat)
        east = normalize(
            np.stack(
                [-np.sin(np.radians(lon)), np.cos(np.radians(lon)), 0 * lon], axis=-1
            )
        )
        north = np.cross(vertical, east)
        zenith = np.radians(group.zenith.to_numpy())[:, None]
        azimuth = np.radians(group.azimuth.to_numpy())[:, None]
        toward = np.cos(zenith) * vertical + np.sin(zenith) * (
            np.sin(azimuth) * east + np.cos(azimuth) * north
        )
        b = np.sum(point * toward, axis=1)
        t = np.full(len(point), np.median(t_e))
        for _ in range(8):
            p, v = interpolate(t, position_e), interpolate(t, velocity_e)
            radius = np.linalg.norm(p, axis=1)
            location = (
                point
                + (-b + np.sqrt(b**2 - np.sum(point**2, axis=1) + radius**2))[:, None]
                * toward
            )
            t = t - np.sum((p - location) * v, axis=1) / np.sum(v * v, axis=1)
        p = interpolate(t, position_e)
        miss = np.linalg.norm(p - location, axis=1)
        frames.append(
            pd.DataFrame(
                {"platform": platform, "lat": lat, "lon": lon, "t": t, "miss": miss}
            )
        )
    return pd.concat(frames, ignore_index=True)


def get_frames(orbit, day, seconds) -> dict[str, np.ndarray]:
    """
    Gets a satellite's Earth-fixed position, horizontal along-track and left directions,
    and its sub-satellite point, altitude, and track azimuth, from its element set.
    """
    track = orbit.to_gp_orbit().get_orbit_track(to_utc(day, seconds))
    position, velocity = track.frame_xyz_and_velocity(itrs)
    position, velocity = position.m.T, velocity.m_per_s.T
    lon, lat, altitude = TO_GEODETIC.transform(*position.T)
    vertical = up(lon, lat)
    along = normalize(
        velocity - np.sum(velocity * vertical, axis=-1)[:, None] * vertical
    )
    east = normalize(
        np.stack(
            [
                -np.sin(np.radians(lon)),
                np.cos(np.radians(lon)),
                np.zeros(np.shape(lon)),
            ],
            axis=-1,
        )
    )
    north = np.cross(vertical, east)
    return {
        "position": position,
        "along": along,
        "left": np.cross(vertical, along),
        "nadir": -vertical,
        "lon": np.asarray(lon),
        "lat": np.asarray(lat),
        "altitude": np.asarray(altitude),
        "azimuth": np.degrees(
            np.arctan2(np.sum(along * east, axis=-1), np.sum(along * north, axis=-1))
        ),
    }


def local_projection(lon: float, lat: float) -> Transformer:
    """Gets an azimuthal equidistant projection centered on a point."""
    return Transformer.from_crs(
        "EPSG:4326", f"+proj=aeqd +lat_0={lat} +lon_0={lon} +units=m", always_xy=True
    )


def project_geometry(geometry, transformer):
    """Projects the polygons of a geometry, as one valid geometry."""
    return unary_union(
        [
            Polygon(
                np.column_stack(
                    transformer.transform(*np.array(part.exterior.coords).T[:2])
                )
            ).buffer(0)
            for part in getattr(geometry, "geoms", [geometry])
        ]
    )


class TileIndex:
    """Finds the tiles of a datatake near a point."""

    def __init__(self, tiles: pd.DataFrame):
        self.geometries = list(tiles.geometry)
        centers = np.array([[g.centroid.x, g.centroid.y] for g in self.geometries])
        self.tree = cKDTree(normalize(ecef(centers[:, 0], centers[:, 1])))

    def near(self, lon: float, lat: float, distance: float) -> list:
        """Gets the tiles whose centers are within a distance (m) of a point."""
        center = normalize(ecef(lon, lat))
        chord = 2 * np.sin(distance / EARTH_RADIUS / 2)
        return [self.geometries[k] for k in self.tree.query_ball_point(center, chord)]


def measure_sections(day, datatakes, tiles, orbits) -> pd.DataFrame:
    """
    Measures the valid data of the tiles of datatakes along lines across the track
    through the sub-satellite point, at times inside the datatakes. Gets the distances
    (m) to the left and right edges, their positions, and the satellite's position.
    """
    rows = []
    for datatake in datatakes.itertuples():
        start, stop = datatake.start + END_MARGIN, datatake.stop - END_MARGIN
        if stop <= start:
            continue
        index = TileIndex(tiles[tiles.datatake == datatake.datatake])
        times = np.arange(start, stop, SECTION_STEP)
        frames = get_frames(orbits[datatake.platform], day, times)
        distances = np.arange(
            -SECTION_HALF_LENGTH, SECTION_HALF_LENGTH + 1, SECTION_SPACING
        )
        for k, seconds in enumerate(times):
            lon0, lat0 = frames["lon"][k], frames["lat"][k]
            nearby = index.near(lon0, lat0, SECTION_HALF_LENGTH + 160e3)
            if not nearby:
                continue
            transformer = local_projection(lon0, lat0)
            area = unary_union([project_geometry(g, transformer) for g in nearby])
            # positive distances to the left of the direction of motion
            azimuth = np.where(
                distances >= 0, frames["azimuth"][k] - 90, frames["azimuth"][k] + 90
            )
            lon, lat, _ = GEOD.fwd(
                np.full(len(distances), lon0),
                np.full(len(distances), lat0),
                azimuth,
                np.abs(distances),
            )
            x, y = transformer.transform(lon, lat)
            line = LineString(np.column_stack([x, y]))
            direction = np.array([x[-1] - x[0], y[-1] - y[0]])
            direction /= np.linalg.norm(direction)
            crossing = line.intersection(area)
            if crossing.is_empty:
                continue
            coordinates = np.concatenate(
                [
                    np.array(part.coords)
                    for part in getattr(crossing, "geoms", [crossing])
                ]
            )
            # distances along the line, from the right end toward the left
            along_line = coordinates @ direction
            low, high = along_line.min(), along_line.max()
            if crossing.length / (high - low) < MIN_SECTION_COMPLETENESS:
                continue
            inverse = Transformer.from_crs(
                transformer.target_crs, "EPSG:4326", always_xy=True
            )
            right_lon, right_lat = inverse.transform(
                *coordinates[np.argmin(along_line)]
            )
            left_lon, left_lat = inverse.transform(*coordinates[np.argmax(along_line)])
            rows.append(
                {
                    "platform": datatake.platform,
                    "datatake": datatake.datatake,
                    "t": seconds,
                    "latitude": lat0,
                    "left": high,
                    "right": -low,
                    "left_lon": left_lon,
                    "left_lat": left_lat,
                    "right_lon": right_lon,
                    "right_lat": right_lat,
                    "altitude": frames["altitude"][k],
                }
            )
    sections = pd.DataFrame(rows)
    print(f"{len(sections)} complete cross sections", flush=True)
    return sections


def fit(day, sections, orbits, reference) -> None:
    """Measures the view angles of the edges of the valid data."""
    angles = []
    for platform, group in sections.groupby("platform"):
        frames = get_frames(orbits[platform], day, group.t.to_numpy())
        for side in ("left", "right"):
            los = normalize(
                ecef(group[f"{side}_lon"].to_numpy(), group[f"{side}_lat"].to_numpy())
                - frames["position"]
            )
            group = group.assign(
                **{
                    f"{side}_angle": np.degrees(
                        np.arctan2(
                            np.sum(los * frames["left"], axis=-1),
                            np.sum(los * frames["nadir"], axis=-1),
                        )
                    )
                }
            )
        angles.append(group)
    angles = pd.concat(angles)
    for platform, group in [("all", angles), *angles.groupby("platform")]:
        left, right = group.left_angle.median(), group.right_angle.median()
        print(
            f"{platform}: {len(group)} sections; left edge {left:.4f} deg, right edge "
            f"{right:.4f} deg, field of view {left - right:.4f} deg, roll "
            f"{(left + right) / 2:.4f} deg; width median "
            f"{(group.left + group.right).median() / 1e3:.2f} km (p5 "
            f"{(group.left + group.right).quantile(0.05) / 1e3:.2f}, p95 "
            f"{(group.left + group.right).quantile(0.95) / 1e3:.2f})"
        )
    print(
        f"imaging times: {len(reference)} grid points, miss distance p95 "
        f"{np.percentile(reference.miss, 95) / 1e3:.3f} km"
    )


def score_swath_width(day, sections, orbit, instrument) -> float:
    """
    Gets the median absolute difference (km) from the widths of the sections, of the
    distance between the rays at the left and right edges of the view.
    """
    observed = (sections.left + sections.right).to_numpy()
    track = orbit.to_gp_orbit().get_orbit_track(to_utc(day, sections.t.to_numpy()))
    left, right = (
        compute_projected_ray_position(
            track,
            instrument.cross_track_field_of_view,
            0,
            instrument.roll_angle,
            instrument.pitch_angle,
            True,
            angle,
        )
        for angle in (0, 180)
    )
    _, _, predicted = GEOD.inv(
        left.longitude.degrees,
        left.latitude.degrees,
        right.longitude.degrees,
        right.latitude.degrees,
    )
    print(
        f"swath width: observed {np.median(observed) / 1e3:.2f} km, predicted "
        f"{np.median(predicted) / 1e3:.2f} km",
        flush=True,
    )
    return float(np.median(np.abs(predicted - observed))) / 1e3


def get_chunks(datatakes: pd.DataFrame) -> list[tuple[int, float, float]]:
    """Divides the interiors of datatakes into chunks of comparable lengths."""
    chunks = []
    for k, datatake in datatakes.iterrows():
        start, stop = datatake.start + END_MARGIN, datatake.stop - END_MARGIN
        count = max(1, int(round((stop - start) / CHUNK_SECONDS)))
        length = (stop - start) / count
        if length < MIN_CHUNK_SECONDS:
            continue
        chunks += [
            (k, start + i * length, start + (i + 1) * length) for i in range(count)
        ]
    return chunks


def score_chunk(chunk, day, datatakes, tiles, orbits, satellite_for) -> dict:
    """Gets the intersection over union of a chunk, by time step."""
    k, start, stop = chunk
    datatake = datatakes.loc[k]
    orbit = orbits[datatake.platform]
    middle = get_frames(orbit, day, [(start + stop) / 2])
    lon0, lat0 = float(middle["lon"][0]), float(middle["lat"][0])
    transformer = local_projection(lon0, lat0)
    # the band between lines across the track at the start and stop of the chunk
    times = np.linspace(start, stop, max(2, int((stop - start) // 10) + 1))
    frames = get_frames(orbit, day, times)
    sides = []
    for offset in (90, -90):
        lon, lat, _ = GEOD.fwd(
            frames["lon"],
            frames["lat"],
            frames["azimuth"] + offset,
            np.full(len(times), BAND_HALF_WIDTH),
        )
        sides.append(np.column_stack(transformer.transform(lon, lat)))
    band = Polygon(np.concatenate([sides[0], sides[1][::-1]])).buffer(0)
    reach = (stop - start) * 7e3 / 2 + BAND_HALF_WIDTH + 160e3
    nearby = TileIndex(tiles[tiles.datatake == datatake.datatake]).near(
        lon0, lat0, reach
    )
    observed = unary_union(
        [project_geometry(g, transformer) for g in nearby]
    ).intersection(band)
    results = {"platform": datatake.platform, "latitude": lat0}
    duration = stop - start
    for time_step in TIME_STEPS:
        n = max(1, int(round(duration / time_step)))
        step = duration / n
        track = compute_ground_track(
            satellite_for(datatake.platform, step),
            to_utc(day, start + step * (np.arange(n) + 0.5)),
        )
        predicted = project_geometry(track.geometry.iloc[0], transformer)
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


def score_area(day, datatakes, tiles, orbits, satellite_for) -> pd.DataFrame:
    """Gets the intersection over union of each chunk, by time step."""
    chunks = get_chunks(datatakes)
    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        table = pd.DataFrame(
            list(
                pool.map(
                    lambda chunk: score_chunk(
                        chunk, day, datatakes, tiles, orbits, satellite_for
                    ),
                    chunks,
                )
            )
        )
    print(f"swath area over {len(table)} chunks:")
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


def observe_passes(points, day, datatakes, tiles, orbits) -> pd.DataFrame:
    """
    Finds the passes of points inside tiles, one per point and datatake, at the time of
    the satellite's nearest approach.
    """
    lon, lat = points.geometry.x.to_numpy(), points.geometry.y.to_numpy()
    tree = STRtree(list(tiles.geometry))
    point_index, tile_index = tree.query(make_points(lon, lat), predicate="within")
    found = pd.DataFrame(
        {"point_id": point_index, "datatake": tiles.datatake.to_numpy()[tile_index]}
    ).drop_duplicates()
    rows = []
    for datatake in datatakes.itertuples():
        ids = found.point_id[found.datatake == datatake.datatake].to_numpy()
        if len(ids) == 0:
            continue
        times = np.arange(datatake.start - 60, datatake.stop + 60, 1.0)
        frames = get_frames(orbits[datatake.platform], day, times)
        _, nearest = cKDTree(normalize(frames["position"])).query(
            normalize(ecef(lon[ids], lat[ids]))
        )
        rows.append(
            pd.DataFrame(
                {
                    "point_id": ids,
                    "platform": datatake.platform,
                    "time": pd.Timestamp(day)
                    + pd.to_timedelta(times[nearest], unit="s"),
                }
            )
        )
    return pd.concat(rows, ignore_index=True)


def score_passes(points, day, datatakes, tiles, orbits, satellite_for) -> dict:
    """Compares the passes predicted with those of lattice points inside tiles."""
    observed = observe_passes(points, day, datatakes, tiles, orbits)
    results = {}
    for platform, group in observed.groupby("platform"):
        group = group.copy()
        ids = np.sort(group.point_id.unique())
        group["point_id"] = np.searchsorted(ids, group.point_id)
        subset = points.iloc[ids].reset_index(drop=True)
        start = group.time.min().to_pydatetime() - SEARCH_WINDOW
        end = group.time.max().to_pydatetime() + SEARCH_WINDOW
        print(
            f"{platform}: predicting passes of {len(subset)} points; "
            f"{len(group)} observed",
            flush=True,
        )
        predicted = predict_passes(
            subset, satellite_for(platform, TIME_STEPS[0]), 0, start, end, []
        )
        results[platform] = match_passes(predicted, group)
        print(json.dumps(results[platform]), flush=True)
    return results


def score_epochs(day, reference, satellite_for) -> pd.DataFrame:
    """Gets the differences (s) of epochs from the imaging times of grid points."""
    sample = reference.iloc[:: max(1, len(reference) // EPOCH_POINTS)].reset_index(
        drop=True
    )

    def observe(row):
        time = day + timedelta(seconds=float(row.t))
        found = collect_observations(
            Point(latitude=float(row.lat), longitude=float(row.lon)),
            satellite_for(row.platform, TIME_STEPS[0]),
            time - SEARCH_WINDOW,
            time + SEARCH_WINDOW,
        )
        if found.empty:
            return np.nan
        differences = (found.epoch - pd.Timestamp(time)).dt.total_seconds().to_numpy()
        return differences[np.argmin(np.abs(differences))]

    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        sample["difference"] = list(pool.map(observe, sample.itertuples()))
    for platform, group in sample.groupby("platform"):
        found = np.isfinite(group.difference)
        print(
            f"epochs, {platform}: {int(found.sum())} of {len(group)} points observed, "
            f"median {group.difference[found].median():+.3f} s, p95 "
            f"{summarize(group.difference[found])['p95']:.3f} s",
            flush=True,
        )
    return sample


def validate(day, tles, entry, orbits, tiles, datatakes, sections, reference, write):
    """Scores the configuration, and writes its results if requested."""

    def satellite_for(platform: str, time_step: float):
        return entry.to_satellite(
            orbits[platform], KEY, time_step=time_step, platform=platform
        )

    metrics = []
    for platform, group in sections.groupby("platform"):
        width = score_swath_width(
            day,
            group,
            orbits[platform],
            satellite_for(platform, TIME_STEPS[0]).instruments[0],
        )
        metrics.append(
            {"name": "swath_width_median", "value": width, "scope": platform}
        )
    epochs = score_epochs(day, reference, satellite_for)
    for platform, group in epochs.groupby("platform"):
        found = group.difference[np.isfinite(group.difference)]
        metrics += [
            {"name": "epoch_p95", "value": summarize(found)["p95"], "scope": platform},
            {
                "name": "epoch_median",
                "value": summarize(found)["median"],
                "scope": platform,
            },
        ]
    lattice = generate_points_fibonacci_lattice(POINT_SPACING)
    passes = score_passes(lattice, day, datatakes, tiles, orbits, satellite_for)
    metrics += [
        {"name": "pass_recall", "value": result["pass_recall"], "scope": platform}
        for platform, result in passes.items()
    ]
    area = score_area(day, datatakes, tiles, orbits, satellite_for)
    for time_step in TIME_STEPS:
        for platform, group in area.groupby("platform"):
            metrics += [
                {
                    "name": "footprint_iou_median",
                    "value": group[f"iou_{time_step}"].median(),
                    "scope": f"{platform}, time step {time_step} s",
                },
                {
                    "name": "footprint_iou_p5",
                    "value": group[f"iou_{time_step}"].quantile(0.05),
                    "scope": f"{platform}, time step {time_step} s",
                },
            ]
    metrics = [{**m, "value": float(m["value"])} for m in metrics]
    print(json.dumps(metrics), flush=True)
    if not write:
        return
    chosen = epochs[
        (epochs.platform == FIXTURE_PLATFORM) & np.isfinite(epochs.difference)
    ]
    rng = np.random.default_rng(1)
    picks = chosen.iloc[
        rng.choice(len(chosen), min(FIXTURE_PASSES, len(chosen)), False)
    ]
    fixture = picks[["lat", "lon", "t"]].to_numpy()
    update_fixture(
        FIXTURE_PATH,
        tle=np.array(tles[FIXTURE_PLATFORM]),
        day=np.array(day.isoformat()),
        passes=fixture,
    )
    update_record(
        RECORD_PATH,
        "msi",
        KEY,
        {
            "digest": entry.configurations[KEY].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": "docs/validation/validate_catalog_msi.py",
            "reference": "Valid data footprints of Sentinel-2 Level-2A tiles; viewing "
            + "angle grids and GPS ephemerides of their metadata",
            "period": f"{day:%Y-%m-%d}",
            "platform": ", ".join(sorted(PLATFORMS)),
            "holdout_gap_days": get_holdout_gap(day, FIT_DAY),
            "metrics": metrics,
            "fixture_metrics": recompute_passes(
                satellite_for(FIXTURE_PLATFORM, TIME_STEPS[0]), 0, day, fixture
            ),
        },
    )
    print(f"wrote {RECORD_PATH}")


def main() -> None:
    """Fits or validates the entry on a day of data."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("day")
    parser.add_argument("--fit", action="store_true")
    parser.add_argument("--write", action="store_true")
    parser.add_argument(
        "--download-only",
        action="store_true",
        help="read and cache the reference data, without scoring",
    )
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)

    day = datetime.fromisoformat(arguments.day).replace(tzinfo=timezone.utc)
    tiles = load_tiles(day, arguments.directory)
    datastrips = load_datastrips(day, arguments.directory, tiles)
    datatakes = get_datatakes(day, tiles, datastrips)
    ephemeris, grids = load_reference_times(day, arguments.directory, tiles)
    print(
        f"{len(datatakes)} datatakes, {len(ephemeris)} ephemeris points, "
        f"{len(grids)} grid points",
        flush=True,
    )
    if arguments.download_only:
        return
    tles = {platform: get_element_set(day, platform) for platform in PLATFORMS}
    orbits = {
        platform: GeneralPerturbationsOrbit.from_tle(tle)
        for platform, tle in tles.items()
    }
    for platform, tle in tles.items():
        print(f"{platform}, {arguments.day}: {tle[0][:32]}", flush=True)
    reference = recover_imaging_times(grids, ephemeris)
    sections = measure_sections(day, datatakes, tiles, orbits)
    if arguments.fit:
        fit(day, sections, orbits, reference)
        return
    validate(
        day,
        tles,
        get_entry("msi"),
        orbits,
        tiles,
        datatakes,
        sections,
        reference,
        arguments.write,
    )


if __name__ == "__main__":
    main()

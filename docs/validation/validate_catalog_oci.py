"""
Validates the OCI catalog entry against the PACE OCI Level 1B product, using the
exact values of the entry through public TAT-C methods. OCI's telescope rotates
across track, tilted fore or aft by about 20 deg to avoid sun glint, and the tilt
changes twice an orbit, at a latitude that changes by day; TAT-C represents it as a
schedule of daily pitch profiles over the argument of latitude (`TiltSchedule`). The
entry's angles are fitted to the positions of the samples on one day (run with
`--fit`), from the spacecraft positions the product records, so the fit needs no
element set. The schedule of the change is fitted with `--switches` to the days of a
period, read with range requests of only the tilts and spacecraft positions, and the
interval around each change in which the product has no valid scans with `--blind`,
at full resolution, on the days given. Running
on a day at least 30 days from the fitting data tests the fitted values on data held
out from the fit.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`. Only the
positions, times, and tilts are read, with range requests (about 1.2 GB a day):

    python docs/validation/validate_catalog_oci.py <data directory> 2026-10-04 --fit
    python docs/validation/validate_catalog_oci.py <data directory> 2026-08-15 --switches 2026-10-08
    python docs/validation/validate_catalog_oci.py <data directory> 2026-08-20 --blind 2026-09-15 2026-10-04 2026-10-07
    python docs/validation/validate_catalog_oci.py <data directory> 2026-11-09
    python docs/validation/validate_catalog_oci.py <data directory> 2026-11-10

A day is reduced to every sixth scan (about 1 s) and every 16th sample (and the last).
The product covers the sunlit part of each orbit, whose dark part is a gap.

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import concurrent.futures
import json
import pathlib
from datetime import datetime, time, timedelta, timezone

import earthaccess
import h5py
import numpy as np
import pandas as pd

from tatc.catalog import get_entry, get_subsolar_latitude
from tatc.schemas import GeneralPerturbationsOrbit
from gesdisc_scanner import Product, Reference, Setup, load_day, validate_swath
from validation_io import (
    TO_GEODETIC,
    RemoteFile,
    ecef,
    get_with_retries,
    normalize,
    up,
)

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
# the angles are fitted to 2026-10-04, and the schedule of the tilt's change to
# 2026-08-15 through 2026-10-08, the last day of fitting data
FIT_DAY = datetime(2026, 10, 8, tzinfo=timezone.utc)
# public element sets from CelesTrak, with epochs after each day of data
ELEMENTS = {
    "2026-11-04": (
        "1 58928U 24025A   26312.25077802  .00000543  00000+0  10874-3 0  9994",
        "2 58928  98.1039 241.4756 0002580 132.8558 227.2862 14.65643827147030",
    ),
    "2026-11-07": (
        "1 58928U 24025A   26313.27483190  .00000501  00000+0  10091-3 0  9993",
        "2 58928  98.1040 242.4857 0002666 129.8006 230.3433 14.65644467147180",
    ),
    "2026-11-08": (
        "1 58928U 24025A   26313.27483190  .00000501  00000+0  10091-3 0  9993",
        "2 58928  98.1040 242.4857 0002666 129.8006 230.3433 14.65644467147180",
    ),
    "2026-11-09": (
        "1 58928U 24025A   26314.09407429  .00000495  00000+0  99766-4 0  9994",
        "2 58928  98.1042 243.2937 0002701 128.3012 231.8434 14.65645200147308",
    ),
    "2026-11-10": (
        "1 58928U 24025A   26315.52774695  .00000494  00000+0  99709-4 0  9991",
        "2 58928  98.1045 244.7079 0002770 125.5620 234.5842 14.65646728147516",
    ),
}
SAMPLES = 1272
SCAN_STRIDE = 6
SCAN_PERIOD = 0.17547  # s
COLUMNS = np.r_[0:SAMPLES:16, SAMPLES - 1]
INCLINATION = 98.0  # deg, approximately, to get the argument of latitude


def read_oci(_product, url: str, session, day: datetime, _directory) -> dict:
    """
    Reads the positions, times, tilts, and spacecraft positions of a granule's scans
    with range requests, and reduces them to every sixth scan and some samples.
    """
    with h5py.File(RemoteFile(url, session, block_size=2**18)) as f:
        start = pd.Timestamp(f.attrs["time_coverage_start"].decode())
        seconds = f["scan_line_attributes/time"][:].astype(float)
        lat = f["geolocation_data/latitude"][:]
        lon = f["geolocation_data/longitude"][:]
        position = f["navigation_data/orb_pos"][:].astype(float) * 1e3
        tilt = f["navigation_data/tilt_angle"][:].astype(float)
        quality = f["scan_line_attributes/scan_quality_flags"][:]
    # the times are seconds of the day on which the granule starts
    seconds = seconds + (start.normalize() - pd.Timestamp(day)).total_seconds()
    seconds = np.where(seconds < seconds[0] - 43200, seconds + 86400, seconds)
    keep = np.arange(0, len(seconds), SCAN_STRIDE)
    sc_lon, sc_lat, sc_alt = TO_GEODETIC.transform(*position[keep].T)
    return {
        "time": seconds[keep],
        "lat": lat[keep][:, COLUMNS],
        "lon": lon[keep][:, COLUMNS],
        "sc_lat": np.asarray(sc_lat),
        "sc_lon": np.asarray(sc_lon),
        "sc_alt": np.asarray(sc_alt),
        "tilt": tilt[keep],
        "quality": quality[keep],
    }


PRODUCT = Product(
    short_name="PACE_OCI_L1B_SCI",
    version="3",
    read=read_oci,
    host="https://obdaac-tea",
)
SETUP = Setup(
    entry_id="oci",
    platform="pace",
    swath="swath",
    beams="",
    size_instrument="OCI",
    published_footprints=(),
    record_path=REPOSITORY / "src" / "tatc" / "resources" / "validation" / "oci.json",
    fixture_path=REPOSITORY / "tests" / "catalog" / "fixtures" / "oci.npz",
    source="docs/validation/validate_catalog_oci.py",
    reference="PACE_OCI_L1B_SCI version 3 (sample positions and times)",
    published_reference="",
    fit_day=FIT_DAY,
    time_steps=(10, 30, 120),
    pass_radius=40,
    point_spacing=200e3,
)


def load_reference(day: datetime, directory: pathlib.Path):
    """Loads a reduced day, with its tilts."""
    data = load_day(PRODUCT, day, directory)
    tilt = data.pop("tilt")
    data.update(pixel=COLUMNS, pixels=SAMPLES, line_stride=SCAN_STRIDE)
    reference = Reference(day, data, SCAN_STRIDE * SCAN_PERIOD)
    valid = np.all(np.abs(data["lat"]) <= 90, axis=1) & (data["quality"] == 0)
    reference.tilt = tilt[valid]
    return reference


def argument_of_latitude(lat, ascending) -> np.ndarray:
    """Gets the approximate argument of latitude (deg) from latitude and direction."""
    u = np.degrees(
        np.arcsin(
            np.clip(np.sin(np.radians(lat)) / np.sin(np.radians(INCLINATION)), -1, 1)
        )
    )
    return np.where(ascending, u % 360, 180 - u)


def fit(reference: Reference) -> None:
    """
    Prints, for each tilt, the pitch (from the along-track tangents), the edges of the
    scan in the tilted view (from the cross-track tangents), and the yaw (the slope of
    the along-track on the cross-track tangents), in both frames, and where the tilt
    changes.
    """
    lon, lat, _ = TO_GEODETIC.transform(*reference.position.T)
    lat = np.asarray(lat)
    zenith = up(lon, lat)
    rays = ecef(reference.lon, reference.lat) - reference.position[:, np.newaxis]
    down = -np.sum(rays * zenith[:, np.newaxis], axis=-1)
    for frame in ("earth_fixed", "inertial"):
        motion = reference.velocity(frame)
        along_axis = normalize(
            motion - np.sum(motion * zenith, axis=-1)[:, None] * zenith
        )
        left_axis = np.cross(zenith, along_axis)
        along = np.sum(rays * along_axis[:, np.newaxis], axis=-1) / down
        cross = np.sum(rays * left_axis[:, np.newaxis], axis=-1) / down
        # orient each scan so that its samples run from right to left
        flipped = cross[:, 0] > cross[:, -1]
        along[flipped], cross[flipped] = along[flipped, ::-1], cross[flipped, ::-1]
        for name, chosen in (("fore", reference.tilt > 0), ("aft", reference.tilt < 0)):
            if not np.any(chosen):
                continue
            pitch = np.degrees(np.arctan(np.median(along[chosen])))
            cosine = np.cos(np.radians(pitch))
            right = -np.degrees(np.arctan(np.median(cross[chosen, 0]) * cosine))
            left = np.degrees(np.arctan(np.median(cross[chosen, -1]) * cosine))
            centered = cross[chosen] - cross[chosen].mean(axis=1, keepdims=True)
            slope = np.sum(centered * along[chosen], axis=1) / np.sum(
                centered**2, axis=1
            )
            print(
                f"{frame}, {name} ({chosen.sum()} scans, tilt "
                f"{np.median(reference.tilt[chosen]):+.3f}): pitch {pitch:+.4f} deg; "
                f"edges {left:.4f} deg left, {right:.4f} deg right, field of view "
                f"{left + right:.4f}, roll {(left - right) / 2:+.4f}; yaw slope "
                f"{np.degrees(np.arctan(np.median(slope))):+.4f} deg (p5 "
                f"{np.degrees(np.arctan(np.percentile(slope, 5))):+.4f}, p95 "
                f"{np.degrees(np.arctan(np.percentile(slope, 95))):+.4f}); "
                f"first sample on the left in {np.mean(flipped[chosen]):.1%}"
            )
    ascending = np.gradient(lat) > 0
    u = argument_of_latitude(lat, ascending)
    sign = np.sign(reference.tilt)
    changes = np.nonzero((np.diff(sign) != 0) & (np.diff(reference.seconds) < 30))[0]
    for i in changes[:40]:
        print(
            f"tilt {reference.tilt[i]:+.2f} -> {reference.tilt[i + 1]:+.2f} at "
            f"{reference.seconds[i] / 3600:.3f} h, latitude {lat[i]:+.3f}, "
            f"{'ascending' if ascending[i] else 'descending'}, argument of latitude "
            f"{u[i]:.3f} deg"
        )
    if len(changes):
        print(
            f"{len(changes)} changes; argument of latitude median "
            f"{np.median(u[changes]):.3f}, p5 {np.percentile(u[changes], 5):.3f}, "
            f"p95 {np.percentile(u[changes], 95):.3f} deg"
        )
    print(
        f"scan period: median {np.median(np.diff(reference.seconds)) / SCAN_STRIDE:.5f} s"
    )


def read_switches(url: str, session) -> list[tuple[str, float]]:
    """
    Reads the times and geodetic latitudes of the spacecraft at which the tilt
    changes from aft to fore in a granule, ascending, with range requests of only the
    scans' times, tilts, and positions.
    """
    with h5py.File(RemoteFile(url, session)) as f:
        start = pd.Timestamp(f.attrs["time_coverage_start"].decode())
        seconds = f["scan_line_attributes/time"][:].astype(float)
        tilt = f["navigation_data/tilt_angle"][:].astype(float)
        position = f["navigation_data/orb_pos"][:].astype(float) * 1e3
    _, lat, _ = TO_GEODETIC.transform(*position.T)
    lat = np.asarray(lat)
    switches = []
    for k in np.nonzero(np.diff(np.sign(tilt)) > 0)[0]:
        if seconds[k] >= 0 and np.isfinite(lat[k]) and lat[k + 1] > lat[k]:
            time_of = start.normalize() + pd.Timedelta(seconds=seconds[k])
            switches.append((time_of.isoformat(), float(lat[k])))
    return switches


def find_switch_granules(day: datetime) -> list[str]:
    """
    Finds the granules of a day that cross latitudes within 12 deg of the subsolar
    latitude, where the tilt changes to fore.
    """
    subsolar = get_subsolar_latitude(day + timedelta(hours=12))
    items = get_with_retries(
        "https://cmr.earthdata.nasa.gov/search/granules.umm_json",
        params={
            "short_name": PRODUCT.short_name,
            "version": PRODUCT.version,
            "temporal": f"{day:%Y-%m-%dT%H:%M:%SZ},"
            f"{day + timedelta(days=1):%Y-%m-%dT%H:%M:%SZ}",
            "bounding_box": f"-180,{subsolar - 12},180,{subsolar + 12}",
            "page_size": 2000,
        },
    ).json()["items"]
    return [
        next(
            u["URL"]
            for u in item["umm"]["RelatedUrls"]
            if u["URL"].startswith(PRODUCT.host) and u["URL"].endswith(".nc")
        )
        for item in items
    ]


def read_blind(url: str, session) -> list[tuple[float, float]]:
    """
    Reads, for each change of tilt from aft to fore in a granule, the times
    (seconds) of the last valid scan before the tilt crosses zero and of the first
    valid scan after, from the middle of the change, at full resolution.
    """
    with h5py.File(RemoteFile(url, session)) as f:
        seconds = f["scan_line_attributes/time"][:].astype(float)
        tilt = f["navigation_data/tilt_angle"][:].astype(float)
        quality = f["scan_line_attributes/scan_quality_flags"][:]
    changes = []
    for k in np.nonzero(np.diff(np.sign(tilt)) > 0)[0]:
        middle = seconds[k] + (seconds[k + 1] - seconds[k]) * (
            -tilt[k] / (tilt[k + 1] - tilt[k])
        )
        offset = seconds - middle
        near = (np.abs(offset) < 90) & (quality == 0)
        before, after = offset[near & (offset < 0)], offset[near & (offset > 0)]
        # only changes with valid scans well before and after, within the granule
        if len(before) and len(after) and before.min() < -30 and after.max() > 30:
            changes.append((float(before.max()), float(after.min())))
    return changes


def fit_blind(days: list[datetime]) -> None:
    """
    Prints the interval around the changes of tilt to fore in which the product has
    no valid scans, over days, at full resolution: from the end of the last valid
    scan (its start plus a scan period) to the start of the first.
    """
    if earthaccess.login(strategy="environment") is None:
        raise RuntimeError("Earthdata login failed: is EARTHDATA_TOKEN set?")
    session = earthaccess.get_requests_https_session()
    changes = []
    for day in days:
        with concurrent.futures.ThreadPoolExecutor(8) as pool:
            for part in pool.map(
                lambda url: read_blind(url, session), find_switch_granules(day)
            ):
                changes += part
    before, after = np.array(changes).T
    print(f"{len(changes)} changes from aft to fore on {len(days)} days")
    print(
        f"last valid scan: median {np.median(before):+.3f} s (from "
        f"{before.min():+.3f} to {before.max():+.3f}); first valid scan: median "
        f"{np.median(after):+.3f} s (from {after.min():+.3f} to {after.max():+.3f})"
    )
    print(
        f"blind from {-(np.median(before) + SCAN_PERIOD):.2f} s before the middle of "
        f"the change to {np.median(after):.2f} s after"
    )


def load_switches(day: datetime, directory: pathlib.Path, session) -> list:
    """
    Reads the changes of tilt of the granules of a day that cross latitudes within
    12 deg of the subsolar latitude, or loads them from a cache.
    """
    cache = directory / f"switch_{day:%Y%m%d}.json"
    if not cache.exists():
        urls = find_switch_granules(day)
        with concurrent.futures.ThreadPoolExecutor(8) as pool:
            parts = list(pool.map(lambda url: read_switches(url, session), urls))
        cache.write_text(
            json.dumps(
                {"granules": len(urls), "switches": sorted(s for p in parts for s in p)}
            )
        )
    return json.loads(cache.read_text())["switches"]


def fit_switches(first: datetime, last: datetime, directory: pathlib.Path) -> None:
    """
    Prints the offsets of the latitude at which the tilt changes to fore from the
    subsolar latitude at noon (UTC) of each day, by the parity of the day of the year,
    over a period of days, and the changes far from their day's median.
    """
    if earthaccess.login(strategy="environment") is None:
        raise RuntimeError("Earthdata login failed: is EARTHDATA_TOKEN set?")
    session = earthaccess.get_requests_https_session()
    rows = []
    for i in range((last - first).days + 1):
        rows += load_switches(first + timedelta(days=i), directory, session)
    switches = pd.DataFrame(rows, columns=["time", "lat"]).drop_duplicates("time")
    switches["time"] = pd.to_datetime(switches.time)
    # granules near midnight are found on both days; each change counts on its day
    switches["day"] = switches.time.dt.date
    switches = switches[
        (switches.day >= first.date()) & (switches.day <= last.date())
    ].copy()
    daily = switches.groupby("day").lat.median()
    subsolar = np.array(
        [
            get_subsolar_latitude(datetime.combine(d, time(12), tzinfo=timezone.utc))
            for d in daily.index
        ]
    )
    offset = daily.to_numpy() - subsolar
    odd = np.array([d.timetuple().tm_yday % 2 == 1 for d in daily.index])
    print(f"{len(daily)} days, {len(switches)} changes from aft to fore")
    for name, chosen in (("odd", odd), ("even", ~odd)):
        print(
            f"{name} days of the year ({chosen.sum()}): offset from the subsolar "
            f"latitude median {np.median(offset[chosen]):+.4f} deg (standard deviation "
            f"{offset[chosen].std():.4f}, from {offset[chosen].min():+.4f} to "
            f"{offset[chosen].max():+.4f})"
        )
    print(
        "every odd day's offset exceeds every even day's: "
        f"{bool(offset[odd].min() > offset[~odd].max())}"
    )
    switches["from_median"] = switches.lat - switches.day.map(daily)
    far = switches[np.abs(switches.from_median) > 0.2]
    for row in far.itertuples():
        print(f"  a change {row.from_median:+.3f} deg from its day's at {row.time}")


def main() -> None:
    """Validates the entry on a day of data, or fits its angles or schedule."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("day")
    parser.add_argument("--fit", action="store_true")
    parser.add_argument("--switches", metavar="LAST_DAY")
    parser.add_argument("--blind", nargs="*", metavar="OTHER_DAY")
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)
    day = datetime.fromisoformat(arguments.day).replace(tzinfo=timezone.utc)
    if arguments.blind is not None:
        fit_blind(
            [day]
            + [
                datetime.fromisoformat(d).replace(tzinfo=timezone.utc)
                for d in arguments.blind
            ]
        )
        return
    if arguments.switches:
        last = datetime.fromisoformat(arguments.switches).replace(tzinfo=timezone.utc)
        fit_switches(day, last, arguments.directory)
        return
    reference = load_reference(day, arguments.directory)
    print(
        f"{arguments.day}: {len(reference.time)} scans from {reference.time[0]}, "
        f"{len(reference.gaps)} gaps over a minute",
        flush=True,
    )
    if arguments.fit:
        fit(reference)
        return
    tle = list(ELEMENTS[arguments.day])
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    entry = get_entry("oci")
    validate_swath(SETUP, day, tle, entry, orbit, reference, arguments.write)
    if arguments.write:
        print(f"wrote {SETUP.record_path} and {SETUP.fixture_path}")


if __name__ == "__main__":
    main()

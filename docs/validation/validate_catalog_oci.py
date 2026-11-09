"""
Validates the OCI catalog entry against the PACE OCI Level 1B product, using the
exact values of the entry through public TAT-C methods. OCI's telescope rotates
across track, tilted fore or aft by about 20 deg to avoid sun glint, and the tilt
changes once an orbit; TAT-C represents it as a pitch profile over the argument of
latitude. The entry's angles are fitted to the positions of the samples on one day
(run with `--fit`), from the spacecraft positions the product records, so the fit
needs no element set; running on a day at least 30 days from it tests the fitted
angles on data held out from the fit.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`. Only the
positions, times, and tilts are read, with range requests (about 1.2 GB a day):

    python docs/validation/validate_catalog_oci.py <data directory> 2026-10-04 --fit
    python docs/validation/validate_catalog_oci.py <data directory> 2026-11-04

A day is reduced to every sixth scan (about 1 s) and every 16th sample (and the last).
The product covers the sunlit part of each orbit, whose dark part is a gap.

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import pathlib
from datetime import datetime, timezone

import h5py
import numpy as np
import pandas as pd

from tatc.catalog import get_entry
from tatc.schemas import GeneralPerturbationsOrbit
from gesdisc_scanner import Product, Reference, Setup, load_day, validate_swath
from validation_io import TO_GEODETIC, RemoteFile, ecef, normalize, up

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
FIT_DAY = datetime(2026, 10, 4, tzinfo=timezone.utc)
# public element sets from CelesTrak, with epochs after each day of data
ELEMENTS = {
    "2026-11-04": (
        "1 58928U 24025A   26312.25077802  .00000543  00000+0  10874-3 0  9994",
        "2 58928  98.1039 241.4756 0002580 132.8558 227.2862 14.65643827147030",
    ),
}
SAMPLES = 1272
SCAN_STRIDE = 6
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
    reference = Reference(day, data, SCAN_STRIDE * 0.17547)
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


def main() -> None:
    """Validates the entry on a day of data, or fits its angles."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("day")
    parser.add_argument("--fit", action="store_true")
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)
    day = datetime.fromisoformat(arguments.day).replace(tzinfo=timezone.utc)
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

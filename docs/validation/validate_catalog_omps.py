"""
Validates the OMPS Nadir Mapper catalog entry against the geolocation of the OMPS
total column product (`OMPS-TC-GEO`) of NOAA-20 and NOAA-21, using the exact values
of the entry through public TAT-C methods. Each platform's scan angles are fitted to
the positions of the samples on one day (run with `--fit`), which needs no element
set: running on a day at least 30 days from it tests the entry on data held out from
the fit, with NOAA's operational element set of the day.

Reads anonymous public data, and needs the `validation` extra:

    python docs/validation/validate_catalog_omps.py <data directory> 2026-10-01 --platform noaa-21 --fit
    python docs/validation/validate_catalog_omps.py <data directory> 2026-11-05 --platform noaa-21

A day is reduced to every second line (2.5 s on NOAA-21, 5 s on NOAA-20) and every
fourth sample (and the last). The Nadir Mapper images the sunlit part of each orbit, whose dark part is a
gap.

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import concurrent.futures
import io
import pathlib
from datetime import datetime, timezone

import h5py
import numpy as np

from tatc.catalog import get_entry
from tatc.schemas import GeneralPerturbationsOrbit
from gesdisc_scanner import Reference, Setup, fit, fit_yaw, validate_swath
from validation_io import (
    TO_GEODETIC,
    fetch,
    iet_to_utc,
    list_keys,
    load_element_set,
    use_archive,
)

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
FIXTURES = REPOSITORY / "tests" / "catalog" / "fixtures"
GROUP = "All_Data/OMPS-TC-GEO_All/"
LINE_STRIDE, SAMPLE_STRIDE = 2, 4
# the day to which each platform's angles were fitted, and its configuration and
# fixture
FIT_DAYS = {
    "noaa-20": datetime(2026, 10, 1, tzinfo=timezone.utc),
    "noaa-21": datetime(2026, 10, 1, tzinfo=timezone.utc),
}
CONFIGURATIONS = {
    "noaa-21": ("swath_noaa_21", FIXTURES / "omps_nm_noaa_21.npz"),
    "noaa-20": ("swath_noaa_20", FIXTURES / "omps_nm_noaa_20.npz"),
}


def read_granule(key: str) -> dict | None:
    """
    Reads a granule's lines: their middle times, the positions of their samples, the
    spacecraft's position, and their quality flags.
    """
    with h5py.File(io.BytesIO(fetch(key))) as f:
        g = f[GROUP]
        samples = int(g["NumberOfIFOVs"][0])
        if samples <= 0 or len(g["MidTime"]) == 0:
            return None
        columns = np.r_[0:samples:SAMPLE_STRIDE, samples - 1]
        return {
            "mid_iet": g["MidTime"][:],
            "lat": g["Latitude"][:, :samples][:, columns],
            "lon": g["Longitude"][:, :samples][:, columns],
            "position": g["SCPosition"][:].astype(float),
            "quality": g["QF1_OMPSTCGEO"][:],
            "samples": np.full(len(g["MidTime"]), samples),
        }


def load_day(day: datetime, directory: pathlib.Path, prefix: str) -> dict:
    """Downloads and reduces the geolocation of a day, or loads it from a cache."""
    path = directory / f"{prefix}_omps_gotco_{day:%Y%m%d}.npz"
    if not path.exists():
        keys = list_keys(f"OMPS-TC-GEO/{day:%Y/%m/%d}/")
        with concurrent.futures.ThreadPoolExecutor(8) as pool:
            parts = [p for p in pool.map(read_granule, keys) if p is not None]
        samples = {int(p["samples"][0]) for p in parts}
        if len(samples) != 1:
            raise ValueError(f"granules differ in their number of samples: {samples}")
        data = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
        order = np.argsort(data["mid_iet"])
        _, first = np.unique(data["mid_iet"][order], return_index=True)
        keep = order[first][::LINE_STRIDE]
        np.savez_compressed(path, **{k: v[keep] for k, v in data.items()})
    return dict(np.load(path))


def load_reference(day: datetime, directory: pathlib.Path, prefix: str) -> Reference:
    """Loads a reduced day as a reference of scans."""
    raw = load_day(day, directory, prefix)
    seconds = (iet_to_utc(raw["mid_iet"]) - day).total_seconds().to_numpy()
    sc_lon, sc_lat, sc_alt = TO_GEODETIC.transform(*raw["position"].T)
    samples = int(raw["samples"][0])
    inside = (seconds >= 0) & (seconds < 86400)
    data = {
        "time": seconds[inside],
        "lat": raw["lat"][inside],
        "lon": raw["lon"][inside],
        "sc_lat": np.asarray(sc_lat)[inside],
        "sc_lon": np.asarray(sc_lon)[inside],
        "sc_alt": np.asarray(sc_alt)[inside],
        "quality": raw["quality"][inside],
        "pixel": np.r_[0:samples:SAMPLE_STRIDE, samples - 1],
        "pixels": samples,
        "line_stride": LINE_STRIDE,
    }
    # the reporting period differs by platform (1.25 s on NOAA-21, 2.5 s on NOAA-20)
    return Reference(day, data, float(np.median(np.diff(data["time"]))))


def make_setup(platform: str) -> Setup:
    """Gets what a run on a platform validates, and where it writes its results."""
    swath, fixture = CONFIGURATIONS[platform]
    return Setup(
        entry_id="omps_nm",
        platform=platform,
        swath=swath,
        beams="",
        size_instrument="OMPS Nadir Mapper",
        published_footprints=(),
        record_path=REPOSITORY
        / "src"
        / "tatc"
        / "resources"
        / "validation"
        / "omps_nm.json",
        fixture_path=fixture,
        source="docs/validation/validate_catalog_omps.py",
        reference="OMPS-TC-GEO (sample positions and times), NOAA Open Data "
        "Dissemination",
        published_reference="",
        fit_day=FIT_DAYS[platform],
        time_steps=(10, 30, 120),
        pass_radius=40,
        point_spacing=200e3,
    )


def main() -> None:
    """Validates the entry on a day of data, or fits its angles."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("day")
    parser.add_argument("--platform", default="noaa-21", choices=sorted(FIT_DAYS))
    parser.add_argument("--fit", action="store_true")
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)
    prefix = use_archive(arguments.platform)
    day = datetime.fromisoformat(arguments.day).replace(tzinfo=timezone.utc)
    reference = load_reference(day, arguments.directory, prefix)
    print(
        f"{arguments.day} ({arguments.platform}): {len(reference.time)} lines from "
        f"{reference.time[0]}, {reference.pixels} samples, {len(reference.gaps)} gaps "
        "over a minute",
        flush=True,
    )
    if arguments.fit:
        for frame in ("inertial", "earth_fixed"):
            fit(reference, frame)
            fit_yaw(reference, frame)
        return
    tle = load_element_set(day, arguments.directory)
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    setup = make_setup(arguments.platform)
    validate_swath(
        setup, day, tle, get_entry("omps_nm"), orbit, reference, arguments.write
    )
    if arguments.write:
        print(f"wrote {setup.record_path} and {setup.fixture_path}")


if __name__ == "__main__":
    main()

"""
Validates the OMPS Nadir Mapper catalog entry against the geolocation of the OMPS
total column product (`OMPS-TC-GEO`) of NOAA-20, NOAA-21, and Suomi NPP (retired on
2026-11-02, whose archive ends that day), using the exact values
of the entry through public TAT-C methods. Each platform's scan angles are fitted to
the positions of the samples on one day (run with `--fit`), which needs no element
set: running on a day at least 30 days from it tests the entry on data held out from
the fit, with NOAA's operational element set of the day.

Reads anonymous public data, and needs the `validation` extra:

    python validation/catalog/validate_catalog_omps.py <data directory> 2026-10-01 --platform noaa-21 --fit
    python validation/catalog/validate_catalog_omps.py <data directory> 2026-11-05 --platform noaa-21
    python validation/catalog/validate_catalog_omps.py <data directory> 2026-10-01 --platform suomi-npp --fit
    python validation/catalog/validate_catalog_omps.py <data directory> 2026-08-25 --platform suomi-npp

A day is reduced to every second line (2.5 s on NOAA-21, 5 s on NOAA-20, 15 s on
Suomi NPP) and every
fourth sample (and the last). The Nadir Mapper images the sunlit part of each orbit, whose dark part is a
gap.

Without a day, it takes that of the platform's record (or, with `--fit`, of the
fit); the platform defaults to NOAA-21.

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import concurrent.futures
import functools
import io
import pathlib
from datetime import datetime

import h5py
import numpy as np

from tatc.schemas import GeneralPerturbationsOrbit
from scanner import Reference, Setup, fit, fit_yaw, validate_swath
from validation_io import (
    ARCHIVES,
    FIXTURES,
    TO_GEODETIC,
    Archive,
    iet_to_utc,
    parse_run,
)

GROUP = "All_Data/OMPS-TC-GEO_All/"
LINE_STRIDE, SAMPLE_STRIDE = 2, 4
# the configuration of each platform's fit and its fixture
CONFIGURATIONS = {
    "noaa-21": ("swath_noaa_21", FIXTURES / "omps_nm_noaa_21.npz"),
    "noaa-20": ("swath_noaa_20", FIXTURES / "omps_nm_noaa_20.npz"),
    "suomi-npp": ("swath_suomi_npp", FIXTURES / "omps_nm_suomi_npp.npz"),
}


def read_granule(archive: Archive, key: str) -> dict | None:
    """
    Reads a granule's lines: their middle times, the positions of their samples, the
    spacecraft's position, and their quality flags.
    """
    with h5py.File(io.BytesIO(archive.fetch(key))) as f:
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


def load_day(day: datetime, directory: pathlib.Path, archive: Archive) -> dict:
    """Downloads and reduces the geolocation of a day, or loads it from a cache."""
    path = directory / f"{archive.prefix}_omps_gotco_{day:%Y%m%d}.npz"
    if not path.exists():
        keys = archive.list_keys(f"OMPS-TC-GEO/{day:%Y/%m/%d}/")
        with concurrent.futures.ThreadPoolExecutor(8) as pool:
            parts = [
                p
                for p in pool.map(functools.partial(read_granule, archive), keys)
                if p is not None
            ]
        samples = {int(p["samples"][0]) for p in parts}
        if len(samples) != 1:
            raise ValueError(f"granules differ in their number of samples: {samples}")
        data = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
        order = np.argsort(data["mid_iet"])
        _, first = np.unique(data["mid_iet"][order], return_index=True)
        keep = order[first][::LINE_STRIDE]
        np.savez_compressed(path, **{k: v[keep] for k, v in data.items()})
    return dict(np.load(path))


def load_reference(
    day: datetime, directory: pathlib.Path, archive: Archive
) -> Reference:
    """Loads a reduced day as a reference of scans."""
    raw = load_day(day, directory, archive)
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
    # the reporting period differs by platform (1.25 s on NOAA-21, 2.5 s on NOAA-20,
    # 7.5 s on Suomi NPP)
    return Reference(day, data, float(np.median(np.diff(data["time"]))))


def make_setup(platform: str) -> Setup:
    """Gets what a run on a platform validates, and where it writes its results."""
    swath, fixture = CONFIGURATIONS[platform]
    return Setup(
        swath=swath,
        beams=None,
        size_instrument=None,
        published_footprints=(),
        fixture_path=fixture,
        reference="OMPS-TC-GEO (sample positions and times), NOAA Open Data "
        "Dissemination",
        published_reference="",
        time_steps=(10, 30, 120),
        pass_radius=40,
        point_spacing=200e3,
    )


def main() -> None:
    """Validates the entry on a day of data, or fits its angles."""
    run = parse_run("omps_nm", __doc__)
    archive = ARCHIVES[run.platform]
    reference = load_reference(run.day, run.directory, archive)
    print(
        f"{run.day:%Y-%m-%d} ({run.platform}): {len(reference.time)} lines from "
        f"{reference.time[0]}, {reference.pixels} samples, {len(reference.gaps)} gaps "
        "over a minute",
        flush=True,
    )
    if run.fit:
        for frame in ("inertial", "earth_fixed"):
            fit(reference, frame)
            fit_yaw(reference, frame)
        return
    tle = run.element_set()
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    # a transfer tests the configuration of the default platform
    entry = run.entry
    default = entry.configurations[entry.default_configuration].platforms[0]
    setup = make_setup(default if run.transfer else run.platform)
    if setup.swath in run.configurations:
        validate_swath(setup, run, reference, orbit, tle)


if __name__ == "__main__":
    main()

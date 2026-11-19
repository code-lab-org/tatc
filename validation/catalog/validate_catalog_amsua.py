"""
Validates the AMSU-A catalog entry against the AMSU-A Level 1B product of Aqua (part
of the AIRS project's products), using the exact values of the entry through public
TAT-C methods. The entry's scan angles are fitted to the positions of the footprints
on one day (run with `--fit`), which needs no element set: running on a day at least
30 days from it tests the entry on data held out from the fit.

Needs the `validation` extra, `pyhdf`, and an Earthdata token in `EARTHDATA_TOKEN`.
The product's files are HDF4 (about 0.5 MB each), so they are downloaded whole:

    python validation/catalog/validate_catalog_amsua.py <data directory> 2026-10-04 --fit
    python validation/catalog/validate_catalog_amsua.py <data directory> 2026-11-05

Without a day, it takes that of the entry's record (or, with `--fit`, of the fit).
Aqua's element sets are those registered in `elements/element_sets.yml`, or those
of a file given with `--elements`.

It also validates the entry's sample size configuration, following
`ValidateSampleSizes`. The gated metrics compare sizes with the footprint published
in the product's documentation, and the projected pixels with the footprints'
positions. `--configuration footprints` (or `swath`) tests one configuration.

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

With `--platform metop-b` or `--platform metop-c`, it tests the same values on Metop,
from the AMSU-A Level 1B product of the EUMETSAT Data Store (EO:EUM:DAT:METOP:AMSUL1,
in EPS native format, whose measurement records are read as versions 3 and 4 of the
record of the ATOVS Level 1b Product Guide, EUM/OPS-EPS/MAN/04/0030, which place the
fields read alike), and adds the results to
the record as a transfer, which does not count toward the tier. Its records give a
time for each scan (the start of the record) but none for each footprint, and no
position of the spacecraft other than its altitude: the spacecraft's ground position
is taken as the midpoint of the two central footprints, as AMSU-A's 30 footprints are
symmetric about nadir. It needs a login to the Data Store (`eumdac auth
--interactive`).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import pathlib
import threading
from datetime import datetime, timedelta

import numpy as np
import pandas as pd

from tatc.schemas import GeneralPerturbationsOrbit
from scanner import (
    Product,
    Reference,
    Setup,
    fit,
    load_day,
    tai93_to_seconds,
    validate_beams,
    validate_swath,
)
from validation_io import (
    FIXTURES,
    TO_GEODETIC,
    download_eumetsat,
    ecef,
    get_with_retries,
    parse_run,
)

SCAN_PERIOD = 8.0
# the HDF4 library is not thread-safe
HDF4_LOCK = threading.Lock()


def read_amsua(_product, url: str, session, day: datetime, directory) -> dict:
    """
    Downloads a granule, and reads its fields: the time at which each scan looks at
    nadir and each footprint's time (seconds from the start of the day), the
    footprints' positions, the spacecraft's position, and whether the spacecraft's
    geolocation or the AMSU-A1 module's data of a scan is flagged. The footprints'
    geolocation flags are not used, as their lowest bit is set for every footprint.
    """
    # pylint: disable=import-outside-toplevel
    import pyhdf.VS  # noqa: F401, needed by HDF.vstart
    from pyhdf.HDF import HDF
    from pyhdf.SD import SD, SDC

    path = pathlib.Path(directory) / "airabrad" / url.rsplit("/", 1)[1]
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        content = get_with_retries(url, session).content
        path.with_suffix(".part").write_bytes(content)
        path.with_suffix(".part").replace(path)
    with HDF4_LOCK:
        sd = SD(str(path), SDC.READ)
        fields = {
            name: sd.select(name)[:] for name in ("Latitude", "Longitude", "Time")
        }
        sd.end()
        hdf = HDF(str(path))
        vs = hdf.vstart()
        for name in (
            "nadirTAI",
            "sat_lat",
            "sat_lon",
            "satheight",
            "satgeoqa",
            "state1",
        ):
            vdata = vs.attach(name)
            fields[name] = np.array(vdata[:]).ravel()
            vdata.detach()
        vs.end()
        hdf.close()
    return {
        "time": tai93_to_seconds(fields["nadirTAI"], day),
        "beam_time": tai93_to_seconds(fields["Time"], day),
        "lat": fields["Latitude"],
        "lon": fields["Longitude"],
        "sc_lat": fields["sat_lat"],
        "sc_lon": fields["sat_lon"],
        "sc_alt": fields["satheight"],
        "quality": (fields["satgeoqa"] != 0) | (fields["state1"] != 0),
    }


PRODUCT = Product(short_name="AIRABRAD", version="005", read=read_amsua)

# the product of Metop's AMSU-A, in the EUMETSAT Data Store, and the satellite's code
# in the names of its products
EPS_COLLECTION = "EO:EUM:DAT:METOP:AMSUL1"
EPS_CODES = {"metop-b": "M01", "metop-c": "M03"}
EPS_EPOCH = pd.Timestamp("2000-01-01", tz="UTC")
# the measurement record (class 8), its versions, and its size (bytes), and the
# offsets (bytes) of its fields, alike in both versions, from the ATOVS Level 1b
# Product Guide (EUM/OPS-EPS/MAN/04/0030, v3B, sections 11.3 and 11.4)
EPS_MDR = (8, (3, 4), 3464)
EPS_DEGRADED = (20, 21)  # DEGRADED_INST_MDR and DEGRADED_PROC_MDR, booleans
EPS_ALTITUDE = 1838  # SPACECRAFT_ALTITUDE, unsigned 32-bit, 0.1 km
EPS_LOCATION = 2082  # EARTH_LOCATION, 30 (latitude, longitude), signed 32-bit, 1e-4 deg
N_FOOTPRINTS = 30


def read_eps(path: pathlib.Path, day: datetime) -> dict:
    """
    Reads the measurement records of a product of Metop's AMSU-A in EPS native
    format: each scan's time (the start of its record, in seconds from the start of
    the day), its footprints' positions, the spacecraft's altitude, its ground
    position (the midpoint of the two central footprints), and whether the record is
    degraded.
    """
    content = path.read_bytes()
    records, offset = [], 0
    while offset < len(content):
        size = int.from_bytes(content[offset + 4 : offset + 8], "big")
        if content[offset] == EPS_MDR[0]:
            if content[offset + 3] not in EPS_MDR[1] or size != EPS_MDR[2]:
                raise ValueError(
                    f"{path.name}: a measurement record of version "
                    f"{content[offset + 3]} and {size} bytes, not of versions "
                    f"{EPS_MDR[1]} of {EPS_MDR[2]} bytes"
                )
            records.append(content[offset : offset + size])
        offset += size
    body = np.frombuffer(b"".join(records), dtype=np.uint8).reshape(len(records), -1)
    days = body[:, 8:10].copy().view(">u2")[:, 0].astype(float)
    milliseconds = body[:, 10:14].copy().view(">u4")[:, 0].astype(float)
    location = (
        body[:, EPS_LOCATION : EPS_LOCATION + 8 * N_FOOTPRINTS]
        .copy()
        .view(">i4")
        .reshape(len(records), N_FOOTPRINTS, 2)
        * 1e-4
    )
    lat, lon = location[..., 0], location[..., 1]
    middle = ecef(lon[:, 14:16], lat[:, 14:16]).mean(axis=1)
    sc_lon, sc_lat, _ = TO_GEODETIC.transform(*middle.T)
    start = (pd.Timestamp(day) - EPS_EPOCH).total_seconds()
    return {
        "time": days * 86400 + milliseconds / 1e3 - start,
        "lat": lat,
        "lon": lon,
        "sc_lat": np.asarray(sc_lat),
        "sc_lon": np.asarray(sc_lon),
        "sc_alt": body[:, EPS_ALTITUDE : EPS_ALTITUDE + 4].copy().view(">u4")[:, 0]
        * 0.1,
        "quality": body[:, EPS_DEGRADED[0]] | body[:, EPS_DEGRADED[1]],
    }


def load_day_eps(day: datetime, directory: pathlib.Path, platform: str) -> dict:
    """
    Downloads and reduces the scans of a day of Metop's AMSU-A, or loads them from a
    cache, keeping each scan once, in time order, within the day.
    """
    cache = directory / f"amsua_{platform}_{day:%Y%m%d}.npz"
    if not cache.exists():
        files = download_eumetsat(
            EPS_COLLECTION,
            day - timedelta(hours=2),
            day + timedelta(days=1),
            directory / f"amsua_{platform}_{day:%Y%m%d}",
            match=f"_{EPS_CODES[platform]}_",
            entry=".nat",
        )
        parts = [read_eps(file, day) for file in files]
        data = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
        _, first = np.unique(np.round(data["time"], 3), return_index=True)
        keep = first[(data["time"][first] >= 0) & (data["time"][first] < 86400)]
        np.savez_compressed(
            cache, granules=len(files), **{k: v[keep] for k, v in data.items()}
        )
        print(f"{len(files)} files, {len(keep)} scans", flush=True)
    return dict(np.load(cache))
# published footprints: instrument, altitude (km), inclination (deg), scan angle
# (deg), and sizes along and across track (km), from the product's documentation,
# which gives 45 km in diameter at nadir, at Aqua's nominal altitude
PUBLISHED_FOOTPRINTS = (("AMSU-A", 705, 98.2, 0, 45.0, 45.0),)
SETUP = Setup(
    swath="swath",
    beams="footprints",
    size_instrument="AMSU-A",
    published_footprints=PUBLISHED_FOOTPRINTS,
    fixture_path=FIXTURES / "amsua.npz",
    reference="AIRABRAD version 005 (footprint positions and times)",
    published_reference="Footprint published in the AIRABRAD documentation; "
    "footprint positions of AIRABRAD version 005",
    time_steps=(SCAN_PERIOD, 32, 128),
    pass_radius=80,
    point_spacing=500e3,
    pixel_stride=20,
    spacing_stride=50,
)


def main() -> None:
    """Validates the entry on a day of data, or fits its angles."""
    run = parse_run("amsua", __doc__)
    data = (
        load_day_eps(run.day, run.directory, run.platform)
        if run.platform in EPS_CODES
        else load_day(PRODUCT, run.day, run.directory)
    )
    reference = Reference(run.day, data, SCAN_PERIOD)
    print(
        f"{run.day:%Y-%m-%d}: {len(reference.time)} scans from {reference.time[0]}, "
        f"{len(reference.gaps)} gaps over a minute",
        flush=True,
    )
    if run.fit:
        for frame in ("inertial", "earth_fixed"):
            fit(reference, frame)
        offsets = reference.beam_seconds - reference.seconds[:, np.newaxis]
        print(f"footprint times from nadir: {np.median(offsets, axis=0)}")
        return
    tle = run.element_set()
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    if SETUP.beams in run.configurations:
        validate_beams(SETUP, run, reference, orbit, tle)
    if SETUP.swath in run.configurations:
        validate_swath(SETUP, run, reference, orbit, tle)


if __name__ == "__main__":
    main()

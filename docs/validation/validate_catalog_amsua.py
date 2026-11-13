"""
Validates the AMSU-A catalog entry against the AMSU-A Level 1B product of Aqua (part
of the AIRS project's products), using the exact values of the entry through public
TAT-C methods. The entry's scan angles are fitted to the positions of the footprints
on one day (run with `--fit`), which needs no element set: running on a day at least
30 days from it tests the entry on data held out from the fit.

Needs the `validation` extra, `pyhdf`, and an Earthdata token in `EARTHDATA_TOKEN`.
The product's files are HDF4 (about 0.5 MB each), so they are downloaded whole:

    python docs/validation/validate_catalog_amsua.py <data directory> 2026-10-04 --fit
    python docs/validation/validate_catalog_amsua.py <data directory> 2026-11-05

It also validates the entry's sample size configuration, following
`ValidateSampleSizes`. The gated metrics compare sizes with the footprint published
in the product's documentation, and the projected pixels with the footprints'
positions.

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import pathlib
import threading
from datetime import datetime, timezone

import numpy as np

from tatc.catalog import get_entry
from tatc.schemas import GeneralPerturbationsOrbit
from gesdisc_scanner import (
    Product,
    Reference,
    Setup,
    fit,
    load_day,
    tai93_to_seconds,
    validate_beams,
    validate_swath,
)
from validation_io import get_with_retries

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
FIT_DAY = datetime(2026, 10, 4, tzinfo=timezone.utc)
# the public element set of `validate_catalog_modis.py` for Aqua, with an epoch after
# the day of data
ELEMENTS = {
    "2026-11-05": (
        "1 27424U 02022A   26311.28865478  .00000831  00000+0  17465-3 0  9993",
        "2 27424  98.4481 283.5818 0000983 104.9545   2.1913 14.62267061304185",
    ),
}
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
# published footprints: instrument, altitude (km), inclination (deg), scan angle
# (deg), and sizes along and across track (km), from the product's documentation,
# which gives 45 km in diameter at nadir, at Aqua's nominal altitude
PUBLISHED_FOOTPRINTS = (("AMSU-A", 705, 98.2, 0, 45.0, 45.0),)
SETUP = Setup(
    entry_id="amsua",
    platform="aqua",
    swath="swath",
    beams="footprints",
    size_instrument="AMSU-A",
    published_footprints=PUBLISHED_FOOTPRINTS,
    record_path=REPOSITORY / "src" / "tatc" / "resources" / "validation" / "amsua.json",
    fixture_path=REPOSITORY / "tests" / "catalog" / "fixtures" / "amsua.npz",
    source="docs/validation/validate_catalog_amsua.py",
    reference="AIRABRAD version 005 (footprint positions and times)",
    published_reference="Footprint published in the AIRABRAD documentation; "
    "footprint positions of AIRABRAD version 005",
    fit_day=FIT_DAY,
    time_steps=(SCAN_PERIOD, 32, 128),
    pass_radius=80,
    point_spacing=500e3,
    pixel_stride=20,
    spacing_stride=50,
)


def main() -> None:
    """Validates the entry on a day of data, or fits its angles."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("day")
    parser.add_argument("--fit", action="store_true")
    parser.add_argument("--sample-sizes-only", action="store_true")
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)

    day = datetime.fromisoformat(arguments.day).replace(tzinfo=timezone.utc)
    reference = Reference(day, load_day(PRODUCT, day, arguments.directory), SCAN_PERIOD)
    print(
        f"{arguments.day}: {len(reference.time)} scans from {reference.time[0]}, "
        f"{len(reference.gaps)} gaps over a minute",
        flush=True,
    )
    if arguments.fit:
        for frame in ("inertial", "earth_fixed"):
            fit(reference, frame)
        offsets = reference.beam_seconds - reference.seconds[:, np.newaxis]
        print(f"footprint times from nadir: {np.median(offsets, axis=0)}")
        return
    tle = list(ELEMENTS[arguments.day])
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    entry = get_entry("amsua")
    validate_beams(SETUP, day, tle, entry, orbit, reference, arguments.write)
    if not arguments.sample_sizes_only:
        validate_swath(SETUP, day, tle, entry, orbit, reference, arguments.write)
    if arguments.write:
        print(f"wrote {SETUP.record_path} and {SETUP.fixture_path}")


if __name__ == "__main__":
    main()

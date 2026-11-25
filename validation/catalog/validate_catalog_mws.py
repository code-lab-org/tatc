"""
Validates the MWS catalog entry against Metop-SG A1's MWS Level 1b product (in BUFR,
EO:EUM:DAT:0892, of the EUMETSAT Data Store), which gives, for each scan, its time
and the positions of its 95 fields of view, using the exact values of the entry
through public TAT-C methods (see `scanner.py` for the scoring).

The product does not record the satellite's position, so positions come from the
element set nearest each day, at the times of the scans. The entry's scan angles are
fitted on one day (run with `--fit`); running on a day at least 30 days from it tests
the fitted angles on data held out from the fit.

Needs the `validation` extra and a login to the EUMETSAT Data Store
(`eumdac auth --interactive`):

    python validation/catalog/validate_catalog_mws.py <data directory> 2026-10-19 --fit
    python validation/catalog/validate_catalog_mws.py <data directory> 2026-11-18

It also validates the entry's sample size configuration (`beams`), following
`ValidateSampleSizes`: sizes are compared with the footprints published in the
EPS-SG MWS Science Plan, which the gate tests, the projected pixels with the
positions of the fields of view, and their spacing with that measured.
`--configuration beams` (or `swath`) tests one configuration.

Without a day, it takes that of the entry's record (or, with `--fit`, of the fit). The
element sets are those of the files in `elements/` (or of a file given with
`--elements`). With `--write`, it writes the entry's validation record and the reduced
fixture from which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import pathlib
from datetime import datetime, timedelta

import numpy as np
from skyfield.framelib import itrs

from tatc.constants import timescale
from tatc.schemas import GeneralPerturbationsOrbit
from scanner import Reference, Setup, fit, validate_beams, validate_swath
from validation_io import FIXTURES, TO_GEODETIC, download_eumetsat, parse_run

COLLECTION = "EO:EUM:DAT:0892"
FIELDS_OF_VIEW = 95
SCAN_PERIOD = 2.254
# published footprints: instrument, altitude (km), inclination (deg), scan angle
# (deg), and sizes along and across track (km), from Figure 8 of the EPS-SG MWS
# Science Plan, which states no altitude or angle: those of the beam widths of the
# AAPP-SG User Manual (827 km) and the outermost field of view (49.312 deg)
PUBLISHED_FOOTPRINTS = (
    ("MWS 23.8-31.4 GHz", 827, 98.7, 0, 39.6, 39.6),
    ("MWS 23.8-31.4 GHz", 827, 98.7, 49.312, 67.6, 132.6),
    ("MWS 50.3-57.3 GHz", 827, 98.7, 0, 20.0, 20.0),
    ("MWS 50.3-57.3 GHz", 827, 98.7, 49.312, 34.2, 66.9),
    ("MWS 89-229 GHz", 827, 98.7, 0, 17.0, 17.0),
    ("MWS 89-229 GHz", 827, 98.7, 49.312, 29.0, 56.7),
)
SETUP = Setup(
    swath="swath",
    beams="beams",
    size_instrument="MWS 50.3-57.3 GHz",
    published_footprints=PUBLISHED_FOOTPRINTS,
    fixture_path=FIXTURES / "mws.npz",
    reference="MWS Level 1b (EO:EUM:DAT:0892, BUFR), positions of the fields of view",
    published_reference="Footprints of the EPS-SG MWS Science Plan; positions of the "
    "fields of view of the MWS Level 1b (EO:EUM:DAT:0892, BUFR)",
    time_steps=(SCAN_PERIOD, 32, 128),
    pass_radius=60,
    point_spacing=500e3,
)


def read_scans(path: pathlib.Path, day: datetime) -> list[tuple]:
    """
    Reads the scans of a file (one in each message): each scan's time (seconds from
    the start of a day) and its fields of view's latitudes and longitudes, in order.
    """
    # imported here, as only the scripts of products in BUFR need it
    import eccodes  # pylint: disable=import-outside-toplevel

    scans = []
    with open(path, "rb") as f:
        while (message := eccodes.codes_bufr_new_from_file(f)) is not None:
            try:
                eccodes.codes_set(message, "unpack", 1)
                fov = eccodes.codes_get_array(message, "fieldOfViewNumber")
                if len(fov) != FIELDS_OF_VIEW:
                    continue
                order = np.argsort(fov)
                date = [
                    int(eccodes.codes_get_array(message, key)[0])
                    for key in ("year", "month", "day", "hour", "minute")
                ]
                second = float(eccodes.codes_get_array(message, "second")[0])
                time = datetime(*date, tzinfo=day.tzinfo) + timedelta(seconds=second)
                scans.append(
                    (
                        (time - day).total_seconds(),
                        eccodes.codes_get_array(message, "latitude")[order],
                        eccodes.codes_get_array(message, "longitude")[order],
                    )
                )
            finally:
                eccodes.codes_release(message)
    return scans


def load_day(run, day: datetime, orbit) -> dict:
    """
    Downloads and reads the scans of a day, with the satellite's geodetic position at
    their times from the orbit, or loads them from a cache.
    """
    path = run.directory / f"mws_{day:%Y%m%d}.npz"
    if path.exists():
        return dict(np.load(path))
    files = download_eumetsat(
        COLLECTION,
        day - timedelta(minutes=5),
        day + timedelta(days=1, minutes=5),
        run.directory / f"mws_{day:%Y%m%d}",
        match="SGA1-MWS-1B-RAD",
        entry=".bin",
    )
    scans = {}
    for file in files:
        for seconds, latitude, longitude in read_scans(file, day):
            if 0 <= seconds < 86400:
                scans[round(seconds, 3)] = (seconds, latitude, longitude)
    keys = sorted(scans)
    seconds = np.array([scans[k][0] for k in keys])
    times = timescale.from_datetimes(
        [day + timedelta(seconds=float(s)) for s in seconds]
    )
    position = orbit.to_gp_orbit().get_orbit_track_at_time(times).frame_xyz(itrs).m
    lon, lat, alt = TO_GEODETIC.transform(*position)
    data = {
        "time": seconds,
        "lat": np.array([scans[k][1] for k in keys]),
        "lon": np.array([scans[k][2] for k in keys]),
        "sc_lat": np.asarray(lat),
        "sc_lon": np.asarray(lon),
        "sc_alt": np.asarray(alt),
    }
    np.savez_compressed(path, granules=len(files), **data)
    print(f"{len(files)} files, {len(keys)} scans", flush=True)
    return data


def main() -> None:
    """Validates the entry on a day of data, or fits its angles."""
    run = parse_run("mws", __doc__)
    tle = run.element_set()
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    reference = Reference(run.day, load_day(run, run.day, orbit), SCAN_PERIOD)
    print(
        f"{run.day:%Y-%m-%d}: {len(reference.time)} scans, "
        f"{len(reference.gaps)} gaps over a minute",
        flush=True,
    )
    if run.fit:
        for frame in ("earth_fixed", "inertial"):
            fit(reference, frame)
        return
    if SETUP.beams in run.configurations:
        validate_beams(SETUP, run, reference, orbit, tle)
    if SETUP.swath in run.configurations:
        validate_swath(SETUP, run, reference, orbit, tle)


if __name__ == "__main__":
    main()

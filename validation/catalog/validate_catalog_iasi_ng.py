"""
Validates the IASI-NG catalog entry against Metop-SG A1's IASI-NG Level 1c product (in
BUFR, EO:EUM:DAT:0892, of the EUMETSAT Data Store), using the exact values of the
entry through public TAT-C methods (see `scanner.py` for the scoring). Each scan holds
14 fields of regard of 16 pixels each, every field of regard at its own time, about
0.82 s apart, in scans every 15.56 s; each field of regard is compared at the center
of its pixels.

The product does not record the satellite's position, so positions come from the
element set nearest each day, at the times of the scans. The entry's scan angles are
fitted on one day (run with `--fit`); running on a day at least 30 days from it tests
the fitted angles on data held out from the fit. As a day of the product is about 40
GB, a run reads the first `HOURS` hours of a day.

Needs the `validation` extra and a login to the EUMETSAT Data Store
(`eumdac auth --interactive`):

    python validation/catalog/validate_catalog_iasi_ng.py <data directory> 2026-10-19 --fit
    python validation/catalog/validate_catalog_iasi_ng.py <data directory> 2026-11-18

It also validates the entry's sample size configuration (`fields_of_regard`),
following `ValidateSampleSizes`: the projected pixels are compared with the centers
of the fields of regard, and their spacing with that measured. The only published
pixel, at nadir in the IASI-NG Science Plan, is compared without the gate.

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
from validation_io import (
    FIXTURES,
    TO_GEODETIC,
    download_eumetsat,
    ecef,
    parse_run,
)

COLLECTION = "EO:EUM:DAT:0892"
FIELDS_OF_REGARD, PIXELS = 14, 16
SCAN_PERIOD = 15.561
# the hours read from the start of a day, about two orbits
HOURS = 3.5
# the published pixel at nadir: instrument, altitude (km), inclination (deg), scan
# angle (deg), and sizes along and across track (km), from Table 2.1 of the IASI-NG
# Science Plan (which states no altitude), compared without the gate, as the pixel
# width comes from it
PUBLISHED_FOOTPRINTS = (("IASI-NG", 835, 98.7, 0, 12, 12),)
SETUP = Setup(
    swath="swath",
    beams="fields_of_regard",
    size_instrument="IASI-NG",
    published_footprints=PUBLISHED_FOOTPRINTS,
    fixture_path=FIXTURES / "iasi_ng.npz",
    reference="IASI-NG Level 1c (EO:EUM:DAT:0892, BUFR), centers of the fields of "
    f"regard, first {HOURS} hours of the day",
    published_reference="Nadir pixel of the IASI-NG Science Plan (not gated); "
    "centers of the fields of regard of the IASI-NG Level 1c (EO:EUM:DAT:0892, "
    f"BUFR), first {HOURS} hours of the day",
    time_steps=(SCAN_PERIOD, 32, 128),
    pass_radius=150,
    point_spacing=500e3,
    # a few hours of scans of 14 fields of regard
    pixel_stride=10,
    spacing_stride=20,
    gate_published=False,
)


def read_scans(path: pathlib.Path, day: datetime) -> list[tuple]:
    """
    Reads the scans of a file (one in each message): the times of each scan's fields
    of regard (seconds from the start of a day) and the latitudes and longitudes of
    the centers of their pixels.
    """
    # imported here, as only the scripts of products in BUFR need it
    import eccodes  # pylint: disable=import-outside-toplevel

    scans = []
    with open(path, "rb") as f:
        while (message := eccodes.codes_bufr_new_from_file(f)) is not None:
            try:
                eccodes.codes_set(message, "unpack", 1)
                fov = eccodes.codes_get_array(message, "fieldOfViewNumber")
                if len(fov) != FIELDS_OF_REGARD * PIXELS:
                    continue
                order = np.argsort(fov)
                date = [
                    int(eccodes.codes_get_array(message, key)[0])
                    for key in ("year", "month", "day", "hour", "minute")
                ]
                base = (datetime(*date, tzinfo=day.tzinfo) - day).total_seconds()
                seconds = eccodes.codes_get_array(message, "second")[order]
                # the message gives one minute, so seconds that pass a minute
                # start again from 0
                seconds = seconds + 60 * np.cumsum(np.r_[0, np.diff(seconds) < 0])
                latitude = eccodes.codes_get_array(message, "latitude")[order]
                longitude = eccodes.codes_get_array(message, "longitude")[order]
                # the center of each field of regard's pixels, on the ellipsoid
                center = (
                    ecef(longitude, latitude)
                    .reshape(FIELDS_OF_REGARD, PIXELS, 3)
                    .mean(axis=1)
                )
                lon, lat, _ = TO_GEODETIC.transform(*center.T)
                scans.append(
                    (
                        base + seconds.reshape(FIELDS_OF_REGARD, PIXELS)[:, 0],
                        np.asarray(lat),
                        np.asarray(lon),
                    )
                )
            finally:
                eccodes.codes_release(message)
    return scans


def load_day(run, day: datetime, orbit) -> dict:
    """
    Downloads and reads the scans of the first hours of a day, with the satellite's
    geodetic position at the times of their middles from the orbit, or loads them
    from a cache.
    """
    path = run.directory / f"iasi_ng_{day:%Y%m%d}.npz"
    if path.exists():
        return dict(np.load(path))
    files = download_eumetsat(
        COLLECTION,
        day,
        day + timedelta(hours=HOURS),
        run.directory / f"iasing_{day:%Y%m%d}",
        match="SGA1-IAS-1C-RAD",
        entry=".bin",
    )
    scans = {}
    for file in files:
        for seconds, latitude, longitude in read_scans(file, day):
            if 0 <= seconds[0] < HOURS * 3600:
                scans[round(seconds[0], 3)] = (seconds, latitude, longitude)
    keys = sorted(scans)
    beam_time = np.array([scans[k][0] for k in keys])
    middle = beam_time.mean(axis=1)
    times = timescale.from_datetimes(
        [day + timedelta(seconds=float(s)) for s in middle]
    )
    position = orbit.to_gp_orbit().get_orbit_track_at_time(times).frame_xyz(itrs).m
    lon, lat, alt = TO_GEODETIC.transform(*position)
    data = {
        "time": middle,
        "beam_time": beam_time,
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
    run = parse_run("iasi_ng", __doc__)
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

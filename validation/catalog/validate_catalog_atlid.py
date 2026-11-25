"""
Validates the ATLID catalog entry against EarthCARE's Level 1B ATLID product
(ATL_NOM_1B), using the exact values of the entry through public TAT-C methods. ATLID
is a single-beam lidar, pointed a few degrees behind nadir. Its pointing is fitted to
the positions of the profiles on one day (run with `--fit`), from the sensor positions
the product records, so the fit needs no element set; running on a day at least 30
days from it tests the fitted pointing on data held out from the fit. A sample of a
day's profiles is scored as for the profiler class: whether the view observes each
within two minutes of its time, the epoch's difference, and its distance across track
from the predicted track of the view's center.

Needs the `validation` extra and a MAAP offline token in `ESA_MAAP_KEY`. Only the
profiles' times and positions are read, with range requests:

    python validation/catalog/validate_catalog_atlid.py <data directory> 2026-09-27 --fit
    python validation/catalog/validate_catalog_atlid.py <data directory> 2026-11-07

It also validates the entry's sample size configuration (`footprints`), following
`validation_io.validate_footprints`: the spacing of profiles along track is compared
with that measured, which the gate tests, and the laser's published footprint,
without the gate.

Without a day, it takes that of the entry's record (or, with `--fit`, of the fit).
The element set is the one nearest the middle of the day in
`elements/earthcare_2026-09-27.3le` (or in a file given with `--elements`). With
`--write`, it writes the entry's validation record and the reduced fixture from which
unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from datetime import datetime

import numpy as np

from tatc.schemas import GeneralPerturbationsOrbit
from earthcare_io import (
    Profiles,
    check_time_scale,
    fit_pointing,
    load_day,
    to_seconds,
    validate_profiler,
    validate_profiler_footprints,
)
from validation_io import FIXTURES, ecef, Footprints, parse_run

CONFIGURATION = "profile"
# published footprints at nadir: instrument, altitude (km), inclination (deg), scan
# angle (deg), and sizes along and across track (km), of the laser on ESA's ATLID
# page, which states no definition, compared without the gate
PUBLISHED_FOOTPRINTS = (("ATLID", 393, 97.05, 0, 0.015, 0.015),)
FOOTPRINTS = Footprints(
    key="footprints",
    published_footprints=PUBLISHED_FOOTPRINTS,
    reference="Laser footprint of ESA's ATLID page (not gated); profile spacing of ATL_NOM_1B",
    fixture_path=FIXTURES / "atlid.npz",
    gate_published=False,
)


def read_atlid(f, day: datetime) -> dict:
    """Reads the times and positions of a frame's profiles."""
    g = f["ScienceData"]
    return {
        "time": to_seconds(g["time"][:], day),
        "lat": g["ellipsoid_latitude"][:],
        "lon": g["ellipsoid_longitude"][:],
        "sc_position": ecef(
            g["sensor_longitude"][:],
            g["sensor_latitude"][:],
            g["sensor_altitude"][:].astype(float),
        ),
    }


def main() -> None:
    """Validates the entry on a day of data, or fits its pointing."""
    run = parse_run("atlid", __doc__)
    data = load_day("ATL_NOM_1B", run.day, run.directory, read_atlid)
    profiles = Profiles(run.day, data)
    print(
        f"{run.day:%Y-%m-%d}: {len(profiles.seconds)} profiles from "
        f"{int(data['frames'])} frames, {len(profiles.gaps)} gaps over a minute",
        flush=True,
    )
    if run.fit:
        check_time_scale("ATL_NOM_1B", run.day, "ScienceData/time")
        fit_pointing(profiles)
        return
    tle = run.element_set()
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    if FOOTPRINTS.key in run.configurations:
        validate_profiler_footprints(run, FOOTPRINTS, profiles, orbit, tle)
    if CONFIGURATION in run.configurations:
        validate_profiler(
            run,
            CONFIGURATION,
            profiles,
            orbit,
            tle,
            "ATL_NOM_1B (profile positions and times)",
        )


if __name__ == "__main__":
    np.seterr(invalid="ignore")
    main()

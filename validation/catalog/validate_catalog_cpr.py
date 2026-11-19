"""
Validates the CPR catalog entry against EarthCARE's Level 1B CPR product (CPR_NOM_1B,
of JAXA), using the exact values of the entry through public TAT-C methods. CPR is a
single-beam Doppler cloud radar, pointed near nadir. Its pointing is fitted to the
positions of the profiles on one day (run with `--fit`), from the spacecraft positions
the product records, so the fit needs no element set; running on a day at least 30
days from it tests the fitted pointing on data held out from the fit. A sample of a
day's profiles is scored as for the profiler class: whether the view observes each
within two minutes of its time, the epoch's difference, and its distance across track
from the predicted track of the view's center.

Needs the `validation` extra and a MAAP offline token in `ESA_MAAP_KEY`. Only the
profiles' times and positions are read, with range requests:

    python validation/catalog/validate_catalog_cpr.py <data directory> 2026-09-27 --fit
    python validation/catalog/validate_catalog_cpr.py <data directory> 2026-11-07

Without a day, it takes that of the entry's record (a fit needs its day). The
element set is the one nearest the middle of the day in
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
)
from validation_io import parse_run

CONFIGURATION = "nadir"


def read_cpr(f, day: datetime) -> dict:
    """Reads the times and positions of a frame's profiles."""
    g = f["ScienceData/Geo"]
    return {
        "time": to_seconds(g["profileTime"][:], day),
        "lat": g["latitude"][:],
        "lon": g["longitude"][:],
        "sc_position": np.stack(
            [g["xPosition"][:], g["yPosition"][:], g["zPosition"][:]], axis=-1
        ),
    }


def main() -> None:
    """Validates the entry on a day of data, or fits its pointing."""
    run = parse_run("cpr", __doc__)
    data = load_day("CPR_NOM_1B", run.day, run.directory, read_cpr)
    profiles = Profiles(run.day, data)
    print(
        f"{run.day:%Y-%m-%d}: {len(profiles.seconds)} profiles from "
        f"{int(data['frames'])} frames, {len(profiles.gaps)} gaps over a minute",
        flush=True,
    )
    if run.fit:
        check_time_scale("CPR_NOM_1B", run.day, "ScienceData/Geo/profileTime")
        fit_pointing(profiles)
        return
    tle = run.element_set()
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    validate_profiler(
        run,
        CONFIGURATION,
        profiles,
        orbit,
        tle,
        "CPR_NOM_1B (profile positions and times)",
    )


if __name__ == "__main__":
    np.seterr(invalid="ignore")
    main()

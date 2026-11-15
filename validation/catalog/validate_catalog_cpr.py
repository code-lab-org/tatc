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

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import pathlib
from datetime import datetime, timezone

import numpy as np

from tatc.catalog import get_entry
from tatc.schemas import GeneralPerturbationsOrbit
from earthcare_io import (
    Profiles,
    check_time_scale,
    fit_pointing,
    load_day,
    to_seconds,
    validate_profiler,
)
from validation_io import (
    read_element_set_file,
    RECORDS,
    FIXTURES,
    source_of,
    ELEMENT_SET_FOLDER,
)

RECORD_PATH = RECORDS / "cpr.json"
FIXTURE_PATH = FIXTURES / "cpr.npz"
NORAD_ID = 59908
FIT_DAY = datetime(2026, 9, 27, tzinfo=timezone.utc)
ELEMENT_SET_FILE = ELEMENT_SET_FOLDER / "earthcare_2026-09-27.3le"
PLATFORM = "earthcare"
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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("day")
    parser.add_argument("--fit", action="store_true")
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)
    day = datetime.fromisoformat(arguments.day).replace(tzinfo=timezone.utc)
    data = load_day("CPR_NOM_1B", day, arguments.directory, read_cpr)
    profiles = Profiles(day, data)
    print(
        f"{arguments.day}: {len(profiles.seconds)} profiles from {int(data['frames'])} "
        f"frames, {len(profiles.gaps)} gaps over a minute",
        flush=True,
    )
    if arguments.fit:
        check_time_scale("CPR_NOM_1B", day, "ScienceData/Geo/profileTime")
        fit_pointing(profiles)
        return
    tle = read_element_set_file(ELEMENT_SET_FILE, NORAD_ID, day)
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    validate_profiler(
        get_entry("cpr"),
        CONFIGURATION,
        PLATFORM,
        profiles,
        orbit,
        tle,
        day,
        FIT_DAY,
        RECORD_PATH,
        FIXTURE_PATH,
        source_of(__file__),
        "CPR_NOM_1B (profile positions and times)",
        arguments.write,
    )


if __name__ == "__main__":
    np.seterr(invalid="ignore")
    main()

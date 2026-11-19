"""
Validates the TMI catalog entry against the GPM-format Level 1C TMI product of TRMM,
using the exact values of the entry through public TAT-C methods. The entry's scan is
fitted to the positions of the samples of each group (S1, 10 GHz; S2, 19 to 37 GHz;
S3, 85.5 GHz) on one day (run with `--fit`), from the spacecraft positions the product
records, so the fit needs no element set; running on a day at least 30 days from it
tests the fitted scan on data held out from the fit. Both days are in TRMM's forward
orientation (the product's `SCorientation` of 0) throughout. Element sets are read
from the file of the history of TRMM's element sets in `elements/` (or one given with
`--elements`). The tests follow those of SSMIS:

* A point is observed by a feed on a pass if it lies in a cell of the lattice of its
  group's samples, at the time it crosses the cone, from the recorded spacecraft positions.
* The swath width is compared with the distance between the first and last samples.
* The swath area over a span of the day (from noon) is compared with the union of
  the cells.

It also validates the entry's sample size configuration, as for AMSR2: sizes at the
center of the scan are compared with the footprints published in the TC PRIMED
documentation, which the gate tests, and the spacing of samples along the arc and
between scans with that measured. To run only these tests, use
`--configuration sample_size`.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`. Only the
needed fields of the product's files are read, with range requests:

    python validation/catalog/validate_catalog_tmi.py <data directory> 2013-06-01 --fit
    python validation/catalog/validate_catalog_tmi.py <data directory> 2013-06-01
    python validation/catalog/validate_catalog_tmi.py <data directory> 2013-07-10

Without a day, it takes that of the entry's record (or, with `--fit`, of the fit).

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import pathlib
from datetime import datetime

import h5py
import numpy as np
import pandas as pd

from tatc.schemas import GeneralPerturbationsOrbit
from conical import (
    compute_center_errors,
    fit_cone_scan,
    measure_arc_spacing,
    validate_cone_swath,
)
from scanner import Product, Reference, load_day
from validation_io import (
    RemoteFile,
    compute_spacing_ratios,
    parse_run,
    score_published,
    write_results,
)

CONFIGURATION = "feeds"
SAMPLE_CONFIGURATION = "sample_size"
GROUPS = ("S1", "S2", "S3")


def read_tmi(_product, url: str, session, day: datetime, _directory) -> dict:
    """
    Reads the times of a granule's scans, the spacecraft's positions, and the
    positions of the samples of each group, with range requests.
    """
    with h5py.File(RemoteFile(url, session)) as f:
        g = f["S1"]
        date = pd.to_datetime(
            pd.DataFrame(
                {
                    "year": g["ScanTime/Year"][:],
                    "month": g["ScanTime/Month"][:],
                    "day": g["ScanTime/DayOfMonth"][:],
                }
            ),
            utc=True,
        )
        part = {
            "time": (date - pd.Timestamp(day)).dt.total_seconds().to_numpy()
            + g["ScanTime/SecondOfDay"][:].astype(float),
            "sc_lat": g["SCstatus/SClatitude"][:],
            "sc_lon": g["SCstatus/SClongitude"][:],
            "sc_alt": g["SCstatus/SCaltitude"][:],
            "orientation": g["SCstatus/SCorientation"][:],
        }
        for group in GROUPS:
            part[f"lat_{group}"] = f[group]["Latitude"][:]
            part[f"lon_{group}"] = f[group]["Longitude"][:]
    return part


PRODUCT = Product(short_name="GPM_1CTRMMTMI", version="08", read=read_tmi)
# the scan period is measured on the day of the fit
SCAN_PERIOD = 1.899
# published footprints: instrument, altitude (km), and sizes along and across track
# (km), from the TC PRIMED documentation, after TRMM's orbit boost of August 2001, at
# the center of the scan, on an orbit of 35 deg inclination
PUBLISHED_FOOTPRINTS = (
    ("TMI 10.65 GHz", 403, 68.1, 41.3),
    ("TMI 19.35 GHz", 403, 34.6, 20.9),
    ("TMI 21.3 GHz", 403, 31.1, 18.9),
    ("TMI 37.0 GHz", 403, 18.3, 11.1),
    ("TMI 85.5 GHz", 403, 7.8, 4.8),
)
PUBLISHED_INCLINATION = 35
# an instrument of each group, whose spacing is measured
SPACING_INSTRUMENTS = {
    "S1": "TMI 10.65 GHz",
    "S2": "TMI 19.35 GHz",
    "S3": "TMI 85.5 GHz",
}
SPACING_SCAN_STRIDE = 200
FIXTURE_SPACING_SAMPLES = 100


def load_references(day: datetime, directory: pathlib.Path) -> dict[str, Reference]:
    """Loads a day, as a reference for each group of samples."""
    data = load_day(PRODUCT, day, directory)
    orientation = data["orientation"]
    print(
        "orientation of the scans: "
        + ", ".join(
            f"{value}: {count}"
            for value, count in zip(*np.unique(orientation, return_counts=True))
        ),
        flush=True,
    )
    references = {}
    for group in GROUPS:
        part = {
            "time": data["time"],
            "sc_lat": data["sc_lat"],
            "sc_lon": data["sc_lon"],
            "sc_alt": data["sc_alt"],
            "lat": data[f"lat_{group}"],
            "lon": data[f"lon_{group}"],
        }
        references[group] = Reference(day, part, SCAN_PERIOD)
    return references


def measure_spacing(reference: Reference) -> tuple[np.ndarray, np.ndarray]:
    """
    Measures the spacing of a group's samples (see `conical.measure_arc_spacing`).
    Returns rows of latitude, longitude, spacing between scans (km), and spacing
    along the arc (km), with the seconds of each row.
    """
    rows = measure_arc_spacing(
        reference.lon,
        reference.lat,
        reference.seconds,
        reference.contiguous_scans()[::SPACING_SCAN_STRIDE],
    )
    seconds = rows[:, 2]
    return np.stack([rows[:, 0], rows[:, 1], rows[:, 4], rows[:, 3]], axis=-1), seconds


def validate_sample_sizes(run, orbit, tle, references) -> None:
    """Scores the sample size configuration, and writes its results if the run writes."""
    entry, day = run.entry, run.day
    within, largest = score_published(
        compute_center_errors(
            entry,
            SAMPLE_CONFIGURATION,
            PUBLISHED_FOOTPRINTS,
            PUBLISHED_INCLINATION,
            time_step=SCAN_PERIOD,
            platform=run.platform,
        )
    )
    print(
        f"published dimensions: {within:.1f}% within 3%, largest error {largest:.2f}%"
    )
    metrics = [
        {"name": "sample_size_within_3pct", "value": within},
        {"name": "sample_size_largest_error", "value": largest},
    ]
    satellite = entry.to_satellite(
        orbit,
        configuration=SAMPLE_CONFIGURATION,
        time_step=SCAN_PERIOD,
        platform=run.platform,
    )
    rng = np.random.default_rng(2)
    fixture_metrics = {"sample_size_within_3pct": within}
    reduced = {}
    for group, name in SPACING_INSTRUMENTS.items():
        rows, seconds = measure_spacing(references[group])
        times = pd.DatetimeIndex(pd.Timestamp(day) + pd.to_timedelta(seconds, "s"))
        # the spacing between scans is along track, and along the arc across track
        between, arc = compute_spacing_ratios(
            satellite, [name] * len(rows), times, rows
        )
        print(
            f"{group} spacing, predicted / measured: {arc:.4f} arc, {between:.4f} scans"
        )
        metrics += [
            {"name": "spacing_arc_ratio_median", "value": arc, "scope": group},
            {"name": "spacing_scan_ratio_median", "value": between, "scope": group},
        ]
        picks = rng.choice(len(rows), FIXTURE_SPACING_SAMPLES, replace=False)
        between, arc = compute_spacing_ratios(
            satellite, [name] * len(picks), times[picks], rows[picks]
        )
        fixture_metrics[f"spacing_arc_ratio_median_{group}"] = arc
        fixture_metrics[f"spacing_scan_ratio_median_{group}"] = between
        reduced[f"spacing_{group}"] = rows[picks]
        reduced[f"spacing_seconds_{group}"] = seconds[picks]
    if not run.write:
        return
    write_results(
        run,
        SAMPLE_CONFIGURATION,
        reference="Footprints of the TC PRIMED documentation; spacing of "
        "GPM_1CTRMMTMI version 08",
        metrics=metrics,
        fixture={"tle": np.array(tle), "day": np.array(day.isoformat()), **reduced},
        fixture_metrics=fixture_metrics,
    )


def main() -> None:
    """Validates the entry on a day of data, or fits its scan."""
    run = parse_run("tmi", __doc__)
    references = load_references(run.day, run.directory)
    reference = references["S1"]
    print(
        f"{run.day:%Y-%m-%d}: {len(reference.time)} scans, {len(reference.gaps)} gaps "
        f"over a minute, median altitude {np.median(reference.altitude) / 1e3:.1f} km",
        flush=True,
    )
    if run.fit:
        for group in GROUPS:
            fit_cone_scan(references[group], range(-45, 45, 15), group)
        return
    tle = run.element_set()
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    if SAMPLE_CONFIGURATION in run.configurations:
        validate_sample_sizes(run, orbit, tle, references)
    if CONFIGURATION in run.configurations:
        validate_cone_swath(
            run,
            references,
            orbit,
            tle,
            key=CONFIGURATION,
            reference="GPM_1CTRMMTMI version 08 (positions of the S1, S2, and S3 "
            "samples)",
            scan_period=SCAN_PERIOD,
            pass_radius=60,
            point_spacing=300e3,
            fixture_passes=80,
        )


if __name__ == "__main__":
    main()

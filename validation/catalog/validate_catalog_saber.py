"""
Validates the SABER catalog entry against TIMED's SABER Level 2A product (version 2.0,
of GATS, Inc.), which records, for every event (a vertical scan), its direction and the
time and the geodetic latitude, longitude, and altitude of the tangent point at each
level of an altitude grid, using the exact values of the entry through public TAT-C
methods (see `sounding.py` for the scoring).

SABER looks to one side of the track and alternates upward and downward scans. TIMED
yaws about every 60 days, which turns SABER's view between the left of the flight
direction (toward the orbit normal, viewing from about 53 deg S to 85 deg N) and the
right (85 deg S to 53 deg N): the `north_viewing` and `south_viewing` configurations,
each tested on a day of its side.

The scan's values come from `--fit`, on a day of data: the cycle period is the median
period from an upward scan to the next, the downward scan's offset the median time
from an upward scan to the next downward scan, the scan's bottom and top those of the
median event, and the duration fitted so that the representative point (at the middle
of the scan) is reached when the median event reaches it. Running on a day at least 30
days from the fit tests the values on data held out from the fit.

Reads public data, and needs the `validation` extra:

    python validation/catalog/validate_catalog_saber.py <data directory> 2026-09-21 --fit
    python validation/catalog/validate_catalog_saber.py <data directory> 2026-11-13 \\
        --configuration south_viewing

Without a day, it takes that of the entry's record (or, with `--fit`, of the fit). The
element sets are registered in `elements/element_sets.yml`. With `--write`, it writes
the entry's validation record and the reduced fixture from which unit tests recompute
the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import json
import pathlib
import re
from datetime import datetime, timedelta, timezone

import netCDF4
import numpy as np

from tatc.analysis.limb_sampling import ScanDirection
from tatc.schemas import GeneralPerturbationsOrbit
from sounding import Scan, fit_duration, validate_limb
from validation_io import get_with_retries, parse_run

URL = "https://data.gats-inc.com/saber/Version2_0/Level2A"
FIELDS = ("time", "tpaltitude", "tplatitude", "tplongitude")
# the fewest valid levels of an event that is compared
MIN_LEVELS = 50


def download_day(day: datetime, directory: pathlib.Path) -> list[pathlib.Path]:
    """
    Downloads the files of the orbits that start on a day or the day before, which
    hold the day's events, if not cached.
    """
    paths = []
    for date_ in (day - timedelta(days=1), day):
        folder = f"{URL}/{date_:%Y}/{date_.timetuple().tm_yday:03d}/"
        # the latest version of each orbit's file, as a folder may hold several
        latest = {}
        for name in sorted(
            set(
                re.findall(
                    r"SABER_L2A_\d{7}_\d+_02\.\d\d\.nc", get_with_retries(folder).text
                )
            )
        ):
            latest[name.split("_")[3]] = name
        for name in latest.values():
            path = directory / name
            if not path.exists():
                path.write_bytes(get_with_retries(folder + name).content)
            paths.append(path)
    return paths


def read_scans(paths: list[pathlib.Path], day: datetime) -> list[Scan]:
    """
    Reads the events that start on a day, with their valid levels in the order of
    their times (`time` counts milliseconds from the start of the orbit's first day).
    """
    scans = []
    for path in paths:
        with netCDF4.Dataset(path) as data:
            data.set_auto_mask(False)
            first_day = datetime.strptime(data.Data_Start_Time[:7], "%Y%j").replace(
                tzinfo=timezone.utc
            )
            values = {name: data[name][:].astype(float) for name in FIELDS}
            modes = data["mode"][:]
        for e, mode in enumerate(modes):
            ok = (
                (values["tpaltitude"][e] > 0)
                & (values["time"][e] >= 0)
                & (np.abs(values["tplatitude"][e]) <= 90)
            )
            if ok.sum() < MIN_LEVELS:
                continue
            order = np.argsort(values["time"][e][ok])
            event = {name: values[name][e][ok][order] for name in FIELDS}
            start = first_day + timedelta(milliseconds=float(event["time"][0]))
            if not day <= start < day + timedelta(days=1):
                continue
            scans.append(
                Scan(
                    start=start,
                    direction=(
                        ScanDirection.UPWARD if mode == 1 else ScanDirection.DOWNWARD
                    ),
                    elapsed=(event["time"] - event["time"][0]) / 1e3,
                    altitude=event["tpaltitude"] * 1e3,
                    latitude=event["tplatitude"],
                    longitude=event["tplongitude"],
                )
            )
    return sorted(scans, key=lambda s: s.start)


def measure_cycle(scans: list[Scan]) -> tuple[float, float]:
    """
    Measures the median period (s) from an upward scan to the next, and the median
    offset (s) of the next downward scan.
    """
    starts = np.array([s.start.timestamp() for s in scans])
    upward = np.array([s.direction == ScanDirection.UPWARD for s in scans])
    spacing = np.diff(starts)
    # consecutive upward, downward, upward scans
    triples = upward[:-2] & ~upward[1:-1] & upward[2:]
    offset = spacing[:-1][triples]
    period = spacing[:-1][triples] + spacing[1:][triples]
    return float(np.median(period)), float(np.median(offset))


def get_side(scans: list[Scan]) -> str:
    """Gets the side of the view, from the extent of the tangent latitudes."""
    latitude = np.concatenate([s.latitude for s in scans])
    return "north_viewing" if latitude.max() > -latitude.min() else "south_viewing"


def main() -> None:
    """Validates the entry on a day of data, or fits its values."""
    run = parse_run("saber", __doc__)
    day = run.day
    scans = read_scans(download_day(day, run.directory), day)
    period, offset = measure_cycle(scans)
    side = get_side(scans)
    print(
        f"{day:%Y-%m-%d}: {len(scans)} events, {side}, cycle period {period:.3f} s, "
        f"downward offset {offset:.3f} s",
        flush=True,
    )
    tle = run.element_set()
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    if run.fit:
        sounding = run.entry.configurations[side].sounding
        bottom = float(np.median([s.altitude.min() for s in scans]))
        top = float(np.median([s.altitude.max() for s in scans]))
        middle = round((bottom + top) / 2)
        elevations = tuple(float(e) for e in np.linspace(bottom, top, 13).round())
        sounding = sounding.model_copy(
            update={"scan_elevations": elevations, "sample_elevation": middle}
        )
        satellite = run.entry.to_satellite(
            orbit, configuration=side, platform=run.platform
        )
        duration = fit_duration(sounding, satellite, scans, middle)
        print(
            json.dumps(
                {
                    "cycle_period": round(period, 3),
                    "downward_offset": round(offset, 3),
                    "scan_elevations": list(elevations),
                    "sample_elevation": float(middle),
                    "scan_duration": round(duration, 3),
                    "event_duration_median": float(
                        np.median([s.elapsed[-1] for s in scans])
                    ),
                }
            )
        )
        return
    for key in run.configurations:
        if key != side:
            print(f"{key}: not tested, as SABER viewed to the other side on this day")
            continue
        satellite = run.entry.to_satellite(
            orbit, configuration=key, platform=run.platform
        )
        validate_limb(
            run,
            key,
            satellite,
            tle,
            scans,
            period,
            "SABER Level 2A version 2.0 (tangent points of the events)",
        )


if __name__ == "__main__":
    main()

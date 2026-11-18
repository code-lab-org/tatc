"""
Validates the TGRS catalog entry (the radio occultation receivers of COSMIC-2) against
UCAR's processing of a day of COSMIC-2 occultations, in the AWS Registry of Open Data
(read with `awsgnssroutils`), using the exact values of the entry through public TAT-C
methods. The database lists every occultation recorded on the day, with its receiver,
transmitter, and direction, and each occultation's calibrated phase (Level 1A) file
records the receiver's and transmitter's Earth-fixed positions.

* Receivers and transmitters are identified by the satellite whose element sets
  reproduce their recorded positions, in a sample of files that holds every receiver
  and transmitter of the day, keeping each satellite's element sets within 20 km of
  them.
* Each recorded occultation's time at the sample elevation is computed from these
  orbits, and the occultation is predicted if the entry's occultations (with every
  identified transmitter of the constellations tracked, over the day) include one of
  the same receiver, transmitter, and direction within 60 s (`occultation_recall`).
  Predicted occultations that were not recorded are not errors, as receivers track
  transmitters selectively; their share is reported but not gated.
* In a sample of occultations, the straight line between the recorded positions gives
  the time and tangent point at the sample elevation, which are compared with the
  predicted occultation's time and position (`epoch_p95`, `tangent_position_median`).

The maximum yaw comes from `--fit`, on a day of data: the largest yaw, from the
receiver's velocity or anti-velocity, of the recorded occultations at the sample
elevation. Running on a day at least 30 days from the fit tests it on data held out
from the fit.

Reads anonymous public data, and needs the `validation` extra and `awsgnssroutils`:

    python validation/catalog/validate_catalog_tgrs.py <data directory> 2026-05-14 --fit
    python validation/catalog/validate_catalog_tgrs.py <data directory> 2026-11-15

Without a day, it takes that of the entry's record (or, with `--fit`, of the fit). The
element sets of the receivers and transmitters of each day are in the files of
`ELEMENT_FILES`. With `--write`, it writes the entry's validation record and the
reduced fixture from which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import collections
import json
import os
import pathlib
from datetime import datetime, timedelta, timezone

import netCDF4
import numpy as np
import pandas as pd
from awsgnssroutils.database import RODatabaseClient
from joblib import Parallel, delayed
from pyproj import Geod, Transformer
from skyfield.api import Distance, wgs84
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric

from tatc.analysis.ro_sampling import _tangent_point_geometry
from tatc.constants import timescale
from tatc.schemas import GeneralPerturbationsOrbit, Satellite
from tatc.utils import compute_tangent_point, compute_vnb_frame
from validation_io import ELEMENT_SET_FOLDER, REPOSITORY, parse_run, write_results

CONFIGURATION = "occultations"
CONSTELLATIONS = {"G": "gps", "R": "glonass", "E": "galileo"}
NOTEBOOK_ELEMENTS = REPOSITORY / "docs" / "validation" / "limb_and_space" / "elements"
# the files of element sets of each day's receivers and transmitters
ELEMENT_FILES = {
    "2026-05-14": [
        NOTEBOOK_ELEMENTS / f"{name}_2026-05-14.3le"
        for name in ("gps", "glonass", "galileo", "cosmic2")
    ],
    "2026-11-15": [
        ELEMENT_SET_FOLDER / f"{name}_2026-11-18.3le"
        for name in ("gps", "glonass", "galileo", "cosmic2")
    ],
}
GPS_EPOCH = datetime(1980, 1, 6, tzinfo=timezone.utc)
# GPS time minus UTC (s), from 2017
GPS_UTC_OFFSET = 18
# the largest error (km) of a kept element set
MAX_ERROR = 20
# the files per transmitter and receiver read to identify it
FILES_PER_SATELLITE = 2
# the occultations per constellation whose tangent points are compared
TRUTH_SAMPLES = 100
FIXTURE_SAMPLES = 60
# the largest difference (s) of the times of paired occultations
PAIR_TOLERANCE = 60
# the window (s) either side of a database time stamp in which an occultation is found
SEARCH_WINDOW = 600
# the window either side of a fixture's occultation in which tests predict it
FIXTURE_WINDOW = timedelta(minutes=10)
N_JOBS = min(16, os.cpu_count() or 1)
GEOD = Geod(ellps="WGS84")
TO_GEODETIC = Transformer.from_crs("EPSG:4978", "EPSG:4979", always_xy=True)


def read_element_sets(paths) -> tuple[dict, dict]:
    """
    Reads files of element sets: the element sets of each NORAD catalog number, by
    epoch, and the satellites' names.
    """
    lines = [
        line.rstrip()
        for path in paths
        for line in pathlib.Path(path).read_text().splitlines()
        if line.strip()
    ]
    element_sets, names = collections.defaultdict(dict), {}
    for i in range(0, len(lines), 3):
        line1, line2 = lines[i + 1], lines[i + 2]
        names[line1[2:7]] = lines[i][2:].strip()
        element_sets[line1[2:7]][float(line1[18:32])] = [line1, line2]
    return element_sets, names


def to_lines(sets: dict) -> list[str]:
    """Gets the lines of element sets, in the order of their epochs."""
    return sum((sets[epoch] for epoch in sorted(sets)), [])


def identify(positions, element_sets) -> tuple:
    """
    Identifies the satellite whose element sets best reproduce recorded Earth-fixed
    positions [(time, m)], keeping its element sets within `MAX_ERROR` of all of
    them: the NORAD catalog number, the kept element sets, and the median error.
    """
    times = timescale.from_datetimes([t for t, _ in positions])
    rotation = itrs.rotation_at(times)
    recorded = np.array([p for _, p in positions]).T

    def errors(sets):
        orbit = GeneralPerturbationsOrbit.from_tle(to_lines(sets))
        p = np.einsum(
            "ij...,j...->i...",
            rotation,
            np.array(orbit.get_orbit_track_at_time(times).position.m),
        )
        return np.linalg.norm(p - recorded, axis=0) / 1e3

    best = {
        norad: min(np.median(errors({e: s})) for e, s in sets.items())
        for norad, sets in element_sets.items()
    }
    norad = min(best, key=best.get)
    kept = {
        e: s
        for e, s in element_sets[norad].items()
        if errors({e: s}).max() <= MAX_ERROR
    }
    return norad, kept, float(np.median(errors(kept))) if kept else np.nan


def read_phase(path: pathlib.Path) -> dict:
    """
    Reads the times (UTC) and Earth-fixed positions (m) of the receiver and
    transmitter in a calibrated phase file, without invalid samples.
    """
    with netCDF4.Dataset(path) as data:
        data.set_auto_mask(False)
        seconds = data["time"][:].astype(float)
        receiver = data["positionLEO"][:].astype(float)
        transmitter = data["positionGNSS"][:].astype(float)
        start = float(data["startTime"][...])
    ok = (
        (np.abs(seconds) < 1e10)
        & np.all(np.abs(receiver) < 1e10, axis=1)
        & np.all(np.abs(transmitter) < 1e10, axis=1)
    )
    base = GPS_EPOCH + timedelta(seconds=start - GPS_UTC_OFFSET)
    return {
        "base": base,
        "seconds": seconds[ok],
        "receiver": receiver[ok],
        "transmitter": transmitter[ok],
    }


def get_truth(path: pathlib.Path, elevation: float):
    """
    Gets the time and the longitude and latitude at which the straight line between
    the recorded positions (at 10 Hz) reaches an elevation, or `None`.
    """
    phase = read_phase(path)
    i = slice(None, None, 10)
    seconds = phase["seconds"][i]
    points = compute_tangent_point(
        phase["receiver"][i].T, phase["transmitter"][i].T - phase["receiver"][i].T
    )
    height = np.asarray(TO_GEODETIC.transform(*points)[2])
    crossings = np.flatnonzero(np.diff(np.sign(height - elevation)))
    if len(crossings) == 0:
        return None
    j = crossings[0]
    f = (elevation - height[j]) / (height[j + 1] - height[j])
    point = points[:, j] + f * (points[:, j + 1] - points[:, j])
    longitude, latitude, _ = TO_GEODETIC.transform(*point)
    return (
        phase["base"]
        + timedelta(seconds=seconds[j] + f * (seconds[j + 1] - seconds[j])),
        float(longitude),
        float(latitude),
    )


def find_crossings(receiver, transmitter, stamps, elevation) -> list[tuple]:
    """
    Gets the time at which the straight line reaches an elevation, and the
    transmitter's yaw in the receiver's frame, within `SEARCH_WINDOW` of each time
    stamp (NaT and NaN if it does not).
    """
    offsets = np.arange(-float(SEARCH_WINDOW), float(SEARCH_WINDOW), 1.0)
    results = []
    for stamp in stamps:
        times = timescale.from_datetimes(
            [stamp + timedelta(seconds=s) for s in offsets]
        )
        rx = receiver.orbit.to_gp_orbit().get_orbit_track_at_time(times)
        tx = transmitter.orbit.to_gp_orbit().get_orbit_track_at_time(times)
        point, sign, _, yaw = _tangent_point_geometry(tx, rx, *compute_vnb_frame(rx))
        height = wgs84.geographic_position_of(
            Geocentric(Distance(m=point).au, None, times)
        ).elevation.m
        height[sign > 0] = np.nan
        s = np.sign(height - elevation)
        crossings = np.flatnonzero(s[:-1] * s[1:] < 0)
        if len(crossings) == 0:
            results.append((pd.NaT, np.nan))
            continue
        j = crossings[0]
        w = (height[j] - elevation) / (height[j] - height[j + 1])
        results.append(
            (
                stamp + timedelta(seconds=offsets[j] + w),
                yaw[j] + w * (((yaw[j + 1] - yaw[j] + 180) % 360) - 180),
            )
        )
    return results


def fold(yaw) -> np.ndarray:
    """Gets the yaw (deg) from the nearer of the velocity and anti-velocity."""
    yaw = np.abs(np.asarray(yaw, dtype=float))
    return np.minimum(yaw, 180 - yaw)


def predict(sounding, receiver, transmitter, start, end) -> pd.DataFrame:
    """Predicts the occultations of a receiver and transmitter in a period."""
    return sounding.collect_observations(
        receiver, transmitter, start, end, time_step=timedelta(seconds=1)
    )


def predict_near(sounding, receiver, transmitter, time, is_rising):
    """
    Predicts an occultation of a receiver and transmitter in the direction given,
    within `FIXTURE_WINDOW` of a time (as the unit tests recompute it): its time and
    position, or `None`.
    """
    found = sounding.collect_observations(
        receiver,
        transmitter,
        time - FIXTURE_WINDOW,
        time + FIXTURE_WINDOW,
        time_step=timedelta(seconds=1),
    )
    found = found[found.is_rising == is_rising]
    if found.empty:
        return None
    k = int(np.argmin(np.abs((found.time - pd.Timestamp(time)).dt.total_seconds())))
    return found.time.iloc[k].to_pydatetime(), found.position.iloc[k]


def pair(predicted: pd.DataFrame, observed: pd.DataFrame) -> pd.DataFrame:
    """
    Pairs predicted and observed occultations one to one, by receiver, transmitter,
    direction, and nearest time within `PAIR_TOLERANCE`.
    """
    candidates = []
    for (rx, tx), p in predicted.groupby(["receiver", "transmitter"]):
        o = observed[(observed.receiver == rx) & (observed.transmitter == tx)]
        if o.empty:
            continue
        dt = np.abs(
            (p.time.values[:, None] - o.time.values[None, :]) / np.timedelta64(1, "s")
        )
        same = p.is_rising.values[:, None] == o.is_rising.values[None, :]
        candidates += [
            (dt[i, j], p.index[i], o.index[j])
            for i, j in zip(*np.nonzero((dt < PAIR_TOLERANCE) & same))
        ]
    used_p, used_o, pairs = set(), set(), []
    for _, i, j in sorted(candidates):
        if i not in used_p and j not in used_o:
            used_p.add(i)
            used_o.add(j)
            pairs.append((i, j))
    return pd.DataFrame(pairs, columns=["predicted", "observed"])


def main() -> None:
    """Validates the entry on a day of data, or fits its maximum yaw."""
    run = parse_run("tgrs", __doc__)
    day = run.day
    key = f"{day:%Y-%m-%d}"
    if key not in ELEMENT_FILES:
        raise SystemExit(f"no element sets of {key}: add their files to ELEMENT_FILES")
    sounding = run.entry.configurations[CONFIGURATION].sounding
    elevation = sounding.sample_elevation
    tracked = {
        code for code, name in CONSTELLATIONS.items() if name in sounding.constellations
    }

    # the day's occultations, and a sample of their files
    database = RODatabaseClient(
        metadata_root=str(run.directory / "metadata"), version="v1.1"
    )
    occultations = database.query(
        missions="cosmic2",
        datetimerange=(key, key + "T23:59:59"),
        silent=True,
    ).filter(availablefiletypes=["ucar_calibratedPhase"])
    records = [r for r in occultations._data if r["transmitter"][0] in tracked]
    rng = np.random.default_rng(int(day.strftime("%Y%m%d")))
    chosen = set()
    by_satellite = collections.defaultdict(list)
    for i, r in enumerate(records):
        by_satellite[r["transmitter"]].append(i)
        by_satellite[r["receiver"]].append(i)
    for indices in by_satellite.values():
        chosen |= set(
            rng.choice(indices, min(FILES_PER_SATELLITE, len(indices)), replace=False)
        )
    truth_indices = set()
    for code in sorted(tracked):
        pool = [i for i, r in enumerate(records) if r["transmitter"][0] == code]
        truth_indices |= set(rng.choice(pool, TRUTH_SAMPLES, replace=False))
    chosen |= truth_indices
    sample = occultations[0:0]
    sample._data = [records[i] for i in sorted(chosen)]
    sample.size = len(sample._data)
    files = run.directory / "files"
    sample.download(
        "ucar_calibratedPhase",
        data_root=str(files),
        keep_aws_structure=False,
        silent=True,
    )

    def local(record) -> pathlib.Path:
        return files / os.path.basename(record["ucar_calibratedPhase"])

    print(
        f"{key}: {len(records)} occultations of {sorted(tracked)}, "
        f"{sample.size} files read",
        flush=True,
    )

    # identify receivers and transmitters
    element_sets, names = read_element_sets(ELEMENT_FILES[key])
    recorded = collections.defaultdict(list)
    for r in sample._data:
        phase = read_phase(local(r))
        for k in (0, len(phase["seconds"]) // 2, len(phase["seconds"]) - 1):
            t = phase["base"] + timedelta(seconds=float(phase["seconds"][k]))
            recorded[r["transmitter"]].append((t, phase["transmitter"][k]))
            recorded[r["receiver"]].append((t, phase["receiver"][k]))
    keys = sorted(recorded)
    identified = Parallel(n_jobs=N_JOBS)(
        delayed(identify)(recorded[k], element_sets) for k in keys
    )
    satellites, lines, errors = {}, {}, {}
    for k, (norad, kept, error) in zip(keys, identified):
        if kept:
            lines[k] = to_lines(kept)
            satellites[k] = Satellite(
                name=k, orbit=GeneralPerturbationsOrbit.from_tle(lines[k])
            )
            errors[k] = (norad, names[norad], error)
    receivers = sorted(k for k in satellites if k.startswith("cosmic2"))
    transmitters = sorted(k for k in satellites if k[0] in tracked)
    print(
        json.dumps(
            {
                "receivers": {k: errors[k][:2] for k in receivers},
                "transmitters": len(transmitters),
                "unidentified": sorted(set(keys) - set(satellites)),
                "transmitter_error_max_median_km": max(
                    errors[k][2] for k in transmitters
                ),
            }
        ),
        flush=True,
    )

    # the recorded occultations, at the sample elevation
    observed = pd.DataFrame(
        [
            {
                "occid": r["occid"],
                "receiver": r["receiver"],
                "transmitter": r["transmitter"],
                "stamp": datetime.fromisoformat(r["time"]).replace(tzinfo=timezone.utc),
                "is_rising": not r["setting"],
            }
            for r in records
            if r["receiver"] in satellites and r["transmitter"] in satellites
        ]
    )
    groups = list(observed.groupby(["receiver", "transmitter"]).groups.items())
    found = Parallel(n_jobs=N_JOBS)(
        delayed(find_crossings)(
            satellites[rx], satellites[tx], list(observed.stamp[index]), elevation
        )
        for (rx, tx), index in groups
    )
    times, yaws = [pd.NaT] * len(observed), [np.nan] * len(observed)
    for (_, index), results in zip(groups, found):
        for k, (t, yaw) in zip(index, results):
            times[k], yaws[k] = t, yaw
    observed["time"] = pd.to_datetime(times, utc=True)
    observed["yaw"] = yaws
    unresolved = int(observed.time.isna().sum())
    observed = observed[
        (observed.time >= day) & (observed.time < day + timedelta(days=1))
    ].reset_index(drop=True)
    folded = fold(observed.yaw)
    print(
        json.dumps(
            {
                "observed": len(observed),
                "unresolved": unresolved,
                "folded_yaw_percentiles_50_99_999_max": np.percentile(
                    folded, [50, 99, 99.9, 100]
                )
                .round(2)
                .tolist(),
            }
        ),
        flush=True,
    )
    if run.fit:
        print(json.dumps({"max_yaw": float(np.ceil(folded.max()))}))
        return

    # the predicted occultations of every receiver and transmitter over the day
    results = Parallel(n_jobs=N_JOBS)(
        delayed(predict)(
            sounding,
            satellites[rx],
            satellites[tx],
            day,
            day + timedelta(days=1),
        )
        for rx in receivers
        for tx in transmitters
    )
    predicted = pd.concat([r for r in results if len(r)], ignore_index=True)
    pairs = pair(predicted, observed)
    recall = 100 * len(pairs) / len(observed)
    precision = {
        CONSTELLATIONS[code]: float(
            100
            * (predicted.loc[pairs.predicted, "transmitter"].str[0] == code).sum()
            / (predicted.transmitter.str[0] == code).sum()
        )
        for code in sorted(tracked)
    }

    # tangent points of the sample, from the recorded positions
    truth_ids = {records[i]["occid"]: records[i] for i in truth_indices}
    rows = []
    for p, o in zip(pairs.predicted, pairs.observed):
        occultation = observed.loc[o]
        if occultation.occid not in truth_ids:
            continue
        truth = get_truth(local(truth_ids[occultation.occid]), elevation)
        if truth is None:
            continue
        prediction = predicted.loc[p]
        rows.append(
            {
                "occid": occultation.occid,
                "receiver": occultation.receiver,
                "transmitter": occultation.transmitter,
                "is_rising": bool(occultation.is_rising),
                "true_time": truth[0],
                "true_longitude": truth[1],
                "true_latitude": truth[2],
                "time_error": (
                    prediction.time.to_pydatetime() - truth[0]
                ).total_seconds(),
                "distance": GEOD.inv(
                    prediction.position.x, prediction.position.y, truth[1], truth[2]
                )[2]
                / 1e3,
            }
        )
    accuracy = pd.DataFrame(rows)
    metrics = {
        "occultation_recall": recall,
        "tangent_position_median": float(accuracy.distance.median()),
        "tangent_position_p95": float(np.percentile(accuracy.distance, 95)),
        "epoch_median": float(accuracy.time_error.median()),
        "epoch_p95": float(np.percentile(np.abs(accuracy.time_error), 95)),
        "observed_occultations": float(len(observed)),
        "predicted_occultations": float(len(predicted)),
        "tangent_points_compared": float(len(accuracy)),
        **{f"recorded_share_{name}": value for name, value in precision.items()},
    }
    print(json.dumps({CONFIGURATION: metrics}), flush=True)
    if not run.write:
        return

    # a sample of the compared occultations, whose predictions tests recompute
    keep = np.sort(rng.choice(len(accuracy), FIXTURE_SAMPLES, replace=False))
    sample_rows = accuracy.iloc[keep]
    used = sorted(set(sample_rows.receiver) | set(sample_rows.transmitter))
    fixture = {
        "day": np.array(day.isoformat()),
        "satellites": np.array(used),
        "element_sets": np.array(["\n".join(lines[k]) for k in used]),
        "receiver": sample_rows.receiver.to_numpy(str),
        "transmitter": sample_rows.transmitter.to_numpy(str),
        "is_rising": sample_rows.is_rising.to_numpy(bool),
        "true_seconds": np.array(
            [(t - day).total_seconds() for t in sample_rows.true_time]
        ),
        "true_longitude": sample_rows.true_longitude.to_numpy(float),
        "true_latitude": sample_rows.true_latitude.to_numpy(float),
    }
    # the results on the fixture, predicted near each occultation as tests do
    distances, time_errors = [], []
    for row in sample_rows.itertuples():
        near = predict_near(
            sounding,
            satellites[row.receiver],
            satellites[row.transmitter],
            row.true_time,
            row.is_rising,
        )
        time_errors.append((near[0] - row.true_time).total_seconds())
        distances.append(
            GEOD.inv(near[1].x, near[1].y, row.true_longitude, row.true_latitude)[2]
            / 1e3
        )
    write_results(
        run,
        CONFIGURATION,
        reference="UCAR COSMIC-2 occultations (AWS GNSS RO database, v1.1) and "
        "calibrated phase files (receiver and transmitter positions)",
        metrics=metrics,
        fixture=fixture,
        fixture_metrics={
            "tangent_position_median": float(np.median(distances)),
            "epoch_p95": float(np.percentile(np.abs(time_errors), 95)),
        },
    )


if __name__ == "__main__":
    np.seterr(invalid="ignore")
    main()

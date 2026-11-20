"""
Scoring of the catalog entries of limb sounders, shared by their validation scripts:
the scans of a day of reference data, with their tangent points; the fit of a scan's
values; and the comparison of the entry's scans with them through public TAT-C
methods.

* The scans are started at the times of the reference's scans, and the representative
  point of each (at the sample elevation) is compared with the reference's tangent
  point at that elevation, interpolated along the scan: the horizontal distance
  (`tangent_position_median`) and the difference of times (`epoch_p95`).
* The period of the cycle of scans is compared with the median period of the
  reference's cycles (`scan_period_error`).
* A day of scans repeating from midnight, as an analysis plans them, is compared with
  the reference's tangent points at the sample elevation: the error of the extent of
  latitudes (`latitude_extent_error`).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import json
from dataclasses import dataclass
from datetime import datetime, timedelta

import numpy as np
from pyproj import Geod

from tatc.analysis import collect_limb_observations
from tatc.analysis.limb_sampling import ScanDirection
from validation_io import write_results

GEOD = Geod(ellps="WGS84")
FIXTURE_SCANS = 80
# the largest difference (seconds) of the start of a scan from the reference's
MAX_START_DIFFERENCE = 0.5


@dataclass
class Scan:
    """
    A reference scan: its start, direction, and samples in the order of their times
    (seconds from the start, and the tangent points' geodetic altitudes in meters,
    latitudes, and longitudes in degrees).
    """

    start: datetime
    direction: ScanDirection
    elapsed: np.ndarray
    altitude: np.ndarray
    latitude: np.ndarray
    longitude: np.ndarray

    def at_elevation(self, elevation: float) -> tuple[float, float, float]:
        """
        Gets the time (seconds from the start), latitude, and longitude at which the
        tangent point first reaches an elevation, interpolated between samples, or
        NaN if the scan does not reach it.
        """
        above = self.altitude >= elevation
        crossings = np.flatnonzero(above[1:] != above[:-1])
        if len(crossings) == 0:
            return np.nan, np.nan, np.nan
        i = crossings[0]
        f = (elevation - self.altitude[i]) / (self.altitude[i + 1] - self.altitude[i])
        longitude = np.unwrap(self.longitude[i : i + 2], period=360)
        return (
            float(self.elapsed[i] + f * (self.elapsed[i + 1] - self.elapsed[i])),
            float(self.latitude[i] + f * (self.latitude[i + 1] - self.latitude[i])),
            float((longitude[0] + f * (longitude[1] - longitude[0]) + 180) % 360 - 180),
        )


def get_reference_points(scans: list[Scan], elevation: float) -> dict[str, np.ndarray]:
    """
    Gets the reference's tangent points of scans at an elevation: their starts
    (seconds from the first scan's day), directions (1 for upward), and the times
    (seconds from the start), latitudes, and longitudes at the elevation.
    """
    day = scans[0].start.replace(hour=0, minute=0, second=0, microsecond=0)
    values = np.array([scan.at_elevation(elevation) for scan in scans]).reshape(-1, 3)
    return {
        "start": np.array([(s.start - day).total_seconds() for s in scans]),
        "upward": np.array([s.direction == ScanDirection.UPWARD for s in scans]),
        "elapsed": values[:, 0],
        "latitude": values[:, 1],
        "longitude": values[:, 2],
    }


def predict_matched(sounding, satellite, day: datetime, starts, upward) -> dict:
    """
    Gets the representative points of the scans of a sounding started at the given
    starts (seconds from a day) in the given directions: the elapsed seconds from
    each start, latitudes, and longitudes (NaN for a scan without one).
    """
    elapsed = np.full(len(starts), np.nan)
    latitude = np.full(len(starts), np.nan)
    longitude = np.full(len(starts), np.nan)
    for direction in (ScanDirection.UPWARD, ScanDirection.DOWNWARD):
        index = np.flatnonzero(upward == (direction == ScanDirection.UPWARD))
        if len(index) == 0:
            continue
        times = [day + timedelta(seconds=float(starts[i])) for i in index]
        observations = collect_limb_observations(
            satellite,
            times,
            sounding.scan_azimuth,
            list(sounding.scan_elevations),
            timedelta(seconds=sounding.scan_duration),
            scan_direction=direction,
            sample_elevation=sounding.sample_elevation,
        )
        seconds = np.array([(t - day).total_seconds() for t in times])
        for row in observations.itertuples():
            row_start = (row.start.to_pydatetime() - day).total_seconds()
            j = int(np.argmin(np.abs(seconds - row_start)))
            if abs(seconds[j] - row_start) > MAX_START_DIFFERENCE:
                continue
            k = index[j]
            elapsed[k] = (row.time.to_pydatetime() - day).total_seconds() - starts[k]
            latitude[k] = row.position.y
            longitude[k] = row.position.x
    return {"elapsed": elapsed, "latitude": latitude, "longitude": longitude}


def compare_matched(predicted: dict, reference: dict) -> dict[str, np.ndarray]:
    """
    Gets the horizontal distances (km) and differences of times (seconds) of
    predicted representative points from the reference's, where both exist.
    """
    ok = np.isfinite(predicted["latitude"]) & np.isfinite(reference["latitude"])
    distance = (
        GEOD.inv(
            predicted["longitude"][ok],
            predicted["latitude"][ok],
            reference["longitude"][ok],
            reference["latitude"][ok],
        )[2]
        / 1e3
    )
    return {
        "distance": distance,
        "time": predicted["elapsed"][ok] - reference["elapsed"][ok],
        "compared": ok,
    }


def predict_latitudes(sounding, satellite, day: datetime) -> np.ndarray:
    """
    Gets the latitudes of the representative points of a day of scans repeating
    from midnight.
    """
    observations = sounding.collect_observations(
        satellite, day, day + timedelta(days=1)
    )
    return np.array([p.y for p in observations.position])


def get_latitude_extent_error(predicted, reference) -> float:
    """
    Gets the larger error (degrees) of the least and greatest latitudes.
    """
    reference = np.asarray(reference)[np.isfinite(reference)]
    return float(
        max(
            abs(np.min(predicted) - np.min(reference)),
            abs(np.max(predicted) - np.max(reference)),
        )
    )


def get_scan_period_error(sounding, measured_period: float) -> float:
    """
    Gets the error (%) of the period of a sounding's cycle of scans.
    """
    return float(100 * abs(sounding.cycle_period - measured_period) / measured_period)


def fit_duration(sounding, satellite, scans: list[Scan], elevation: float) -> float:
    """
    Fits a scan duration so that the scans' representative points (at an
    elevation) are reached at the median time at which the reference's scans reach
    it: as a scan's samples are timed in proportion to the change of the viewing
    angle, the elapsed time at an elevation is proportional to the duration, so a
    scan of 1 s gives the proportion.
    """
    unit = sounding.model_copy(
        update={"scan_duration": 1.0, "sample_elevation": elevation}
    )
    reference = get_reference_points(scans, elevation)
    day = scans[0].start.replace(hour=0, minute=0, second=0, microsecond=0)
    predicted = predict_matched(
        unit, satellite, day, reference["start"], reference["upward"]
    )
    proportion = np.nanmedian(predicted["elapsed"])
    return float(np.nanmedian(reference["elapsed"]) / proportion)


def validate_limb(
    run,
    key: str,
    satellite,
    tle: list[str],
    scans: list[Scan],
    measured_period: float,
    reference_name: str,
    seed: int = 1,
) -> None:
    """
    Scores a configuration of a limb sounder on a day of reference scans, prints
    its metrics, and writes its results if the run writes.
    """
    sounding = run.entry.configurations[key].sounding
    day = run.day
    scans = [s for s in scans if s.direction in {x.direction for x in sounding.scans}]
    reference = get_reference_points(scans, sounding.sample_elevation)
    predicted = predict_matched(
        sounding, satellite, day, reference["start"], reference["upward"]
    )
    matched = compare_matched(predicted, reference)
    latitudes = predict_latitudes(sounding, satellite, day)
    metrics = {
        "tangent_position_median": float(np.median(matched["distance"])),
        "tangent_position_p95": float(np.percentile(matched["distance"], 95)),
        "epoch_median": float(np.median(matched["time"])),
        "epoch_p95": float(np.percentile(np.abs(matched["time"]), 95)),
        "scan_period_error": get_scan_period_error(sounding, measured_period),
        "latitude_extent_error": get_latitude_extent_error(
            latitudes, reference["latitude"]
        ),
        "scans_compared": float(matched["compared"].sum()),
        "reference_scans": float(len(scans)),
        "scans_per_day": float(len(latitudes)),
    }
    print(json.dumps({key: metrics}), flush=True)
    if not run.write:
        return
    # a sample of the compared scans, from which tests recompute the results
    rng = np.random.default_rng(seed)
    keep = np.sort(
        rng.choice(
            np.flatnonzero(matched["compared"]),
            min(FIXTURE_SCANS, int(matched["compared"].sum())),
            replace=False,
        )
    )
    sample = {k: v[keep] for k, v in reference.items()}
    sample_matched = compare_matched({k: v[keep] for k, v in predicted.items()}, sample)
    # keyed by configuration, whose days may differ (`day_<key>`, as fixtures name
    # the days of their parts)
    fixture = {
        f"tle_{key}": np.array(tle),
        f"day_{key}": np.array(day.isoformat()),
    }
    fixture |= {
        f"{key}_{name}": values
        for name, values in {
            "start": sample["start"],
            "upward": sample["upward"],
            "elapsed": sample["elapsed"],
            "latitude": sample["latitude"],
            "longitude": sample["longitude"],
            "latitude_extent": np.array(
                [np.nanmin(reference["latitude"]), np.nanmax(reference["latitude"])]
            ),
            "measured_period": np.array(measured_period),
        }.items()
    }
    write_results(
        run,
        key,
        reference=reference_name,
        metrics=metrics,
        fixture=fixture,
        fixture_metrics={
            "tangent_position_median": float(np.median(sample_matched["distance"])),
            "epoch_p95": float(np.percentile(np.abs(sample_matched["time"]), 95)),
            "scan_period_error": metrics["scan_period_error"],
            "latitude_extent_error": metrics["latitude_extent_error"],
        },
    )

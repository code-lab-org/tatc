"""
Validates the MLS catalog entry against Aura's MLS Level 1B orbit and attitude product
(ML1OA, version 6, of NASA GES DISC), which records, for every minor frame (about 1/6
s) of every major frame (a vertical scan, about 24.6 s), the geodetic latitude,
longitude, and altitude of the GHz radiometer's tangent point, using the exact values
of the entry through public TAT-C methods (see `sounding.py` for the scoring). The
upward sweep of each major frame is compared; major frames whose line of sight leaves
the orbit plane by more than 20 km (attitude maneuvers) are excluded.

The scan's values come from `--fit`, on a day of data: the cycle period is the median
period of major frames, the scan's bottom and top are those of the median sweep, and
the sample elevation is the median altitude of the Level 1B tangent points nearest the
positions of the Level 2 temperature profiles (ML2T), at which the scan's duration is
fitted so that the representative point is reached when the median sweep reaches
it. TAT-C's scans sweep at a constant angular rate, whereas MLS dwells in the lower
atmosphere, so only the representative point is fitted. Running on a day at least 30
days from the fit tests the values on data held out from the fit.

Needs the `validation` extra and NASA Earthdata credentials (`EARTHDATA_TOKEN` or
~/.netrc):

    python validation/catalog/validate_catalog_mls.py <data directory> 2026-08-14 --fit
    python validation/catalog/validate_catalog_mls.py <data directory> 2026-11-15

Without a day, it takes that of the entry's record (or, with `--fit`, of the fit). The
element sets are registered in `elements/element_sets.yml`. With `--write`, it writes
the entry's validation record and the reduced fixture from which unit tests recompute
the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import json
import pathlib
from datetime import datetime, timedelta, timezone

import earthaccess
import netCDF4
import numpy as np
from pyproj import Geod

from tatc.analysis.limb_sampling import ScanDirection
from tatc.schemas import GeneralPerturbationsOrbit
from sounding import Scan, fit_duration, validate_limb
from validation_io import login_earthdata, parse_run

CONFIGURATION = "scans"
TAI93 = datetime(1993, 1, 1, tzinfo=timezone.utc)
# the minor frames of the GHz radiometer in a major frame
GHZ_FRAMES = 125
# the largest offset (m) of the line of sight from the orbit plane of a nominal scan
MAX_ORBIT_OFFSET = 20e3
GEOD = Geod(ellps="WGS84")


def download(short_name: str, day: datetime, directory: pathlib.Path) -> pathlib.Path:
    """Downloads a day's file of a version 6 MLS product, if not cached."""
    pattern = {"ML1OA": "L1BOA", "ML2T": "L2GP-Temperature"}[short_name]
    name = f"*{pattern}_v06-*_{day:%Y}d{day.timetuple().tm_yday:03d}.*"
    cached = sorted(directory.glob(name))
    if cached:
        return cached[0]
    login_earthdata()
    granules = earthaccess.search_data(
        short_name=short_name,
        version="006",
        temporal=(f"{day:%Y-%m-%d}", f"{day:%Y-%m-%d}"),
    )
    earthaccess.download(granules, str(directory))
    return sorted(directory.glob(name))[0]


class LevelOne:
    """
    The scans of a day of the Level 1B orbit and attitude product: the upward sweep
    (the first minor frames, to the top of the median sweep) of each nominal major
    frame that starts in the day.
    """

    def __init__(self, path: pathlib.Path, day: datetime):
        with netCDF4.Dataset(path) as data:
            data.set_auto_mask(False)
            start_tai = data["MAFStartTimeTAI"][:].astype(float)
            start_utc = [
                datetime.strptime(s, "%Y-%m-%dT%H:%M:%S.%fZ").replace(
                    tzinfo=timezone.utc
                )
                for s in data["MAFStartTimeUTC"][:]
            ]
            frame_tai = data["sc/MIF_TAI"][:, :GHZ_FRAMES].astype(float)
            latitude, longitude, altitude = (
                data["GHz/" + v][:, :GHZ_FRAMES].astype(float)
                for v in ("GeodLat", "Lon", "GeodAlt")
            )
            offset = data["GHz/OrbY"][:, :GHZ_FRAMES].astype(float)
        # the offset of TAI from UTC, implied by the file's own times
        self.tai_offset = float(
            np.median(start_tai - [(t - TAI93).total_seconds() for t in start_utc])
        )
        valid = (np.abs(latitude) <= 90) & (altitude > -10e3) & (altitude < 200e3)
        nominal = np.all(np.abs(offset) < MAX_ORBIT_OFFSET, axis=1)
        in_day = np.array([day <= t < day + timedelta(days=1) for t in start_utc])
        median = np.nanmedian(np.where(valid, altitude, np.nan)[nominal], axis=0)
        self.top = int(np.nanargmax(median))
        self.bottom_altitude = float(median[0])
        self.top_altitude = float(median[self.top])
        # the median period of consecutive nominal major frames
        consecutive = nominal[1:] & nominal[:-1]
        self.period = float(np.median(np.diff(start_tai)[consecutive]))
        self.scans = []
        sweep = slice(0, self.top + 1)
        for i in np.flatnonzero(nominal & in_day):
            ok = valid[i, sweep]
            if ok.sum() < 2:
                continue
            self.scans.append(
                Scan(
                    start=start_utc[i],
                    direction=ScanDirection.UPWARD,
                    elapsed=(frame_tai[i, sweep] - start_tai[i])[ok],
                    altitude=altitude[i, sweep][ok],
                    latitude=latitude[i, sweep][ok],
                    longitude=longitude[i, sweep][ok],
                )
            )
        self.all_frames = {
            "seconds": (frame_tai - self.tai_offset) - (day - TAI93).total_seconds(),
            "latitude": np.where(valid, latitude, np.nan),
            "longitude": np.where(valid, longitude, np.nan),
            "altitude": np.where(valid, altitude, np.nan),
        }


def measure_profile_elevation(level_one: LevelOne, path: pathlib.Path, day) -> float:
    """
    Measures the elevation (m) of the positions of the Level 2 temperature profiles:
    the median altitude of the Level 1B tangent points nearest them, within a minute.
    """
    with netCDF4.Dataset(path) as data:
        data.set_auto_mask(False)
        geolocation = data["HDFEOS/SWATHS/Temperature/Geolocation Fields"]
        seconds = geolocation["Time"][:].astype(float) - level_one.tai_offset
        seconds -= (day - TAI93).total_seconds()
        latitude = geolocation["Latitude"][:].astype(float)
        longitude = geolocation["Longitude"][:].astype(float)
    frames = level_one.all_frames
    keep = (np.abs(latitude) <= 90) & (seconds >= 0) & (seconds < 86400)
    altitudes, distances = [], []
    for s, la, lo in zip(seconds[keep], latitude[keep], longitude[keep]):
        near = (np.abs(frames["seconds"] - s) < 60) & np.isfinite(frames["altitude"])
        if not near.any():
            continue
        distance = GEOD.inv(
            np.full(near.sum(), lo),
            np.full(near.sum(), la),
            frames["longitude"][near],
            frames["latitude"][near],
        )[2]
        k = int(np.argmin(distance))
        altitudes.append(frames["altitude"][near][k])
        distances.append(distance[k] / 1e3)
    print(
        json.dumps(
            {
                "profiles": len(altitudes),
                "profile_elevation_median": float(np.median(altitudes)),
                "profile_elevation_quartiles": np.percentile(
                    altitudes, [25, 75]
                ).tolist(),
                "nearest_distance_median_km": float(np.median(distances)),
            }
        ),
        flush=True,
    )
    return float(np.median(altitudes))


def main() -> None:
    """Validates the entry on a day of data, or fits its values."""
    run = parse_run("mls", __doc__)
    day = run.day
    level_one = LevelOne(download("ML1OA", day, run.directory), day)
    print(
        f"{day:%Y-%m-%d}: {len(level_one.scans)} nominal scans, period "
        f"{level_one.period:.6f} s, sweep {level_one.bottom_altitude:.0f} to "
        f"{level_one.top_altitude:.0f} m",
        flush=True,
    )
    tle = run.element_set()
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    satellite = run.entry.to_satellite(
        orbit, configuration=CONFIGURATION, platform=run.platform
    )
    if run.fit:
        elevation = measure_profile_elevation(
            level_one, download("ML2T", day, run.directory), day
        )
        sounding = run.entry.configurations[CONFIGURATION].sounding
        bottom, top = level_one.bottom_altitude, level_one.top_altitude
        elevations = sorted(
            set(np.linspace(bottom, top, 11).round()) | {round(elevation)}
        )
        sounding = sounding.model_copy(update={"scan_elevations": tuple(elevations)})
        duration = fit_duration(sounding, satellite, level_one.scans, round(elevation))
        print(
            json.dumps(
                {
                    "cycle_period": round(level_one.period, 4),
                    "scan_elevations": [float(e) for e in elevations],
                    "sample_elevation": float(round(elevation)),
                    "scan_duration": round(duration, 3),
                }
            )
        )
        return
    validate_limb(
        run,
        CONFIGURATION,
        satellite,
        tle,
        level_one.scans,
        level_one.period,
        "ML1OA version 6 (GHz tangent points of the upward sweeps of nominal major "
        "frames)",
    )


if __name__ == "__main__":
    main()

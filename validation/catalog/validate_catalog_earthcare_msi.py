"""
Validates the EarthCARE MSI catalog entry against EarthCARE's Level 1C MSI product
(MSI_RGR_1C, regridded to the pixels of the thermal infrared bands), using the exact
values of the entry through public TAT-C methods. MSI is a pushbroom imager of 150 km,
tilted away from the sun. The product gives each pixel's position and each line's time
but no spacecraft positions, so they are interpolated to the lines' times from those
of the ATLID product of the same day (ATL_NOM_1B), and the fit needs no element set.
The view's angles are fitted to one day (run with `--fit`, then `--fit-view` to center
the view in TAT-C's frame): running on a day at least 30 days from it tests the fitted
angles on data held out from the fit. The swath is scored as for other pushbroom
swaths: passes of points within the lattice of its pixels, the swath width, and the
swath area over the span of the day compared.

Needs the `validation` extra and a MAAP offline token in `ESA_MAAP_KEY`. Only the
pixels' positions and the lines' times are read, of every fourth orbit, and reduced
to every eighth line and some of the 357 columns of imaged pixels (of the grid's 384),
with range requests:

    python validation/catalog/validate_catalog_earthcare_msi.py <data directory> 2026-09-27 --fit
    python validation/catalog/validate_catalog_earthcare_msi.py <data directory> 2026-09-27 --fit-view
    python validation/catalog/validate_catalog_earthcare_msi.py <data directory> 2026-11-07

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
from earthcare_io import check_time_scale, load_day, to_seconds
from gesdisc_scanner import Reference, Setup, fit, fit_view, fit_yaw, validate_swath
from validation_io import (
    TO_GEODETIC,
    read_element_set_file,
    RECORDS,
    FIXTURES,
    source_of,
    ELEMENT_SET_FOLDER,
)
from validate_catalog_atlid import read_atlid

NORAD_ID = 59908
ELEMENT_SET_FILE = ELEMENT_SET_FOLDER / "earthcare_2026-09-27.3le"
PRODUCT = "MSI_RGR_1C"
# the product's grid has 384 columns, of which 12 to 368 hold the imaged pixels and
# the others are always fill
FIRST_COLUMN = 12
PIXELS = 357
LINE_STRIDE = 8
COLUMNS = np.r_[0:PIXELS:8, PIXELS - 1]
ORBIT_STRIDE = 4
# the body's yaw, and the along-track angle at nadir from `--fit` on the day of the
# fit (Earth-fixed frame), from which `--fit-view` centers the view with the
# instrument's pitch; `--fit` gives a slope of the line of -0.106 deg from the
# geodetic nadir, but in TAT-C's frame, rolled and pitched, an unyawed view slopes
# only 0.04 deg, which the entry does not model
FIT_YAW = 0.0
FIT_PITCH = -0.6198
SETUP = Setup(
    entry_id="earthcare_msi",
    platform="earthcare",
    swath="swath",
    beams="",
    size_instrument="",
    published_footprints=(),
    record_path=RECORDS / "earthcare_msi.json",
    fixture_path=FIXTURES / "earthcare_msi.npz",
    source=source_of(__file__),
    reference=f"{PRODUCT} (pixel positions of every fourth orbit, every eighth line), "
    "with spacecraft positions of ATL_NOM_1B",
    published_reference="",
    time_steps=(10, 30, 120),
    pass_radius=20,
    point_spacing=100e3,
)


def read_msi(f, day: datetime) -> dict:
    """Reads a frame's line times and pixel positions, reduced."""
    g = f["ScienceData"]
    lines = np.arange(0, g["time"].shape[0], LINE_STRIDE)
    return {
        "time": to_seconds(g["time"][:], day)[lines],
        "lat": g["latitude"][:][lines][:, FIRST_COLUMN + COLUMNS],
        "lon": g["longitude"][:][lines][:, FIRST_COLUMN + COLUMNS],
    }


def load_reference(day: datetime, directory: pathlib.Path) -> Reference:
    """
    Loads a reduced day, with the spacecraft's position interpolated, in Earth-fixed
    coordinates, from the sensor positions the ATLID product records.
    """
    data = load_day(PRODUCT, day, directory, read_msi, orbit_stride=ORBIT_STRIDE)
    lidar = load_day("ATL_NOM_1B", day, directory, read_atlid)
    order = np.argsort(lidar["time"])
    seconds, position = lidar["time"][order], lidar["sc_position"][order]
    near = np.abs(
        seconds[np.clip(np.searchsorted(seconds, data["time"]), 0, len(seconds) - 1)]
        - data["time"]
    )
    within = (data["time"] >= seconds[0]) & (data["time"] <= seconds[-1]) & (near < 5)
    data = {k: v[within] for k, v in data.items() if k != "frames"}
    interpolated = np.stack(
        [np.interp(data["time"], seconds, position[:, axis]) for axis in range(3)],
        axis=-1,
    )
    lon, lat, height = TO_GEODETIC.transform(*interpolated.T)
    data.update(
        sc_lat=lat,
        sc_lon=lon,
        sc_alt=height,
        pixel=COLUMNS,
        pixels=PIXELS,
        line_stride=LINE_STRIDE,
    )
    period = float(np.median(np.diff(data["time"])))
    return Reference(day, data, period)


def main() -> None:
    """Validates the entry on a day of data, or fits its angles."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("day")
    parser.add_argument("--fit", action="store_true")
    parser.add_argument(
        "--fit-view",
        action="store_true",
        help="center the view in TAT-C's frame, with the fitted yaw",
    )
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)
    day = datetime.fromisoformat(arguments.day).replace(tzinfo=timezone.utc)
    reference = load_reference(day, arguments.directory)
    print(
        f"{arguments.day}: {len(reference.time)} lines (every {LINE_STRIDE}), "
        f"{len(reference.gaps)} gaps over a minute, line period "
        f"{reference.scan_period:.4f} s",
        flush=True,
    )
    if arguments.fit:
        check_time_scale(PRODUCT, day, "ScienceData/time")
        for frame in ("inertial", "earth_fixed"):
            fit(reference, frame)
            fit_yaw(reference, frame)
        return
    if arguments.fit_view:
        for frame in ("inertial", "earth_fixed"):
            fit_view(reference, frame, FIT_YAW, FIT_PITCH, instrument_pitch=True)
        return
    tle = read_element_set_file(ELEMENT_SET_FILE, NORAD_ID, day)
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    entry = get_entry("earthcare_msi")
    validate_swath(SETUP, day, tle, entry, orbit, reference, arguments.write)
    if arguments.write:
        print(f"wrote {SETUP.record_path} and {SETUP.fixture_path}")


if __name__ == "__main__":
    main()

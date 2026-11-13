"""
Compares the swaths of the catalog's swept configurations at different time
steps, without reference data: TAT-C's swath against itself. For each
configuration whose instruments sweep along track (see
`tatc.schemas.PointedInstrument.is_swept`), over segments of ten minutes:

* the swept swath at 120 s against the swept swath at 5 s, which shows that it
  does not depend on the time step;
* the swept swath against the union of views tailored to 2 s, as TAT-C modeled
  swaths before swept footprints (views that tile the ground track near nadir,
  but leave thin gaps between them where wide views turn from one time step to
  the next);
* the union of views tailored to 120 s against that at 2 s, which shows how far
  views tailored to long time steps reach beyond the swath.

The orbit of each configuration is the element set of its fixture.

    python docs/validation/compare_swath_time_steps.py [entry ...]

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import pathlib
import warnings
from datetime import datetime, timedelta

import numpy as np
from pyproj import Transformer
from shapely.ops import transform as transform_geometry

from tatc.analysis import compute_ground_track
from tatc.catalog import get_entry, get_entry_ids
from tatc.schemas import GeneralPerturbationsOrbit

FIXTURES = (
    pathlib.Path(__file__).resolve().parents[2] / "tests" / "catalog" / "fixtures"
)
DURATION = 600  # seconds
OFFSETS = (3, 9, 15)  # hours after the start of the fixture's day


def get_fixture(
    entry_id: str, key: str, platform: str
) -> tuple[pathlib.Path, str] | None:
    """
    Gets the fixture whose element set orbits a configuration's platform, and the
    suffix of its element set and day: the configuration's own, if it has them.
    """
    for name in (
        f"{entry_id}_{key}",
        f"{entry_id}_{platform.replace('-', '_')}",
        entry_id,
    ):
        path = FIXTURES / f"{name}.npz"
        if path.exists():
            with np.load(path, allow_pickle=False) as data:
                for suffix in (f"_{key}", ""):
                    if f"tle{suffix}" in data and f"day{suffix}" in data:
                        return path, suffix
    return None


def get_swath(entry, key, orbit, platform, day, start, step, swept):
    """Gets the swath of a configuration over a segment at a time step."""
    satellite = entry.to_satellite(
        orbit,
        key,
        time_step=step,
        platform=platform,
        start=day,
        end=day + timedelta(days=1),
    )
    if not swept:
        satellite = satellite.model_copy(
            update={
                "instruments": [
                    instrument.model_copy(update={"is_swept": False})
                    for instrument in satellite.instruments
                ]
            }
        )
    times = [
        start + timedelta(seconds=step * (k + 0.5))
        for k in range(int(round(DURATION / step)))
    ]
    grouping = str(entry.configurations[key].grouping)
    index = None if grouping.endswith("union") else 0
    track = compute_ground_track(satellite, times, instrument_index=index)
    return None if track.empty else track.geometry.iloc[0]


def compare(a, b) -> float:
    """Gets the intersection over union of two swaths, in a local projection."""
    center = a.centroid
    forward = Transformer.from_crs(
        "EPSG:4326",
        f"+proj=aeqd +lat_0={center.y} +lon_0={center.x} +units=m",
        always_xy=True,
    ).transform
    a, b = (transform_geometry(forward, g).buffer(0) for g in (a, b))
    return a.intersection(b).area / a.union(b).area


def compare_configuration(
    entry, key: str, fixture: tuple[pathlib.Path, str]
) -> list[tuple]:
    """Gets the comparisons of a configuration's swaths over each segment."""
    platform = entry.configurations[key].platforms[0]
    path, suffix = fixture
    with np.load(path, allow_pickle=False) as data:
        orbit = GeneralPerturbationsOrbit.from_tle(list(data[f"tle{suffix}"]))
        day = datetime.fromisoformat(str(data[f"day{suffix}"]))
    results = []
    for hours in OFFSETS:
        start = day + timedelta(hours=hours)
        swaths = {
            (step, swept): get_swath(
                entry, key, orbit, platform, day, start, step, swept
            )
            for step, swept in ((2, False), (5, True), (120, True), (120, False))
        }
        views = swaths[(2, False)]
        if views is None or views.bounds[2] - views.bounds[0] > 180:
            continue
        results.append(
            (
                compare(swaths[(5, True)], swaths[(120, True)]),
                compare(views, swaths[(120, True)]),
                compare(views, swaths[(120, False)]),
            )
        )
    return results


def main() -> None:
    """Compares the swaths of each swept configuration."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("entries", nargs="*")
    arguments = parser.parse_args()
    warnings.simplefilter("ignore")
    print("configuration                    swept 120 s   swept 120 s    views 120 s")
    print("                                 vs 5 s        vs views 2 s   vs views 2 s")
    for entry_id in arguments.entries or get_entry_ids():
        entry = get_entry(entry_id)
        for key, configuration in entry.configurations.items():
            if not any(s.parameters.get("is_swept") for s in configuration.instruments):
                continue
            fixture = get_fixture(entry_id, key, configuration.platforms[0])
            if fixture is None:
                print(f"{entry_id}/{key}: no element set")
                continue
            results = compare_configuration(entry, key, fixture)
            if not results:
                print(f"{entry_id}/{key}: no segment away from the anti-meridian")
                continue
            print(
                f"{entry_id + '/' + key:32s} "
                + "  ".join(f"{value:12.4f}" for value in np.median(results, axis=0)),
                flush=True,
            )


if __name__ == "__main__":
    main()

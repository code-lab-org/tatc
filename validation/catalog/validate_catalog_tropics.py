"""
Validates the TROPICS catalog entry against the Level 1B brightness temperature
products of the TROPICS constellation's CubeSats (TROPICS0nBRTTL1B, of NASA GES DISC),
using the exact values of the entry through public TAT-C methods. Each satellite's
radiometer spins once every 2 s, sampling 81 spots across track in about 0.67 s, and
each product gives every spot's time and position (for each of five bands of
channels) and the spacecraft's Earth-fixed position by scan. The spots of the 118 GHz
band (band 3, whose positions bands 1 and 2 share) are compared.

The entry's scan angles are fitted on one day (run with `--fit`), from the spacecraft
positions the product records, so the fit needs no element set; running on a day at
least 30 days from it tests the fitted angles on data held out from the fit. Each
satellite has its own configuration, fitted and tested on its own data.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`. Only the
needed fields of the product's files are read, with range requests:

    python validation/catalog/validate_catalog_tropics.py <data directory> 2024-08-15 --fit --platform tropics-03
    python validation/catalog/validate_catalog_tropics.py <data directory> 2024-10-01 --platform tropics-03

Without a day, it takes that of the entry's record (or, with `--fit`, of the fit).
The satellites' element sets are read from the file of their history in `elements/`
(or from a file given with `--elements`). With `--write`, it writes the entry's
validation record and the reduced fixture from which unit tests recompute the
results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from datetime import datetime

import h5py
import numpy as np
import pandas as pd

from pyproj import Geod
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric

from tatc.constants import timescale
from tatc.schemas import GeneralPerturbationsOrbit
from tatc.utils import compute_projected_ray_position
from scanner import Product, Reference, Setup, fit, load_day, validate_swath
from validation_io import FIXTURES, TO_GEODETIC, RemoteFile, parse_run

SCAN_PERIOD = 2.0
GEOD = Geod(ellps="WGS84")
# the scans whose edges are fitted
FIT_SCANS = 3000
# the band whose spots are compared (band 3, 118 GHz, whose positions bands 1 and 2
# share)
BAND = 2
# the start of the TROPICS epoch time (TAI seconds since 2000), as UTC, and TAI - UTC
# since 2017
TET_EPOCH = pd.Timestamp("2000-01-01", tz="UTC")
TAI_UTC = 37


def tet_to_seconds(tet, day: datetime) -> np.ndarray:
    """Converts TROPICS epoch times to UTC seconds from the start of a day."""
    return (
        np.asarray(tet, dtype=float)
        - TAI_UTC
        - (pd.Timestamp(day) - TET_EPOCH).total_seconds()
    )


def read_tropics(_product, url: str, session, day: datetime, _directory) -> dict:
    """
    Reads the fields of a granule with range requests: the middle spot's time of each
    scan and each spot's time (seconds from the start of the day), the spots'
    positions, and the spacecraft's geodetic position.
    """
    with h5py.File(RemoteFile(url, session)) as f:
        beam_time = tet_to_seconds(f["time"][:], day)
        position = f["scPosECEF"][:].astype(float) * 1e3
        lon, lat, alt = TO_GEODETIC.transform(*position)
        return {
            "time": beam_time[:, beam_time.shape[1] // 2],
            "beam_time": beam_time,
            "lat": f["latitude"][BAND],
            "lon": f["longitude"][BAND],
            "sc_lat": np.asarray(lat),
            "sc_lon": np.asarray(lon),
            "sc_alt": np.asarray(alt),
        }


def fit_edges(reference: Reference, iterations: int = 6) -> None:
    """
    Prints the view whose edges TAT-C projects onto the outermost spots: for each
    edge, the roll (deg, positive to the left of the inertial velocity) of the ray
    that `compute_projected_ray_position` projects, from the spacecraft's position
    and velocity of the product, at the median distance across track of zero from
    the spot, by secant steps; then the field of view between the edges and its
    roll. Angles measured from the geodetic nadir (`fit`) do not reproduce TAT-C's
    projection of a view this wide. Scans flown with the scan reversed are left out.
    """
    angles = reference.cross_track_angles("inertial")
    usual = np.sign(np.median(angles[:, -1] - angles[:, 0]))
    contiguous = (
        np.r_[False, reference.is_next_contiguous]
        & np.r_[reference.is_next_contiguous, False]
        & (np.sign(angles[:, -1] - angles[:, 0]) == usual)
    )
    candidates = np.flatnonzero(contiguous)
    pick = candidates[np.linspace(0, len(candidates) - 1, FIT_SCANS).astype(int)]
    edges = {}
    for spot in (0, reference.beams - 1):
        seconds = reference.beam_seconds[pick, spot]
        t = timescale.from_datetimes(
            [reference.day + pd.Timedelta(seconds=float(s)) for s in seconds]
        )
        # the spacecraft's position and velocity at the spot's time
        position = np.array(
            [
                np.interp(seconds, reference.seconds, reference.position[:, i])
                for i in range(3)
            ]
        ).T
        velocity = np.array(
            [
                np.interp(
                    seconds, reference.seconds, reference.velocity("inertial")[:, i]
                )
                for i in range(3)
            ]
        ).T
        rotation = itrs.rotation_at(t)
        au_m = 149597870700.0
        track = Geocentric(
            np.einsum("jin,nj->in", rotation, position) / au_m,
            np.einsum("jin,nj->in", rotation, velocity) / au_m * 86400,
            t,
            center=399,
        )
        heading = GEOD.inv(
            reference.lon[pick, reference.beams // 2] * 0
            + TO_GEODETIC.transform(*position.T)[0],
            TO_GEODETIC.transform(*position.T)[1],
            TO_GEODETIC.transform(*(position + velocity).T)[0],
            TO_GEODETIC.transform(*(position + velocity).T)[1],
        )[0]

        def across(angle):
            ray = compute_projected_ray_position(
                track, 0, 0, roll_angle=angle, velocity_frame="inertial"
            )
            azimuth, _, distance = GEOD.inv(
                ray.longitude.degrees,
                ray.latitude.degrees,
                reference.lon[pick, spot],
                reference.lat[pick, spot],
            )
            # positive to the left of the motion
            return -float(np.median(distance * np.sin(np.radians(azimuth - heading))))

        a0 = float(np.median(angles[pick, spot]))
        a1 = a0 + 0.1
        e0, e1 = across(a0), across(a1)
        for _ in range(iterations):
            if e1 == e0:
                break
            a0, a1, e0 = a1, a1 - e1 * (a1 - a0) / (e1 - e0), e1
            e1 = across(a1)
        edges[spot] = a1
        print(f"spot {spot}: roll {a1:.4f} deg, median miss across {e1:.3f} m")
    left, right = max(edges.values()), min(edges.values())
    print(
        f"TAT-C view: cross-track field of view {left - right:.4f} deg, roll "
        f"{(left + right) / 2:.4f} deg ({len(pick)} scans)"
    )


def product(number: str) -> Product:
    """Gets the Level 1B product of a satellite of the constellation."""
    return Product(
        short_name=f"TROPICS{number}BRTTL1B",
        version="0.2" if number == "05" else "1.0",
        read=read_tropics,
    )


PLATFORMS = {f"tropics-{n}": n for n in ("03", "05", "06")}


def setup(platform: str) -> Setup:
    """Gets what is validated on a satellite."""
    number = PLATFORMS[platform]
    name = product(number).short_name
    return Setup(
        swath=f"swath_tropics_{number}",
        beams=None,
        size_instrument=None,
        published_footprints=(),
        fixture_path=FIXTURES / f"tropics_{number}.npz",
        reference=f"{name} (118 GHz spot positions and times)",
        published_reference="",
        time_steps=(SCAN_PERIOD, 32, 128),
        pass_radius=60,
        point_spacing=500e3,
    )


def main() -> None:
    """Validates the entry on a day of data, or fits its angles."""
    run = parse_run("tropics", __doc__)
    number = PLATFORMS[run.platform]
    reference = Reference(
        run.day, load_day(product(number), run.day, run.directory), SCAN_PERIOD
    )
    print(
        f"{run.day:%Y-%m-%d} ({run.platform}): {len(reference.time)} scans from "
        f"{reference.time[0]}, {len(reference.gaps)} gaps over a minute",
        flush=True,
    )
    if run.fit:
        for frame in ("inertial", "earth_fixed"):
            fit(reference, frame)
        fit_edges(reference)
        offsets = reference.beam_seconds - reference.seconds[:, np.newaxis]
        print(f"spot times from the middle spot: {np.median(offsets, axis=0)[[0, -1]]}")
        return
    tle = run.element_set()
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    validate_swath(
        setup(run.platform),
        run,
        reference,
        orbit,
        tle,
        fixture_extras=None,
    )


if __name__ == "__main__":
    main()

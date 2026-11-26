"""
Builds the table of global-mean thermospheric densities of
`tatc.utils.compute_atmospheric_density`
(`src/tatc/resources/atmosphere/nrlmsis_density.csv`) from the NRLMSIS 2.1
empirical model (Emmert et al., 2022, https://doi.org/10.1029/2022SW003219),
as implemented by pymsis (https://swxtrec.github.io/pymsis/).

The table gives the mass density (kg/m^3) at heights from 100 to 1500 km
above the Earth's mean radius (the altitude of `tatc.schemas` orbits: the
semimajor axis minus `tatc.constants.EARTH_MEAN_RADIUS`), every 5 km. Each
value is the mean over the sphere at that height, weighted by area, over local
times (24 longitudes at 00 and 12 UT) and seasons (the 15th day of each month),
for the three long-term reference levels of solar and geomagnetic activity of
ECSS-E-ST-10-04C (Space environment, 15 November 2008), Table 6-3: low
(F10.7 = 65, Ap = 0), moderate (F10.7 = 140, Ap = 15), and high (F10.7 = 250,
Ap = 45), with the 81-day mean of F10.7 equal to its daily value, which
clause 6.2.2 b of the standard says to use for periods longer than 27 days.

It reports the errors of interpolating the table log-linearly in height (at
the heights halfway between its rows), and the spread of the local density
about the global mean (which an orbit averages over its revolutions).

Usage:

    python validation/atmosphere/build_density_table.py

Requires `pymsis` (`pip install tatc[validation]`).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import os
from datetime import datetime, timezone

import numpy as np
import pymsis

# WGS 84, as tatc.constants
EQUATORIAL_RADIUS = 6378137.0  # m
POLAR_RADIUS = 6356752.314245179  # m
MEAN_RADIUS = (2 * EQUATORIAL_RADIUS + POLAR_RADIUS) / 3  # m

VERSION = 2.1
LEVELS = {  # ECSS-E-ST-10-04C, Table 6-3, long-term values
    "low": (65.0, 0.0),
    "moderate": (140.0, 15.0),
    "high": (250.0, 45.0),
}
HEIGHTS = np.arange(100, 1500 + 1, 2.5) * 1e3  # m, every 2.5 km (5 km written)
LATITUDES = np.arange(-87.5, 90, 5.0)  # deg
LONGITUDES = np.arange(0, 360, 15.0)  # deg
DATES = np.array(
    [
        np.datetime64(datetime(2021, month, 15, hour), "s")
        for month in range(1, 13)
        for hour in (0, 12)
    ]
)
OUTPUT = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "src",
    "tatc",
    "resources",
    "atmosphere",
    "nrlmsis_density.csv",
)


def get_surface_radius(latitude):
    """
    Gets the distance (m) from the Earth's center to the WGS 84 ellipsoid at a
    latitude (degrees; geodetic and geocentric latitudes differ by at most 0.2
    degrees, which changes the radius by less than 25 m).
    """
    phi = np.radians(latitude)
    a, b = EQUATORIAL_RADIUS, POLAR_RADIUS
    return np.sqrt(
        ((a**2 * np.cos(phi)) ** 2 + (b**2 * np.sin(phi)) ** 2)
        / ((a * np.cos(phi)) ** 2 + (b * np.sin(phi)) ** 2)
    )


def compute_local_densities(f107, ap):
    """
    Computes the densities (kg/m^3) at every date, longitude, latitude, and
    height above the mean radius: an array (date, longitude, latitude, height).
    """
    n_dates = len(DATES)
    densities = np.empty((n_dates, len(LONGITUDES), len(LATITUDES), len(HEIGHTS)))
    for k, latitude in enumerate(LATITUDES):
        # heights above the ellipsoid of the points at heights above the
        # mean radius
        altitudes = (MEAN_RADIUS + HEIGHTS - get_surface_radius(latitude)) / 1e3
        output = pymsis.calculate(
            DATES,
            LONGITUDES,
            [latitude],
            altitudes,
            f107s=np.full(n_dates, f107),
            f107as=np.full(n_dates, f107),
            aps=np.full((n_dates, 7), ap),
            version=VERSION,
        )
        densities[:, :, k, :] = output[:, :, 0, :, pymsis.Variable.MASS_DENSITY]
    return densities


def main():
    """
    Builds the table and reports its errors.
    """
    weights = np.cos(np.radians(LATITUDES))
    weights = weights / weights.sum()
    means = {}
    lines = []
    for level, (f107, ap) in LEVELS.items():
        local = compute_local_densities(f107, ap)
        # mean over dates and longitudes, then over latitudes by area
        mean = np.einsum("dolh,l->h", local, weights) / (
            local.shape[0] * local.shape[1]
        )
        means[level] = mean
        # errors of log-linear interpolation of the 5 km rows at 2.5 km
        rows, middles = mean[::2], mean[1::2]
        interpolated = np.sqrt(rows[:-1] * rows[1:])
        error = np.abs(interpolated / middles - 1)
        lines.append(
            f"{level}: log-linear interpolation error at 2.5 km offsets: "
            + f"max {error.max():.2%}, mean {error.mean():.2%}"
        )
        for height in (400e3, 800e3):
            j = int(np.argmin(np.abs(HEIGHTS - height)))
            ratio = (local[..., j] / mean[j]).ravel()
            ratio_weights = np.broadcast_to(
                weights[None, None, :], local.shape[:3]
            ).ravel()
            order = np.argsort(ratio)
            cumulative = np.cumsum(ratio_weights[order]) / ratio_weights.sum()
            p05, p95 = np.interp([0.05, 0.95], cumulative, ratio[order])
            lines.append(
                f"{level}: local/global-mean density at {height / 1e3:.0f} km: "
                + f"5th percentile {p05:.2f}, 95th percentile {p95:.2f}"
            )
    for line in lines:
        print(line)

    os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
    with open(OUTPUT, "w", encoding="utf-8", newline="\n") as file:
        file.write(
            "# Global-mean mass density (kg/m^3) of NRLMSIS 2.1 (Emmert et al., 2022,\n"
            + "# https://doi.org/10.1029/2022SW003219; pymsis "
            + f"{pymsis.__version__}) at heights (km) above\n"
            + "# the Earth's mean radius, averaged over the sphere by area, local\n"
            + "# times, and seasons, for the long-term reference levels of solar and\n"
            + "# geomagnetic activity of ECSS-E-ST-10-04C, Table 6-3: low (F10.7 = 65,\n"
            + "# Ap = 0), moderate (F10.7 = 140, Ap = 15), high (F10.7 = 250, Ap = 45),\n"
            + "# with the 81-day mean F10.7 equal to the daily value.\n"
            + "# Built by validation/atmosphere/build_density_table.py on "
            + f"{datetime.now(timezone.utc):%Y-%m-%d}:\n"
        )
        for line in lines:
            file.write(f"#   {line}\n")
        file.write("altitude_km,low,moderate,high\n")
        for i in range(0, len(HEIGHTS), 2):
            values = ",".join(f"{means[level][i]:.5e}" for level in LEVELS)
            file.write(f"{HEIGHTS[i] / 1e3:.0f},{values}\n")
    print(f"wrote {os.path.normpath(OUTPUT)}")


if __name__ == "__main__":
    main()

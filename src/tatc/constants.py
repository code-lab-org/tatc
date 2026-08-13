"""
Numerical constants.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import os

import numpy as np
from skyfield.api import Loader, load

# load ephemeris file
resources_dir = os.path.join(os.path.dirname(__file__), "resources")
de421_loader = Loader(resources_dir)
de421 = de421_loader("de421.bsp")

# load timescale
timescale = load.timescale()

# time properties
EARTH_SOLAR_DAY_S = 86400
EARTH_SIDEREAL_DAY_S = 86164.0905

# wgs84 oblate spheroid parameters
EARTH_FLATTENING = 1 / 298.257223563
EARTH_EQUATORIAL_RADIUS = 6378137.0
EARTH_POLAR_RADIUS = 6356752.314245179
EARTH_MU = 3.986004418e14
# derived from the WGS84/EGM96 fully-normalized zonal coefficient C_bar_20 =
# -0.484165143790e-3 (NIMA TR8350.2, "Department of Defense World Geodetic
# System 1984"), converted to unnormalized form via J2 = -C_bar_20 * sqrt(5)
EARTH_J2 = 1.0826261738504e-3

# derived properties based on wgs84 oblate spheroid
EARTH_ECCENTRICITY = np.sqrt(2 * EARTH_FLATTENING - EARTH_FLATTENING**2)
EARTH_EQUATORIAL_CIRCUMFERENCE = 2 * np.pi * EARTH_EQUATORIAL_RADIUS
EARTH_POLAR_CIRCUMFERENCE = 2 * np.pi * EARTH_POLAR_RADIUS
EARTH_SURFACE_AREA = (
    2 * np.pi * EARTH_EQUATORIAL_RADIUS**2
    + np.pi
    * EARTH_POLAR_RADIUS**2
    / EARTH_ECCENTRICITY
    * np.log((1 + EARTH_ECCENTRICITY) / (1 - EARTH_ECCENTRICITY))
)
EARTH_MEAN_RADIUS = (2 * EARTH_EQUATORIAL_RADIUS + EARTH_POLAR_RADIUS) / 3
EARTH_J2_CRITICAL_INCLINATION = np.degrees(np.arccos(np.sqrt(1 / 5)))

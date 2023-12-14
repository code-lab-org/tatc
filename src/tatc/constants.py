# -*- coding: utf-8 -*-
"""
Numerical constants.

@author: Paul T. Grogan <pgrogan@stevens.edu>
"""

import numpy as np
from skyfield.api import load


# load ephemeris file
de421 = load("de421.bsp")

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

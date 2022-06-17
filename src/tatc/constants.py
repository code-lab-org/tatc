# -*- coding: utf-8 -*-
"""
Numerical constants.

@author: Paul T. Grogan <pgrogan@stevens.edu>
"""

import numpy as np
from skyfield.api import load


# load ephemeris file
de421 = load("de421.bsp")

# wgs84 oblate spheroid parameters
earth_flattening = 1 / 298.257223563
earth_equatorial_radius = 6378137.0
earth_polar_radius = 6356752.314245179
earth_mu = 3.986004418e14

# derived properties based on wgs84 oblate spheroid
earth_eccentricity = np.sqrt(2 * earth_flattening - earth_flattening ** 2)
earth_surface_area = (
    2 * np.pi * earth_equatorial_radius ** 2
    + np.pi
    * earth_polar_radius ** 2
    / earth_eccentricity
    * np.log((1 + earth_eccentricity) / (1 - earth_eccentricity))
)
earth_mean_radius = (2 * earth_equatorial_radius + earth_polar_radius) / 3

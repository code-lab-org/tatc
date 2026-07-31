"""
Surface utility functions.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import numpy as np
from numba import njit

from .. import constants


@njit
def compute_number_samples(distance: float) -> int:
    """
    Compute the number of global samples required to achieve a typical
    sample distance (meters) assuming equal spacing, by dividing the
    Earth's surface area (assuming a mean sphere) by the area of a
    spherical cap whose angular radius corresponds to half the sample
    distance. The result is truncated (not rounded), so the achieved
    average sample spacing is always at least the requested distance,
    never less. Requires distance > 0.

    Args:
        distance (float): The typical distance between samples (meters).

    Returns:
        int: The number of global samples.
    """
    # compute the angular distance of each sample (assuming mean sphere)
    theta = distance / constants.EARTH_MEAN_RADIUS
    # compute the distance from the center of earth to conic plane (assuming sphere)
    radius = constants.EARTH_MEAN_RADIUS * np.cos(theta / 2)
    # compute the distance from the conic plane to the surface (assuming sphere)
    height = constants.EARTH_MEAN_RADIUS - radius
    # compute the sperical cap area covered by the sample (assuming sphere)
    # https://en.wikipedia.org/wiki/Spherical_cap
    sample_area = 2 * np.pi * constants.EARTH_MEAN_RADIUS * height
    # return the fraction of earth-to-sample area
    return int(constants.EARTH_SURFACE_AREA / sample_area)

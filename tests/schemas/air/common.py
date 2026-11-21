"""
Trajectories shared by the tests of aircraft.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import numpy as np
import pandas as pd
from pyproj import Geod

from tatc.constants import EARTH_ECCENTRICITY, EARTH_EQUATORIAL_RADIUS

GEOD = Geod(ellps="WGS84")
START = pd.Timestamp("2026-03-01T18:00:00Z")


def surface_distance(distance, latitude, azimuth, altitude):
    """
    Gets the distance (meters) over the WGS 84 ellipsoid below a distance
    flown horizontally at an altitude (meters), toward an azimuth (degrees)
    from a latitude (degrees): scaled by the radius of curvature of the
    ellipsoid in that azimuth.
    """
    lat, azi = np.radians(latitude), np.radians(azimuth)
    e2 = EARTH_ECCENTRICITY**2
    w = 1 - e2 * np.sin(lat) ** 2
    meridian = EARTH_EQUATORIAL_RADIUS * (1 - e2) / w**1.5
    prime = EARTH_EQUATORIAL_RADIUS / np.sqrt(w)
    radius = 1 / (np.cos(azi) ** 2 / meridian + np.sin(azi) ** 2 / prime)
    return distance * radius / (radius + altitude)


def straight_trajectory(
    n=60,
    altitude=10000.0,
    speed=200.0,
    track=30.0,
    step=1.0,
    latitude=35.0,
    longitude=-110.0,
    **attitude,
):
    """
    Gets the states of an aircraft flying along a geodesic at a constant
    ground speed at its altitude (a data frame), with optional constant
    `heading`, `roll`, and `pitch` angles (degrees). The track is the
    geodesic's azimuth at each state.
    """
    distance = surface_distance(speed * step * np.arange(n), latitude, track, altitude)
    lon, lat, back = GEOD.fwd(
        np.full(n, longitude), np.full(n, latitude), np.full(n, track), distance
    )
    data = pd.DataFrame(
        {
            "time": START + pd.to_timedelta(step * np.arange(n), unit="s"),
            "latitude": lat,
            "longitude": lon,
            "altitude": altitude,
            "track": np.mod(np.asarray(back) + 180, 360),
            "ground_speed": speed,
        }
    )
    for name, value in attitude.items():
        data[name] = value
    return data

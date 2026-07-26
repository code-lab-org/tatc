"""
Methods to analyze dilusion of precision.

@author: Michael P. Jones <mpj@mit.edu>
@author: Paul T. Grogan <paul.grogan@asu.edu>
"""
from __future__ import annotations

import warnings
from datetime import datetime
from enum import Enum

import geopandas as gpd
import numpy as np
import pandas as pd
from skyfield.api import EarthSatellite, wgs84

from ..constants import timescale
from ..schemas import Point, Satellite


class DopMethod(str, Enum):
    """
    Enumeration of Dilusion of Precision (DOP) calculation methods.
    """

    GDOP = "gdop"  # Geometric Dilusion of Precision
    PDOP = "pdop"  # Position (3D) Dilusion of Precision
    HDOP = "hdop"  # Horizontal Dilusion of Precision
    VDOP = "vdop"  # Vertical Dilusion of Precision
    TDOP = "tdop"  # Time Dilusion of Precision


def compute_dop(
    times: list[datetime],
    point: Point,
    satellites: list[Satellite],
    min_elevation: float,
    dop_method: DopMethod,
    min_count_visible: int = 3,
) -> gpd.GeoDataFrame:
    """
    Calculate the specified dilusion of precision value based on inputs.

    Args:
        times: a vector of datetimes at which to measure dilusion of precision
        point: a ground point intended to view satellites
        satellites: the list of satellites to be viewed by the ground point
        min_elevation: the minimum elevation angle (deg) to consider a satellite visible
        dop_method: dilusion of precision calculation method
        min_count_visible: minimum number of visible satellites for a valid measurement

    Outputs:
    - geopandas.GeoDataFrame: the dop for the given user location and satellite.

    """
    # construct skyfield satellites for each satellite
    sk_sats = [
        satellite.orbit.to_gp_orbit().get_closest_element(times[0]).to_skyfield()
        for satellite in satellites
    ]

    # construct skyfield times for each datetime
    sk_times = timescale.from_datetimes(times)

    # construct skyfield geodetic position for user
    sk_position = wgs84.latlon(point.latitude, point.longitude)

    # compute elevation/azimuth angles and range
    altazs = [(sk_sat - sk_position).at(sk_times).altaz() for sk_sat in sk_sats]
    el = np.array([i[0].radians for i in altazs])
    az = np.array([i[1].radians for i in altazs])
    r = np.array([i[2].m for i in altazs])

    # compute number of visible satellites
    n = np.sum(el >= np.deg2rad(min_elevation), axis=0)

    # compute position of visible satellites (in local frame)
    x = r * np.cos(el) * np.cos(az)
    y = r * np.cos(el) * np.sin(az)
    z = r * np.sin(el)

    def _dop(i):
        """
        Compute the dilution of precision value for time index i.
        """
        if n[i] <= min_count_visible:
            return np.nan
        mask = el[:, i] >= np.deg2rad(min_elevation)
        # H is a nx4 matrix where n is the number of visible satellites
        H = np.column_stack(
            (
                x[mask, i] / r[mask, i],
                y[mask, i] / r[mask, i],
                z[mask, i] / r[mask, i],
                np.ones((sum(mask), 1)),
            )
        )
        # calculate the pseudoinverse of H
        try:
            H_inv = np.linalg.solve(H.T @ H, np.eye(4))
        except np.linalg.LinAlgError:
            # If the H matrix nearly singular, linalg will not be able to solve, return NaN and a warning
            warnings.warn("H matrix could not be inverted, NaN DOP value returned.")
            return np.nan
        # compute and return the dop value
        if dop_method == DopMethod.GDOP:
            return np.sqrt(np.trace(H_inv))
        if dop_method == DopMethod.PDOP:
            return np.sqrt(H_inv[0, 0] + H_inv[1, 1] + H_inv[2, 2])
        if dop_method == DopMethod.HDOP:
            return np.sqrt(H_inv[0, 0] + H_inv[1, 1])
        if dop_method == DopMethod.VDOP:
            return np.sqrt(H_inv[2, 2])
        if dop_method == DopMethod.TDOP:
            return np.sqrt(H_inv[3, 3])
        raise ValueError("Invalid DOP method")

    dop = np.array([_dop(i) for i in range(len(times))])

    columns = {
        "dop": pd.Series(dop, dtype="float", index=times),
        "geometry": pd.Series(
            gpd.points_from_xy(
                [point.latitude] * len(dop), [point.longitude] * len(dop)
            ),
            dtype="object",
            index=times,
        ),
    }
    dop_df = gpd.GeoDataFrame(columns, crs="EPSG:4326")

    return dop_df

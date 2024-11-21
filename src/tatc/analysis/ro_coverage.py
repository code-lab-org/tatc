# -*- coding: utf-8 -*-
"""
Methods to perform radio occultation (RO) coverage analysis.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from typing import List, Tuple, Union
from datetime import datetime
from itertools import chain

import geopandas as gpd
import numpy as np
from shapely.geometry import MultiPoint, Point
from skyfield.api import EarthSatellite, wgs84, Distance, Velocity, Time
from skyfield.positionlib import ICRF

from ..schemas.satellite import Satellite

from ..constants import timescale


def _collect_ro_series(
    transmitter: Satellite,
    ts: Time,
    rx_pv: ICRF,
    rx_n_u: List[float],
    rx_t_u: List[float],
    max_yaw: float,
    range_elevation: Tuple[float],
):
    # construct transmitter
    tx_tle = transmitter.orbit.to_tle()
    tx = EarthSatellite(tx_tle.tle[0], tx_tle.tle[1], transmitter.name)
    # transmitter position (x_tx), velocity (v_tx)
    tx_pv = tx.at(ts)
    # relative position, velocity of transmitter from receiver
    # x_(rx,tx) = x_tx - x_rx; v_(rx,tx) = v_tx - v_rx
    rx_tx_pv = tx_pv - rx_pv
    # tangent point position (m)
    # x_tp = x_tx - x_(rx,tx) . [ x_tx . x_(rx,tx) ] / || x_(rx,tx) ||
    tp_p = tx_pv.position.m - np.einsum(
        "ij,j->ij",
        rx_tx_pv.position.m,
        np.divide(
            np.einsum("ij,ij->j", tx_pv.position.m, rx_tx_pv.position.m),
            np.einsum("ij,ij->j", rx_tx_pv.position.m, rx_tx_pv.position.m),
        ),
    )
    # tangent point velocity (m/s) - derived using chain rule
    # v_tp = v_tx - v_(rx,tx) . [ x_tx . x_(rx,tx) ] / || x_(rx,tx) ||
    #        - x_(rx,tx) . [
    #           [ v_tx . x_(rx,tx) ] + [ x_tx . v_(rx,tx) ] ] / || x_(rx,tx) || ]
    #           - 2 * [ v_(rx,tx) . x_(rx,tx) ] * [ x_tx . x_(rx,tx) ] / || x_(rx,tx) ||^2
    #        ]
    tp_v = (
        tx_pv.velocity.m_per_s
        - np.einsum(
            "ij,j->ij",
            rx_tx_pv.velocity.m_per_s,
            np.divide(
                np.einsum("ij,ij->j", tx_pv.position.m, rx_tx_pv.position.m),
                np.einsum("ij,ij->j", rx_tx_pv.position.m, rx_tx_pv.position.m),
            ),
        )
        - np.einsum(
            "ij,j->ij",
            rx_tx_pv.position.m,
            (
                np.divide(
                    (
                        np.einsum(
                            "ij,ij->j", tx_pv.velocity.m_per_s, rx_tx_pv.position.m
                        )
                        + np.einsum(
                            "ij,ij->j", tx_pv.position.m, rx_tx_pv.velocity.m_per_s
                        )
                    ),
                    np.einsum("ij,ij->j", rx_tx_pv.position.m, rx_tx_pv.position.m),
                )
                - 2
                * np.divide(
                    np.multiply(
                        np.einsum(
                            "ij,ij->j", rx_tx_pv.velocity.m_per_s, rx_tx_pv.position.m
                        ),
                        np.einsum("ij,ij->j", tx_pv.position.m, rx_tx_pv.position.m),
                    ),
                    np.power(
                        np.einsum("ij,ij->j", rx_tx_pv.position.m, rx_tx_pv.position.m),
                        2,
                    ),
                )
            ),
        )
    )
    # intersecting (-1) or parallel (+1) view of tangent point
    tp_sign = np.sign(
        np.einsum("ij,ij->j", tp_p - tx_pv.position.m, tp_p - rx_pv.position.m)
    )
    # relative transmitter position from receiver in plane normal to receiver orbit
    rx_tx_p_rx_n_plane = rx_tx_pv.position.m - np.einsum(
        "ij,j->ij", rx_n_u, np.einsum("ij,ij->j", rx_n_u, rx_tx_pv.position.m)
    )
    # relative transmitter position from receiver in plane tangent to receiver orbit
    rx_tx_p_rx_t_plane = rx_tx_pv.position.m - np.einsum(
        "ij,j->ij", rx_t_u, np.einsum("ij,ij->j", rx_t_u, rx_tx_pv.position.m)
    )
    # transmitter pitch angle in receiver body-fixed frame
    rx_tx_pitch = np.degrees(
        np.arccos(
            np.divide(
                np.einsum("ij,ij->j", rx_tx_p_rx_n_plane, rx_pv.velocity.m_per_s),
                np.multiply(
                    np.linalg.norm(rx_tx_p_rx_n_plane, axis=0),
                    np.linalg.norm(rx_pv.velocity.m_per_s, axis=0),
                ),
            )
        )
    )
    # transmitter yaw angle in receiver body-fixed frame
    rx_tx_yaw = np.degrees(
        np.arccos(
            np.divide(
                np.einsum("ij,ij->j", rx_tx_p_rx_t_plane, rx_pv.velocity.m_per_s),
                np.multiply(
                    np.linalg.norm(rx_tx_p_rx_t_plane, axis=0),
                    np.linalg.norm(rx_pv.velocity.m_per_s, axis=0),
                ),
            )
        )
    )
    # occultation observations
    occ_obs = []
    # occultation arc
    occ_arc = None
    # valid if tangent point intersects and transmitter view angle below maximum
    valid = np.logical_and(
        tp_sign < 0, rx_tx_yaw % (180 - max_yaw) < max_yaw
    )
    # events occur when validity changes value
    is_event = np.diff(valid)
    # loop over valid times
    for j in np.nonzero(valid)[0]:
        # tangent point inertial position
        tpp_icrf = ICRF(
            Distance(m=tp_p[:, j]).au,
            Velocity(km_per_s=tp_v[:, j] / 1000).au_per_d,
            ts[j],
            399,
        )
        # tangent point geodetic position
        tpp_geo = wgs84.geographic_position_of(tpp_icrf)
        # check if the tangent point height is within elevation range
        if (
            tpp_geo.elevation.m > range_elevation[0]
            and tpp_geo.elevation.m < range_elevation[1]
        ):
            if occ_arc is None:
                # start of new RO observation
                occ_arc = {
                    "tx": tx.name,
                    "is_rising": rx_tx_pitch[j] < 90,
                    "points": [],
                }

            # azimuth of transmitter from geodetic tangent point (clockwise from North)
            tp_tx_azmimuth = (tx - tpp_geo).at(ts[j]).altaz()[1].degrees

            occ_arc["points"].append(
                {
                    "time": ts[j].utc_datetime(),
                    "tangent_point": tpp_geo,
                    "rx_tx_pitch": rx_tx_pitch[j],
                    "rx_tx_yaw": rx_tx_yaw[j],
                    "tp_tx_azimuth": tp_tx_azmimuth,
                }
            )
            if j + 1 >= len(ts) or is_event[j]:
                # end of RO observation due to validity or boundary constraint
                occ_obs.append(occ_arc)
                occ_arc = None
        elif occ_arc is not None:
            # end of RO observation due to elevation constraints
            occ_obs.append(occ_arc)
            occ_arc = None
    return occ_obs


def collect_ro_observations(
    receiver: Satellite,
    transmitters: Union[Satellite, List[Satellite]],
    times: List[datetime],
    sample_elevation: float = 0,
    max_yaw: float = 65,
    range_elevation: Tuple[float] = (-200e3, 60e3),
) -> gpd.GeoDataFrame:
    """
    Collects Radio Occultation (RO) observations.

    Args:
        receiver (Satellite): the satellite with a RO receiver.
        transmitters (typing.Union[Satellite,typing.List[Satellite]]): the satellite(s) with a RO transmitter.
        times (typing.List[datetime.datetime]): The list of datetimes to sample.
        sample_elevation: (float): the elevation (m) at which to sample observation attributes.
        max_yaw (float): the maximum transmitter yaw angle (from receiver body-fixed frame) for a valid obsevation.
        range_elevation: (typing.Tuple[float]): the lower and upper bound on tangent point elevation (m) for a valid observation.
    """
    # construct receiver
    rx_tle = receiver.orbit.to_tle()
    rx = EarthSatellite(rx_tle.tle[0], rx_tle.tle[1], receiver.name)
    # define time steps
    ts = timescale.from_datetimes(times)
    # receiver position, velocity
    rx_pv = rx.at(ts)
    # unit vector normal to receiver orbit plane
    rx_n_u = np.einsum(
        "iik->ik",
        np.cross(rx_pv.position.m.T[:, None, :], rx_pv.velocity.m_per_s.T[None, :, :]),
    ).T
    rx_n_u = np.divide(rx_n_u, np.linalg.norm(rx_n_u, axis=0))
    # unit vector tangent to receiver orbit plane
    rx_t_u = np.divide(rx_pv.position.m, np.linalg.norm(rx_pv.position.m, axis=0))
    # generate observations
    obs = list(
        chain.from_iterable(
            [
                _collect_ro_series(
                    transmitter, ts, rx_pv, rx_n_u, rx_t_u, max_yaw, range_elevation
                )
                for transmitter in (
                    transmitters if isinstance(transmitters, list) else [transmitters]
                )
            ]
        )
    )
    # format results
    return gpd.GeoDataFrame(
        [
            {
                "receiver": receiver.name,
                "transmitter": o["tx"],
                "is_rising": o["is_rising"],
                "geometry": MultiPoint(
                    [
                        [
                            point["tangent_point"].longitude.degrees,
                            point["tangent_point"].latitude.degrees,
                            point["tangent_point"].elevation.m,
                        ]
                        for point in o["points"]
                    ]
                ),
                "position": Point(
                    o["points"][sample_index]["tangent_point"].longitude.degrees,
                    o["points"][sample_index]["tangent_point"].latitude.degrees,
                    o["points"][sample_index]["tangent_point"].elevation.m,
                ),
                "rx_tx_pitch": o["points"][sample_index]["rx_tx_pitch"],
                "rx_tx_yaw": o["points"][sample_index]["rx_tx_yaw"],
                "tp_tx_azimuth": o["points"][sample_index]["tp_tx_azimuth"],
                "start": o["points"][0]["time"],
                "end": o["points"][-1]["time"],
                "time": o["points"][sample_index]["time"],
            }
            for o in obs
            for sample_index in [
                min(
                    range(len(o["points"])),
                    key=lambda i: abs(
                        o["points"][i]["tangent_point"].elevation.m - sample_elevation
                    ),
                )
            ]
        ]
    ).sort_values("time", ignore_index=True)

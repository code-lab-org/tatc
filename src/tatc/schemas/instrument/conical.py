"""
Object schemas for conically scanning instruments.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt
from pydantic import Field, model_validator
from shapely import MultiPolygon, Polygon, unary_union
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric
from skyfield.toposlib import GeographicPosition

from ... import config, constants
from ...utils.geometry import _build_split_polygons, project_polygon_to_elevation
from ...utils.projection import (
    VelocityFrame,
    _compute_projected_rays,
    _shift_along_track,
    compute_cone_and_azimuth,
    compute_ground_sample_size,
    compute_projected_ray_position,
)
from .simple import Instrument


def _cone_ray(cone_angle: float, azimuth: float) -> tuple[float, float]:
    """
    Gets the roll and pitch angles (degrees) of a rigidly rotated view
    centered on a ray at a cone angle and scan azimuth (degrees).

    Args:
        cone_angle (float): The angle (degrees) of the ray from nadir.
        azimuth (float): The scan azimuth (degrees) of the ray, from the
            along-track direction, positive to the left.

    Returns:
        tuple[float, float]: The roll and pitch angles (degrees).
    """
    cone, azimuth = np.radians(cone_angle), np.radians(azimuth)
    return (
        float(np.degrees(np.arctan2(np.sin(cone) * np.sin(azimuth), np.cos(cone)))),
        float(np.degrees(np.arcsin(np.sin(cone) * np.cos(azimuth)))),
    )


class ConicalInstrument(Instrument):
    """
    Remote sensing instrument that scans a cone about the nadir: its view
    sweeps rays at a constant angle from nadir (the cone angle) through a
    sector of scan azimuths, tracing an arc on the ground (for example, ahead
    of the satellite). A point is observed when it crosses the cone within
    the scan sector. The instantaneous footprint is the region the arc sweeps
    as the satellite moves along track by the distance subtended at nadir by
    the along-track field of view, which stands in for the motion of the arc
    over an integration time.
    """

    cone_angle: float = Field(
        ...,
        description="Angle (degrees) of the scanned cone from nadir.",
        gt=0,
        lt=90,
    )
    along_track_field_of_view: float = Field(
        ...,
        description="Angle (degrees) subtended at nadir by the along-track distance "
        + "swept by the instantaneous footprint (standing in for the motion of the "
        + "arc over an integration time).",
        gt=0,
        lt=180,
    )
    scan_center_azimuth: float = Field(
        default=0,
        description="Scan azimuth (degrees) at the center of the scan sector, "
        + "from the along-track direction, positive to the left of the direction "
        + "of motion (0: forward; 180: aft).",
        ge=-180,
        le=180,
    )
    scan_half_width: float = Field(
        default=180,
        description="Half width (degrees) of the scan sector on either side of "
        + "its center (180: a full rotation).",
        gt=0,
        le=180,
    )
    beam_width: float | None = Field(
        default=None,
        description="Angular width (degrees) of the beam (for example, its "
        + "half-power width), which sets the size of a sample on the ground "
        + "(see `compute_sample_size`), or `None` if not modeled.",
        gt=0,
        lt=90,
    )
    scan_sample_spacing: float | None = Field(
        default=None,
        description="Scan azimuth (degrees) between the centers of neighboring "
        + "samples along the scanned arc (see `compute_sample_spacing`), or "
        + "`None` if not modeled.",
        gt=0,
        le=180,
    )
    velocity_frame: VelocityFrame = Field(
        default=VelocityFrame.EARTH_FIXED,
        description="Reference frame of the velocity vector that defines the "
        + "along-track direction (`earth_fixed`: aligned with the ground track, "
        + "as for a yaw-steered spacecraft; `inertial`: aligned with the orbit "
        + "plane, as for a spacecraft without yaw steering).",
    )

    @model_validator(mode="before")
    @classmethod
    def default_field_of_regard(cls, data: Any) -> Any:
        """
        Defaults the field of regard to contain the instantaneous footprint,
        with a margin of 1 degree.
        """
        if isinstance(data, dict) and data.get("field_of_regard") is None:
            if "cone_angle" in data and "along_track_field_of_view" in data:
                data = {
                    **data,
                    "field_of_regard": min(
                        180,
                        2
                        * (
                            float(data["cone_angle"])
                            + float(data["along_track_field_of_view"]) / 2
                            + 1
                        ),
                    ),
                }
        return data

    def get_swath_width(self, height: float) -> float:
        """
        Gets the instrument swath width projected to the Earth's surface: the
        cross-track extent of the scanned arc, assuming a spherical Earth.

        Args:
            height (float): Height (meters) above surface of the observation.

        Returns:
            float: The swath width (meters).
        """
        radius = constants.EARTH_MEAN_RADIUS
        cone = np.radians(self.cone_angle)
        # earth central angle from nadir to the cone on the surface (limited
        # to the horizon)
        central = np.arcsin(min(1.0, (radius + height) / radius * np.sin(cone))) - cone
        # scan azimuths at the ends of the sector and, if within it, at the
        # extreme left and right of the cone
        azimuths = [
            self.scan_center_azimuth - self.scan_half_width,
            self.scan_center_azimuth + self.scan_half_width,
        ] + [azimuth for azimuth in (-90, 90) if self._in_sector(azimuth)]
        cross = np.arcsin(np.sin(central) * np.sin(np.radians(azimuths)))
        return float(radius * (np.max(cross) - np.min(cross)))

    def _in_sector(self, azimuth: npt.ArrayLike) -> npt.NDArray[np.bool_]:
        """
        Determines if scan azimuths (degrees) lie within the scan sector.
        """
        offset = (np.asarray(azimuth) - self.scan_center_azimuth + 180) % 360 - 180
        return np.abs(offset) <= self.scan_half_width

    def _sub_sectors(self) -> list[tuple[float, float]]:
        """
        Divides the scan sector at the cross-track directions (scan azimuths
        of -90 and 90 degrees), where moving along track moves the arc along
        itself, so that each part sweeps a simple polygon.
        """
        if self.scan_half_width >= 180:
            # a full rotation, cut only at the cross-track directions
            return [(-90, 90), (90, 270)]
        lower = self.scan_center_azimuth - self.scan_half_width
        upper = self.scan_center_azimuth + self.scan_half_width
        cuts = [
            a
            for a in np.arange(-450, 451, 180)
            if lower < a < upper
            and not np.isclose(a, lower)
            and not np.isclose(a, upper)
        ]
        bounds = [lower, *cuts, upper]
        return list(zip(bounds[:-1], bounds[1:]))

    def _sweep_time(self, orbit_track: Geocentric) -> npt.NDArray:
        """
        Gets the time (seconds) over which the satellite's ground track
        advances by half the along-track distance swept by the footprint.
        """
        position, velocity = orbit_track.frame_xyz_and_velocity(itrs)
        radius = np.linalg.norm(np.array(position.m), axis=0)
        altitude = radius - constants.EARTH_MEAN_RADIUS
        ground_speed = (
            np.linalg.norm(np.array(velocity.m_per_s), axis=0)
            * constants.EARTH_MEAN_RADIUS
            / radius
        )
        return (
            altitude
            * np.tan(np.radians(self.along_track_field_of_view / 2))
            / ground_speed
        )

    def compute_footprint(
        self,
        orbit_track: Geocentric,
        number_points: int | None = None,
        elevation: float = 0,
        ground_elevation: float = 0,
    ) -> list[Polygon | MultiPolygon]:
        """
        Compute the instanteous instrument footprint: the region swept by the
        scanned arc as the satellite moves along track by the distance
        subtended at nadir by the along-track field of view.

        Args:
            orbit_track (skyfield.positionlib.Geocentric): The satellite position/velocity.
            number_points (int | None): The number of points along a full
                rotation of each edge of the arc. Defaults to the runtime
                configuration.
            elevation (float): The elevation (meters) at which project the footprint.
            ground_elevation (float): The elevation (meters) of the ground,
                which blocks rays toward an elevation above the instrument.

        Returns:
            list[shapely.geometry.Polygon | shapely.geometry.MultiPolygon]: The instrument footprint(s).
        """
        if number_points is None:
            number_points = config.get_rc().footprint_points_elliptical
        sweep = self._sweep_time(orbit_track)
        fractions = np.linspace(1, -1, 5)
        tracks = {f: _shift_along_track(orbit_track, f * sweep) for f in fractions}

        def project(rays):
            # rays as (shift fraction, cone azimuth) pairs, projected together
            # for each shift fraction
            shifts = np.array([f for f, _ in rays])
            angles = np.array([_cone_ray(self.cone_angle, a) for _, a in rays])
            geos = np.empty((2, len(rays)) + np.shape(orbit_track.t))  # type: ignore
            for f in np.unique(shifts):
                index = np.flatnonzero(shifts == f)
                geos[:, index] = _compute_projected_rays(
                    tracks[f],
                    0,
                    0,
                    angles[index, :1],
                    angles[index, 1:],
                    False,
                    np.zeros(len(index)),
                    elevation,
                    self.velocity_frame,
                    self.nadir_reference,
                    ground_elevation=ground_elevation,
                    attitude=self._attitude,
                )[:2]
            return np.stack([geos[0], geos[1]], axis=-1)

        # boundary of each sub-sector: the forward arc, the end of the arc
        # moving aft, the aft arc (reversed), and the start of the arc moving
        # forward
        # (with the rays of all sub-sectors projected together)
        rays = []
        for lower, upper in self._sub_sectors():
            n = max(4, int(np.ceil(number_points * (upper - lower) / 360)) + 1)
            azimuths = np.linspace(lower, upper, n)
            rays.append(
                [(fractions[0], a) for a in azimuths]
                + [(f, upper) for f in fractions[1:-1]]
                + [(fractions[-1], a) for a in azimuths[::-1]]
                + [(f, lower) for f in fractions[::-1][1:-1]]
            )
        rings = np.split(
            project([ray for sector in rays for ray in sector]),
            np.cumsum([len(sector) for sector in rays])[:-1],
        )
        is_vectorized = len(np.shape(orbit_track.t)) > 0  # type: ignore
        # polygons of each sub-sector (one per time), built at once and split
        # along the anti-meridian and poles, and repaired if invalid
        parts = []
        for ring in rings:
            coords = np.swapaxes(ring, 0, 1) if is_vectorized else ring[np.newaxis]
            parts.append(_build_split_polygons(coords[..., 0], coords[..., 1]))
        footprints = []
        for i in range(len(parts[0])):
            footprint = (
                parts[0][i] if len(parts) == 1 else unary_union([p[i] for p in parts])
            )
            footprints.append(project_polygon_to_elevation(footprint, elevation))
        return footprints

    def compute_sample_size(
        self, orbit_track: Geocentric, target: GeographicPosition
    ) -> tuple[npt.NDArray, npt.NDArray]:
        """
        Computes the size on the ground (meters) of the beam viewing a
        target, along track and across track (see
        `tatc.utils.compute_ground_sample_size`): the beam is circular, but
        its footprint is longer in the direction of the view (along track,
        for a view ahead of the satellite), as the view is oblique. The sizes
        are `NaN` if the beam width is not modeled.

        Args:
            orbit_track (skyfield.positionlib.Geocentric): The satellite position/velocity.
            target (skyfield.toposlib.GeographicPosition): The target position
                (a single target, or one target for each time of the orbit track).

        Returns:
            tuple[numpy.typing.NDArray, numpy.typing.NDArray]: the along-track
                and cross-track sizes (meters).
        """
        if self.beam_width is None:
            nan = np.full(
                np.broadcast(
                    np.asarray(orbit_track.t.tt), np.asarray(target.latitude.degrees)
                ).shape,
                np.nan,
            )
            return nan, nan.copy()
        return compute_ground_sample_size(
            orbit_track,
            target,
            self.beam_width,
            self.beam_width,
            self.velocity_frame,
            nadir_reference=self.nadir_reference,
            attitude=self._attitude,
        )

    def compute_sample_spacing(
        self, orbit_track: Geocentric, target: GeographicPosition
    ) -> tuple[npt.NDArray, npt.NDArray]:
        """
        Computes the spacing on the ground (meters) between the centers of
        neighboring samples viewing a target: along track, the distance the
        scanned arc advances between scans (subtended at nadir by the
        along-track field of view), and across track, the distance between
        samples along the scanned arc (the scan sample spacing times the
        arc's radius on the ground, from the sub-satellite point to the
        target). The cross-track spacing is `NaN` if the scan sample spacing
        is not modeled.

        Args:
            orbit_track (skyfield.positionlib.Geocentric): The satellite position/velocity.
            target (skyfield.toposlib.GeographicPosition): The target position
                (a single target, or one target for each time of the orbit track).

        Returns:
            tuple[numpy.typing.NDArray, numpy.typing.NDArray]: the along-track
                and cross-track spacings (meters).
        """
        position = np.reshape(np.array(orbit_track.frame_xyz(itrs).m), (3, -1))
        radius = np.linalg.norm(position, axis=0)
        # the along-track distance swept per scan (see _sweep_time)
        along = (
            2
            * (radius - constants.EARTH_MEAN_RADIUS)
            * np.tan(np.radians(self.along_track_field_of_view / 2))
        )
        target_xyz = np.reshape(np.array(target.itrs_xyz.m), (3, -1))
        # the Earth central angle from the sub-satellite point to the target
        central = np.arccos(
            np.clip(
                np.sum(position * target_xyz, axis=0)
                / radius
                / np.linalg.norm(target_xyz, axis=0),
                -1,
                1,
            )
        )
        along, central = np.broadcast_arrays(along, central)
        cross = np.full(central.shape, np.nan)
        if self.scan_sample_spacing is not None:
            cross = (
                central
                * np.linalg.norm(target_xyz, axis=0)
                * np.radians(self.scan_sample_spacing)
            )
        if (
            np.ndim(np.asarray(orbit_track.t.tt)) == 0
            and np.ndim(target.latitude.degrees) == 0
        ):
            return along[0], cross[0]
        return along, cross

    def compute_footprint_center(
        self,
        orbit_track: Geocentric,
        elevation: float = 0,
        ground_elevation: float = 0,
    ) -> GeographicPosition:
        """
        Compute the center of an instaneous instrument footprint: the point
        on the cone at the center of the scan sector.

        Args:
            orbit_track (skyfield.positionlib.Geocentric): The satellite position/velocity.
            elevation (float): The elevation (meters) at which project the footprint.
            ground_elevation (float): The elevation (meters) of the ground,
                which blocks rays toward an elevation above the instrument.

        Returns:
            skyfield.toposlib.GeographicPosition: The instrument footprint center.
        """
        roll, pitch = _cone_ray(self.cone_angle, self.scan_center_azimuth)
        return compute_projected_ray_position(
            orbit_track,
            0,
            0,
            roll,
            pitch,
            elevation=elevation,
            velocity_frame=self.velocity_frame,
            nadir_reference=self.nadir_reference,
            ground_elevation=ground_elevation,
            attitude=self._attitude,
        )

    def is_in_field_of_view(
        self,
        orbit_track: Geocentric,
        target: GeographicPosition,
    ) -> npt.NDArray[np.bool_]:
        """
        Determines if a target lies within the instantaneous footprint: if it
        crosses the cone within the scan sector as the satellite moves along
        track by the distance subtended at nadir by the along-track field of
        view (linearly interpolating its cone angle and scan azimuth between
        the start, middle, and end of the motion). Does not check whether
        the target is above the satellite's horizon.

        Args:
            orbit_track (skyfield.positionlib.Geocentric): The satellite position/velocity.
            target (skyfield.toposlib.GeographicPosition): The target position.

        Returns:
            numpy.typing.NDArray: Array of indicators: `True` if the target is in the field of view.
        """
        sweep = self._sweep_time(orbit_track)
        angles = [
            compute_cone_and_azimuth(
                _shift_along_track(orbit_track, f * sweep),
                target,
                self.velocity_frame,
                self.nadir_reference,
                attitude=self._attitude,
            )
            for f in (-1, 0, 1)
        ]
        inside = np.zeros(np.size(orbit_track.t.tt), dtype=bool)  # type: ignore
        for (cone_a, azimuth_a), (cone_b, azimuth_b) in zip(angles[:-1], angles[1:]):
            cone_a, cone_b = np.atleast_1d(cone_a), np.atleast_1d(cone_b)
            crosses = (cone_a - self.cone_angle) * (cone_b - self.cone_angle) <= 0
            with np.errstate(divide="ignore", invalid="ignore"):
                fraction = np.where(
                    cone_a == cone_b, 0, (self.cone_angle - cone_a) / (cone_b - cone_a)
                )
            # interpolate the scan azimuth along the shorter way around
            change = (np.atleast_1d(azimuth_b) - azimuth_a + 180) % 360 - 180
            inside |= crosses & self._in_sector(azimuth_a + fraction * change)
        return inside

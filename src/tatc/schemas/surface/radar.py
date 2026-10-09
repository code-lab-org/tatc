"""
Object schemas for ground-based radar stations.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from enum import Enum

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, Field, model_validator
from shapely.geometry import MultiPolygon, Polygon

from ... import config
from ...utils.radar import compute_radar_footprint, compute_radar_footprint_profile
from ...utils.radar import compute_radar_ground_range_bounds
from .point import Point


class TerrainMask(BaseModel):
    """
    Azimuthal terrain blockage mask for a ground-based radar station: the
    minimum usable elevation angle at each azimuth, due to terrain
    blocking the station's line of sight at low grazing angles. Typically
    derived by ray-tracing a digital elevation model (DEM) outward from
    the station to find the highest terrain obstruction angle along each
    azimuthal direction.

    This represents one blocking angle per azimuth (the common case for a
    DEM-derived terrain horizon, where the nearest/highest ridge along a
    ray dominates): any elevation angle at or above the masked value is
    assumed clear, and any below it is blocked. It does not model
    multiple, separately-blocked elevation bands along the same azimuth
    (e.g. a nearby building blocking only a narrow low-elevation slice
    while a distant mountain blocks a separate, higher slice).
    """

    azimuth: list[float] = Field(
        ...,
        description="Azimuth samples (decimal degrees, clockwise from "
        + "north, in [0, 360)), strictly increasing, at which the minimum "
        + "elevation angle is specified. Values are linearly interpolated "
        + "between samples, wrapping around 0/360 degrees.",
        min_length=2,
    )
    min_elevation_angle: list[float] = Field(
        ...,
        description="Minimum usable elevation angle (decimal degrees), "
        + "due to terrain blockage, at each corresponding azimuth.",
        min_length=2,
    )

    @model_validator(mode="after")
    def validate_profile(self) -> TerrainMask:
        """
        Validates that `azimuth` and `min_elevation_angle` are the same
        length, and that `azimuth` values fall in [0, 360) and are
        strictly increasing.
        """
        if len(self.azimuth) != len(self.min_elevation_angle):
            raise ValueError(
                "azimuth and min_elevation_angle must have the same length"
            )
        if any(a < 0 or a >= 360 for a in self.azimuth):
            raise ValueError("azimuth values must be in [0, 360)")
        if any(b <= a for a, b in zip(self.azimuth, self.azimuth[1:])):
            raise ValueError("azimuth values must be strictly increasing")
        return self

    def get_min_elevation_angle(self, azimuth: float) -> float:
        """
        Interpolates (periodically, wrapping at 0/360 degrees) the minimum
        usable elevation angle at a specified azimuth.

        Args:
            azimuth (float): Azimuth (decimal degrees, clockwise from north).

        Returns:
            float: The interpolated minimum usable elevation angle (degrees).
        """
        return float(self.get_min_elevation_angles(azimuth))

    def get_min_elevation_angles(self, azimuths: npt.ArrayLike) -> npt.NDArray:
        """
        Interpolates (periodically, wrapping at 0/360 degrees) the minimum
        usable elevation angles at specified azimuths, vectorized across
        azimuths (see `get_min_elevation_angle`).

        Args:
            azimuths (numpy.typing.ArrayLike): Azimuths (decimal degrees,
                clockwise from north).

        Returns:
            numpy.typing.NDArray: The interpolated minimum usable elevation
            angles (degrees).
        """
        samples = np.concatenate(
            ([self.azimuth[-1] - 360], self.azimuth, [self.azimuth[0] + 360])
        )
        angles = np.concatenate(
            (
                [self.min_elevation_angle[-1]],
                self.min_elevation_angle,
                [self.min_elevation_angle[0]],
            )
        )
        return np.interp(np.asarray(azimuths) % 360, samples, angles)


class RadarBand(str, Enum):
    """
    Common weather radar frequency bands. This is an informational tag
    only: TAT-C does not model transient, weather-dependent propagation
    effects (e.g. rain attenuation), which is the dominant practical
    difference between bands in reality -- a heavy storm can attenuate an
    X-band signal to the point of near-total extinction a short distance
    behind it, a C-band signal moderately, and an S-band signal barely at
    all. The static effect this schema does represent is a shorter
    typical maximum range at shorter wavelengths, reflecting that overall
    attenuation/sensitivity tradeoff in clear-to-moderate conditions; see
    `RadarStation.from_band`.
    """

    S = "S"
    C = "C"
    X = "X"


# illustrative, order-of-magnitude nominal maximum range (meters) by band,
# NOT a model of any specific radar: S-band mirrors the NEXRAD WSR-88D Level II
# reflectivity range (negligible rain attenuation, long range); C-band and X-band are shortened to loosely
# reflect their greater susceptibility to rain attenuation, consistent with
# the shorter ranges typical of real C-band national networks and X-band
# short-range/gap-filling/mobile radars, respectively
_RADAR_BAND_NOMINAL_MAX_RANGE = {
    RadarBand.S: 460000,
    RadarBand.C: 150000,
    RadarBand.X: 60000,
}


class RadarStation(Point):
    """
    Ground-based radar station (e.g. a NOAA NEXRAD WSR-88D weather radar)
    in the WGS 84 coordinate system.

    Absent a `terrain_mask`, the coverage footprint computed by this
    schema assumes an idealized, azimuthally symmetric sensor with clear
    line of sight in every direction. Supplying a `terrain_mask` (e.g.
    derived from a digital elevation model) relaxes that assumption by
    raising the effective minimum elevation angle in blocked directions.
    Coverage extends half the `beam_width` beyond the lowest and highest
    scanned beam centers, treating the elevation angles in between as
    continuously sampled (gaps between widely spaced tilts are not
    represented). In all cases, this schema uses the standard-atmosphere "4/3 Earth
    radius" refraction approximation, which does not capture anomalous
    propagation (ducting/sub-refraction) under non-standard weather
    conditions, and reports a 2-D ground-range profile bounding where a
    target at a given height is observable, not a full 3-D reconstruction
    of the scanned volume.
    """

    name: str = Field(..., description="Radar station name", examples=["KOUN"])
    elevation: float = Field(
        default=0,
        description="Elevation (meters) of the radar antenna, including "
        + "tower height. Coverage depends only on target elevations "
        + "relative to this value, so it must use the same vertical "
        + "reference as target elevations and any terrain mask: heights "
        + "above mean sea level (as in NOAA's NEXRAD station tables and "
        + "digital elevation models such as Copernicus DEM GLO-30) are "
        + "recommended, noting that they differ from heights above the "
        + "WGS 84 ellipsoid by the local geoid undulation.",
    )
    max_range: float = Field(
        default=460000,
        description="Maximum unambiguous slant range (meters) for radar "
        + "detection, limited by the pulse repetition frequency (PRF) and "
        + "receiver hardware. Defaults to 460,000 m (460 km), the range of "
        + "NEXRAD WSR-88D Level II reflectivity on its lowest (surveillance) "
        + "tilts; the conventional base reflectivity product is displayed "
        + "to 230,000 m.",
        gt=0,
        examples=[460000],
    )
    min_elevation_angle: float = Field(
        default=0.5,
        description="Nominal (hardware/operational) lowest scanned "
        + "elevation angle (decimal degrees) above local horizontal. "
        + "Defaults to 0.5 degrees, the lowest elevation cut in most "
        + "NEXRAD volume coverage patterns. The effective minimum "
        + "elevation angle at a given azimuth is never below this value, "
        + "but may be raised locally by `terrain_mask`. Negative values "
        + "represent scanning below local horizontal, as at some elevated "
        + "sites (e.g. NEXRAD KFSX scans at about -0.2 degrees); "
        + "intersection of such a beam with the terrain is represented "
        + "only through `terrain_mask`.",
        ge=-90,
        le=90,
    )
    max_elevation_angle: float = Field(
        default=19.5,
        description="Highest scanned elevation angle (decimal degrees) "
        + "above local horizontal. Defaults to 19.5 degrees, a common "
        + "highest elevation cut in NEXRAD volume coverage patterns.",
        ge=0,
        le=90,
    )
    terrain_mask: TerrainMask | None = Field(
        default=None,
        description="Optional azimuthal terrain blockage mask (e.g. "
        + "derived from a digital elevation model), raising the effective "
        + "minimum elevation angle in specific directions. When omitted, "
        + "coverage is assumed azimuthally symmetric. A raised scan angle "
        + "shortens the outer range at which a target above the station is "
        + "observable, and removes coverage of a target below the station "
        + "wherever it lifts the lowest usable angle to or above local "
        + "horizontal.",
    )
    beam_width: float = Field(
        default=0.95,
        description="Half-power (-3 dB) beam width (decimal degrees). "
        + "Defaults to 0.95 degrees, the NEXRAD WSR-88D beam width. Half "
        + "of this width extends the observed elevation angles below the "
        + "lowest scanned beam center (`min_elevation_angle`) and above "
        + "the highest (`max_elevation_angle`); set to 0 to bound coverage "
        + "by the beam centers alone.",
        ge=0,
        le=180,
    )
    band: RadarBand | None = Field(
        default=None,
        description="Optional informational tag for this station's "
        + "nominal operating frequency band (S, C, or X). It does not "
        + "affect any coverage computation directly -- set `max_range` "
        + "(and other fields) explicitly to represent a band's practical "
        + "effect; see `from_band` for illustrative nominal defaults.",
    )

    @model_validator(mode="after")
    def max_elevation_angle_ge_min_elevation_angle(self) -> RadarStation:
        """
        Validates that the maximum elevation angle is not less than the
        minimum elevation angle.
        """
        if self.max_elevation_angle < self.min_elevation_angle:
            raise ValueError(
                "max_elevation_angle must be greater than or equal to min_elevation_angle"
            )
        return self

    @classmethod
    def from_band(cls, band: RadarBand, **kwargs) -> RadarStation:
        """
        Constructs a `RadarStation` with an illustrative nominal
        `max_range` typical of a given weather radar frequency band, for
        a reasonable starting point without needing to know specific
        hardware numbers.

        These are representative orders of magnitude, not a model of any
        particular radar, and the real driver of a band's practical range
        -- precipitation attenuation, which TAT-C does not model as a
        transient, weather-dependent effect -- is only indirectly
        represented by this static shortening: negligible (effectively no
        shortening) at S-band, moderate at C-band, and severe (so a much
        shorter nominal range) at X-band. Pass `max_range` as a keyword
        argument to override the band default.

        Args:
            band (RadarBand): The nominal operating frequency band.
            **kwargs: Other `RadarStation` fields (e.g. `name`,
                `latitude`, `longitude`, `elevation`), including an
                optional `max_range` override.

        Returns:
            RadarStation: A station with the band tagged and its nominal
            `max_range` applied (unless overridden).
        """
        defaults = {"band": band, "max_range": _RADAR_BAND_NOMINAL_MAX_RANGE[band]}
        return cls(**{**defaults, **kwargs})

    def get_effective_max_elevation_angle(self) -> float:
        """
        Gets the highest observed elevation angle: the upper half-power
        edge of the highest scanned beam.

        Returns:
            float: The effective maximum elevation angle (degrees),
            `max_elevation_angle` plus half the `beam_width`, capped at 90.
        """
        return min(self.max_elevation_angle + self.beam_width / 2, 90)

    def get_effective_min_elevation_angle(self, azimuth: float) -> float:
        """
        Gets the lowest observed elevation angle at a specified azimuth,
        accounting for this station's `beam_width` and `terrain_mask` (if
        any).

        Args:
            azimuth (float): Azimuth (decimal degrees, clockwise from north).

        Returns:
            float: The effective minimum elevation angle (degrees): the
            lower half-power edge of the lowest scanned beam
            (`min_elevation_angle` minus half the `beam_width`), raised to
            any terrain blockage at that azimuth (terrain hides the part
            of the beam below it), and capped at the effective maximum
            elevation angle (an azimuth blocked beyond it is simply
            unobservable, not scannable at a steeper angle still).
        """
        effective = max(self.min_elevation_angle - self.beam_width / 2, -90)
        if self.terrain_mask is not None:
            effective = max(
                effective, self.terrain_mask.get_min_elevation_angle(azimuth)
            )
        return min(effective, self.get_effective_max_elevation_angle())

    def compute_ground_ranges(self, elevation: float) -> tuple[float, float] | None:
        """
        Computes the idealized, azimuthally symmetric ground-range annulus
        (meters) within which this station can detect a target at a
        specified elevation, between the lower half-power edge of the
        lowest scanned beam and the upper edge of the highest (this does
        not account for `terrain_mask`; see `compute_ground_range_profile`
        for the azimuth-resolved bounds).

        Args:
            elevation (float): The elevation (meters) of the observed
                target, in the same vertical reference as this station's
                `elevation`.

        Returns:
            tuple[float, float] | None: The `(inner_ground_range,
            outer_ground_range)` bounds (meters), or `None` if the target
            is not observable at any ground range.
        """
        return compute_radar_ground_range_bounds(
            max(self.min_elevation_angle - self.beam_width / 2, -90),
            self.get_effective_max_elevation_angle(),
            self.max_range,
            elevation,
            self.elevation,
        )

    def compute_ground_range_profile(
        self, elevation: float, number_points: int | None = None
    ) -> list[tuple[float, float, float]]:
        """
        Computes the ground-range annulus bounds at a sampled set of
        azimuths, accounting for this station's `terrain_mask` (if any).

        Args:
            elevation (float): The elevation (meters) of the observed
                target, in the same vertical reference as this station's
                `elevation`.
            number_points (int | None): The number of azimuth samples to
                generate (evenly spaced over one full revolution).
                Defaults to the runtime configuration.

        Returns:
            list[tuple[float, float, float]]: A list of `(azimuth,
            inner_ground_range, outer_ground_range)` tuples; a fully
            blocked azimuth reports `(azimuth, 0, 0)`.
        """
        if number_points is None:
            number_points = config.get_rc().footprint_points_radar_azimuthal
        azimuths = np.linspace(0, 360, number_points, endpoint=False)
        # effective minimum elevation angle at every azimuth at once (see
        # get_effective_min_elevation_angle)
        min_elevation_angles = np.full(
            len(azimuths), max(self.min_elevation_angle - self.beam_width / 2, -90)
        )
        if self.terrain_mask is not None:
            min_elevation_angles = np.maximum(
                min_elevation_angles,
                self.terrain_mask.get_min_elevation_angles(azimuths),
            )
        min_elevation_angles = np.minimum(
            min_elevation_angles, self.get_effective_max_elevation_angle()
        )
        # ground range bounds, computed once per distinct elevation angle
        bounds = {}
        profile = []
        for azimuth, min_elevation_angle in zip(azimuths, min_elevation_angles):
            min_elevation_angle = float(min_elevation_angle)
            if min_elevation_angle not in bounds:
                bounds[min_elevation_angle] = compute_radar_ground_range_bounds(
                    min_elevation_angle,
                    self.get_effective_max_elevation_angle(),
                    self.max_range,
                    elevation,
                    self.elevation,
                ) or (0.0, 0.0)
            inner_ground_range, outer_ground_range = bounds[min_elevation_angle]
            profile.append((float(azimuth), inner_ground_range, outer_ground_range))
        return profile

    def compute_footprint(
        self, elevation: float, number_points: int | None = None
    ) -> Polygon | MultiPolygon:
        """
        Computes this station's static ground coverage footprint (a disk,
        or an annulus if an overhead "cone of silence" applies or, for a
        target below the station, the lowest beam only descends to the
        target some distance out) at a specified target elevation. When `terrain_mask` is set, the
        footprint is azimuthally irregular, reflecting blocked directions.

        Args:
            elevation (float): The elevation (meters) of the observed
                target, in the same vertical reference as this station's
                `elevation`.
            number_points (int | None): The number of azimuth samples used
                when this station has a `terrain_mask` (ignored otherwise,
                since the symmetric case uses a faster closed-form circle/
                annulus). Defaults to the runtime configuration.

        Returns:
            shapely.geometry.Polygon | shapely.geometry.MultiPolygon: The
            radar coverage footprint, or an empty `Polygon` if the target
            is not observable at any ground range.
        """
        if self.terrain_mask is None:
            ranges = self.compute_ground_ranges(elevation)
            if ranges is None:
                return Polygon()
            inner_ground_range, outer_ground_range = ranges
            return compute_radar_footprint(
                self.longitude,
                self.latitude,
                inner_ground_range,
                outer_ground_range,
                elevation,
            )
        profile = self.compute_ground_range_profile(elevation, number_points)
        if all(outer_ground_range <= 0 for _, _, outer_ground_range in profile):
            return Polygon()
        # for a target above the station, the inner ground range is governed
        # only by the effective maximum elevation angle (not terrain), so it
        # is uniform across non-blocked samples and subtracted as a circle;
        # for a target below the station, it depends on the (terrain-raised)
        # lowest usable angle, so it is subtracted as a sampled profile
        inner_ground_ranges = [inner for _, inner, outer in profile if outer > 0]
        if max(inner_ground_ranges) - min(inner_ground_ranges) < 1e-6:
            inner_ground_range = inner_ground_ranges[0]
        else:
            inner_ground_range = [inner for _, inner, _ in profile]
        return compute_radar_footprint_profile(
            self.longitude,
            self.latitude,
            [azimuth for azimuth, _, _ in profile],
            [outer for _, _, outer in profile],
            inner_ground_range,
            elevation,
        )

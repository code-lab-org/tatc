"""
Object schemas for limb sounders and radio occultation receivers, which profile
the atmosphere through the limb rather than viewing a scene.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timedelta
from enum import Enum
from typing import Any, ClassVar, Literal, NoReturn

from pydantic import BaseModel, ConfigDict, Field, model_validator

from .simple import Instrument


class ScanDirection(str, Enum):
    """
    Enumeration of the two broad classes of limb sounder vertical scan:
    sweeping from the lowest to the highest requested tangent point
    elevation (`UPWARD`), or the reverse (`DOWNWARD`).
    """

    UPWARD = "upward"
    DOWNWARD = "downward"


def _default_scan_direction(scan_azimuth: float) -> ScanDirection:
    """
    Chooses a default scan direction from the sensor's viewing azimuth: a
    forward-looking sensor (`scan_azimuth` closer to 0 deg than to 180 deg)
    defaults to an upward scan, and a rearward-looking one (closer to 180
    deg than to 0 deg) defaults to a downward scan. Exactly sideways
    (90 or 270 deg, equidistant from both) defaults to upward. Compared as
    exact angular distances (not trigonometrically), so this tie is exact
    rather than subject to floating-point rounding noise.
    """
    normalized = scan_azimuth % 360.0
    distance_to_forward = min(normalized, 360.0 - normalized)
    distance_to_rearward = abs(normalized - 180.0)
    return (
        ScanDirection.UPWARD
        if distance_to_forward <= distance_to_rearward
        else ScanDirection.DOWNWARD
    )


# the fields of the base instrument that a sounder does not use, as its
# observations are profiles through the limb rather than views of targets
_UNUSED_FIELDS = (
    "field_of_regard",
    "min_access_time",
    "access_time_fixed",
    "nadir_reference",
    "req_self_sunlit",
    "req_target_sunlit",
    "min_target_solar_elevation",
    "max_target_solar_elevation",
    "min_target_glint_angle",
    "max_target_glint_angle",
)


class _Sounder(Instrument):
    """
    Base class of instruments that profile the atmosphere through the limb, whose
    observations are collected with their own analyses rather than as views of
    targets. Only a satellite carries one, and its view is measured from the
    satellite's inertial velocity, without the attitude of its body.
    """

    _analysis: ClassVar[str] = ""

    @model_validator(mode="after")
    def _check_unused_fields(self) -> "_Sounder":
        defaults = Instrument.model_fields
        unused = [
            field
            for field in _UNUSED_FIELDS
            # pylint: disable-next=unsubscriptable-object
            if getattr(self, field) != defaults[field].get_default()
        ]
        if unused:
            raise ValueError(
                f"{type(self).__name__} {self.name!r} does not use {unused}: its "
                f"observations are collected with {self._analysis}, not as views "
                "of targets"
            )
        return self

    def _not_a_view(self) -> NoReturn:
        raise TypeError(
            f"{type(self).__name__} {self.name!r} has no view of targets: collect "
            f"its observations with tatc.analysis.{self._analysis}"
        )

    def get_swath_width(self, height: float) -> float:
        """
        Not defined for a sounder, which raises a `TypeError`.
        """
        del height
        self._not_a_view()

    def get_min_elevation_angle(self, height: float) -> float:
        """
        Not defined for a sounder, which raises a `TypeError`.
        """
        del height
        self._not_a_view()

    def compute_footprint(self, *args: Any, **kwargs: Any) -> Any:
        """
        Not defined for a sounder, which raises a `TypeError`.
        """
        del args, kwargs
        self._not_a_view()

    def compute_footprint_center(self, *args: Any, **kwargs: Any) -> Any:
        """
        Not defined for a sounder, which raises a `TypeError`.
        """
        del args, kwargs
        self._not_a_view()

    def is_valid_observation(self, *args: Any, **kwargs: Any) -> Any:
        """
        Not defined for a sounder, which raises a `TypeError`.
        """
        del args, kwargs
        self._not_a_view()


class ScanStart(BaseModel):
    """
    The start of a vertical scan in a limb sounder's cycle of scans.
    """

    model_config = ConfigDict(extra="forbid")

    offset: timedelta = Field(
        default=timedelta(0),
        ge=timedelta(0),
        description="The time from the start of the cycle.",
    )
    direction: ScanDirection = Field(
        ..., description="Whether the scan sweeps `upward` or `downward`."
    )


class LimbSounder(_Sounder):
    """
    Limb sounder, which views the limb at a fixed azimuth from the satellite's
    inertial velocity and sweeps its tangent point through a vertical scan of
    elevations (see `tatc.analysis.collect_limb_soundings`). Its scans repeat in
    a cycle (one or more scans, such as alternating upward and downward scans)
    from the start of an analysis: the phase of the cycle, which depends on the
    instrument's clock, is not modeled. Only a satellite carries it.
    """

    _analysis: ClassVar[str] = "collect_limb_soundings"

    type: Literal["limb"] = Field(  # type: ignore[assignment]
        default="limb", description="Instrument type discriminator."
    )
    scan_azimuth: float = Field(
        ...,
        ge=0,
        lt=360,
        description="The viewing azimuth (degrees) from the inertial velocity, in "
        + "the local horizontal plane (0 forward, 90 to the left, toward the orbit "
        + "normal).",
        examples=[0],
    )
    scan_elevations: list[float] = Field(
        ...,
        min_length=2,
        description="The tangent point elevations (meters) of the samples of a "
        + "scan, in any order (the direction of each scan sets the order).",
        examples=[[10000, 30000, 50000, 70000, 90000]],
    )
    scan_duration: timedelta = Field(
        ...,
        gt=timedelta(0),
        description="The time to sweep the elevations, at a constant angular rate.",
        examples=[timedelta(seconds=20)],
    )
    scans: list[ScanStart] | None = Field(
        default=None,
        min_length=1,
        description="The starts of the scans of a cycle, at increasing offsets, "
        + "or `None` for one scan at the start of the cycle in the default "
        + "direction (upward for a view closer to forward than to rearward, "
        + "else downward).",
    )
    cycle_period: timedelta | None = Field(
        default=None,
        gt=timedelta(0),
        description="The time between starts of the cycle of scans, or `None` "
        + "for the scan duration.",
        examples=[timedelta(seconds=25)],
    )
    sample_elevation: float | None = Field(
        default=None,
        description="The tangent point elevation (meters) of the representative "
        + "position and time of each scan, or `None` for the middle of the scan's "
        + "elevations.",
        examples=[25000],
    )

    @model_validator(mode="after")
    def _check_scans(self) -> "LimbSounder":
        offsets = [scan.offset for scan in self.get_scans()]
        if offsets != sorted(set(offsets)) or offsets[-1] >= self.get_cycle_period():
            raise ValueError(
                "scans start at increasing offsets within the cycle period"
            )
        if (
            not min(self.scan_elevations)
            <= self.get_sample_elevation()
            <= max(self.scan_elevations)
        ):
            raise ValueError("the sample elevation is within the scan's elevations")
        return self

    def get_scans(self) -> list[ScanStart]:
        """
        Gets the starts of the scans of a cycle: those given, or one scan at the
        start of the cycle in the default direction.

        Returns:
            list[ScanStart]: The starts of the scans.
        """
        if self.scans is not None:
            return list(self.scans)
        return [ScanStart(direction=_default_scan_direction(self.scan_azimuth))]

    def get_cycle_period(self) -> timedelta:
        """
        Gets the time between starts of the cycle of scans: that given, or the
        scan duration.

        Returns:
            datetime.timedelta: The period.
        """
        return self.scan_duration if self.cycle_period is None else self.cycle_period

    def get_sample_elevation(self) -> float:
        """
        Gets the tangent point elevation (meters) of the representative position
        and time of each scan: that given, or the middle of the scan's elevations.

        Returns:
            float: The elevation (meters).
        """
        if self.sample_elevation is not None:
            return self.sample_elevation
        return (min(self.scan_elevations) + max(self.scan_elevations)) / 2

    def get_scan_times(
        self, start: datetime, end: datetime
    ) -> dict[ScanDirection, list[datetime]]:
        """
        Gets the start times of the scans of each direction, in cycles that repeat
        from the start of an analysis until its end.

        Args:
            start (datetime.datetime): The start of the analysis.
            end (datetime.datetime): The end of the analysis.

        Returns:
            dict[ScanDirection, list[datetime.datetime]]: The start times, by
            direction.
        """
        scans = self.get_scans()
        period = self.get_cycle_period()
        times: dict[ScanDirection, list[datetime]] = {
            scan.direction: [] for scan in scans
        }
        if end <= start:
            return times
        for cycle in range((end - start) // period + 1):
            for scan in scans:
                time_ = start + cycle * period + scan.offset
                if time_ < end:
                    times[scan.direction].append(time_)
        return times


class OccultationReceiver(_Sounder):
    """
    Radio occultation receiver, which tracks the signals of navigation satellites
    (the transmitters of an analysis) setting and rising behind the Earth's limb
    (see `tatc.analysis.collect_occultations`). Refraction is not modeled: the
    tangent point is that of the straight line between the satellites. Only a
    satellite carries it.
    """

    _analysis: ClassVar[str] = "collect_occultations"

    type: Literal["occultation"] = Field(  # type: ignore[assignment]
        default="occultation", description="Instrument type discriminator."
    )
    constellations: list[str] = Field(
        default=[],
        description="The constellations of transmitters tracked (such as `gps`, "
        + "`glonass`, `galileo`, or `beidou`, the identifiers of the catalog's "
        + "navigation systems; see `tatc.catalog.get_transmitters`). An analysis "
        + "is given its transmitters, so these do not affect it.",
        examples=[["gps", "galileo"]],
    )
    max_yaw: float = Field(
        default=65,
        gt=0,
        le=90,
        description="The maximum yaw (degrees) of a transmitter from the "
        + "receiver's velocity or anti-velocity, for a valid occultation.",
    )
    range_elevation: tuple[float, float] = Field(
        default=(-200e3, 60e3),
        description="The lower and upper tangent point elevations (meters, of the "
        + "straight line) of a valid occultation.",
    )
    sample_elevation: float = Field(
        default=-80e3,
        description="The tangent point elevation (meters, of the straight line) of "
        + "the representative position and time of each occultation.",
    )
    min_profile_duration: timedelta = Field(
        default=timedelta(seconds=30),
        ge=timedelta(0),
        description="The minimum duration of a valid occultation.",
    )

    @model_validator(mode="after")
    def _check_elevations(self) -> "OccultationReceiver":
        lower, upper = self.range_elevation
        if not lower < self.sample_elevation < upper:
            raise ValueError("the sample elevation is within the range of elevations")
        return self

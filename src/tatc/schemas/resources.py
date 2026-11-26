"""
Object schemas for the performance and resources of instruments and
spacecraft, which feed sizing, data volume, and observation operators but not
the geometry of observations.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Iterable
from datetime import timedelta

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, ConfigDict, Field, model_validator

from ..utils.decay import DEFAULT_DRAG_COEFFICIENT


class InstrumentPerformance(BaseModel):
    """
    The nominal performance of an instrument at a reference altitude, scaled to
    other altitudes by power laws: the ground sample distance as
    `(altitude / reference_altitude) ** ground_sample_distance_exponent` and the
    signal-to-noise ratio as
    `(reference_altitude / altitude) ** signal_to_noise_exponent`.

    Typical exponents of the ground sample distance are 1 for a fixed angular
    resolution (an imager's or radiometer's instantaneous field of view, a
    real-aperture radar's beam, or a lidar's divergence), 0.5 for a
    pulse-limited radar altimeter, and 0 for a synthetic aperture radar, whose
    resolution is set by its bandwidth and antenna rather than its range. Typical
    exponents of the signal-to-noise ratio are 0 for a passive instrument whose
    instantaneous field of view is fixed (the radiance it collects does not
    depend on range), 2 for a lidar or a radar viewing a target that fills its
    beam, 3 for a synthetic aperture radar viewing a distributed target, and 4
    for a radar viewing a point target.
    """

    model_config = ConfigDict(extra="forbid")

    reference_altitude: float = Field(
        ...,
        gt=0,
        description="The altitude (meters) at which the performance is stated.",
        examples=[824000],
    )
    ground_sample_distance: float | None = Field(
        default=None,
        gt=0,
        description="The nominal ground sample distance (meters) at nadir at the "
        + "reference altitude: the spatial resolution published for the "
        + "instrument's principal product, or `None` if not given.",
        examples=[750],
    )
    ground_sample_distance_exponent: float = Field(
        default=1,
        ge=0,
        le=4,
        description="The exponent of the ratio of the altitude to the reference "
        + "altitude by which the ground sample distance scales (1 for a fixed "
        + "angular resolution, 0 for a synthetic aperture radar).",
    )
    signal_to_noise_ratio: float | None = Field(
        default=None,
        gt=0,
        description="The nominal signal-to-noise ratio at the reference altitude, "
        + "for the reference scene of its source, or `None` if not given.",
        examples=[200],
    )
    signal_to_noise_exponent: float = Field(
        default=0,
        ge=0,
        le=4,
        description="The exponent of the ratio of the reference altitude to the "
        + "altitude by which the signal-to-noise ratio scales (0 for a passive "
        + "instrument with a fixed instantaneous field of view, 2 for a lidar).",
    )

    def get_noise_scale(self, altitude: npt.ArrayLike) -> float | npt.NDArray:
        """
        Gets the factor by which noise relative to the signal grows from the
        reference altitude: the inverse of the scaling of the signal-to-noise
        ratio, with which to scale an observation error (such as a noise
        equivalent temperature difference) stated at the reference altitude.

        Args:
            altitude (numpy.typing.ArrayLike): The altitude (meters).

        Returns:
            float | numpy.typing.NDArray: The factor.
        """
        return _scale(altitude, self.reference_altitude, self.signal_to_noise_exponent)

    def get_signal_to_noise_ratio(self, altitude: npt.ArrayLike) -> float | npt.NDArray:
        """
        Gets the signal-to-noise ratio at an altitude.

        Args:
            altitude (numpy.typing.ArrayLike): The altitude (meters).

        Returns:
            float | numpy.typing.NDArray: The signal-to-noise ratio.

        Raises:
            ValueError: If the signal-to-noise ratio is not given.
        """
        if self.signal_to_noise_ratio is None:
            raise ValueError("the signal-to-noise ratio is not given")
        return self.signal_to_noise_ratio / self.get_noise_scale(altitude)

    def get_ground_sample_distance(
        self, altitude: npt.ArrayLike
    ) -> float | npt.NDArray:
        """
        Gets the ground sample distance (meters) at nadir at an altitude.

        Args:
            altitude (numpy.typing.ArrayLike): The altitude (meters).

        Returns:
            float | numpy.typing.NDArray: The ground sample distance (meters).

        Raises:
            ValueError: If the ground sample distance is not given.
        """
        if self.ground_sample_distance is None:
            raise ValueError("the ground sample distance is not given")
        return self.ground_sample_distance * _scale(
            altitude, self.reference_altitude, self.ground_sample_distance_exponent
        )


def _scale(
    altitude: npt.ArrayLike, reference_altitude: float, exponent: float
) -> float | npt.NDArray:
    """
    Gets `(altitude / reference_altitude) ** exponent`, raising a `ValueError`
    for an altitude that is not positive.
    """
    ratio = np.asarray(altitude, dtype=float) / reference_altitude
    if np.any(ratio <= 0):
        raise ValueError("the altitude must be positive")
    scale = ratio**exponent
    return float(scale) if np.ndim(scale) == 0 else scale


class InstrumentResources(BaseModel):
    """
    The resources an instrument uses on its platform. A value that is not
    given is `None`. Where several instruments model the parts or alternative
    views of one physical instrument, one of them carries its resources (see
    `tatc.catalog.InstrumentEntry.build`), so that totals (see `total`) count
    each physical instrument once.
    """

    model_config = ConfigDict(extra="forbid")

    mass: float | None = Field(
        default=None, ge=0, description="The mass (kg).", examples=[275]
    )
    power: float | None = Field(
        default=None,
        ge=0,
        description="The average power (watts) while it operates.",
        examples=[200],
    )
    data_rate: float | None = Field(
        default=None,
        ge=0,
        description="The average rate (bits per second) at which it produces "
        + "data while it observes.",
        examples=[10.5e6],
    )

    def get_data_volume(self, duration: timedelta) -> float:
        """
        Gets the volume of data (bits) produced while observing for a duration.

        Args:
            duration (datetime.timedelta): The duration of observation.

        Returns:
            float: The volume (bits).

        Raises:
            ValueError: If the data rate is not given.
        """
        if self.data_rate is None:
            raise ValueError("the data rate is not given")
        return self.data_rate * duration.total_seconds()

    @classmethod
    def total(
        cls, resources: Iterable[InstrumentResources | None]
    ) -> InstrumentResources:
        """
        Totals the resources of instruments: each value is the sum of the
        values given, or `None` if none is given. An instrument without a value
        adds nothing, as for the parts of a physical instrument that do not
        carry its resources.

        Args:
            resources (Iterable[InstrumentResources | None]): The resources.

        Returns:
            InstrumentResources: The totals.
        """
        given = [r for r in resources if r is not None]
        totals = {}
        for field in cls.model_fields:  # pylint: disable=not-an-iterable
            values = [getattr(r, field) for r in given if getattr(r, field) is not None]
            totals[field] = float(sum(values)) if values else None
        return cls(**totals)


class SpacecraftResources(BaseModel):
    """
    The mass, power, life, and data handling of a spacecraft, from which
    to size systems and their data paths. A value that is not given is `None`.
    """

    model_config = ConfigDict(extra="forbid")

    launch_mass: float | None = Field(
        default=None,
        gt=0,
        description="The mass (kg) at launch, with propellant.",
        examples=[2294],
    )
    dry_mass: float | None = Field(
        default=None,
        gt=0,
        description="The mass (kg) without propellant.",
        examples=[1929],
    )
    power: float | None = Field(
        default=None,
        gt=0,
        description="The power (watts) of the spacecraft, as published (such as "
        + "the end-of-life generation of its solar array).",
        examples=[1932],
    )
    design_life: float | None = Field(
        default=None,
        gt=0,
        description="The design life (years).",
        examples=[7],
    )
    storage: float | None = Field(
        default=None,
        gt=0,
        description="The capacity (bits) of its onboard data storage.",
        examples=[3.4e12],
    )
    downlink_rate: float | None = Field(
        default=None,
        gt=0,
        description="The rate (bits per second) of its principal downlink of "
        + "mission data.",
        examples=[300e6],
    )
    drag_area: float | None = Field(
        default=None,
        gt=0,
        description="The mean area (square meters) it presents to the flow of "
        + "air along its orbit, by which drag decays it (see "
        + "`tatc.utils.compute_orbit_lifetime`).",
        examples=[10.0],
    )
    drag_coefficient: float | None = Field(
        default=None,
        gt=0,
        description="The drag coefficient of its drag area, or `None` for "
        + "`tatc.utils.DEFAULT_DRAG_COEFFICIENT` (2.2).",
        examples=[2.2],
    )
    specific_impulse: float | None = Field(
        default=None,
        gt=0,
        description="The specific impulse (seconds) of its propulsion, by which "
        + "propellant gives velocity (see `tatc.sizing.compute_propellant_mass`), "
        + "or `None` if not given.",
        examples=[220],
    )

    @model_validator(mode="after")
    def _check_masses(self) -> SpacecraftResources:
        if (
            self.launch_mass is not None
            and self.dry_mass is not None
            and self.dry_mass > self.launch_mass
        ):
            raise ValueError("dry_mass must not exceed launch_mass")
        return self

    def get_ballistic_coefficient(self, mass: float | None = None) -> float:
        """
        Gets the ballistic coefficient (kg/m^2) by which drag decays the
        spacecraft's orbit: its mass over its drag coefficient times its drag
        area. By default, its mass is its dry mass (as at the end of its
        life, when its orbit decays without propulsion), or else its launch
        mass.

        Args:
            mass (float | None): The mass (kg), or `None` for the dry mass,
                or else the launch mass.

        Returns:
            float: The ballistic coefficient (kg/m^2).

        Raises:
            ValueError: If the drag area, or a mass, is not given.
        """
        if self.drag_area is None:
            raise ValueError("the drag area is not given")
        if mass is None:
            mass = self.dry_mass if self.dry_mass is not None else self.launch_mass
        if mass is None:
            raise ValueError("the mass is not given")
        coefficient = (
            DEFAULT_DRAG_COEFFICIENT
            if self.drag_coefficient is None
            else self.drag_coefficient
        )
        return mass / (coefficient * self.drag_area)

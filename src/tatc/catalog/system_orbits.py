"""
The orbits of the catalog's satellite systems (see `tatc.catalog.system_schemas`):
the planes of orbits at the plane grade, and the geostationary slots and repeat
ground tracks (`tatc.schemas.RepeatGroundTrackOrbit`) at the ground track grade.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import math
from datetime import datetime, time, timezone
from typing import Annotated, Any, Literal

from pydantic import AwareDatetime, Field, model_validator

from ..schemas import (
    CircularOrbit,
    GeosynchronousOrbit,
    RepeatGroundTrackOrbit,
    SunSynchronousOrbit,
)
from .records import _CatalogModel

SLOT_EPOCH = datetime(2026, 1, 1, tzinfo=timezone.utc)
"""
The epoch of a geostationary orbit built without one: as it is maintained at its
longitude, any epoch gives the same subpoint.
"""

# the optional values of orbits at the plane grade, which a configuration gives
# for every orbit or for none, relative to its other orbits
RELATIVE_VALUES = ("phase", "right_ascension_ascending_node")


class SunSynchronousPlane(_CatalogModel):
    """
    The plane of a sun-synchronous orbit, at the plane grade: its mean altitude
    and the local time of its ascending node, and, for an orbit whose local time
    is no longer maintained (and drifts), its inclination and the time at which
    its local time is cataloged. Built at an epoch as a `SunSynchronousOrbit`
    (whose inclination is sun-synchronous, so that its local time stays
    constant), or, with an inclination, as a `CircularOrbit` whose node has
    drifted from the cataloged local time at the rate of its inclination from
    the reference time to the epoch, and drifts on
    during an analysis (at SGP4's secular rate, with which TAT-C propagates it,
    so that building it at any epoch gives the same orbit). Its phase is that of
    a constellation whose operator
    maintains the satellites' relative phases, if given, or else spaced evenly
    among the satellites built (see `OrbitConfiguration.build_orbits`).
    """

    type: Literal["sso"] = Field(default="sso", description="Orbit discriminator.")
    mean_altitude: float = Field(
        ..., gt=0, description="Mean altitude (meters), as TAT-C defines it."
    )
    equator_crossing_time: time = Field(
        ...,
        description="Local mean solar time of the ascending node (at the "
        + "`reference_time` of an orbit whose local time drifts).",
    )
    inclination: float | None = Field(
        default=None,
        gt=0,
        lt=180,
        description="Inclination (degrees) of an orbit whose local time drifts, or "
        + "None (the default) for the sun-synchronous inclination.",
    )
    reference_time: AwareDatetime | None = Field(
        default=None,
        description="The time at which the local time of an orbit whose local time "
        + "drifts is cataloged (the epoch of the element set to which it was "
        + "fitted), given with its inclination.",
    )
    phase: float | None = Field(
        default=None,
        ge=0,
        lt=360,
        description="Argument of latitude (degrees) at the epoch, relative to the "
        + "other satellites of the configuration, for a constellation whose "
        + "operator maintains their relative phases, or None (the default) for an "
        + "arbitrary phase.",
    )

    @model_validator(mode="after")
    def _check_drift(self) -> SunSynchronousPlane:
        if (self.inclination is None) != (self.reference_time is None):
            raise ValueError(
                "an orbit whose local time drifts has an inclination and a "
                "reference_time, and a sun-synchronous orbit neither"
            )
        return self

    def build(
        self,
        epoch: datetime | None = None,
        phase: float = 0,
        right_ascension_ascending_node: float = 0,
    ) -> SunSynchronousOrbit | CircularOrbit:
        """
        Builds the orbit at an epoch.

        Args:
            epoch (datetime | None): The epoch.
            phase (float): The argument of latitude (degrees) at the epoch.
            right_ascension_ascending_node (float): Not used, as the local time
                sets the node.

        Returns:
            SunSynchronousOrbit | CircularOrbit: The orbit.

        Raises:
            ValueError: If no epoch is given.
        """
        del right_ascension_ascending_node
        epoch = _require_epoch(epoch)
        if self.inclination is None:
            return SunSynchronousOrbit(
                mean_altitude=self.mean_altitude,
                equator_crossing_time=self.equator_crossing_time,
                epoch=epoch,
                true_anomaly=phase,
            )
        cataloged = SunSynchronousOrbit(
            mean_altitude=self.mean_altitude,
            equator_crossing_time=self.equator_crossing_time,
            epoch=self.reference_time,
        ).get_right_ascension_ascending_node()
        drifting = CircularOrbit(
            mean_altitude=self.mean_altitude,
            inclination=self.inclination,
            right_ascension_ascending_node=cataloged,
            epoch=self.reference_time,
        )
        # SGP4's secular rate of the node (radians per minute)
        rate = drifting.to_gp_orbit().elements[0].to_satrec().nodedot
        minutes = (epoch - self.reference_time).total_seconds() / 60
        node = cataloged + math.degrees(rate * minutes)
        return drifting.model_copy(
            update={
                "right_ascension_ascending_node": float(node % 360),
                "epoch": epoch,
                "true_anomaly": phase,
            }
        )

    def describe(self) -> str:
        """
        Describes the orbit.

        Returns:
            str: The description.
        """
        text = (
            f"Sun-synchronous, {self.mean_altitude / 1e3:.2f} km, ascending node "
            f"at {self.equator_crossing_time:%H:%M:%S}"
        )
        if self.inclination is not None:
            text += (
                f" on {self.reference_time:%Y-%m-%d}, drifting at an inclination "
                f"of {self.inclination} deg"
            )
        if self.phase is not None:
            text += f", phase {self.phase:g} deg"
        return text


class InclinedPlane(_CatalogModel):
    """
    The plane of an inclined circular orbit, at the plane grade: its mean
    altitude and inclination. Built as a `CircularOrbit` at an epoch. Its right
    ascension of ascending node and phase are those of a constellation whose
    operator maintains the satellites' relative planes and phases, if given, or
    else spaced evenly among the satellites built (see
    `OrbitConfiguration.build_orbits`).
    """

    type: Literal["circular"] = Field(
        default="circular", description="Orbit discriminator."
    )
    mean_altitude: float = Field(
        ..., gt=0, description="Mean altitude (meters), as TAT-C defines it."
    )
    inclination: float = Field(..., gt=0, lt=180, description="Inclination (degrees).")
    right_ascension_ascending_node: float | None = Field(
        default=None,
        ge=0,
        lt=360,
        description="Right ascension of ascending node (degrees) at the epoch, "
        + "relative to the other satellites of the configuration, or None (the "
        + "default) for an arbitrary one.",
    )
    phase: float | None = Field(
        default=None,
        ge=0,
        lt=360,
        description="Argument of latitude (degrees) at the epoch, relative to the "
        + "other satellites of the configuration, or None (the default) for an "
        + "arbitrary phase.",
    )

    @model_validator(mode="after")
    def _check_relative(self) -> InclinedPlane:
        if self.phase is not None and self.right_ascension_ascending_node is None:
            raise ValueError(
                "a relative phase is of a plane relative to the others: give the "
                "right_ascension_ascending_node too"
            )
        return self

    def build(
        self,
        epoch: datetime | None = None,
        phase: float = 0,
        right_ascension_ascending_node: float = 0,
    ) -> CircularOrbit:
        """
        Builds the orbit at an epoch.

        Args:
            epoch (datetime | None): The epoch.
            phase (float): The argument of latitude (degrees) at the epoch.
            right_ascension_ascending_node (float): The right ascension of
                ascending node (degrees) at the epoch.

        Returns:
            CircularOrbit: The orbit.

        Raises:
            ValueError: If no epoch is given.
        """
        epoch = _require_epoch(epoch)
        return CircularOrbit(
            mean_altitude=self.mean_altitude,
            inclination=self.inclination,
            right_ascension_ascending_node=right_ascension_ascending_node,
            epoch=epoch,
            true_anomaly=phase,
        )

    def describe(self) -> str:
        """
        Describes the orbit.

        Returns:
            str: The description.
        """
        text = f"{self.mean_altitude / 1e3:.2f} km, inclination {self.inclination} deg"
        if self.right_ascension_ascending_node is not None:
            text += f", node {self.right_ascension_ascending_node:g} deg"
        if self.phase is not None:
            text += f", phase {self.phase:g} deg"
        return text


class GeostationarySlot(_CatalogModel):
    """
    The longitude of a geostationary orbit, at the ground track grade. Built as
    a `GeosynchronousOrbit` at an epoch, maintained at its longitude.
    """

    type: Literal["geosynchronous"] = Field(
        default="geosynchronous", description="Orbit discriminator."
    )
    longitude: float = Field(
        ..., ge=-180, le=180, description="Longitude (degrees) of the subpoint."
    )
    inclination: float = Field(
        default=0, ge=0, lt=180, description="Inclination (degrees)."
    )

    def build(
        self,
        epoch: datetime | None = None,
        phase: float = 0,
        right_ascension_ascending_node: float = 0,
    ) -> GeosynchronousOrbit:
        """
        Builds the orbit at an epoch.

        Args:
            epoch (datetime | None): The epoch (by default, `SLOT_EPOCH`: as the
                orbit is maintained at its longitude, any epoch gives the same
                subpoint).
            phase (float): Not used, as the longitude sets the position.
            right_ascension_ascending_node (float): Not used.

        Returns:
            GeosynchronousOrbit: The orbit.
        """
        del phase, right_ascension_ascending_node
        return GeosynchronousOrbit(
            longitude=self.longitude,
            inclination=self.inclination,
            epoch=SLOT_EPOCH if epoch is None else epoch,
        )

    def describe(self) -> str:
        """
        Describes the orbit.

        Returns:
            str: The description.
        """
        return (
            f"Geostationary at longitude {self.longitude} deg, inclination "
            f"{self.inclination} deg"
        )


class RepeatGroundTrack(RepeatGroundTrackOrbit):
    """
    An orbit maintained on a repeat ground track, at the ground track grade: a
    `RepeatGroundTrackOrbit` pinned to an ascending node, built for any epoch.
    """

    def build(
        self,
        epoch: datetime | None = None,
        phase: float = 0,
        right_ascension_ascending_node: float = 0,
    ) -> RepeatGroundTrackOrbit:
        """
        Builds the orbit, which is the same at any epoch.

        Args:
            epoch (datetime | None): Not used, as the node pins the orbit.
            phase (float): Not used.
            right_ascension_ascending_node (float): Not used.

        Returns:
            RepeatGroundTrackOrbit: The orbit.
        """
        del epoch, phase, right_ascension_ascending_node
        return RepeatGroundTrackOrbit(**self.model_dump())

    def describe(self) -> str:
        """
        Describes the orbit.

        Returns:
            str: The description.
        """
        node = self.ascending_node
        return (
            f"Repeat of {self.revolutions} orbits in {self.days} days "
            f"({self.get_mean_altitude() / 1e3:.2f} km), ascending node at "
            f"longitude {node.longitude} deg on {node.time:%Y-%m-%d %H:%M:%S} UTC"
        )


CatalogOrbit = Annotated[
    SunSynchronousPlane | InclinedPlane | GeostationarySlot | RepeatGroundTrack,
    Field(discriminator="type"),
]
"""An orbit of a catalog system, by its `type`."""


def describe_orbit(orbit: CatalogOrbit) -> str:
    """
    Describes an orbit of a configuration, such as "Sun-synchronous, 834.10 km,
    ascending node at 13:30:03".

    Args:
        orbit (CatalogOrbit): The orbit.

    Returns:
        str: The description.
    """
    return orbit.describe()


def _require_epoch(epoch: datetime | None) -> datetime:
    """
    Gets the epoch of an orbit at the plane grade, which needs one.
    """
    if epoch is None:
        raise ValueError("an orbit at the plane grade is built at an epoch")
    return epoch


def get_orbit_values(orbit: Any) -> set[str]:
    """
    Gets the values of an orbit that need provenance: its fields, without its
    `type` or unset optional values (such as a derived inclination).

    Args:
        orbit (CatalogOrbit): The orbit.

    Returns:
        set[str]: The names of the fields.
    """
    return {field for field, value in orbit if field != "type" and value is not None}

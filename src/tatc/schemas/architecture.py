"""
Object schemas for mission architectures: observing platforms of every kind,
the ground stations and relays that receive their data, and the windows of
time in which each is active.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from datetime import datetime, timedelta, timezone
from typing import TYPE_CHECKING, Any

from pydantic import AwareDatetime, BaseModel, Field, field_validator, model_validator

from .._warnings import warn
from .air import Aircraft
from .communication import RelayCommunications
from .space import Satellite
from .space.base_constellation import BaseConstellation
from .surface import GroundStation, RadarStation, SurfacePlatform

if TYPE_CHECKING:  # pragma: no cover
    from ..analysis.architecture import ArchitectureEvaluation


def _as_utc(time: datetime) -> datetime:
    """Gets a time with a time zone (UTC if it has none), for comparisons."""
    return time if time.tzinfo is not None else time.replace(tzinfo=timezone.utc)


def _expand_constellations(values: Any) -> Any:
    """
    Expands the constellations of a list of satellites into their members
    (see `BaseConstellation.generate_members`).
    """
    if not isinstance(values, (list, tuple)):
        return values
    expanded = []
    for value in values:
        if isinstance(value, BaseConstellation):
            expanded.extend(value.generate_members())
        else:
            expanded.append(value)
    return expanded


class ActiveWindow(BaseModel):
    """
    A window of time in which a member of an architecture is active (for
    example, between a satellite's launch and its retirement, or after a
    ground station is built): open-ended without a `start` or an `end`.
    """

    start: AwareDatetime | None = Field(
        None, description="The start of the window, or `None` for none."
    )
    end: AwareDatetime | None = Field(
        None, description="The end of the window, or `None` for none."
    )

    @model_validator(mode="after")
    def _check_order(self) -> ActiveWindow:
        if self.start is not None and self.end is not None and self.end <= self.start:
            raise ValueError(
                f"the window ends ({self.end}) before it starts ({self.start})"
            )
        return self

    def overlap(
        self, start: datetime, end: datetime
    ) -> tuple[datetime, datetime] | None:
        """
        Gets the part of a period within the window.

        Args:
            start (datetime): The start of the period.
            end (datetime): The end of the period.

        Returns:
            tuple[datetime, datetime] | None: the start and end of the part
                within the window, or `None` if none is (or only an instant)
        """
        first = start if self.start is None else max(_as_utc(start), self.start)
        last = end if self.end is None else min(_as_utc(end), self.end)
        if _as_utc(last) <= _as_utc(first):
            return None
        return first, last

    def contains(self, time: datetime) -> bool:
        """
        Checks whether a time is within the window (from its start, inclusive,
        to its end, exclusive).

        Args:
            time (datetime): The time.

        Returns:
            bool: `True` if the time is within the window
        """
        time = _as_utc(time)
        return (self.start is None or self.start <= time) and (
            self.end is None or time < self.end
        )


class Relay(BaseModel):
    """
    A data relay (such as a constellation of tracking and data relay
    satellites in geostationary orbits): a satellite of an architecture
    contacts the relay while it has a line of sight to one of the relay's
    satellites that grazes the Earth no lower than `min_grazing_altitude`
    (and within `max_range`, if given). With `ground_stations`, the relay's
    satellite must also be in view of one of them at the same time (its
    elevation angle above the station's minimum); without them, each relay
    satellite is assumed to be always in contact with the ground (as one in
    a geostationary orbit is with a station that it serves). Contacts act as
    downlinks (see `tatc.analysis.collect_relay_contacts`).

    The relay's satellites carry no instruments of interest: they are not
    observing platforms of the architecture.
    """

    name: str = Field(..., description="Relay name.", examples=["relay 1"])
    satellites: list[Satellite] = Field(
        ..., min_length=1, description="The relay satellites."
    )
    ground_stations: list[GroundStation] = Field(
        [],
        description="The ground stations that receive the relay's data, or none "
        + "if each relay satellite is always in contact with the ground.",
    )
    min_grazing_altitude: float = Field(
        0,
        description="The minimum altitude (meters) above the WGS 84 ellipsoid of "
        + "the line of sight between a satellite and a relay satellite.",
    )
    max_range: float | None = Field(
        None,
        gt=0,
        description="The maximum range (meters) between a satellite and a relay "
        + "satellite, or `None` for none.",
    )
    min_access_time: timedelta = Field(
        timedelta(0),
        description="The minimum duration of a contact.",
        examples=[timedelta(seconds=60)],
    )
    communications: RelayCommunications | None = Field(
        None,
        description="The antennas, receivers, and feeder link of each relay "
        + "satellite, or `None` for any number of antennas without a link "
        + "budget (see `tatc.analysis.compute_data_flow`).",
    )

    _expand = field_validator("satellites", mode="before")(_expand_constellations)


ObservingPlatform = Satellite | Aircraft | SurfacePlatform | RadarStation
"""An observing platform of an architecture."""

Member = Satellite | Aircraft | SurfacePlatform | RadarStation | GroundStation | Relay
"""A member of an architecture: an observing platform, a ground station, or
a relay."""

_FIELDS: dict[type, str] = {
    Satellite: "satellites",
    Aircraft: "aircraft",
    SurfacePlatform: "surface_platforms",
    RadarStation: "radar_stations",
    GroundStation: "ground_stations",
    Relay: "relays",
}
"""The field of each type of member, in the order of the fields (and of the
kinds of platforms of the analyses, see `tatc.analysis.check._PLATFORM_KINDS`)."""

_PLATFORM_FIELDS = ("satellites", "aircraft", "surface_platforms", "radar_stations")
"""The fields of observing platforms."""


def _field_of(member: object) -> str:
    """
    Gets the field of an architecture that holds a member, raising a
    `TypeError` if it is not a member.
    """
    # (radar stations and surface platforms are points, as ground stations
    # are, so the most specific type is found first by exact type, then by
    # subclass)
    for kind, field in _FIELDS.items():
        if type(member) is kind:  # pylint: disable=unidiomatic-typecheck
            return field
    for kind, field in _FIELDS.items():
        if isinstance(member, kind):
            return field
    names = ", ".join(kind.__name__ for kind in _FIELDS)
    raise TypeError(
        f"members of an architecture are a {names}, not a {type(member).__name__}"
    )


class Architecture(BaseModel):
    """
    Mission architecture: observing platforms of every kind (satellites,
    aircraft, fixed surface platforms, and radar stations), the ground
    stations and relays that receive the data of its satellites, and the
    windows in which its members are active.

    Every member has a distinct name, by which `windows` refers to it: a
    member without a window is always active. The analyses of an
    architecture (see `evaluate`) analyze each member only over the part of
    the period within its window: a platform observes, and a satellite
    contacts ground stations and relays, only while active (and a contact
    with a ground station or a relay only while it is active, too).
    """

    name: str = Field("architecture", description="Name of this architecture.")
    satellites: list[Satellite] = Field(
        [],
        description="Observing satellites (the members of constellations given "
        + "are expanded, see `BaseConstellation.generate_members`).",
    )
    aircraft: list[Aircraft] = Field([], description="Observing aircraft.")
    surface_platforms: list[SurfacePlatform] = Field(
        [], description="Observing fixed surface platforms."
    )
    radar_stations: list[RadarStation] = Field([], description="Radar stations.")
    ground_stations: list[GroundStation] = Field(
        [], description="Ground stations that receive the data of the satellites."
    )
    relays: list[Relay] = Field(
        [], description="Relays that receive the data of the satellites."
    )
    windows: dict[str, ActiveWindow] = Field(
        {},
        description="The window in which each member (by name) is active, if "
        + "not always.",
    )

    _expand = field_validator("satellites", mode="before")(_expand_constellations)

    @model_validator(mode="before")
    @classmethod
    def _rename_stations(cls, data: Any) -> Any:
        if isinstance(data, dict) and "stations" in data:
            if "ground_stations" in data:
                raise ValueError("give ground_stations, not also stations")
            warn(
                "the stations of an architecture are now ground_stations; "
                + "stations will be removed",
                DeprecationWarning,
            )
            data = {
                ("ground_stations" if key == "stations" else key): value
                for key, value in data.items()
            }
        return data

    @model_validator(mode="after")
    def _check_names(self) -> Architecture:
        names = [member.name for member in self.members]
        repeated = sorted({name for name in names if names.count(name) > 1})
        if repeated:
            raise ValueError(
                f"the members of an architecture have distinct names, not {repeated}"
            )
        unknown = sorted(set(self.windows) - set(names))
        if unknown:
            raise ValueError(
                f"windows of members that are not in the architecture: {unknown}"
            )
        return self

    @property
    def platforms(self) -> list[ObservingPlatform]:
        """list: The observing platforms: the satellites, aircraft, surface
        platforms, and radar stations, in that order."""
        return [member for field in _PLATFORM_FIELDS for member in getattr(self, field)]

    @property
    def members(self) -> list[Member]:
        """list: The members: the observing platforms (see `platforms`),
        ground stations, and relays, in that order."""
        return [*self.platforms, *self.ground_stations, *self.relays]

    def get_member(self, name: str) -> Member:
        """
        Gets a member by name.

        Args:
            name (str): The name.

        Returns:
            Member: the member

        Raises:
            KeyError: if no member has the name.
        """
        for member in self.members:
            if member.name == name:
                return member
        raise KeyError(f"no member of the architecture is named {name!r}")

    def get_window(self, name: str) -> ActiveWindow:
        """
        Gets the window in which a member is active (unbounded if it has none).

        Args:
            name (str): The name of the member.

        Returns:
            ActiveWindow: the window
        """
        return self.windows.get(name, ActiveWindow())

    def get_active_period(
        self, name: str, start: datetime, end: datetime
    ) -> tuple[datetime, datetime] | None:
        """
        Gets the part of a period in which a member is active (see
        `ActiveWindow.overlap`).

        Args:
            name (str): The name of the member.
            start (datetime): The start of the period.
            end (datetime): The end of the period.

        Returns:
            tuple[datetime, datetime] | None: the part of the period, or
                `None` if the member is not active in it
        """
        return self.get_window(name).overlap(start, end)

    def is_active(self, name: str, time: datetime) -> bool:
        """
        Checks whether a member is active at a time.

        Args:
            name (str): The name of the member.
            time (datetime): The time.

        Returns:
            bool: `True` if the member is active
        """
        return self.get_window(name).contains(time)

    def with_members(
        self,
        *members: Member | Iterable[Member],
        windows: Mapping[str, ActiveWindow] | None = None,
        name: str | None = None,
    ) -> Architecture:
        """
        Gets a copy of the architecture with other members (each in the field
        of its kind), such as to measure the value of a platform it lacks.

        Args:
            *members (Member | Iterable[Member]): The members (or lists of
                them, such as a constellation's satellites).
            windows (Mapping[str, ActiveWindow] | None): The windows of the
                new members, if any.
            name (str | None): The name of the copy, by default the same.

        Returns:
            Architecture: the copy
        """
        fields = {field: list(getattr(self, field)) for field in _FIELDS.values()}
        for member in _flatten(members):
            fields[_field_of(member)].append(member)
        return Architecture(
            name=self.name if name is None else name,
            **fields,
            windows={**self.windows, **(windows or {})},
        )

    def without_members(
        self, *names: str | Iterable[str], name: str | None = None
    ) -> Architecture:
        """
        Gets a copy of the architecture without members (and their windows),
        such as to measure the value of a member it has.

        Args:
            *names (str | Iterable[str]): The names of the members (or lists
                of them).
            name (str | None): The name of the copy, by default the same.

        Returns:
            Architecture: the copy

        Raises:
            KeyError: if a name is not of a member.
        """
        removed = set(_flatten_names(names))
        unknown = sorted(removed - {member.name for member in self.members})
        if unknown:
            raise KeyError(f"members that are not in the architecture: {unknown}")
        return Architecture(
            name=self.name if name is None else name,
            **{
                field: [m for m in getattr(self, field) if m.name not in removed]
                for field in _FIELDS.values()
            },
            windows={k: v for k, v in self.windows.items() if k not in removed},
        )

    def merge(self, other: Architecture, name: str | None = None) -> Architecture:
        """
        Gets an architecture of the members (and windows) of this one and
        another, such as a design's and a baseline's.

        Args:
            other (Architecture): The other architecture.
            name (str | None): The name of the result, by default this one's.

        Returns:
            Architecture: the merged architecture
        """
        return self.with_members(other.members, windows=other.windows, name=name)

    def evaluate(
        self,
        targets: Any,
        start: datetime,
        end: datetime,
        metrics: Callable[[Any], Mapping[str, Any]] | None = None,
        **options: Any,
    ) -> ArchitectureEvaluation:
        """
        Evaluates the architecture: the observations of the targets by its
        platforms, the contacts of its satellites with its ground stations
        and relays, the latency of each observation, and its metrics (see
        `tatc.analysis.evaluate_architecture`).

        Args:
            targets (Any): The points and regions of interest.
            start (datetime): The start of the period.
            end (datetime): The end of the period.
            metrics (Callable | None): A function of the observations that
                gives their metrics (a mapping), if any.
            **options (Any): Other options of
                `tatc.analysis.evaluate_architecture`, such as `time_step`.

        Returns:
            tatc.analysis.ArchitectureEvaluation: the evaluation
        """
        # pylint: disable-next=import-outside-toplevel
        from ..analysis.architecture import evaluate_architecture

        return evaluate_architecture(self, targets, start, end, metrics, **options)


def _flatten(values: Iterable[Any]) -> list[Any]:
    """Flattens members given alone or in lists (but not as models)."""
    flat = []
    for value in values:
        if isinstance(value, BaseModel):
            flat.append(value)
        else:
            flat.extend(value)
    return flat


def _flatten_names(values: Iterable[Any]) -> list[str]:
    """Flattens names given alone or in lists."""
    flat = []
    for value in values:
        if isinstance(value, str):
            flat.append(value)
        else:
            flat.extend(value)
    return flat

"""
Object schema for general perturbations orbital elements.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import csv
import json
from datetime import datetime, timedelta, timezone

import numpy as np
from pydantic import BaseModel, Field
from sgp4 import exporter, omm
from sgp4.api import WGS72, Satrec
from sgp4.conveniences import sat_epoch_datetime
from skyfield.api import EarthSatellite

from ... import constants, utils
from ...utils.cache import get_cached
from ...utils.propagation import (
    RepeatCycleSearch,
    _compute_repeat_element,
    _search_repeat_cycle,
)

SUN_SYNCHRONOUS_NODAL_DAY_TOLERANCE_S = 10
"""Maximum difference (seconds) between the nodal day of a sun-synchronous orbit and a mean solar day."""


class GeneralPerturbationsElements(BaseModel):
    """General perturbations orbital elements for a satellite."""

    object_name: str | None = Field(default=None, description="Object name.")
    epoch: datetime = Field(..., description="Epoch.")
    mean_motion: float = Field(..., description="Mean motion (degrees/second).", gt=0)
    eccentricity: float = Field(..., description="Eccentricity.", ge=0, le=1)
    inclination: float = Field(..., description="Inclination (degrees).", ge=0, le=180)
    ra_of_asc_node: float = Field(
        ..., description="Right ascension of ascending node (degrees).", ge=0, lt=360
    )
    arg_of_pericenter: float = Field(
        ..., description="Argument of pericenter (degrees).", ge=0, lt=360
    )
    mean_anomaly: float = Field(
        ..., description="Mean anomaly (degrees).", ge=0, lt=360
    )
    norad_cat_id: int = Field(default=0, description="NORAD catalog identifier.", ge=0)
    bstar: float = Field(default=0, description="Starred ballistic coefficient.")
    mean_motion_dot: float = Field(
        default=0, description="First derivative of mean motion (degrees/second^2)."
    )
    mean_motion_ddot: float = Field(
        default=0, description="Second derivative of mean motion (degrees/second^3)."
    )
    classification: str = Field(default="U", description="Classification type.")
    international_designator: str = Field(
        default="00000A", description="International designator."
    )
    ephemeris_type: int = Field(default=0, description="Ephemeris type.")
    element_set_num: int = Field(default=0, description="Element set number.")
    revolution_num: int = Field(default=0, description="Revolution number at epoch.")

    def _get_key(self) -> tuple:
        """
        Gets this element's field values, which key the values cached on it.

        Returns:
            tuple: the field values
        """
        # pylint: disable-next=not-an-iterable
        return tuple(getattr(self, name) for name in type(self).model_fields)

    @property
    def has_drag(self) -> bool:
        """True, if any drag term (B* or a derivative of mean motion) is nonzero."""
        return (
            self.bstar != 0 or self.mean_motion_dot != 0 or self.mean_motion_ddot != 0
        )

    def without_drag(self) -> GeneralPerturbationsElements:
        """
        Gets this element without drag, as for an orbit maintained against
        drag: its drag terms (B* and the derivatives of mean motion) are set
        to zero, while the Earth's oblateness (and, for deep-space orbits,
        the Moon and Sun) still perturb the orbit.

        Returns:
            GeneralPerturbationsElements: a copy without drag (or this
                element, if it has none)
        """
        if not self.has_drag:
            return self
        return self.model_copy(
            update={"bstar": 0, "mean_motion_dot": 0, "mean_motion_ddot": 0}
        )

    @classmethod
    def from_satrec(cls, satrec: Satrec) -> GeneralPerturbationsElements:
        """
        Creates a GP elements object from a Satrec object.

        Args:
            satrec (Satrec): The Satrec object.

        Returns:
            GeneralPerturbationsElements: the GP elements
        """
        return GeneralPerturbationsElements(
            epoch=sat_epoch_datetime(satrec),
            mean_motion=np.degrees(satrec.no_kozai) / 60,
            eccentricity=satrec.ecco,
            inclination=np.degrees(satrec.inclo),
            ra_of_asc_node=np.degrees(satrec.nodeo),
            arg_of_pericenter=np.degrees(satrec.argpo),
            mean_anomaly=np.degrees(satrec.mo),
            norad_cat_id=satrec.satnum,
            bstar=satrec.bstar,
            mean_motion_dot=np.degrees(satrec.ndot) / 60**2,
            mean_motion_ddot=np.degrees(satrec.nddot) / 60**3,
            classification=satrec.classification,
            international_designator=satrec.intldesg,
            ephemeris_type=satrec.ephtype,
            element_set_num=satrec.elnum,
            revolution_num=satrec.revnum,
        )

    def to_satrec(self) -> Satrec:
        """
        Converts this GP elements object to a Satrec object.

        Returns:
            Satrec: the Satrec object
        """
        satrec = Satrec()
        satrec.classification = self.classification
        satrec.intldesg = self.international_designator
        satrec.ephtype = self.ephemeris_type
        satrec.elnum = self.element_set_num
        satrec.revnum = self.revolution_num
        satrec.sgp4init(
            WGS72,
            "i",
            self.norad_cat_id,
            (self.epoch - datetime(1949, 12, 31, tzinfo=timezone.utc))
            / timedelta(days=1),
            self.bstar,
            np.radians(self.mean_motion_dot) * 60**2,
            np.radians(self.mean_motion_ddot) * 60**3,
            self.eccentricity,
            np.radians(self.arg_of_pericenter),
            np.radians(self.inclination),
            np.radians(self.mean_anomaly),
            np.radians(self.mean_motion) * 60,
            np.radians(self.ra_of_asc_node),
        )
        return satrec

    def get_orbit_period(self) -> timedelta:
        """
        Gets the approximate orbit period.

        Returns:
            timedelta: the orbit period
        """
        return timedelta(
            seconds=utils.orbital.mean_motion_to_orbit_period(self.mean_motion)
        )

    def get_semimajor_axis(self) -> float:
        """
        Gets the semimajor axis.

        Returns:
            float: the semimajor axis (meters)
        """

        return utils.orbital.mean_motion_to_semimajor_axis(self.mean_motion)

    def get_mean_altitude(self) -> float:
        """
        Gets the mean altitude.

        Returns:
            float: the mean altitude (meters)
        """
        return self.get_semimajor_axis() - constants.EARTH_MEAN_RADIUS

    def get_true_anomaly(self) -> float:
        """
        Gets the true anomaly.

        Returns:
            float: the true anomaly (degrees)
        """
        return utils.orbital.mean_anomaly_to_true_anomaly(
            self.mean_anomaly, self.eccentricity
        )

    @classmethod
    def from_tle(cls, tle_lines: tuple[str, str]) -> GeneralPerturbationsElements:
        """
        Creates a GP elements object from two line element (TLE) lines.

        Args:
            tle_lines (tuple[str, str]): The two TLE lines.

        Returns:
            GeneralPerturbationsElements: the GP elements
        """
        return GeneralPerturbationsElements.from_satrec(
            Satrec.twoline2rv(tle_lines[0], tle_lines[1])
        )

    def to_tle(self) -> tuple[str, str]:
        """
        Converts this GP elements object to a two line element (TLE) representation.

        Returns:
            tuple[str, str]: the two line elements
        """
        return exporter.export_tle(self.to_satrec())

    @classmethod
    def from_omm_dict(cls, omm_dict: dict) -> GeneralPerturbationsElements:
        """
        Creates a GP elements object from an OMM dictionary.

        Args:
            omm_dict (dict): The OMM dictionary.

        Returns:
            GeneralPerturbationsElements: the GP elements
        """
        satrec = Satrec()
        omm.initialize(satrec, omm_dict)
        elements = GeneralPerturbationsElements.from_satrec(satrec)
        # object_name has no equivalent on Satrec, so from_satrec can never
        # recover it; restore it directly from the OMM dictionary
        return elements.model_copy(update={"object_name": omm_dict.get("OBJECT_NAME")})

    def to_omm_dict(self) -> dict:
        """
        Converts this GP elements object to an OMM dictionary.

        Returns:
            dict: the OMM dictionary
        """
        return exporter.export_omm(self.to_satrec(), self.object_name)

    @classmethod
    def from_omm_csv(cls, omm_csv: list[str]) -> GeneralPerturbationsElements:
        """
        Creates a GP elements object from OMM CSV lines. Only the first
        data row is used; all subsequent rows are ignored.

        Args:
            omm_csv (list[str]): The OMM CSV lines, including a header row.

        Returns:
            GeneralPerturbationsElements: the GP elements
        """
        for fields in csv.DictReader(omm_csv):
            return GeneralPerturbationsElements.from_omm_dict(fields)
        raise ValueError("No OMM CSV lines found.")

    @classmethod
    def from_omm_json(cls, omm_json: str) -> GeneralPerturbationsElements:
        """
        Creates a GP elements object from an OMM JSON string. Only the
        first entry in the JSON array is used; all subsequent entries
        are ignored.

        Args:
            omm_json (str): The OMM JSON string, encoding a list of OMM
                records.

        Returns:
            GeneralPerturbationsElements: the GP elements
        """
        for fields in json.loads(omm_json):
            return GeneralPerturbationsElements.from_omm_dict(fields)
        raise ValueError("No OMM JSON lines found.")

    def to_skyfield(self, remove_drag: bool = False) -> EarthSatellite:
        """
        Converts this GP elements object to a Skyfield `EarthSatellite`,
        which can be used to propagate this orbital state via SGP4. The
        satellite is cached.

        Args:
            remove_drag (bool): True, to propagate without drag (see `without_drag`).

        Returns:
            skyfield.api.EarthSatellite: the Skyfield EarthSatellite
        """
        element = self.without_drag() if remove_drag else self
        return get_cached(
            self,
            "skyfield_without_drag" if remove_drag else "skyfield",
            self._get_key(),
            lambda: EarthSatellite.from_omm(constants.timescale, element.to_omm_dict()),
        )

    def get_nodal_period_and_day(self) -> tuple[float, float]:
        """
        Gets the nodal period (between ascending node crossings) and the
        nodal day (the period of the Earth's rotation relative to the
        precessing ascending node), from the SGP4 model's own secular rates
        for mean anomaly, argument of perigee, and right ascension of
        ascending node (mdot, argpdot, and nodedot, in radians/minute).

        Returns:
            tuple[float, float]: the nodal period and nodal day (seconds)
        """
        model = self.to_satrec()
        nodal_period = 2 * np.pi / (model.mdot + model.argpdot) * 60
        # Earth's rotation rate (radians/minute), as SGP4's rates
        earth_rotation_rate = constants.EARTH_ROTATION_RATE * 60
        nodal_day = 2 * np.pi / (earth_rotation_rate - model.nodedot) * 60
        return nodal_period, nodal_day

    def is_sun_synchronous(self) -> bool:
        """
        Checks whether this element is sun-synchronous: whether its nodal day
        is within `SUN_SYNCHRONOUS_NODAL_DAY_TOLERANCE_S` of a mean solar day.

        Returns:
            bool: True, if this element is sun-synchronous
        """
        _, nodal_day = self.get_nodal_period_and_day()
        return (
            abs(nodal_day - constants.EARTH_SOLAR_DAY_S)
            < SUN_SYNCHRONOUS_NODAL_DAY_TOLERANCE_S
        )

    def get_repeat_element(
        self, repeat_cycle: timedelta
    ) -> GeneralPerturbationsElements:
        """
        Gets a copy of this element maintained on the repeat ground track
        with the approximate repeat cycle: its drag terms (B* and the
        derivatives of mean motion) are set to zero and its mean motion is
        adjusted so that a whole number of nodal periods spans the refined
        repeat cycle (see `refine_repeat_cycle`) exactly. For a
        sun-synchronous orbit, its inclination is also adjusted so that its
        nodal day is exactly a mean solar day, as for an orbit maintained at
        a constant local time of ascending node. A satellite propagated with
        the copy returns to its initial Earth-fixed position after each
        repeat cycle, so that repeated cycles have no discontinuity at their
        ends. The copy is cached.

        A single element's mean motion reflects the satellite's position
        within its maintenance band (between maneuvers), so it differs
        slightly from the exact repeat: by 67 m of semimajor axis for an
        ICESat-2 element set, which, over its 91-day repeat cycle,
        accumulates to a 113 s difference of the time of ascending node
        crossing. Near the element's epoch, the copy departs from the
        element's own propagation at the rate of this difference (for
        ICESat-2, about 1.2 s per day).

        Args:
            repeat_cycle (timedelta): The approximate repeat cycle.

        Returns:
            GeneralPerturbationsElements: the maintained element
        """
        return get_cached(
            self,
            "repeat_element",
            (self._get_key(), repeat_cycle),
            lambda: _compute_repeat_element(self, repeat_cycle),
        )

    def get_repeat_cycle(
        self,
        max_delta_position: float | None = None,
        max_delta_velocity: float | None = None,
        max_search_duration: timedelta | None = None,
        max_delta_semimajor_axis: float | None = None,
    ) -> timedelta | None:
        """
        Compute this element's repeat cycle. Reuses a previously-computed
        repeat cycle if available.

        Uses the classical repeat-ground-track condition: the orbit
        repeats once a whole number of orbits (paced by the nodal period)
        fits a whole number of nodal days (paced by the node precession
        rate) -- the same rational-commensurability principle behind
        published repeat cycles like Landsat-8's 233 orbits/16 days. The
        nodal period and nodal day are read directly from the underlying
        SGP4 model's own secular rates (mdot, argpdot, nodedot), rather
        than re-derived independently, since SGP4 initialization applies a
        Kozai-to-Brouwer mean element correction that a from-scratch J2
        calculation (starting from the TLE's mean motion converted to a
        semimajor axis via plain Kepler's third law) would otherwise miss
        -- a small (~0.1%) but real discrepancy that is enough to make a
        genuine multi-week repeat cycle miss its tolerance entirely.

        A whole number of nodal days is confirmed as the repeat cycle in
        either of two ways, both checking that the propagated satellite
        returns within `max_delta_position` and `max_delta_velocity` of its
        initial position and velocity (which rejects candidates that the
        secular rates alone do not rule out, for example for an orbit whose
        argument of perigee precesses):

        1. The element itself, propagated without drag (as for an orbit
           maintained against drag), returns to its initial state.
        2. The element's semimajor axis is within `max_delta_semimajor_axis`
           of the exact repeat, and the element maintained on that repeat
           ground track (see `get_repeat_element`) returns to its initial
           state. A maintained orbit's semimajor axis varies within a band
           (for example, about 200 m for Landsat 9) as drag lowers it and
           maneuvers raise it, so an element set taken anywhere within the
           band can be far enough from the exact repeat for its ground
           track to drift beyond `max_delta_position` over a multi-week
           repeat cycle. Exact repeats of up to D nodal days are spaced by
           about 1/D^2 orbits per day, so this applies only to repeat cycles
           whose exact repeats are spaced by at least three times the
           tolerance (in low Earth orbit with the default 100 m, up to about
           32 days), beyond which it would admit chance near-repeats.

        The first (shortest) confirmed candidate is the reported repeat
        cycle: the whole number of nodal days (or, for a sun-synchronous
        orbit, mean solar days) of the element maintained on its repeat
        ground track (see `refine_repeat_cycle` and `get_repeat_element`).

        Long repeat cycles cannot be identified reliably from a single
        element set: among the many whole numbers of nodal days within a
        long search, some happen to be closer to an exact repeat than the
        true one. For an ICESat-2 element set, whose semimajor axis is 67 m
        from its 91-day repeat, chance repeats after 62 and 95 days are
        within 14 and 37 m. Searches longer than about 40 days risk such
        false repeat cycles; declare a long repeat cycle instead (see
        `GeneralPerturbationsOrbit.repeat_cycle`).

        Args:
            max_delta_position (float | None): the maximum difference in position (m) allowed for a repeat.
            max_delta_velocity (float | None): the maximum difference in velocity (m/s) allowed for a repeat.
            max_search_duration (timedelta | None): the maximum period of time to search for repeats.
            max_delta_semimajor_axis (float | None): the maximum difference (m) between the
                semimajor axis and that of an exact repeat for a candidate repeat.

        Returns:
            timedelta: the repeat cycle duration (if it exists)
        """
        search = RepeatCycleSearch.resolve(
            max_delta_position,
            max_delta_velocity,
            max_search_duration,
            max_delta_semimajor_axis,
        )
        return get_cached(
            self,
            "repeat_cycle",
            (self._get_key(), search),
            lambda: _search_repeat_cycle(self, search),
        )

    def refine_repeat_cycle(self, repeat_cycle: timedelta) -> timedelta:
        """
        Refines the approximate repeat cycle of an orbit maintained on a
        repeat ground track (for example, a nominal number of days) to the
        nearest whole number of nodal days (the period of the Earth's
        rotation relative to the orbit's precessing ascending node, read from
        the SGP4 model's secular rates, as in `get_repeat_cycle`), after which
        a repeat ground track orbit returns over the same ground track. For a
        sun-synchronous orbit (see `is_sun_synchronous`), the repeat cycle is
        refined to the nearest whole number of mean solar days instead: a
        maintained sun-synchronous orbit keeps its local time of ascending
        node (so that its nodal day is, on average, exactly a mean solar day),
        whereas the elements' nodal precession at their epoch typically
        differs slightly (for Landsat 9, by about 0.3 seconds per day), which
        would accumulate over many cycles.

        Unlike `get_repeat_cycle`, the satellite's return to its initial
        Earth-fixed position is not verified: over a long repeat cycle, a
        small difference between the elements' mean motion and the exact
        repeat (well within the orbit's maintenance) accumulates to tens of
        kilometers, and chance near-repeats at other durations can be closer.

        Args:
            repeat_cycle (timedelta): The approximate repeat cycle.

        Returns:
            timedelta: the refined repeat cycle
        """
        _, nodal_day = self.get_nodal_period_and_day()
        if (
            abs(nodal_day - constants.EARTH_SOLAR_DAY_S)
            < SUN_SYNCHRONOUS_NODAL_DAY_TOLERANCE_S
        ):
            nodal_day = constants.EARTH_SOLAR_DAY_S
        days = max(1, round(repeat_cycle.total_seconds() / nodal_day))
        return timedelta(seconds=days * nodal_day)

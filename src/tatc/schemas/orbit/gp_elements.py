"""
Object schema for general perturbations orbital elements.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import csv
import json
from datetime import datetime, timedelta, timezone

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, Field
from sgp4 import exporter, omm
from sgp4.api import WGS72, Satrec
from sgp4.conveniences import sat_epoch_datetime
from skyfield.api import EarthSatellite, Time
from skyfield.framelib import itrs
from skyfield.searchlib import find_minima

from ... import config, constants, utils

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

    def to_skyfield(self) -> EarthSatellite:
        """
        Converts this GP elements object to a Skyfield `EarthSatellite`,
        which can be used to propagate this orbital state via SGP4.

        Returns:
            skyfield.api.EarthSatellite: the Skyfield EarthSatellite
        """
        return EarthSatellite.from_omm(constants.timescale, self.to_omm_dict())

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
        earth_rotation_rate = 2 * np.pi / constants.EARTH_SIDEREAL_DAY_S * 60
        nodal_day = 2 * np.pi / (earth_rotation_rate - model.nodedot) * 60
        return nodal_period, nodal_day

    def get_repeat_element(
        self, repeat_cycle: timedelta
    ) -> GeneralPerturbationsElements:
        """
        Gets a copy of this element maintained on the repeat ground track
        with the approximate repeat cycle: its drag terms (B* and the
        derivatives of mean motion) are set to zero and its mean motion is
        adjusted so that a whole number of nodal periods spans the refined
        repeat cycle (see `refine_repeat_cycle`) exactly. A satellite
        propagated with the copy returns to its initial position and time of
        ascending node crossing after each repeat cycle, so that repeating
        its first cycle has no discontinuity at the end of each cycle.

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
        element = self.model_copy(
            update={"bstar": 0, "mean_motion_dot": 0, "mean_motion_ddot": 0}
        )
        nodal_period, _ = element.get_nodal_period_and_day()
        orbits = max(
            1,
            round(
                element.refine_repeat_cycle(repeat_cycle).total_seconds() / nodal_period
            ),
        )
        for _ in range(10):
            nodal_period, _ = element.get_nodal_period_and_day()
            refined = element.refine_repeat_cycle(repeat_cycle).total_seconds()
            residual = orbits * nodal_period - refined
            if abs(residual) < 1e-6:
                break
            # the nodal period varies (nearly) inversely with mean motion
            element = element.model_copy(
                update={"mean_motion": element.mean_motion * (1 + residual / refined)}
            )
        return element

    def get_repeat_cycle(
        self,
        max_delta_position: float | None = None,
        max_delta_velocity: float | None = None,
        max_search_duration: timedelta | None = None,
        lazy_load: bool | None = None,
        max_delta_semimajor_axis: float | None = None,
    ) -> timedelta | None:
        """
        Compute this element's repeat cycle. Lazy-loads a previously-computed
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
           ground track (see `get_repeat_element`: without drag, and with
           the exact repeat's mean motion) returns to its initial state. A
           maintained orbit's semimajor axis varies within a band (for
           example, about 200 m for Landsat 9) as drag lowers it and
           maneuvers raise it, so an element set taken anywhere within the
           band can be far enough from the exact repeat for its ground
           track to drift beyond `max_delta_position` over a multi-week
           repeat cycle. Exact repeats of up to D nodal days are spaced by
           about 1/D^2 orbits per day, so this applies only to repeat cycles
           whose exact repeats are spaced by at least three times the
           tolerance (in low Earth orbit with the default 100 m, up to about
           32 days), beyond which it would admit chance near-repeats.

        The first (shortest) confirmed candidate is the reported repeat
        cycle, rounded to a whole number of nodal days (or, for a
        sun-synchronous orbit, mean solar days; see `refine_repeat_cycle`).

        Long repeat cycles cannot be identified reliably from a single
        element set: among the many whole numbers of nodal days within a
        long search, some happen to be closer to an exact repeat than the
        true one. For an ICESat-2 element set, whose semimajor axis is 67 m
        from its 91-day repeat, chance repeats after 62 and 95 days are
        within 14 and 37 m. Searches longer than about 40 days risk such
        false repeat cycles; declare a long repeat cycle instead (see
        `GeneralPerturbationsOrbit.repeat_cycle`).

        This is scoped to a single element on purpose: a
        `GeneralPerturbationsOrbit` with multiple elements may span a
        significant maneuver (altitude change, plane change, etc.)
        partway through its history, after which this element's repeat
        cycle (if any) may no longer apply. See
        `GeneralPerturbationsOrbit.get_repeat_cycle`, which checks every
        element's own repeat cycle for mutual consistency before
        reporting one for the whole orbit.

        Args:
            max_delta_position (float | None): the maximum difference in position (m) allowed for a repeat.
            max_delta_velocity (float | None): the maximum difference in velocity (m/s) allowed for a repeat.
            max_search_duration (timedelta | None): the maximum period of time to search for repeats.
            lazy_load (bool | None): True, if the previously-computed repeat cycle should be loaded.
            max_delta_semimajor_axis (float | None): the maximum difference (m) between the
                semimajor axis and that of an exact repeat for a candidate repeat.

        Returns:
            timedelta: the repeat cycle duration (if it exists)
        """
        # load defaults
        if max_delta_position is None:
            max_delta_position = config.get_rc().repeat_cycle_delta_position_m
        if max_delta_velocity is None:
            max_delta_velocity = config.get_rc().repeat_cycle_delta_velocity_m_per_s
        if max_search_duration is None:
            max_search_duration = timedelta(
                days=config.get_rc().repeat_cycle_search_duration_days
            )
        if lazy_load is None:
            lazy_load = config.get_rc().repeat_cycle_lazy_load
        if max_delta_semimajor_axis is None:
            max_delta_semimajor_axis = (
                config.get_rc().repeat_cycle_delta_semimajor_axis_m
            )

        # keyed by the element's field values and the search options, so
        # that a copy with changed fields (e.g. from `model_copy(update=...)`)
        # or a search with other options is recomputed
        key = (
            tuple(self.model_dump().values()),
            max_delta_position,
            max_delta_velocity,
            max_search_duration,
            max_delta_semimajor_axis,
        )
        cached = self.__dict__.get("repeat_cycle") if lazy_load else None
        repeat_cycle = cached[1] if cached is not None and cached[0] == key else None
        if repeat_cycle is None:
            epoch = self.epoch
            # analytic repeat ground track candidates: how many nodal days
            # (D) are needed for a whole number of orbits (C) to elapse
            element = self.model_copy(
                update={"bstar": 0, "mean_motion_dot": 0, "mean_motion_ddot": 0}
            )
            nodal_period, nodal_day = element.get_nodal_period_and_day()
            orbits_per_day = nodal_day / nodal_period
            max_days = int(max_search_duration.total_seconds() / nodal_day)
            days_range = np.arange(1, max_days + 1)
            orbit_counts = np.round(orbits_per_day * days_range)
            residual_orbits = orbits_per_day * days_range - orbit_counts
            # ground-track drift (m) at the equator implied by missing a whole
            # orbit count by residual_orbits
            ground_track_spacing = (
                2 * np.pi * constants.EARTH_MEAN_RADIUS / orbits_per_day
            )
            drift = np.abs(residual_orbits) * ground_track_spacing
            # difference (m) between the semimajor axis and that of the exact
            # repeat, from the sensitivity of orbits per nodal day to the
            # semimajor axis (by a finite difference of mean motion)
            perturbed = element.model_copy(
                update={"mean_motion": element.mean_motion * (1 + 1e-6)}
            )
            perturbed_period, perturbed_day = perturbed.get_nodal_period_and_day()
            sensitivity = (perturbed_day / perturbed_period - orbits_per_day) / (
                perturbed.get_semimajor_axis() - element.get_semimajor_axis()
            )
            delta_semimajor_axis = (
                orbit_counts / days_range - orbits_per_day
            ) / sensitivity
            # exact repeats of up to D nodal days are spaced by about 1/D^2
            # orbits per day: apply the semimajor axis tolerance only where
            # they are spaced by at least three times the tolerance, beyond
            # which it would admit chance near-repeats
            short = 3 * max_delta_semimajor_axis * abs(sensitivity) * days_range**2 <= 1
            # generous margin: the drift is estimated from secular rates,
            # while the verification below also includes periodic terms
            near_drift = drift < 3 * max_delta_position
            near_semimajor_axis = short & (
                np.abs(delta_semimajor_axis) < max_delta_semimajor_axis
            )
            # verify candidates without drag (as for an orbit maintained
            # against drag), which would otherwise move the propagated
            # satellite away from its repeat ground track
            satellite = element.to_skyfield()

            def find_closest_approach(
                satellite: EarthSatellite, center: datetime, period: float
            ) -> tuple[float, float] | None:
                # initial position and velocity in the Earth-fixed frame
                position_0, velocity_0 = satellite.at(
                    constants.timescale.from_datetime(epoch)
                ).frame_xyz_and_velocity(itrs)
                p_0_m = np.array(position_0.m)
                v_0_m_per_s = np.array(velocity_0.m_per_s)

                def position_error(t: Time) -> npt.NDArray[np.float64]:
                    position, _ = satellite.at(t).frame_xyz_and_velocity(itrs)
                    return np.linalg.norm((np.array(position.m).T - p_0_m.T).T, axis=0)

                position_error.rough_period = period / 86400  # type: ignore
                window = timedelta(seconds=period / 2)
                times, errors = find_minima(
                    constants.timescale.from_datetime(center - window),
                    constants.timescale.from_datetime(center + window),
                    position_error,
                )
                if len(times) == 0:
                    return None
                t_min = times[np.argmin(errors)]
                position, velocity = satellite.at(t_min).frame_xyz_and_velocity(itrs)
                return (
                    float(np.linalg.norm(np.array(position.m) - p_0_m)),
                    float(np.linalg.norm(np.array(velocity.m_per_s) - v_0_m_per_s)),
                )

            # assign zero repeat cycle value to avoid recalculation, unless
            # a candidate below is confirmed
            def is_confirmed(result: tuple[float, float] | None) -> bool:
                return (
                    result is not None
                    and result[0] < max_delta_position
                    and result[1] < max_delta_velocity
                )

            repeat_cycle = timedelta(0)
            for i in np.flatnonzero(near_drift | near_semimajor_axis):
                candidate = timedelta(seconds=int(days_range[i]) * nodal_day)
                # verify a candidate near the element's mean motion with the
                # element itself (without drag), and one near the semimajor
                # axis of the exact repeat with the element maintained on its
                # repeat ground track (with the exact repeat's mean motion)
                confirmed = near_drift[i] and is_confirmed(
                    find_closest_approach(satellite, epoch + candidate, nodal_period)
                )
                if not confirmed and near_semimajor_axis[i]:
                    maintained = self.get_repeat_element(candidate)
                    confirmed = is_confirmed(
                        find_closest_approach(
                            maintained.to_skyfield(),
                            epoch + maintained.refine_repeat_cycle(candidate),
                            nodal_period,
                        )
                    )
                if confirmed:
                    # report a whole number of nodal (or, for a
                    # sun-synchronous orbit, mean solar) days, as for a
                    # maintained orbit
                    repeat_cycle = self.refine_repeat_cycle(candidate)
                    break
            self.__dict__["repeat_cycle"] = (key, repeat_cycle)  # type: ignore
        if repeat_cycle is not None and repeat_cycle > timedelta(0):
            return repeat_cycle
        return None

    def refine_repeat_cycle(self, repeat_cycle: timedelta) -> timedelta:
        """
        Refines the approximate repeat cycle of an orbit maintained on a
        repeat ground track (for example, a nominal number of days) to the
        nearest whole number of nodal days (the period of the Earth's
        rotation relative to the orbit's precessing ascending node, read from
        the SGP4 model's secular rates, as in `get_repeat_cycle`), after which
        a repeat ground track orbit returns over the same ground track. For a
        sun-synchronous orbit (whose nodal day is within
        `SUN_SYNCHRONOUS_NODAL_DAY_TOLERANCE_S` of a mean solar day), the
        repeat cycle is refined to the nearest whole number of mean solar days
        instead: a maintained sun-synchronous orbit keeps its local time of
        ascending node (so that its nodal day is, on average, exactly a mean
        solar day), whereas the elements' nodal precession at their epoch
        typically differs slightly (for Landsat 9, by about 0.3 seconds per
        day), which would accumulate over many cycles.

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
        if abs(nodal_day - 86400) < SUN_SYNCHRONOUS_NODAL_DAY_TOLERANCE_S:
            nodal_day = 86400
        days = max(1, round(repeat_cycle.total_seconds() / nodal_day))
        return timedelta(seconds=days * nodal_day)

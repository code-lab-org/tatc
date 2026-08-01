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

    def get_repeat_cycle(
        self,
        max_delta_position: float | None = None,
        max_delta_velocity: float | None = None,
        max_search_duration: timedelta | None = None,
        lazy_load: bool | None = None,
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
        genuine multi-week repeat cycle miss its tolerance entirely. Each
        analytically-predicted candidate is confirmed by directly
        propagating the real orbit and checking that both position and
        velocity match the initial state within tolerance; the first
        (shortest) candidate that does so is the reported repeat cycle.

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

        if lazy_load:
            repeat_cycle = self.__dict__.get("repeat_cycle")
        else:
            repeat_cycle = None
        if repeat_cycle is None:
            epoch = self.epoch
            satellite = self.to_skyfield()
            # record the initial position and velocity in Earth-centered Earth-fixed frame
            position_0, velocity_0 = satellite.at(
                constants.timescale.from_datetime(epoch)
            ).frame_xyz_and_velocity(itrs)
            p_0_m = np.array(position_0.m)
            v_0_m_per_s = np.array(velocity_0.m_per_s)

            # analytic repeat ground track candidates: how many nodal days
            # (D) are needed for a whole number of orbits (C) to elapse.
            # mdot/argpdot/nodedot (rad/minute) are SGP4's own secular
            # rates for mean anomaly, argument of perigee, and RAAN.
            model = satellite.model
            nodal_period = 2 * np.pi / (model.mdot + model.argpdot) * 60
            earth_rotation_rate = 2 * np.pi / constants.EARTH_SIDEREAL_DAY_S * 60
            nodal_day = 2 * np.pi / (earth_rotation_rate - model.nodedot) * 60
            orbits_per_day = nodal_day / nodal_period
            max_days = int(max_search_duration.total_seconds() / nodal_day)
            days_range = np.arange(1, max_days + 1)
            orbit_counts = np.round(orbits_per_day * days_range)
            residual_orbits = np.abs(orbits_per_day * days_range - orbit_counts)
            # approximate ground-track drift (m) implied by missing a whole
            # orbit count by residual_orbits, at the equator; only a coarse
            # heuristic to shortlist candidates worth verifying by direct
            # propagation below, not itself a pass/fail criterion
            ground_track_spacing = (
                2 * np.pi * constants.EARTH_MEAN_RADIUS / orbits_per_day
            )
            approx_drift = residual_orbits * ground_track_spacing
            # generous margin: this estimate is J2-only, while the actual
            # verification below propagates the real (e.g. SGP4) orbit
            candidate_days = days_range[approx_drift < 3 * max_delta_position]

            def position_error(t: Time) -> npt.NDArray[np.float64]:
                position, _ = satellite.at(t).frame_xyz_and_velocity(itrs)
                return np.linalg.norm((np.array(position.m).T - p_0_m.T).T, axis=0)

            position_error.rough_period = nodal_period / 86400  # type: ignore

            def find_closest_approach(
                center: datetime,
            ) -> tuple[datetime, float, float] | None:
                window = timedelta(seconds=nodal_period / 2)
                times, errors = find_minima(
                    constants.timescale.from_datetime(center - window),
                    constants.timescale.from_datetime(center + window),
                    position_error,
                )
                if len(times) == 0:
                    return None
                t_min = times[np.argmin(errors)]
                position, velocity = satellite.at(t_min).frame_xyz_and_velocity(itrs)
                delta_position = float(np.linalg.norm(np.array(position.m) - p_0_m))
                delta_velocity = float(
                    np.linalg.norm(np.array(velocity.m_per_s) - v_0_m_per_s)
                )
                return t_min.utc_datetime(), delta_position, delta_velocity

            # assign zero repeat cycle value to avoid recalculation, unless
            # a candidate below is confirmed
            repeat_cycle = timedelta(0)
            for days in candidate_days:
                predicted = epoch + timedelta(seconds=int(days) * nodal_day)
                result = find_closest_approach(predicted)
                if result is None:
                    continue
                t_min, delta_position, delta_velocity = result
                if (
                    delta_position < max_delta_position
                    and delta_velocity < max_delta_velocity
                ):
                    repeat_cycle = t_min - epoch
                    break
            self.__dict__["repeat_cycle"] = repeat_cycle  # type: ignore
        if repeat_cycle is not None and repeat_cycle > timedelta(0):
            return repeat_cycle
        return None

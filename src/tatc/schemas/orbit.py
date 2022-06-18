# -*- coding: utf-8 -*-
"""
Object schemas for satellite orbits.

@author: Paul T. Grogan <pgrogan@stevens.edu>
"""

from datetime import datetime, time, timedelta, timezone
import numpy as np
from pydantic import BaseModel, Field, validator
from sgp4.api import Satrec, WGS72
from sgp4 import exporter
from sgp4.conveniences import sat_epoch_datetime
from typing import List, Optional
from typing_extensions import Literal
import re

from .. import constants, utils


class TwoLineElements(BaseModel):
    """
    Orbit defined with standard two line elements.

    :param type: The type of orbit
    :type type: literal["tle"], restricted: literal["tle"]
    :param tle: The two line elements.
    :type tle: list[str]
    """

    type: Literal["tle"] = Field("tle")
    tle: List[str] = Field(
        ...,
        description="Two line elements.",
        min_items=2,
        max_items=2,
        example=[
            "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
            "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
        ],
    )

    def get_catalog_number(self) -> int:
        """
        Gets the TLE catalog number.
        """
        return int(self.tle[0][2:7])

    def get_classification(self) -> str:
        """
        Gets the TLE classification type (U: unclassified; C: classified).
        """
        return self.tle[0][7]

    def get_international_designator(self) -> str:
        """
        Gets the TLE international designator.
        """
        year = self.tle[0][9:11]
        launch = self.tle[0][11:14]
        piece = self.tle[0][14:17]
        return str(
            ("" if not year.isdigit() else "20" if int(year) < 57 else "19")
            + year
            + "-"
            + launch
            + piece
        ).strip()

    def get_epoch(self) -> datetime:
        """
        Gets the TLE epoch time.
        """
        year = int(self.tle[0][18:20])
        days = float(self.tle[0][20:32])
        return datetime(
            year + (2000 if year < 57 else 1900), 1, 1, tzinfo=timezone.utc
        ) + timedelta(days=days)

    def get_first_derivative_mean_motion(self) -> float:
        """
        Gets the first derivative of mean motion (ballistic coefficient).
        """
        return float(self.tle[0][33:43])

    def get_second_derivative_mean_motion(self) -> float:
        """
        Gets the second derivative of mean motion.
        """
        return float("0." + self.tle[0][44:50].strip()) * 10 ** (
            int(self.tle[0][50:52])
        )

    def get_b_star(self) -> float:
        """
        Gets the b-star term (drag or radiation pressure coefficient).
        """
        return float("0." + self.tle[0][53:59].strip()) * 10 ** (
            int(self.tle[0][59:61])
        )

    def get_ephemeris_type(self) -> int:
        """
        Gets the TLE ephemeris type.
        """
        return int(self.tle[0][62])

    def get_element_set_number(self) -> int:
        """
        Gets the TLE element set number.
        """
        return int(self.tle[0][64:68])

    def get_inclination(self) -> float:
        """
        Gets the orbit inclination (decimal degrees).
        """
        return float(self.tle[1][8:16])

    def get_right_ascension_ascending_node(self) -> float:
        """
        Gets the right ascension of ascending node (decimal degrees).
        """
        return float(self.tle[1][17:25])

    def get_eccentricity(self) -> float:
        """
        Gets the eccentricity.
        """
        return float("0." + self.tle[1][26:33].strip())

    def get_perigee_argument(self) -> float:
        """
        Gets the argument of perigee (decimal degrees).
        """
        return float(self.tle[1][34:42])

    def get_mean_anomaly(self) -> float:
        """
        Gets the mean anomaly (decimal degrees).
        """
        return float(self.tle[1][43:51])

    def get_mean_motion(self) -> float:
        """
        Gets the mean motion (revolutions per day).
        """
        return float(self.tle[1][52:63])

    def get_revolution_number_at_epoch(self) -> int:
        """
        Gets the revolution number at epoch.
        """
        return int(self.tle[1][63:68])

    def get_semimajor_axis(self) -> float:
        """
        Gets the semimajor axis (meters).
        """
        mean_motion_rad_s = self.get_mean_motion() * 2 * np.pi / 86400
        return np.power(
            constants.earth_mu / mean_motion_rad_s**2,
            1 / 3,
        )

    def get_altitude(self) -> float:
        """
        Gets the altitude (meters).
        """
        return self.get_semimajor_axis() - constants.earth_mean_radius

    def get_true_anomaly(self) -> float:
        """
        Gets the true anomaly (decimal degrees).
        """
        return utils.mean_anomaly_to_true_anomaly(
            self.get_mean_anomaly(), self.get_eccentricity()
        )

    @validator("tle")
    def valid_tle(cls, v):
        """
        Validate the two line element set
        """
        # based on orekit's TLE.isFormatOK function
        if len(v[0]) != 69:
            raise ValueError("Invalid tle: line 1 incorrect length.")
        if len(v[1]) != 69:
            raise ValueError("Invalid tle: line 2 incorrect length.")

        line_1_pattern = r"1 [ 0-9A-HJ-NP-Z][ 0-9]{4}[A-Z] [ 0-9]{5}[ A-Z]{3} [ 0-9]{5}[.][ 0-9]{8} (?:(?:[ 0+-][.][ 0-9]{8})|(?: [ +-][.][ 0-9]{7})) [ +-][ 0-9]{5}[+-][ 0-9] [ +-][ 0-9]{5}[+-][ 0-9] [ 0-9] [ 0-9]{4}[ 0-9]"
        if re.match(line_1_pattern, v[0]) is None:
            raise ValueError("Invalid tle: line 1 does not match pattern.")
        line_2_pattern = r"2 [ 0-9A-HJ-NP-Z][ 0-9]{4} [ 0-9]{3}[.][ 0-9]{4} [ 0-9]{3}[.][ 0-9]{4} [ 0-9]{7} [ 0-9]{3}[.][ 0-9]{4} [ 0-9]{3}[.][ 0-9]{4} [ 0-9]{2}[.][ 0-9]{13}[ 0-9]"
        if re.match(line_2_pattern, v[1]) is None:
            raise ValueError("Invalid tle: line 2 does not match pattern.")

        def checksum(line):
            sum = 0
            for i in range(68):
                if line[i].isdigit():
                    sum += int(line[i])
                elif line[i] == "-":
                    sum += 1
            return sum % 10

        if int(v[0][68]) != checksum(v[0]):
            raise ValueError("Invalid tle: line 1 checksum failed.")
        if int(v[1][68]) != checksum(v[1]):
            raise ValueError("Invalid tle: line 2 checksum failed.")
        return v

    def get_derived_orbit(self, delta_mean_anomaly, delta_raan):
        """
        Gets a derived orbit.

        Args:
            :param delta_mean_anomaly:  Delta mean anomaly (degrees).
            :type delta_mean_anomaly: float
            :param delta_raan:  Delta right ascension of ascending node (degrees).
            :type delta_raan: float

        Returns:
            TwoLineElements
        """
        lead_tle = Satrec.twoline2rv(self.tle[0], self.tle[1])
        epoch = sat_epoch_datetime(lead_tle)
        satrec = Satrec()
        satrec.sgp4init(
            WGS72,
            "i",
            0,
            (epoch - datetime(1950, 1, 1, tzinfo=timezone.utc)) / timedelta(days=1),
            lead_tle.bstar,
            lead_tle.ndot,
            lead_tle.nddot,
            lead_tle.ecco,
            lead_tle.argpo,
            lead_tle.inclo,
            np.mod(lead_tle.mo + np.radians(delta_mean_anomaly), 2 * np.pi),
            lead_tle.no_kozai,
            np.mod(lead_tle.nodeo + np.radians(delta_raan), 2 * np.pi),
        )
        tle1, tle2 = exporter.export_tle(satrec)
        return TwoLineElements(tle=[tle1.replace("\x00", "U"), tle2])

    def to_tle(self):
        """
        Return the two line element set for this orbit.
        """

        return self


class OrbitBase(BaseModel):
    """
    Base class for orbit definition"

    :param altitude:  Mean altitude (meters).
    :type altitude: float
    :param true_anomaly: True anomaly (degrees).
    :type true_anomaly: float, default: 0
    :param epoch: Timestamp (epoch) of the  initial orbital state.
    :type epoch: :class:`datetime.datetime`, optional, default: None
    """

    altitude: float = Field(..., description="Mean altitude (meters).")
    true_anomaly: float = Field(0, description="True anomaly (degrees).", ge=0, lt=360)
    epoch: Optional[datetime] = Field(
        datetime.now(tz=timezone.utc),
        description="Timestamp (epoch) of the initial orbital state.",
    )

    def get_mean_anomaly(self) -> float:
        """
        Gets the mean anomaly (decimal degrees).
        """
        return utils.true_anomaly_to_mean_anomaly(self.true_anomaly)

    def get_semimajor_axis(self) -> float:
        """
        Gets the semimajor axis (meters).
        """
        return constants.earth_mean_radius + self.altitude

    def get_mean_motion(self) -> float:
        """
        Gets the mean motion (revolutions per day).
        """
        return 1 / (self.get_orbit_period()/timedelta(days=1))

    def get_orbit_period(self) -> timedelta:
        """
        Gets the approximate orbit period.
        """
        return timedelta(seconds=utils.compute_orbit_period(self.altitude))


class CircularOrbit(OrbitBase):
    """
    Orbit specification using Keplerian elements for elliptical motion -- circular motion case.

    :param type: The type of orbit.
    :type type: literal["circular"], restricted: literal["circular"]
    :param inclination: Orbit inclination (degrees).
    :type inclination: float
    :param right_ascension_ascending_node: Right ascension of ascending node (degrees).
    :type right_ascension_ascending_node: float (0 to 360), default: 0
    """

    type: Literal["circular"] = Field("circular")
    inclination: float = Field(0, description="Inclination (degrees).", ge=0, lt=180)
    right_ascension_ascending_node: float = Field(
        0, description="Right ascension of ascending node (degrees).", ge=0, lt=360
    )

    def get_derived_orbit(self, delta_mean_anomaly, delta_raan):
        """
        Gets a derived orbit.

        Args:
            :param delta_mean_anomaly:  Delta mean anomaly (degrees).
            :type delta_mean_anomaly: float
            :param delta_raan:  Delta right ascension of ascending node (degrees).
            :type delta_raan: float

        Returns:
            CircularOrbit
        """
        true_anomaly = utils.mean_anomaly_to_true_anomaly(
            self.get_mean_anomaly() + delta_mean_anomaly
        )
        raan = self.right_ascension_ascending_node + delta_raan
        return CircularOrbit(
            altitude=self.altitude,
            true_anomaly=true_anomaly,
            epoch=self.epoch,
            inclination=self.inclination,
            right_ascension_ascending_node=raan,
        )

    def to_tle(self) -> TwoLineElements:
        """
        Create a two line element set representation of the orbit.
        """
        return KeplerianOrbit(
            altitude=self.altitude,
            true_anomaly=self.true_anomaly,
            epoch=self.epoch,
            inclination=self.inclination,
            right_ascension_ascending_node=self.right_ascension_ascending_node,
        ).to_tle()


class SunSynchronousOrbit(OrbitBase):
    """
    Orbit defined by sun synchronous parameters.
    """

    type: Literal["sso"] = Field("sso")
    altitude: float = Field(
        ...,
        description="Mean altitude (meters).",
        ge=0,
        lt=12352000 - constants.earth_mean_radius,
    )
    equator_crossing_time: time = Field(
        ..., description="Equator crossing time (local solar time)."
    )
    equator_crossing_ascending: bool = Field(
        True,
        description="True, if the equator crossing time is ascending (south-to-north).",
    )

    def get_inclination(self) -> float:
        """
        Gets the inclination (decimal degrees).
        """
        return np.degrees(
            np.arccos(-np.power(self.get_semimajor_axis() / 12352000, 7 / 2))
        )

    def get_right_ascension_ascending_node(self) -> float:
        """
        Gets the right ascension of ascending node (decimal degrees).
        """
        ect_day = timedelta(
            hours=self.equator_crossing_time.hour,
            minutes=self.equator_crossing_time.minute,
            seconds=self.equator_crossing_time.second,
            microseconds=self.equator_crossing_time.microsecond,
        ) / timedelta(days=1)
        t = constants.timescale.from_datetime(self.epoch)
        sun = constants.de421["sun"]
        earth = constants.de421["earth"]
        ra, _, _ = earth.at(t).observe(sun).radec()
        return (
            ra._degrees + 360 * ect_day + 180 * self.equator_crossing_ascending
        ) % 360

    def get_derived_orbit(self, delta_mean_anomaly, delta_raan):
        """
        Gets a derived orbit.

        Args:
            :param delta_mean_anomaly:  Delta mean anomaly (degrees).
            :type delta_mean_anomaly: float
            :param delta_raan:  Delta right ascension of ascending node (degrees).
            :type delta_raan: float

        Returns:
            CircularOrbit
        """
        true_anomaly = utils.mean_anomaly_to_true_anomaly(
            self.get_mean_anomaly() + delta_mean_anomaly
        )
        raan = self.get_right_ascension_ascending_node() + delta_raan
        # TODO the resulting orbit *is* still sun-synchronous; however, with a
        # different equator-crossing time. Need to do the math to determine.
        # For now, simply return a circular orbit with the correct parameters.
        return CircularOrbit(
            altitude=self.altitude,
            true_anomaly=true_anomaly,
            epoch=self.epoch,
            inclination=self.get_inclination(),
            right_ascension_ascending_node=raan,
        )

    def to_tle(self) -> TwoLineElements:
        """
        Create a two line element set representation of the orbit.
        """
        semimajor_axis = constants.earth_mean_radius + self.altitude
        return KeplerianOrbit(
            altitude=self.altitude,
            inclination=self.get_inclination(),
            right_ascension_ascending_node=self.get_right_ascension_ascending_node(),
            true_anomaly=self.true_anomaly,
            epoch=self.epoch,
        ).to_tle()


class KeplerianOrbit(CircularOrbit):
    """
    Orbit specification using Keplerian elements for elliptical motion.
    """

    type: Literal["keplerian"] = Field("keplerian")
    eccentricity: float = Field(0, description="Eccentricity.", ge=0)
    perigee_argument: float = Field(
        0, description="Perigee argument (degrees).", ge=0, lt=360
    )

    def get_mean_anomaly(self) -> float:
        """
        Gets the mean anomaly (decimal degrees).
        """
        return utils.true_anomaly_to_mean_anomaly(self.true_anomaly, self.eccentricity)

    def get_derived_orbit(self, delta_mean_anomaly, delta_raan):
        """
        Gets a derived orbit.

        Args:
            :param delta_mean_anomaly:  Delta mean anomaly (degrees).
            :type delta_mean_anomaly: float
            :param delta_raan:  Delta right ascension of ascending node (degrees).
            :type delta_raan: float

        Returns:
            KeplerianOrbit
        """
        true_anomaly = utils.mean_anomaly_to_true_anomaly(
            self.get_mean_anomaly() + delta_mean_anomaly, eccentricity=self.eccentricity
        )
        raan = self.right_ascension_ascending_node + delta_raan
        return KeplerianOrbit(
            altitude=self.altitude,
            true_anomaly=true_anomaly,
            epoch=self.epoch,
            inclination=self.inclination,
            right_ascension_ascending_node=raan,
            eccentricity=self.eccentricity,
            perigee_argument=self.perigee_argument,
        )

    def to_tle(self) -> TwoLineElements:
        satrec = Satrec()
        satrec.sgp4init(
            WGS72,
            "i",
            0,
            (self.epoch - datetime(1950, 1, 1, tzinfo=timezone.utc))
            / timedelta(days=1),
            0,
            0.0,
            0.0,
            self.eccentricity,
            np.radians(self.perigee_argument),
            np.radians(self.inclination),
            np.radians(self.get_mean_anomaly()),
            self.get_mean_motion() * 2 * np.pi / 1440,
            np.radians(self.right_ascension_ascending_node),
        )
        tle1, tle2 = exporter.export_tle(satrec)
        return TwoLineElements(tle=[tle1.replace("\x00", "U"), tle2])
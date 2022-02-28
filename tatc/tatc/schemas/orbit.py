# -*- coding: utf-8 -*-
"""
Object schemas for satellite orbits.

@author: Paul T. Grogan <pgrogan@stevens.edu>
"""

from datetime import datetime, time, timedelta, timezone
import numpy as np
from pydantic import BaseModel, Field, validator
from skyfield.api import load
from sgp4.api import Satrec, WGS72
from sgp4 import exporter
from typing import List, Optional
from typing_extensions import Literal
import re

from .. import constants


class Orbit(BaseModel):
    type: str = Field(..., description="Type of orbit.")


class TwoLineElements(Orbit):
    """
    Orbit defined with standard two line elements.
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

    @validator("tle")
    def valid_tle(cls, v):
        """

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

    def to_tle(self):
        return self


class OrbitBase(Orbit):
    """
    Base class for orbit definition
    """
    altitude: float = Field(..., description="Mean altitude (meters).")
    true_anomaly: float = Field(0, description="True anomaly (degrees).", ge=0, lt=360)
    epoch: Optional[datetime] = Field(
        None, description="Timestamp (epoch) of the initial orbital state."
    )


class CircularOrbit(OrbitBase):
    """
    Orbit specification using Keplerian elements for elliptical motion -- circular motion case.
    """

    type: Literal["circular"] = Field("circular")
    inclination: float = Field(0, description="Inclination (degrees).", ge=0, lt=180)
    right_ascension_ascending_node: float = Field(
        0, description="Right ascension of ascending node (degrees).", ge=0, lt=360
    )

    def to_tle(self) -> TwoLineElements:
        return KeplerianOrbit(
            altitude=self.altitude,
            inclination=self.inclination,
            right_ascension_ascending_node=self.right_ascension_ascending_node,
            true_anomaly=self.true_anomaly,
            epoch=self.epoch,
        ).to_tle()


class SunSynchronousOrbit(OrbitBase):
    """
    Orbit defined by sun synchronous parameters.
    """

    type: Literal["sso"] = Field("sso")
    equator_crossing_time: time = Field(
        ..., description="Equator crossing time (local solar time)."
    )
    equator_crossing_ascending: bool = Field(
        True,
        description="True, if the equator crossing time is ascending (south-to-north).",
    )

    def to_tle(self) -> TwoLineElements:
        semimajor_axis = constants.earth_mean_radius + self.altitude
        inclination = np.arccos(-np.power(semimajor_axis / 12352000, 7 / 2))
        ect_day = (
            timedelta(
                hours=self.equator_crossing_time.hour,
                minutes=self.equator_crossing_time.minute,
                seconds=self.equator_crossing_time.second,
                microseconds=self.equator_crossing_time.microsecond,
            )
            / timedelta(days=1)
        )
        if self.epoch is None:
            epoch = datetime.now(tz=timezone.utc)
        else:
            epoch = self.epoch
        ts = load.timescale()
        t = ts.from_datetime(epoch)
        eph = load("de421.bsp")
        sun = eph["sun"]
        earth = eph["earth"]
        ra, _, _ = earth.at(t).observe(sun).radec()
        right_ascension_ascending_node = (
            ra.radians + 2 * np.pi * ect_day - np.pi * self.equator_crossing_ascending
        ) % (2 * np.pi)
        if right_ascension_ascending_node < 0:
            right_ascension_ascending_node += 2 * np.pi
        return KeplerianOrbit(
            altitude=self.altitude,
            inclination=np.degrees(inclination),
            right_ascension_ascending_node=np.degrees(right_ascension_ascending_node),
            true_anomaly=self.true_anomaly,
            epoch=epoch,
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

    def to_tle(self) -> TwoLineElements:
        semimajor_axis = constants.earth_mean_radius + self.altitude
        mean_motion = np.sqrt(constants.earth_mu / semimajor_axis ** 3)
        eccentricity = self.eccentricity
        inclination = np.radians(self.inclination)
        perigee_argument = np.radians(self.perigee_argument)
        right_ascension_ascending_node = np.radians(self.right_ascension_ascending_node)
        true_anomaly = np.radians(self.true_anomaly)
        mean_anomaly = (
            true_anomaly
            - 2 * eccentricity * np.sin(true_anomaly)
            + (3 / 4 * eccentricity ** 2 + 1 / 8 * eccentricity ** 4)
            * np.sin(2 * true_anomaly)
            - 1 / 3 * eccentricity ** 3 * np.sin(3 * true_anomaly)
            + 5 / 32 * eccentricity ** 4 * np.sin(4 * true_anomaly)
        )
        if self.epoch is None:
            epoch = datetime.now(tz=timezone.utc)
        else:
            epoch = self.epoch
        satrec = Satrec()
        satrec.sgp4init(
            WGS72,
            "i",
            0,
            (epoch - datetime(1949, 12, 31, tzinfo=timezone.utc)) / timedelta(days=1),
            0,
            0.0,
            0.0,
            eccentricity,
            perigee_argument,
            inclination,
            mean_anomaly,
            mean_motion * 60,
            right_ascension_ascending_node,
        )
        tle1, tle2 = exporter.export_tle(satrec)
        return TwoLineElements(tle=[tle1.replace("\x00", "U"), tle2])

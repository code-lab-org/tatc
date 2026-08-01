"""
Object schemas for orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timedelta, timezone

from pydantic import BaseModel, Field

from ... import config, constants, utils
from .gp import GeneralPerturbationsOrbit


class OrbitBase(BaseModel):
    """
    Base class for orbits.
    """

    true_anomaly: float = Field(default=0, description="True anomaly (degrees).", ge=0, lt=360)
    epoch: datetime = Field(
        default=datetime(2020, 1, 1, tzinfo=timezone.utc),
        description="Timestamp (epoch) of the initial orbital state.",
    )

    def get_true_anomaly(self) -> float:
        """
        Gets the true anomaly (decimal degrees).

        Returns:
            float: the true anomaly
        """
        return self.true_anomaly

    def get_epoch(self) -> datetime:
        """
        Gets the timestamp (epoch) of the initial orbital state.

        Returns:
            datetime: the epoch
        """
        return self.epoch

    def get_mean_anomaly(self) -> float:
        """
        Gets the mean anomaly (decimal degrees).

        Returns:
            float: the mean anomaly
        """
        return utils.orbital.true_anomaly_to_mean_anomaly(self.true_anomaly)

    def get_semimajor_axis(self) -> float:
        """
        Gets the semimajor axis (meters). Must be implemented by
        subclasses, since OrbitBase has no universal way to derive it.

        Returns:
            float: the semimajor axis
        """
        raise NotImplementedError(
            "get_semimajor_axis() must be implemented in subclasses."
        )

    def get_mean_altitude(self) -> float:
        """
        Gets the mean altitude (meters) above the WGS 84 mean radius.
        Derived from the semimajor axis by default; subclasses that store
        altitude directly (e.g. circular orbits) may override this for
        efficiency.

        Returns:
            float: the mean altitude
        """
        return self.get_semimajor_axis() - constants.EARTH_MEAN_RADIUS

    def get_mean_motion(self) -> float:
        """
        Gets the mean motion (revolutions per day). Derived from the
        semimajor axis by default.

        Returns:
            float: the mean motion
        """
        return utils.orbital.semimajor_axis_to_mean_motion(self.get_semimajor_axis())

    def get_orbit_period(self) -> timedelta:
        """
        Gets the approximate orbit period. Derived from the semimajor axis
        by default. Subclasses for which the orbit period is instead the
        defining (independent) quantity -- with semimajor axis derived
        from it, rather than the other way around -- must override this
        method (and cannot rely on this default, which would be
        circular).

        Returns:
            timedelta: the orbit period
        """
        return timedelta(
            seconds=utils.orbital.semimajor_axis_to_orbit_period(
                self.get_semimajor_axis()
            )
        )

    def get_inclination(self) -> float:
        """
        Gets the inclination (degrees). Must be implemented by subclasses,
        since OrbitBase has no universal way to derive it.

        Returns:
            float: the inclination
        """
        raise NotImplementedError(
            "get_inclination() must be implemented in subclasses."
        )

    def get_right_ascension_ascending_node(self) -> float:
        """
        Gets the right ascension of ascending node (degrees). Must be
        implemented by subclasses, since OrbitBase has no universal way to
        derive it.

        Returns:
            float: the right ascension of ascending node
        """
        raise NotImplementedError(
            "get_right_ascension_ascending_node() must be implemented in subclasses."
        )

    def get_eccentricity(self) -> float:
        """
        Gets the eccentricity (float between 0 and 1). Must be
        implemented by subclasses, since OrbitBase has no universal way to
        derive it.

        Returns:
            float: the eccentricity
        """
        raise NotImplementedError(
            "get_eccentricity() must be implemented in subclasses."
        )

    def get_perigee_argument(self) -> float:
        """
        Gets the perigee argument (degrees). Must be implemented by
        subclasses, since OrbitBase has no universal way to derive it.

        Returns:
            float: the perigee argument
        """
        raise NotImplementedError(
            "get_perigee_argument() must be implemented in subclasses."
        )

    def to_gp_orbit(self, lazy_load: bool | None = None) -> GeneralPerturbationsOrbit:
        """
        Converts this orbit to a general perturbations orbit representation.
        Lazy-loads a previously-computed conversion if available, since
        `_compute_gp_orbit()` can be expensive (e.g. requires SGP4 fitting).
        Subclasses must implement `_compute_gp_orbit()` rather than
        overriding this method directly.

        Args:
            lazy_load (bool | None): True, if this gp orbit should be lazy-loaded.

        Returns:
            GeneralPerturbationsOrbit: the general perturbations orbit
        """
        if lazy_load is None:
            lazy_load = config.rc.gp_orbit_lazy_load
        if lazy_load:
            gp_orbit = self.__dict__.get("gp_orbit")
        else:
            gp_orbit = None
        if gp_orbit is None:
            gp_orbit = self._compute_gp_orbit()
            self.__dict__["gp_orbit"] = gp_orbit  # type: ignore
        return gp_orbit

    def _compute_gp_orbit(self) -> GeneralPerturbationsOrbit:
        """
        Computes a general perturbations orbit representation of this
        orbit, without caching, by constructing an intermediate
        KeplerianOrbit from this orbit's getter methods. Subclasses for
        which this generic Keplerian-element construction does not apply
        (e.g. KeplerianOrbit itself, which would otherwise recurse
        infinitely) must override this method.

        Returns:
            GeneralPerturbationsOrbit: the general perturbations orbit
        """
        # deferred import to avoid a circular import (keplerian.py imports
        # OrbitBase from this module); the cycle is real but harmless since
        # this import only runs after both modules have finished loading
        # pylint: disable-next=import-outside-toplevel,cyclic-import
        from .keplerian import KeplerianOrbit

        return KeplerianOrbit(
            semimajor_axis=self.get_semimajor_axis(),
            inclination=self.get_inclination(),
            right_ascension_ascending_node=self.get_right_ascension_ascending_node(),
            eccentricity=self.get_eccentricity(),
            perigee_argument=self.get_perigee_argument(),
            true_anomaly=self.get_true_anomaly(),
            epoch=self.get_epoch(),
        ).to_gp_orbit()

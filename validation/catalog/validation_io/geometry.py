"""
Geometry shared by the validation scripts: coordinates, times of the JPSS
products, the cone of conical scanners, and the lattices of samples that scanners
image, in which observed passes are found.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import numpy as np
import pandas as pd
import shapely
from pyproj import Transformer
from scipy.spatial import cKDTree

from tatc.constants import EARTH_ROTATION_RATE

TO_ECEF = Transformer.from_crs("EPSG:4979", "EPSG:4978", always_xy=True)
TO_GEODETIC = Transformer.from_crs("EPSG:4978", "EPSG:4979", always_xy=True)
EARTH_ROTATION = 7.2921150e-5
"""
The Earth's rotation rate (rad/s) with which the scripts derive inertial
velocities from products' Earth-fixed ones (the IERS conventional value; TAT-C's
`EARTH_ROTATION_RATE` differs from it in the ninth digit).
"""
IET_EPOCH = pd.Timestamp("1958-01-01", tz="UTC")
LEAP_SECONDS = 37


def transform_geometry(transform, geometry):
    """
    Transforms the coordinates of a geometry by a function of their x and y
    arrays (such as a pyproj transform), as the deprecated
    `shapely.ops.transform` did.
    """
    return shapely.transform(
        geometry, lambda xy: np.column_stack(transform(xy[:, 0], xy[:, 1]))
    )


def iet_to_utc(iet) -> pd.DatetimeIndex:
    """Converts IDPS Epoch Times (microseconds since 1958, with leap seconds) to UTC."""
    return IET_EPOCH + pd.to_timedelta(
        np.asarray(iet) - LEAP_SECONDS * 1_000_000, unit="us"
    )


def ecef(lon, lat, height=0.0) -> np.ndarray:
    """Gets the Earth-fixed position (m) of geodetic coordinates."""
    return np.stack(
        TO_ECEF.transform(lon, lat, np.broadcast_to(height, np.shape(lon))), axis=-1
    )


def up(lon, lat) -> np.ndarray:
    """Gets the geodetic up unit vector."""
    lon, lat = np.radians(lon), np.radians(lat)
    return np.stack(
        [np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=-1
    )


def normalize(x: np.ndarray) -> np.ndarray:
    """Gets unit vectors along the last axis."""
    return x / np.linalg.norm(x, axis=-1, keepdims=True)


def inertial_velocity(position: np.ndarray, velocity: np.ndarray) -> np.ndarray:
    """
    Gets the inertial velocity (m/s) of Earth-fixed positions and velocities, with
    TAT-C's rotation rate.
    """
    return velocity + EARTH_ROTATION_RATE * np.stack(
        [-position[..., 1], position[..., 0], np.zeros(position.shape[:-1])], axis=-1
    )


def fit_cone(position, velocity, lon, lat) -> dict[str, float]:
    """
    Measures the cone of a conical scanner from the Earth-fixed positions (m) and
    velocities (m/s) of its scans and the geodetic positions of their footprints
    (scans by samples): the median angle of the lines of sight from geodetic nadir,
    and the scan azimuths (deg, positive to the left) of the first and last samples
    from the horizontal direction of the inertial velocity, with the center and half
    width of the sector they span (medians over scans), and the median azimuth step
    between samples. The center from the Earth-fixed velocity is also given, for
    comparison.
    """
    lines = normalize(ecef(lon, lat) - position[:, np.newaxis])
    sub_lon, sub_lat, _ = TO_GEODETIC.transform(*position.T)
    nadir = -up(sub_lon, sub_lat)
    horizontal = (
        lines
        - np.sum(lines * nadir[:, np.newaxis], axis=-1, keepdims=True)
        * nadir[:, np.newaxis]
    )

    def azimuth(v):
        forward = normalize(v - np.sum(v * nadir, axis=-1, keepdims=True) * nadir)
        left = np.cross(forward, nadir)
        return np.degrees(
            np.arctan2(
                np.sum(horizontal * left[:, np.newaxis], axis=-1),
                np.sum(horizontal * forward[:, np.newaxis], axis=-1),
            )
        )

    inertial = azimuth(inertial_velocity(position, velocity))
    earth = azimuth(velocity)
    return {
        "cone_angle": float(
            np.median(
                np.degrees(np.arccos(np.sum(lines * nadir[:, np.newaxis], axis=-1)))
            )
        ),
        "first_azimuth": float(np.median(inertial[:, 0])),
        "last_azimuth": float(np.median(inertial[:, -1])),
        "scan_center_azimuth": float(np.median(inertial[:, -1] + inertial[:, 0]) / 2),
        "scan_half_width": float(
            np.median(np.abs(inertial[:, -1] - inertial[:, 0])) / 2
        ),
        "scan_sample_spacing": float(np.median(np.diff(inertial, axis=1))),
        "scan_center_azimuth_earth_fixed": float(
            np.median(earth[:, -1] + earth[:, 0]) / 2
        ),
    }


def in_lattice_cells(cells, is_next_contiguous, points, neighbors=2) -> np.ndarray:
    """
    Gets whether unit vectors lie in the cells of a lattice of unit vectors, with rows
    of scans and columns of samples, given whether each scan is followed by its next.
    The cells around a number of the nearest vertices are tested, which needs to be
    larger for elongated cells.
    """
    rows, columns = cells.shape[:2]
    _, nearest = cKDTree(cells.reshape(-1, 3)).query(points, k=neighbors, workers=-1)
    inside = np.zeros(len(points), dtype=bool)
    for column in range(neighbors):
        k, j = np.divmod(nearest[:, column], columns)
        for dk in (-1, 0):
            for dj in (-1, 0):
                kk, jj = k + dk, j + dj
                valid = (kk >= 0) & (kk < rows - 1) & (jj >= 0) & (jj < columns - 1)
                kk = np.clip(kk, 0, rows - 2)
                jj = np.clip(jj, 0, columns - 2)
                valid &= is_next_contiguous[kk]
                corners = [
                    cells[kk, jj],
                    cells[kk, jj + 1],
                    cells[kk + 1, jj + 1],
                    cells[kk + 1, jj],
                ]
                signs = np.stack(
                    [
                        np.sum(np.cross(a, b) * points, axis=-1)
                        for a, b in zip(corners, corners[1:] + corners[:1])
                    ]
                )
                inside |= valid & (
                    np.all(signs >= 0, axis=0) | np.all(signs <= 0, axis=0)
                )
    return inside


class ScanLattice:
    """
    A day of a scanner's scans, whose samples (beams, fields of regard, or frames)
    form a lattice of geodetic up unit vectors: `unit` (scans by samples), with
    each scan's `time` and whether it is followed by the next (`is_next_contiguous`).
    A product gives the times of its samples (`sample_seconds`, scans by samples,
    from any origin) and the time of a pass between two scans (`pass_time`).
    """

    unit: np.ndarray
    time: pd.DatetimeIndex
    is_next_contiguous: np.ndarray

    @property
    def sample_seconds(self) -> np.ndarray:
        """The time of each sample (s, from any origin), scans by samples."""
        raise NotImplementedError

    def pass_time(self, before: int, sample: int, fraction: float) -> pd.Timestamp:
        """
        Gets the time of a pass a fraction of the way from a sample of a scan to the
        same sample of the next scan.
        """
        raise NotImplementedError

    def in_cells(self, points: np.ndarray, first: int, last: int) -> np.ndarray:
        """Gets whether unit vectors lie in cells of samples on a range of scans."""
        return in_lattice_cells(
            self.unit[first : last + 1], self.is_next_contiguous[first:last], points
        )

    def observe_passes(self, points, radius: float) -> pd.DataFrame:
        """
        Finds the passes on which the scanner observed points, and their times: a
        point is observed on a pass if it lies in a cell of the lattice near the
        sample nearest it (within `radius`, in km), at the time interpolated between
        the scans before and after it.
        """
        point_unit = up(points.geometry.x.to_numpy(), points.geometry.y.to_numpy())
        flat = self.unit.reshape(-1, 3)
        flat_seconds = self.sample_seconds.ravel()
        samples = self.unit.shape[1]
        tree = cKDTree(flat)
        rows = []
        for p, hits in enumerate(
            tree.query_ball_point(
                point_unit, r=2 * np.sin(radius / 6371 / 2), workers=-1
            )
        ):
            if not hits:
                continue
            hits = np.array(hits)
            hits = hits[np.argsort(flat_seconds[hits])]
            for group in np.split(
                hits, np.nonzero(np.diff(flat_seconds[hits]) > 1200)[0] + 1
            ):
                k, j = np.divmod(
                    group[
                        np.argmin(np.linalg.norm(flat[group] - point_unit[p], axis=-1))
                    ],
                    samples,
                )
                first, last = max(k - 1, 0), min(k + 1, len(self.time) - 1)
                if not self.in_cells(point_unit[[p]], first, last)[0]:
                    continue
                # interpolate between the scans before and after the point
                before = (
                    k - 1
                    if np.sum(
                        (point_unit[p] - self.unit[k, j])
                        * (self.unit[k, j] - self.unit[max(k - 1, 0), j])
                    )
                    < 0
                    else k
                )
                before = int(np.clip(before, 0, len(self.time) - 2))
                a, b = self.unit[before, j], self.unit[before + 1, j]
                fraction = np.clip(
                    np.sum((point_unit[p] - a) * (b - a)) / np.sum((b - a) ** 2), 0, 1
                )
                rows.append(
                    {"point_id": p, "time": self.pass_time(before, j, fraction)}
                )
        return pd.DataFrame(rows).astype({"point_id": "int64"})

    def observe_orbit(
        self,
        orbit,
        grid,
        margin: float,
        begin: pd.Timestamp | None = None,
        times: pd.DatetimeIndex | None = None,
        floor: bool = False,
    ):
        """
        Gets the cells of a grid compared over one orbit from noon (or from `begin`),
        within 300 km of the scans and more than a `margin` (s) from the orbit's ends,
        and which of them were observed. The orbit's scans are selected by `time`;
        its start, end, and the times of its scans are those of `times` (by default,
        `time`), floored to the microsecond if `floor`.

        Returns:
            tuple: The start and end of the orbit, and whether each cell is compared
            and was observed.
        """
        lon, lat, _ = grid
        unit = up(lon, lat)
        times = self.time if times is None else times
        if begin is None:
            begin = self.time[0].floor("D") + pd.Timedelta(hours=12)
        period = pd.Timedelta(days=1 / orbit.elements[0].mean_motion * 360 / 86400)
        scans = np.nonzero((self.time >= begin) & (self.time < begin + period))[0]
        first, last = times[scans[0]], times[scans[-1]]
        distance, nearest = cKDTree(self.unit[scans].reshape(-1, 3)).query(
            unit, workers=-1
        )
        elapsed = (times[scans] - first).total_seconds().to_numpy()
        nearest_time = elapsed[nearest // self.unit.shape[1]]
        compared = (
            (distance < 2 * np.sin(300 / 6371 / 2))
            & (nearest_time > margin)
            & (nearest_time < (last - first).total_seconds() - margin)
        )
        observed = np.zeros(len(unit), dtype=bool)
        observed[compared] = self.in_cells(unit[compared], scans[0], scans[-1])
        if floor:
            first, last = first.floor("us"), last.floor("us")
        return first.to_pydatetime(), last.to_pydatetime(), compared, observed

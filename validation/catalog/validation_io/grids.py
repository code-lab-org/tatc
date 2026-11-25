"""
Scoring the sample size configurations of geostationary imagers, whose pixels are a
fixed grid of scan angles: the projected pixels against the positions of the
product's grid (seen from the nominal position, to which the products are
navigated), the spacing of pixels against that of the grid, and the published sizes.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import pathlib
from dataclasses import dataclass
from datetime import datetime

import numpy as np
import pandas as pd

from tatc.catalog.metrics import summarize
from tatc.schemas import GeosynchronousOrbit

from .geometry import ecef
from .records import write_results
from .scoring import compute_published_errors, compute_spacing_ratios, score_published


@dataclass(frozen=True)
class FixedGrid:
    """
    The sample size configuration of a geostationary imager: its key, the instrument
    whose pixels are compared, the published pixels (see
    `scoring.compute_published_errors`), whether they test the gate (sizes at the
    sub-satellite point, which the grid's angle defines, are compared without it),
    the reference, the fixture, and the samples of pixels.
    """

    key: str
    instrument: str
    published_footprints: tuple
    reference: str
    fixture_path: pathlib.Path
    gate_published: bool = False
    samples: int = 2000
    fixture_samples: int = 200


def build_satellite(entry, key: str, longitude: float, time: datetime):
    """
    Builds the satellite of a configuration at the nominal position of a fixed grid
    (a geosynchronous orbit at its longitude, at the time of the data).
    """
    orbit = GeosynchronousOrbit(longitude=longitude, epoch=time, inclination=0)
    return entry.build(configuration=key, orbit=orbit).to_satellite(
        orbit, name=entry.name
    )


def project_pixels(satellite, instrument: str, time: datetime, rows, columns, size):
    """
    Projects pixels of a grid of `size` rows and columns (rows from the north,
    columns from the west), as TAT-C numbers them from the right of the motion
    (east, in the inertial frame of a geostationary platform): from the south and
    from the east. Gets their latitudes and longitudes (deg).
    """
    names = [i.name for i in satellite.instruments]
    view = satellite.get_instrument(names.index(instrument))
    track = satellite.orbit.to_gp_orbit().get_orbit_track(time)
    lat, lon = np.empty(len(rows)), np.empty(len(rows))
    for k, (row, column) in enumerate(zip(rows, columns)):
        position = view.compute_projected_pixel_position(
            track, size - 1 - int(row), size - 1 - int(column)
        )
        lat[k], lon[k] = position.latitude.degrees, position.longitude.degrees
    return lat, lon


def measure_grid_spacing(lat: np.ndarray, lon: np.ndarray, rows, columns):
    """
    Measures the spacing (km) of a grid's pixels, east-west and north-south, from
    half the distance between the pixels on either side. Gets rows of latitude,
    longitude, east-west (along track), and north-south (across track) spacing.
    """
    east_west = (
        np.linalg.norm(
            ecef(lon[rows, columns + 1], lat[rows, columns + 1])
            - ecef(lon[rows, columns - 1], lat[rows, columns - 1]),
            axis=-1,
        )
        / 2e3
    )
    north_south = (
        np.linalg.norm(
            ecef(lon[rows + 1, columns], lat[rows + 1, columns])
            - ecef(lon[rows - 1, columns], lat[rows - 1, columns]),
            axis=-1,
        )
        / 2e3
    )
    return np.stack(
        [lat[rows, columns], lon[rows, columns], east_west, north_south], axis=-1
    )


def validate_fixed_grid(  # pylint: disable=too-many-arguments,too-many-locals
    grid: FixedGrid,
    run,
    lat: np.ndarray,
    lon: np.ndarray,
    longitude: float,
    time: datetime,
) -> None:
    """
    Scores the sample size configuration of a geostationary imager in a run (see
    `cli.parse_run`) against a product's grid of pixel latitudes and longitudes
    (rows from the north, columns from the west; `NaN` off the Earth's disk) of a
    scan at a time: the published sizes, the projected pixels of a sample of the
    disk from the nominal position (the grid's longitude), and the spacing of the
    pixels. Writes its results if the run writes.
    """
    entry = run.entry
    within, largest = score_published(
        compute_published_errors(entry, grid.key, grid.published_footprints)
    )
    print(
        f"published dimensions: {within:.1f}% within 3%, largest error {largest:.2f}%"
        + ("" if grid.gate_published else " (not gated)")
    )
    published = {"sample_size_within_3pct": within} if grid.gate_published else {}
    size = lat.shape[0]
    # pixels of the disk whose neighbors are on it
    on_disk = np.isfinite(lat) & np.isfinite(lon)
    inner = np.zeros_like(on_disk)
    inner[1:-1, 1:-1] = (
        on_disk[1:-1, 1:-1]
        & on_disk[:-2, 1:-1]
        & on_disk[2:, 1:-1]
        & on_disk[1:-1, :-2]
        & on_disk[1:-1, 2:]
    )
    candidates = np.argwhere(inner)
    rng = np.random.default_rng(5)
    picks = candidates[
        np.sort(rng.choice(len(candidates), grid.samples, replace=False))
    ]
    rows, columns = picks[:, 0], picks[:, 1]
    satellite = build_satellite(entry, grid.key, longitude, time)
    predicted_lat, predicted_lon = project_pixels(
        satellite, grid.instrument, time, rows, columns, size
    )
    difference = (
        np.linalg.norm(
            ecef(predicted_lon, predicted_lat)
            - ecef(lon[rows, columns], lat[rows, columns]),
            axis=-1,
        )
        / 1e3
    )
    stats = summarize(difference)
    print(
        f"pixel position: median {stats['median']:.3f} km, "
        f"95th percentile {stats['p95']:.3f} km"
    )
    spacing = measure_grid_spacing(lat, lon, rows, columns)
    times = pd.DatetimeIndex([pd.Timestamp(time)] * len(spacing))
    along, cross = compute_spacing_ratios(
        satellite, [grid.instrument] * len(spacing), times, spacing
    )
    print(
        f"spacing, predicted / measured: {along:.4f} east-west, {cross:.4f} "
        "north-south"
    )
    if not run.write:
        return
    reduced = rng.choice(len(rows), grid.fixture_samples, replace=False)
    fixture_along, fixture_cross = compute_spacing_ratios(
        satellite, [grid.instrument] * len(reduced), times[reduced], spacing[reduced]
    )
    write_results(
        run,
        grid.key,
        reference=grid.reference,
        metrics=[
            *({"name": name, "value": value} for name, value in published.items()),
            {"name": "sample_size_largest_error", "value": largest},
            {"name": "pixel_position_median", "value": stats["median"]},
            {"name": "pixel_position_p95", "value": stats["p95"]},
            {"name": "spacing_along_ratio_median", "value": along},
            {"name": "spacing_cross_ratio_median", "value": cross},
        ],
        fixture={
            "time": np.array(time.isoformat()),
            "longitude": np.array(longitude),
            "size": np.array(size),
            "pixels": np.stack(
                [
                    rows[reduced],
                    columns[reduced],
                    lat[rows[reduced], columns[reduced]],
                    lon[rows[reduced], columns[reduced]],
                ],
                axis=-1,
            ),
            "spacing": spacing[reduced],
        },
        fixture_path=grid.fixture_path,
        fixture_metrics={
            **(published or {"sample_size_largest_error": largest}),
            "pixel_position_median": summarize(difference[reduced])["median"],
            "spacing_along_ratio_median": fixture_along,
            "spacing_cross_ratio_median": fixture_cross,
        },
        period=f"{time:%Y-%m-%dT%H:%MZ}",
    )

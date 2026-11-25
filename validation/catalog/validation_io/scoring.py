"""
Scoring through public TAT-C methods: observations of points, sample sizes against
published footprints, predicted and observed passes, and swath areas on a grid.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import concurrent.futures
from datetime import timedelta

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import contains_xy, prepare
from shapely.ops import unary_union

from tatc.analysis import (
    collect_observations,
    compute_ground_track,
    compute_sample_sizes,
)
from tatc.catalog.metrics import compute_mask_iou, compute_pass_scores, summarize
from tatc.schemas import CircularOrbit, Point
from tatc.utils import VelocityFrame, compute_projected_ray_position

PUBLISHED_EPOCH = pd.Timestamp("2026-10-01", tz="UTC").to_pydatetime()
SIZE_TOLERANCE = 0.03
GRID_STEP = 0.1


def to_observations(satellite, instruments, epochs, lon, lat) -> gpd.GeoDataFrame:
    """Builds observations of points by named instruments of a satellite."""
    return gpd.GeoDataFrame(
        {
            "satellite": satellite.name,
            "instrument": instruments,
            "epoch": pd.DatetimeIndex(epochs),
        },
        geometry=gpd.points_from_xy(np.ravel(lon), np.ravel(lat)),
        crs="EPSG:4326",
    )


def compute_published_errors(entry, key: str, footprints) -> np.ndarray:
    """
    Gets the relative errors of sample sizes against published footprints, given as
    rows of instrument name, altitude (km), inclination (deg), scan angle (deg), and
    the sizes (km) along and across track. Gives two errors for each row.
    """
    errors = []
    for name, altitude, inclination, angle, along, cross in footprints:
        orbit = CircularOrbit(
            altitude=altitude * 1e3, inclination=inclination, epoch=PUBLISHED_EPOCH
        )
        # a day, for a configuration that schedules its attitude
        satellite = entry.build(
            configuration=key,
            orbit=orbit,
            start=PUBLISHED_EPOCH,
            end=PUBLISHED_EPOCH + timedelta(days=1),
        ).to_satellite(orbit, name=entry.name)
        time = PUBLISHED_EPOCH + timedelta(minutes=10)
        target = compute_projected_ray_position(
            orbit.to_gp_orbit().get_orbit_track(time),
            0,
            0,
            angle,
            0,
            velocity_frame=VelocityFrame.INERTIAL,
        )
        sizes = compute_sample_sizes(
            to_observations(
                satellite,
                [name],
                [time],
                [target.longitude.degrees],
                [target.latitude.degrees],
            ),
            satellite,
        )
        errors += [
            sizes.sample_along.iloc[0] / 1e3 / along - 1,
            sizes.sample_cross.iloc[0] / 1e3 / cross - 1,
        ]
        print(f"{name}, {angle} deg: {errors[-2]:+.2%} along, {errors[-1]:+.2%} across")
    return np.array(errors)


def score_published(errors: np.ndarray) -> tuple[float, float]:
    """Gets the percent of dimensions within tolerance and the largest error (%)."""
    return (
        100 * float(np.mean(np.abs(errors) <= SIZE_TOLERANCE)),
        100 * float(np.abs(errors).max()),
    )


def compute_spacing_ratios(satellite, names, epochs, rows) -> tuple[float, float]:
    """
    Gets the median ratios of predicted to measured spacing along and across track,
    for rows of latitude, longitude, along-track (km), and cross-track spacing (km).
    """
    sizes = compute_sample_sizes(
        to_observations(satellite, names, epochs, rows[:, 1], rows[:, 0]), satellite
    )
    return (
        float(np.nanmedian(sizes.spacing_along.to_numpy() / 1e3 / rows[:, 2])),
        float(np.nanmedian(sizes.spacing_cross.to_numpy() / 1e3 / rows[:, 3])),
    )


def predict_passes(points, satellite, index, start, end, gaps) -> pd.DataFrame:
    """Predicts the observations of points, outside periods without data."""

    def predict(item):
        point_id, geometry = item
        return collect_observations(
            Point(latitude=geometry.y, longitude=geometry.x),
            satellite,
            start,
            end,
            index,
        ).assign(point_id=point_id)

    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        predicted = pd.concat(list(pool.map(predict, enumerate(points.geometry))))
    predicted = predicted.reset_index(drop=True)
    predicted["point_id"] = predicted.point_id.astype("int64")
    in_gap = np.zeros(len(predicted), dtype=bool)
    for first, last in gaps:
        in_gap |= ((predicted.epoch >= first) & (predicted.epoch <= last)).to_numpy()
    return predicted[~in_gap]


def match_passes(predicted, observed) -> dict:
    """Compares predicted passes with observed ones, with `point_id` and `time`."""
    merged = pd.merge_asof(
        observed.assign(time=observed.time.dt.as_unit("ns")).sort_values("time"),
        predicted[["point_id", "epoch"]]
        .assign(epoch=predicted.epoch.dt.as_unit("ns"))
        .sort_values("epoch"),
        left_on="time",
        right_on="epoch",
        by="point_id",
        direction="nearest",
        tolerance=pd.Timedelta(minutes=10),
    )
    hit = merged.epoch.notna()
    differences = (merged.epoch - merged.time).dt.total_seconds()[hit]
    return {
        **compute_pass_scores(int(hit.sum()), len(observed), len(predicted)),
        "epoch_p95": summarize(differences)["p95"],
        "epoch_median": summarize(differences)["median"],
        "observed": len(observed),
        "predicted": len(predicted),
        "matched": int(hit.sum()),
    }


def score_passes(points, satellite, index, observed, start, end, gaps) -> dict:
    """
    Compares predicted and observed passes of points, given observations with
    `point_id` and `time`, and the periods without data as pairs of times.
    """
    return match_passes(
        predict_passes(points, satellite, index, start, end, gaps), observed
    )


def recompute_passes(satellite, index, day, passes) -> dict[str, float]:
    """
    Recomputes observed passes, given as rows of latitude, longitude, and seconds from
    the start of a day. Gets the percent found and the epochs' 95th percentile difference.
    """
    differences = []
    for latitude, longitude, seconds in passes:
        time = day + timedelta(seconds=float(seconds))
        found = collect_observations(
            Point(latitude=float(latitude), longitude=float(longitude)),
            satellite,
            time - timedelta(minutes=10),
            time + timedelta(minutes=10),
            index,
        )
        differences.append(
            min(
                ((e - pd.Timestamp(time)).total_seconds() for e in found.epoch), key=abs
            )
            if not found.empty
            else np.nan
        )
    differences = np.array(differences, dtype=float)
    return {
        "passes_found": 100 * float(np.isfinite(differences).mean()),
        "epoch_p95": float(summarize(differences)["p95"]),
    }


def build_grid() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Gets the longitudes, latitudes, and areas of the cells of a global grid."""
    lon, lat = np.meshgrid(
        np.arange(-180 + GRID_STEP / 2, 180, GRID_STEP),
        np.arange(-90 + GRID_STEP / 2, 90, GRID_STEP),
    )
    weights = np.cos(np.radians(lat)) * (GRID_STEP * 111.32) ** 2
    return lon.ravel(), lat.ravel(), weights.ravel()


def score_swath_area(
    satellite_for, index, first, last, time_steps, grid, compared, observed
) -> float:
    """
    Gets the lowest IoU of a swath with the observed area, over time steps, for an
    instrument's index, or for the union of a satellite's instruments if `None`.
    """
    lon, lat, weights = grid
    scores = []
    for time_step in time_steps:
        count = int((last - first).total_seconds() // time_step)
        times = [first + timedelta(seconds=time_step * (i + 0.5)) for i in range(count)]
        satellite = satellite_for(time_step)
        indices = range(len(satellite.instruments)) if index is None else [index]
        geometry = unary_union(
            [
                compute_ground_track(
                    satellite,
                    times,
                    instrument_index=i,
                    time_step=timedelta(seconds=time_step),
                ).geometry.iloc[0]
                for i in indices
            ]
        )
        prepare(geometry)
        predicted = np.zeros(len(lon), dtype=bool)
        predicted[compared] = contains_xy(geometry, lon[compared], lat[compared])
        scores.append(compute_mask_iou(observed, predicted, weights))
        print(f"  swath IoU at {time_step:.4g} s: {scores[-1]:.5f}", flush=True)
    return min(scores)


def pick_passes(points, observed, day, count: int, seed=1) -> np.ndarray:
    """
    Picks observed passes for a test fixture, at random (from a seed or a random
    generator): rows of the point's latitude and longitude, and the seconds of the
    pass from the start of the day.
    """
    rng = np.random.default_rng(seed)
    picks = observed.iloc[rng.choice(len(observed), count, replace=False)]
    return np.stack(
        [
            points.geometry.y.to_numpy()[picks.point_id],
            points.geometry.x.to_numpy()[picks.point_id],
            (picks.time - pd.Timestamp(day)).dt.total_seconds().to_numpy(),
        ],
        axis=-1,
    )

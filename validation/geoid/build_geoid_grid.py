"""
Builds the grid of EGM2008 geoid heights of `tatc.utils.compute_geoid_height`
(`src/tatc/resources/egm2008_15.npz`) from the 2.5-minute EGM2008 grid of the
National Geospatial-Intelligence Agency (public domain), as distributed by PROJ
(https://cdn.proj.org/us_nga_egm08_25.tif), and reports the errors of grids of
other spacings, interpolated bilinearly, at the nodes of the 2.5-minute grid.
With PROJ's grid in its data directory, it also compares the bundled grid with
PROJ's transformation of WGS 84 heights to EGM2008 heights.

Usage:

    python validation/geoid/build_geoid_grid.py path/to/us_nga_egm08_25.tif

Requires `rasterio` (`pip install tatc[preprocess]`).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import os

import numpy as np
import rasterio
from pyproj import Transformer, datadir

SOURCE_SPACING = 2.5 / 60  # degrees
STEP = 6  # 15 minutes (0.25 degrees)
OUTPUT = os.path.join(
    os.path.dirname(__file__), "..", "..", "src", "tatc", "resources", "egm2008_15.npz"
)


def read_source(path):
    """
    Reads the 2.5-minute grid: rows from 90 to -90 degrees latitude, columns
    from -180 degrees longitude (with the column at 180 degrees appended).
    """
    with rasterio.open(path) as dataset:
        # the grid's nodes are at its pixels' points
        assert dataset.tags()["AREA_OR_POINT"] == "Point"
        assert np.isclose(dataset.transform.c + dataset.transform.a / 2, -180)
        assert np.isclose(dataset.transform.f + dataset.transform.e / 2, 90)
        heights = dataset.read(1).astype(float)
    return np.concatenate([heights, heights[:, :1]], axis=1)


def report_errors(heights, step):
    """
    Reports the errors (meters) of the grid subsampled every `step` nodes,
    interpolated bilinearly at every node, weighted by area.
    """
    coarse = heights[::step, ::step]
    rows, columns = heights.shape
    i = np.arange(rows) / step
    j = np.arange(columns) / step
    i0 = np.minimum(np.floor(i).astype(int), coarse.shape[0] - 2)
    j0 = np.minimum(np.floor(j).astype(int), coarse.shape[1] - 2)
    u, v = (i - i0)[:, None], (j - j0)[None, :]
    interpolated = (
        (1 - u) * (1 - v) * coarse[i0][:, j0]
        + (1 - u) * v * coarse[i0][:, j0 + 1]
        + u * (1 - v) * coarse[i0 + 1][:, j0]
        + u * v * coarse[i0 + 1][:, j0 + 1]
    )
    error = np.abs(interpolated - heights).ravel()
    latitude = 90 - np.arange(rows) * SOURCE_SPACING
    weight = np.repeat(np.cos(np.radians(latitude)), columns)
    order = np.argsort(error)
    cumulative = np.cumsum(weight[order]) / weight.sum()
    rms = np.sqrt(np.sum(weight * error**2) / weight.sum())
    p99 = error[order][np.searchsorted(cumulative, 0.99)]
    print(
        f"{step * SOURCE_SPACING:.4f} degrees: rms {rms:.3f} m, 99th percentile "
        f"{p99:.2f} m, max {error.max():.2f} m, {coarse[:, :-1].size * 2 / 1024:.0f} KB"
    )


def compare_with_proj(path):
    """
    Compares the bundled grid with PROJ's transformation of WGS 84 heights to
    EGM2008 heights, with the 2.5-minute grid, at random points.
    """
    from tatc.utils import compute_geoid_height  # pylint: disable=C0415

    datadir.append_data_dir(os.path.dirname(os.path.abspath(path)))
    transformer = Transformer.from_crs("EPSG:4979", "EPSG:4326+3855", always_xy=True)
    rng = np.random.default_rng(0)
    longitude = rng.uniform(-180, 180, 10000)
    latitude = np.degrees(np.arcsin(rng.uniform(-1, 1, 10000)))
    _, _, height = transformer.transform(longitude, latitude, np.zeros(10000))
    difference = compute_geoid_height(longitude, latitude) + np.asarray(height)
    print(
        f"bundled grid - PROJ: rms {np.sqrt(np.mean(difference**2)):.3f} m, "
        f"max {np.abs(difference).max():.2f} m"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", help="path to us_nga_egm08_25.tif")
    args = parser.parse_args()
    heights = read_source(args.source)
    for step in (3, 6, 12, 24):
        report_errors(heights, step)
    coarse = heights[::STEP, :-1:STEP]
    np.savez_compressed(
        OUTPUT,
        height_cm=np.round(coarse * 100).astype(np.int16),
        spacing_deg=np.array(STEP * SOURCE_SPACING),
    )
    print(f"wrote {os.path.abspath(OUTPUT)} ({os.path.getsize(OUTPUT) / 1024:.0f} KB)")
    compare_with_proj(args.source)


if __name__ == "__main__":
    main()

"""
Scoring the sample size configurations of beam instruments (lidars, profilers, and
altimeters), whose footprints point near nadir: their sizes against published
footprints, and the spacing of their samples along track against that measured in
a product.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import pathlib
import warnings
from dataclasses import dataclass
from datetime import timedelta

import numpy as np
import pandas as pd

from .geometry import ecef
from .records import write_results
from .scoring import compute_published_errors, compute_spacing_ratios, score_published


@dataclass(frozen=True)
class Footprints:
    """
    The sample size configuration of a beam instrument: its key, the published
    footprints (see `scoring.compute_published_errors`), the reference, the
    fixture, and whether the published footprints test the gate (those too few, or
    from which the footprint's angle comes, are compared without it).
    """

    key: str
    published_footprints: tuple
    reference: str
    fixture_path: pathlib.Path
    gate_published: bool = True
    fixture_spacing: int = 200


def measure_along_spacing(
    seconds: np.ndarray,
    lat: np.ndarray,
    lon: np.ndarray,
    period: float,
    samples: int = 1,
) -> np.ndarray:
    """
    Measures the spacing (km) of a beam's samples along its track, from half the
    distance between the kept samples before and after each, where both are one
    period (s) of kept samples away, divided by the instrument's samples between
    kept ones (of a product reduced to every so many, or of averages of so many).
    Gets rows of latitude, longitude, along-track spacing (km), a cross-track
    spacing (`NaN`, as a beam has one sample across track), and seconds.
    """
    xyz = ecef(lon, lat)
    gaps = np.diff(seconds)
    regular = np.abs(gaps - period) < 0.25 * period
    middle = np.nonzero(regular[:-1] & regular[1:])[0] + 1
    along = np.linalg.norm(xyz[middle + 1] - xyz[middle - 1], axis=-1) / (2e3 * samples)
    return np.stack(
        [
            lat[middle],
            lon[middle],
            along,
            np.full(len(middle), np.nan),
            seconds[middle],
        ],
        axis=-1,
    )


def validate_footprints(  # pylint: disable=too-many-arguments
    footprints: Footprints,
    run,
    satellite,
    names: list[str],
    rows: np.ndarray,
    tle,
) -> None:
    """
    Scores the sample size configuration of a beam instrument in a run (see
    `cli.parse_run`): the sizes of its footprints against the published ones, and
    the median ratio of predicted to measured spacing along track, for rows of
    latitude, longitude, along-track spacing (km), cross-track spacing (`NaN`), and
    seconds from the start of the day (see `measure_along_spacing`), each of the
    instrument named in `names`. Writes its results if the run writes.
    """
    entry = run.entry
    within, largest = score_published(
        compute_published_errors(entry, footprints.key, footprints.published_footprints)
    )
    print(
        f"published dimensions: {within:.1f}% within 3%, largest error {largest:.2f}%"
        + ("" if footprints.gate_published else " (not gated)")
    )
    published = {"sample_size_within_3pct": within} if footprints.gate_published else {}
    times = pd.DatetimeIndex(
        [pd.Timestamp(run.day) + timedelta(seconds=float(s)) for s in rows[:, 4]]
    )
    along, _ = _ratios(satellite, names, times, rows)
    error = 100 * abs(along - 1)
    print(f"spacing along track, predicted / measured: {along:.4f} ({error:.2f}%)")
    if not run.write:
        return
    rng = np.random.default_rng(4)
    picks = np.sort(
        rng.choice(len(rows), min(footprints.fixture_spacing, len(rows)), replace=False)
    )
    fixture_along, _ = _ratios(
        satellite, [names[i] for i in picks], times[picks], rows[picks]
    )
    write_results(
        run,
        footprints.key,
        reference=footprints.reference,
        metrics=[
            *({"name": name, "value": value} for name, value in published.items()),
            {"name": "sample_size_largest_error", "value": largest},
            {"name": "spacing_along_ratio_median", "value": along},
            {"name": "spacing_along_error", "value": error},
        ],
        transfer_metrics={
            "spacing_along_ratio_median": along,
            "spacing_along_error": error,
        },
        fixture={
            "tle": np.array(tle),
            "day": np.array(run.day.isoformat()),
            "spacing": rows[picks, :4],
            "spacing_seconds": rows[picks, 4],
            "spacing_instruments": np.array([names[i] for i in picks]),
        },
        fixture_path=footprints.fixture_path,
        fixture_metrics={
            **(published or {"sample_size_largest_error": largest}),
            "spacing_along_ratio_median": fixture_along,
        },
    )


def _ratios(satellite, names, times, rows) -> tuple[float, float]:
    """The median ratios of predicted to measured spacing, ignoring cross track."""
    with warnings.catch_warnings():
        # a beam has no spacing across track, whose median is NaN
        warnings.simplefilter("ignore", RuntimeWarning)
        return compute_spacing_ratios(satellite, list(names), times, rows[:, :4])

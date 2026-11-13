"""
Metrics that compare analysis results with reference data.

These depend only on `numpy`, so they can run in unit tests as well as in
validation notebooks.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt


def summarize(values: npt.ArrayLike) -> dict[str, float]:
    """
    Summarizes a sample of residuals, ignoring values that are not finite.

    Args:
        values (npt.ArrayLike): The residuals.

    Returns:
        dict[str, float]: The number of values `n`, `median`, and the 95th
        percentile `p95` and maximum `max` of their absolute values.

    Raises:
        ValueError: If there are no finite values.
    """
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        raise ValueError("no finite values to summarize")
    return {
        "n": int(x.size),
        "median": float(np.median(x)),
        "p95": float(np.percentile(np.abs(x), 95)),
        "max": float(np.abs(x).max()),
    }


def compute_pass_scores(
    matched: int, observed: int, predicted: int
) -> dict[str, float]:
    """
    Computes the agreement between observed and predicted passes.

    Args:
        matched (int): The number of observed passes that were predicted.
        observed (int): The number of observed passes.
        predicted (int): The number of predicted passes.

    Returns:
        dict[str, float]: The percentages `pass_recall`, of observed passes
        that were predicted, and `pass_precision`, of predicted passes that
        were observed.

    Raises:
        ValueError: If there are no passes, or more matches than passes.
    """
    if observed <= 0 or predicted <= 0:
        raise ValueError("observed and predicted passes must be positive")
    if not 0 <= matched <= min(observed, predicted):
        raise ValueError("matched passes must be between 0 and the number of passes")
    return {
        "pass_recall": 100 * matched / observed,
        "pass_precision": 100 * matched / predicted,
    }


def compute_mask_iou(
    first: npt.ArrayLike,
    second: npt.ArrayLike,
    weights: npt.ArrayLike | None = None,
) -> float:
    """
    Computes the intersection over union of two masks on the same grid.

    Args:
        first (npt.ArrayLike): The first boolean mask.
        second (npt.ArrayLike): The second boolean mask.
        weights (npt.ArrayLike | None): The weight of each cell (for example,
            its area), or `None` for equal weights.

    Returns:
        float: The intersection over union.

    Raises:
        ValueError: If the union of the masks is empty.
    """
    a, b = np.asarray(first, dtype=bool), np.asarray(second, dtype=bool)
    w = np.ones(a.shape) if weights is None else np.asarray(weights, dtype=float)
    union = float(np.sum(w[a | b]))
    if union <= 0:
        raise ValueError("the union of the masks is empty")
    return float(np.sum(w[a & b])) / union


def compute_product_limit_error(
    zenith: npt.ArrayLike, observed: npt.ArrayLike, bound: float
) -> float:
    """
    Computes how far the boundary of predicted observations lies from a
    product's declared bound of local zenith angle.

    Args:
        zenith (npt.ArrayLike): The local zenith angle (degrees) of each pixel.
        observed (npt.ArrayLike): Whether each pixel is predicted to be observed.
        bound (float): The declared bound (degrees).

    Returns:
        float: The largest distance (degrees) from the bound of the greatest
        angle observed and the smallest angle not observed.

    Raises:
        ValueError: If the pixels are not on both sides of the boundary.
    """
    angle, seen = np.asarray(zenith, dtype=float), np.asarray(observed, dtype=bool)
    if seen.all() or not seen.any():
        raise ValueError("pixels must include both observed and not observed")
    return float(max(abs(angle[seen].max() - bound), abs(angle[~seen].min() - bound)))

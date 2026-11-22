"""
Schemas of the validation records of catalog entries, which validation scripts
write, and of the periods of data they cite.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Iterable
from datetime import date, datetime, timedelta, timezone
from typing import Annotated

from pydantic import AfterValidator, BaseModel, ConfigDict, Field, model_validator


class _CatalogModel(BaseModel):
    """
    Base class that rejects unknown fields, as catalog files are written by hand.
    """

    model_config = ConfigDict(extra="forbid", frozen=True)


def _parse_day(text: str) -> date:
    """
    Parses an ISO 8601 date, or the UTC date of an ISO 8601 time with a zone.
    """
    try:
        return date.fromisoformat(text)
    except ValueError:
        moment = datetime.fromisoformat(text.replace("Z", "+00:00"))
        if moment.tzinfo is None:
            raise ValueError(f"the time {text!r} needs a time zone") from None
        return moment.astimezone(timezone.utc).date()


def get_period_days(period: str) -> tuple[tuple[date, date], ...]:
    """
    Gets the days of a period: one or more items separated by commas, each a
    date or a time with a zone in ISO 8601 (such as `2026-10-02` or
    `2026-10-03T12:00Z`), or a range of dates (`2026-08-15 to 2026-10-08`).

    Args:
        period (str): The period.

    Returns:
        tuple[tuple[datetime.date, datetime.date], ...]: The first and last day
        (UTC) of each item.

    Raises:
        ValueError: If the period is not of this form.
    """
    intervals = []
    for item in period.split(","):
        first, separator, last = item.strip().partition(" to ")
        start = _parse_day(first)
        end = _parse_day(last) if separator else start
        if end < start:
            raise ValueError(f"the range {item.strip()!r} ends before it starts")
        intervals.append((start, end))
    return tuple(intervals)


def format_period(days: Iterable[date | str]) -> str:
    """
    Formats days as a period (the inverse of `get_period_days`), with runs of
    consecutive days as ranges, such as `2026-10-02, 2026-10-05 to 2026-10-07`.

    Args:
        days (Iterable[datetime.date | str]): The days, as dates or in ISO 8601,
            in any order and with any repeats.

    Returns:
        str: The period.

    Raises:
        ValueError: If there are no days.
    """
    runs: list[list[date]] = []
    for day in sorted(
        {d if isinstance(d, date) else date.fromisoformat(d) for d in days}
    ):
        if runs and day - runs[-1][-1] == timedelta(days=1):
            runs[-1].append(day)
        else:
            runs.append([day])
    if not runs:
        raise ValueError("a period has at least one day")
    return ", ".join(
        f"{run[0]}" if len(run) == 1 else f"{run[0]} to {run[-1]}" for run in runs
    )


def _check_period(period: str) -> str:
    get_period_days(period)
    return period


Period = Annotated[str, AfterValidator(_check_period)]
"""
A period of dates (see `get_period_days`), such as `2026-10-02`,
`2026-10-03T12:00Z`, `2026-08-15 to 2026-10-08`, or `2026-11-01, 2026-11-03`.
"""


def _get_days_apart(first: str, second: str) -> int:
    """
    Gets the days between the nearest days of two periods (0 if they share one).
    """
    return min(
        max(0, (b_start - a_end).days, (a_start - b_end).days)
        for a_start, a_end in get_period_days(first)
        for b_start, b_end in get_period_days(second)
    )


class Metric(_CatalogModel):
    """
    A value of a metric from a comparison with reference data.
    """

    name: str = Field(..., description="The metric (see `get_gates`).")
    value: float = Field(
        ..., allow_inf_nan=False, description="The value, which is finite."
    )
    scope: str | None = Field(
        default=None,
        description="The part or condition of the test it applies to (such as an "
        + "instrument or a time step), if not all of it. A gate takes the worst "
        + "value of its metric over scopes.",
    )


class Transfer(_CatalogModel):
    """
    The results of a configuration on a platform on which it was not tested,
    which do not count toward its validation.
    """

    period: Period = Field(..., description="The period of the reference data.")
    metrics: tuple[Metric, ...] = Field(..., min_length=1, description="The results.")


class ConfigurationRecord(_CatalogModel):
    """
    The results of validating a configuration, written by a validation script.
    """

    digest: str = Field(..., description="The `Configuration.digest` tested.")
    standard_version: int = Field(
        ..., description="The `STANDARD_VERSION` of the test."
    )
    tatc_version: str = Field(..., description="The version of TAT-C tested.")
    source: str = Field(..., description="The script of the results.")
    reference: str = Field(..., description="The reference data.")
    period: Period = Field(..., description="The period of the reference data.")
    platforms: tuple[str, ...] = Field(
        ..., min_length=1, description="The platforms of the reference data."
    )
    platform_periods: dict[str, Period] = Field(
        default_factory=dict,
        description="The period of each platform's reference data, if they differ "
        + "(as for the spacecraft of a system, whose values are fitted to each "
        + "spacecraft's own data), for every platform tested or for none.",
    )
    metrics: tuple[Metric, ...] = Field(..., description="The results.")
    transfers: dict[str, Transfer] = Field(
        default_factory=dict,
        description="Results on platforms on which the configuration was not "
        + "tested, by platform.",
    )
    fixture_metrics: dict[str, float] = Field(
        default_factory=dict,
        description="Results from the reduced test fixture, which detect "
        + "regressions but do not themselves establish agreement with the "
        + "reference data.",
    )

    @model_validator(mode="after")
    def _check_platform_periods(self) -> "ConfigurationRecord":
        if self.platform_periods and set(self.platform_periods) != set(self.platforms):
            raise ValueError(
                "platform_periods gives the period of every platform tested, or of "
                f"none: {sorted(self.platform_periods)} of {list(self.platforms)}"
            )
        return self


class ValidationRecord(_CatalogModel):
    """
    The results of validating the configurations of an entry.
    """

    entry_id: str = Field(..., description="The identifier of the entry.")
    configurations: dict[str, ConfigurationRecord] = Field(
        ..., description="The results of each configuration."
    )

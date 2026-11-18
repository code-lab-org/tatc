"""
The command line shared by the validation scripts, with defaults from the catalog:
the days of a run default to those of the entry's record (or, with `--fit`, to
those of the data to which its values were fitted), and the platforms to those of
the entry.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import pathlib
import sys
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from datetime import date, datetime, timedelta, timezone

from tatc.catalog import InstrumentEntry, get_entry
from tatc.catalog.standard import MIN_HOLDOUT_GAP_DAYS

from .archives import get_element_set
from .records import source_of


def parse_day(text: str) -> datetime:
    """
    Parses a day, as an ISO 8601 date (its start) or a time with a zone (such as
    `2026-10-03T12:00Z`), in UTC.
    """
    text = text.strip()
    try:
        day = date.fromisoformat(text)
    except ValueError:
        moment = datetime.fromisoformat(text.replace("Z", "+00:00"))
        if moment.tzinfo is None:
            raise ValueError(f"the time {text!r} needs a time zone") from None
        return moment.astimezone(timezone.utc)
    return datetime(day.year, day.month, day.day, tzinfo=timezone.utc)


def _format_day(day: datetime) -> str:
    """Formats a day, with its time if it has one."""
    return f"{day:%Y-%m-%d}" if day.time() == datetime.min.time() else f"{day:%Y-%m-%dT%H:%MZ}"


@dataclass
class Run:
    """
    What a run of a validation script tests: the entry's `configurations` on a
    `platform` (as a `transfer`, if they were not tested on it), on the data of
    `days`; or, with `fit`, the values it fits. Options that a script adds are
    its attributes too.
    """

    entry: InstrumentEntry
    directory: pathlib.Path
    days: list[datetime]
    platform: str
    transfer: bool
    configurations: tuple[str, ...]
    fit: bool
    write: bool
    elements: pathlib.Path | None
    source: str
    options: argparse.Namespace = field(default_factory=argparse.Namespace)

    def __getattr__(self, name):
        # the options a script adds
        options = self.__dict__.get("options")
        if options is not None and hasattr(options, name):
            return getattr(options, name)
        raise AttributeError(name)

    @property
    def day(self) -> datetime:
        """The day of the data (the first, of several)."""
        return self.days[0]

    @property
    def period(self) -> str:
        """The days of the data, as a record states them."""
        return ", ".join(f"{day:%Y-%m-%d}" for day in self.days)

    def element_set(
        self,
        day: datetime | None = None,
        *,
        within: timedelta | None = None,
        directory: pathlib.Path | None = None,
    ) -> list[str]:
        """
        Gets the element set of the run's platform on a day (by default, the run's),
        from the file given with `--elements`, if any (see
        `archives.get_element_set`).
        """
        return get_element_set(
            self.platform,
            self.day if day is None else day,
            file=self.elements,
            within=within,
            directory=self.directory if directory is None else directory,
        )


def _measured_periods(entry, keys, platform) -> set[str]:
    """Gets the periods of the data to which values were fitted on a platform."""
    periods = set()
    for key in keys:
        configuration = entry.configurations[key]
        provenances = [
            p for i in configuration.instruments for p in i.provenance.values()
        ]
        if configuration.tilt_schedule is not None:
            provenances += list(configuration.tilt_schedule.provenance.values())
        periods |= {
            p.period
            for p in provenances
            if p.kind == "measured" and platform in p.platforms
        }
    return periods


def _recorded_periods(entry, keys, platform, transfer) -> dict[str, str]:
    """Gets the periods of the data of the records (or transfers) of configurations."""
    periods = {}
    for key in keys:
        record = None if entry.record is None else entry.record.configurations.get(key)
        if record is not None and transfer:
            record = record.transfers.get(platform)
        if record is not None:
            periods[key] = record.period
    return periods


def _default_days(parser, entry, keys, platform, fit, transfer) -> list[datetime]:
    """
    Gets the days of a run from the catalog: those of the fit of the
    configurations, or of their records.
    """
    if fit:
        periods = _measured_periods(entry, keys, platform)
        what = f"the values fitted on {platform}"
    else:
        recorded = _recorded_periods(entry, keys, platform, transfer)
        periods = set(recorded.values())
        what = f"the records of {', '.join(keys)}" + (
            f" on {platform}" if transfer else ""
        )
        if len(periods) > 1:
            parser.error(
                f"{what} have different periods ({recorded}); give the days, or "
                "one --configuration"
            )
    if len(periods) != 1:
        parser.error(f"give the days: {what} have no single period ({periods or 'none'})")
    period = periods.pop()
    if " to " in period:
        parser.error(f"give the days: the period of {what} is a range ({period})")
    return [parse_day(item) for item in period.split(",")]


def parse_run(
    entry_id: str,
    description: str,
    *,
    configurations: Sequence[str] | None = None,
    default_platform: str | None = None,
    many_days: bool = False,
    add_arguments: Callable[[argparse.ArgumentParser], None] | None = None,
    arguments: Sequence[str] | None = None,
) -> Run:
    """
    Parses the command line of a validation script:

        <script> DIRECTORY [DAY ...] [--platform P] [--transfer]
            [--configuration KEY ...] [--fit] [--write] [--elements FILE]

    Args:
        entry_id (str): The identifier of the entry the script validates.
        description (str): The script's usage, its docstring.
        configurations (Sequence[str] | None): The configurations the script can
            validate, by default all of the entry's.
        default_platform (str | None): The platform, if not given, by default the
            first of the entry's default configuration.
        many_days (bool): Whether a run may pool the data of several days.
        add_arguments (Callable | None): Adds the script's own options to the
            parser.
        arguments (Sequence[str] | None): The arguments, by default those of the
            command line.

    Returns:
        Run: The run.
    """
    entry = get_entry(entry_id)
    available = tuple(entry.configurations if configurations is None else configurations)
    if default_platform is None:
        default_platform = entry.configurations[entry.default_configuration].platforms[0]
    platforms = sorted(platform.id for platform in entry.platforms)
    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "directory",
        type=pathlib.Path,
        help="where downloaded data are cached",
    )
    parser.add_argument(
        "days",
        nargs="*",
        metavar="DAY",
        help="the day of the data (UTC, as 2026-10-02), "
        + ("or several days to pool, " if many_days else "")
        + "by default that of the entry's record, or with --fit that of the fit",
    )
    parser.add_argument(
        "--platform",
        choices=platforms,
        default=default_platform,
        help=f"the platform whose data are compared (default: {default_platform})",
    )
    parser.add_argument(
        "--transfer",
        action="store_true",
        help=f"test the configurations of {default_platform} on --platform, as a "
        "transfer that does not count toward the tier (implied for a platform "
        "without configurations of its own)",
    )
    parser.add_argument(
        "--configuration",
        choices=available,
        action="append",
        help="a configuration to test (repeatable; default: all of the platform's)",
    )
    parser.add_argument(
        "--fit",
        action="store_true",
        help="measure the fitted values and print them, without writing",
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help="write the validation record and the reduced test fixture",
    )
    parser.add_argument(
        "--elements",
        type=pathlib.Path,
        help="a file of element sets, from which to take the platform's",
    )
    if add_arguments is not None:
        add_arguments(parser)
    options = parser.parse_args(arguments)

    platform = options.platform
    own = [k for k in available if platform in entry.configurations[k].platforms]
    transfer = options.transfer or not own
    if transfer and platform == default_platform:
        parser.error(f"a transfer is to a platform other than {default_platform}")
    if transfer and options.fit:
        parser.error("--fit measures a platform's own configurations, not a transfer")
    tested_on = default_platform if transfer else platform
    keys = [k for k in available if tested_on in entry.configurations[k].platforms]
    if options.configuration:
        unknown = sorted(set(options.configuration) - set(keys))
        if unknown:
            parser.error(f"{unknown} are not configurations tested on {tested_on}")
        keys = [k for k in keys if k in options.configuration]
    if options.days:
        try:
            days = [parse_day(text) for text in options.days]
        except ValueError as error:
            parser.error(str(error))
    else:
        days = _default_days(parser, entry, keys, platform, options.fit, transfer)
        print(f"days from the catalog: {', '.join(map(_format_day, days))}")
    if len(days) > 1 and not many_days:
        parser.error("give one day")
    if options.write and not options.fit and not transfer:
        period = ", ".join(f"{d:%Y-%m-%d}" for d in days)
        for key in keys:
            gap = entry.configurations[key].get_holdout_gap(period)
            if gap is not None and gap < MIN_HOLDOUT_GAP_DAYS:
                print(
                    f"warning: {key}'s fitted values are {gap} days from these data "
                    f"(fewer than {MIN_HOLDOUT_GAP_DAYS}), which its validation "
                    "reports as a gap"
                )
    options.directory.mkdir(parents=True, exist_ok=True)
    common = ("directory", "days", "platform", "transfer", "configuration", "fit")
    extra = argparse.Namespace(
        **{
            k: v
            for k, v in vars(options).items()
            if k not in (*common, "write", "elements")
        }
    )
    return Run(
        entry=entry,
        directory=options.directory,
        days=days,
        platform=platform,
        transfer=transfer,
        configurations=tuple(keys),
        fit=options.fit,
        write=options.write,
        elements=options.elements,
        source=source_of(sys.modules["__main__"].__file__),
        options=extra,
    )

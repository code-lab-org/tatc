"""
The validation records of the catalog's entries and the reduced test fixtures from
which unit tests recompute their results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import json
import pathlib

import numpy as np
import pandas as pd

import tatc
from tatc.catalog import STANDARD_VERSION, InstrumentEntry
from tatc.catalog.records import ConfigurationRecord, Transfer, ValidationRecord
from tatc.constants import timescale

# where the scripts write validation records and test fixtures
REPOSITORY = pathlib.Path(__file__).resolve().parents[3]
RECORDS = REPOSITORY / "src" / "tatc" / "resources" / "validation"
FIXTURES = REPOSITORY / "tests" / "catalog" / "fixtures"


def source_of(script):
    """
    The path of a script relative to the repository, as validation records
    cite their source.
    """
    return pathlib.Path(script).resolve().relative_to(REPOSITORY).as_posix()


def _utc_date(value) -> pd.Timestamp:
    """Returns the start of the UTC day of a time, read as UTC if naive."""
    time = pd.Timestamp(value)
    time = time.tz_localize("UTC") if time.tzinfo is None else time.tz_convert("UTC")
    return time.normalize()


def get_fixture_days(arrays: dict) -> list:
    """
    Returns the UTC dates a fixture's arrays record, from its `time` and `day`
    arrays and its `day_<mode>` arrays.
    """
    keys = [k for k in arrays if k in ("day", "time") or k.startswith("day_")]
    return sorted({_utc_date(str(arrays[k])) for k in keys})


def get_earth_orientation_samples(day) -> tuple:
    """
    Returns the delta T (TT - UT1, seconds) the installed skyfield gives from a
    day before to two days after a day, as TT Julian dates and values, sampled
    hourly and at the knots of its daily table, so that interpolating the samples
    linearly gives the values skyfield gives.
    """
    start = _utc_date(day) - pd.Timedelta(days=1)
    tt = timescale.from_datetimes(
        [(start + pd.Timedelta(hours=h)).to_pydatetime() for h in range(72 + 1)]
    ).tt
    table = getattr(timescale, "delta_t_table", None)
    if table is not None:
        knots = np.asarray(table[0])
        tt = np.union1d(tt, knots[(knots > tt[0]) & (knots < tt[-1])])
    return tt, np.asarray(timescale.delta_t_function(tt), dtype=float)


def add_earth_orientation(fixture: dict, written=None) -> dict:
    """
    Returns a fixture's arrays with the delta T (TT - UT1) samples of its days
    (`delta_t_tt`, `delta_t`), which the unit tests compute with instead of the
    installed skyfield's table, so that recorded values stay reproducible across
    skyfield releases. Samples are computed for the days of the arrays `written`
    (by default, all of them); samples of the fixture's other days are kept, as
    its records of those days were made with them.
    """
    written = fixture if written is None else written
    samples = [get_earth_orientation_samples(d) for d in get_fixture_days(written)]
    tt = np.asarray(fixture.get("delta_t_tt", []), dtype=float)
    delta_t = np.asarray(fixture.get("delta_t", []), dtype=float)
    keep = np.ones(len(tt), dtype=bool)
    for day_tt, _ in samples:
        keep &= (tt < day_tt[0]) | (tt > day_tt[-1])
    tt = np.concatenate([tt[keep]] + [day_tt for day_tt, _ in samples])
    delta_t = np.concatenate([delta_t[keep]] + [values for _, values in samples])
    tt, first = np.unique(tt, return_index=True)
    return {**fixture, "delta_t_tt": tt, "delta_t": delta_t[first]}


def update_fixture(path: pathlib.Path, **arrays) -> None:
    """
    Writes arrays to a fixture, keeping the others, with the delta T samples of
    the days written (see `add_earth_orientation`).
    """
    existing = {}
    if path.exists():
        with np.load(path, allow_pickle=False) as data:
            existing = {key: data[key] for key in data.files}
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **add_earth_orientation({**existing, **arrays}, arrays))
    print(f"wrote {path}")


def _to_metrics(metrics: list[dict] | dict[str, float]) -> list[dict]:
    """Gets metrics as records state them, from a mapping of names to values."""
    if isinstance(metrics, dict):
        return [{"name": name, "value": float(v)} for name, v in metrics.items()]
    return metrics


def _record_path(entry: InstrumentEntry) -> pathlib.Path:
    """Gets the path of an entry's record."""
    return RECORDS / f"{entry.id}.json"


def _read_record(entry: InstrumentEntry) -> dict:
    """Reads an entry's record, or starts one."""
    path = _record_path(entry)
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    return {"entry_id": entry.id, "configurations": {}}


def _write_record(entry: InstrumentEntry, record: dict) -> None:
    """Checks an entry's record against its schema and writes it."""
    path = _record_path(entry)
    checked = ValidationRecord(**record).model_dump(mode="json", exclude_none=True)
    path.write_text(json.dumps(checked, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {path}")


def write_record(
    entry: InstrumentEntry,
    key: str,
    *,
    source: str,
    reference: str,
    period: str,
    platforms: list[str],
    metrics: list[dict] | dict[str, float],
    fixture_metrics: dict[str, float] | None = None,
) -> None:
    """
    Writes the results of a configuration to its entry's record, keeping those of
    the others, with the configuration's digest and the versions of TAT-C and of
    the standard. Its transfers are kept if the configuration is unchanged.
    Metrics are records of `name`, `value`, and optionally `scope`, or a mapping
    of names to values.
    """
    record = _read_record(entry)
    digest = entry.configurations[key].digest
    previous = record["configurations"].get(key, {})
    record["configurations"][key] = ConfigurationRecord(
        digest=digest,
        standard_version=STANDARD_VERSION,
        tatc_version=tatc.__version__,
        source=source,
        reference=reference,
        period=period,
        platforms=platforms,
        metrics=_to_metrics(metrics),
        transfers=(
            previous.get("transfers", {}) if previous.get("digest") == digest else {}
        ),
        fixture_metrics=fixture_metrics or {},
    ).model_dump(mode="json", exclude_none=True)
    _write_record(entry, record)


def write_transfer(
    entry: InstrumentEntry,
    key: str,
    platform: str,
    *,
    period: str,
    metrics: list[dict] | dict[str, float],
) -> None:
    """
    Writes the results of a configuration on a platform on which it was not
    tested, which do not count toward its validation, to its record (with
    metrics as for `write_record`).

    Raises:
        ValueError: If the configuration has no record, or changed after it was
            tested, as its next record would discard the transfer.
    """
    record = _read_record(entry)
    tested = record["configurations"].get(key)
    if tested is None or tested["digest"] != entry.configurations[key].digest:
        raise ValueError(
            f"{entry.id}:{key} has no record of the configuration as cataloged; "
            "write its record before its transfers"
        )
    tested.setdefault("transfers", {})[platform] = Transfer(
        period=period, metrics=_to_metrics(metrics)
    ).model_dump(mode="json", exclude_none=True)
    _write_record(entry, record)


def write_results(
    run,
    key: str,
    *,
    reference: str,
    metrics: list[dict] | dict[str, float],
    fixture: dict | None = None,
    fixture_path: pathlib.Path | None = None,
    fixture_metrics: dict[str, float] | None = None,
    transfer_metrics: list[dict] | dict[str, float] | None = None,
    period: str | None = None,
) -> None:
    """
    Writes the results of a configuration in a run (see `cli.parse_run`). On a
    platform on which it was not tested, its `transfer_metrics` (by default, its
    metrics) are written to its record as a transfer. Otherwise its record is
    written, for the platforms on which it was tested, after the arrays of a
    `fixture`, if any, are written to the test fixture at `fixture_path`.

    Args:
        run (cli.Run): The run.
        key (str): The identifier of the configuration.
        reference (str): The reference data.
        metrics (list[dict] | dict[str, float]): The results (see `write_record`).
        fixture (dict | None): Arrays to write to the test fixture, if any.
        fixture_path (pathlib.Path | None): The test fixture, by default that of
            the entry (`tests/catalog/fixtures/<entry>.npz`).
        fixture_metrics (dict[str, float] | None): The results on the fixture.
        transfer_metrics (list[dict] | dict[str, float] | None): The results of a
            transfer, if not the metrics.
        period (str | None): The period of the reference data, by default the
            run's days.
    """
    period = run.period if period is None else period
    if run.transfer:
        write_transfer(
            run.entry,
            key,
            run.platform,
            period=period,
            metrics=metrics if transfer_metrics is None else transfer_metrics,
        )
        return
    if fixture is not None:
        update_fixture(
            FIXTURES / f"{run.entry.id}.npz" if fixture_path is None else fixture_path,
            **fixture,
        )
    write_record(
        run.entry,
        key,
        source=run.source,
        reference=reference,
        period=period,
        platforms=list(run.entry.configurations[key].platforms),
        metrics=metrics,
        fixture_metrics=fixture_metrics,
    )

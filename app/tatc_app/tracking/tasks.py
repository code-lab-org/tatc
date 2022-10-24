from datetime import datetime
from shapely.geometry import shape
import time
import json
from geojson_pydantic import FeatureCollection

from tatc.analysis.track import collect_orbit_track, collect_ground_track
from tatc import schemas

from ..worker import app


@app.task
def collect_orbit_track_task(
    satellite: str, instrument: str, times: list, elevation: float, mask: str
) -> str:
    """
    Task to collect orbit track.

    Args:
        satellites (str): JSON serialized :class:`tatc.schemas.Satellite` object.
        instrument (str): JSON serialized :class:`tatc.schemas.Instrument` object.
        times (str): List of ISO 8601 serialized times for which to compute orbit track.
        elevation (float): Elevation (meters) above datum in the WGS 84 coordinate system.
        mask (str): Optional GeoJSON serialized mask to constrain points.

    Returns:
        str: GeoJSON serialized orbit track.
    """
    results = collect_orbit_track(
        schemas.satellite.Satellite.parse_raw(satellite),
        schemas.instrument.Instrument.parse_raw(instrument),
        [datetime.fromisoformat(time) for time in times],
        elevation,
        shape(mask) if mask is not None else None,
    )
    # serialize Timestamp
    results["time"] = results["time"].apply(lambda t: t.isoformat())
    return results.to_json(show_bbox=False, drop_id=True)


@app.task
def collect_ground_track_task(
    satellite: str, instrument: str, times: list, elevation: float, mask: str, crs: str
) -> str:
    """
    Task to collect ground track.

    Args:
        satellites (str): JSON serialized :class:`tatc.schemas.Satellite` object.
        instrument (str): JSON serialized :class:`tatc.schemas.Instrument` object.
        times (list): List of ISO 8601 serialized times for which to compute orbit track.
        elevation (float): Elevation (meters) above datum in the WGS 84 coordinate system.
        mask (str): Optional GeoJSON serialized mask to constrain points.
        crs (str): Coordinate reference system (CRS) in which to project ground track.

    Returns:
        str: GeoJSON serialized ground track.
    """
    results = collect_ground_track(
        schemas.satellite.Satellite.parse_raw(satellite),
        schemas.instrument.Instrument.parse_raw(instrument),
        [datetime.fromisoformat(time) for time in times],
        elevation,
        shape(mask) if mask is not None else None,
        crs,
    )
    # serialize Timestamp
    results["time"] = results["time"].apply(lambda t: t.isoformat())
    return results.to_json(show_bbox=False, drop_id=True)

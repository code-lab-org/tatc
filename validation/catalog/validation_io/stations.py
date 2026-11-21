"""
The stations of the catalog's networks compared with independent sources: the
NEXRAD Level II archive on AWS (the radars' own metadata and the volumes they
scanned), NCEI's station tables, and the stations that report data on days.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import io
import logging
import pathlib
import re
import struct
from datetime import date

import numpy as np
import pandas as pd
from pyproj import Geod

from .archives import get_with_retries

LEVEL_2 = "https://unidata-nexrad-level2.s3.amazonaws.com/"
NCEI_NEXRAD = "https://www.ncei.noaa.gov/access/homr/file/nexrad-stations.txt"
GEOD = Geod(ellps="WGS84")
FEET = 0.3048


def list_level_2_stations(day: date) -> list[str]:
    """Lists the stations with Level II volumes on a day (UTC)."""
    listing = get_with_retries(
        LEVEL_2,
        params={"list-type": "2", "prefix": f"{day:%Y/%m/%d}/", "delimiter": "/"},
    ).text
    return re.findall(rf"<Prefix>{day:%Y/%m/%d}/([A-Z0-9]+)/</Prefix>", listing)


def list_level_2_volumes(day: date, station: str) -> list[str]:
    """Lists the keys of a station's Level II volumes on a day (UTC)."""
    keys, token = [], None
    while True:
        params = {"list-type": "2", "prefix": f"{day:%Y/%m/%d}/{station}/"}
        if token:
            params["continuation-token"] = token
        listing = get_with_retries(LEVEL_2, params=params).text
        keys += re.findall(r"<Key>(.*?)</Key>", listing)
        match = re.search(
            r"<NextContinuationToken>(.*?)</NextContinuationToken>", listing
        )
        if match is None:
            return [key for key in keys if key.endswith("_V06")]
        token = match.group(1)


def fetch_level_2(key: str, directory: pathlib.Path, records: int | None = 3) -> bytes:
    """
    Fetches a Level II volume (cached in a directory): its first records (each a
    compressed block of radials, after a 24-byte header), or all of it.
    """
    path = directory / (key.replace("/", "_") + (f".{records}" if records else ""))
    if path.exists():
        return path.read_bytes()
    if records is None:
        data = get_with_retries(LEVEL_2 + key).content
    else:
        data, chunk = b"", 2**20
        while True:
            response = get_with_retries(
                LEVEL_2 + key,
                headers={"Range": f"bytes={len(data)}-{len(data) + chunk - 1}"},
            )
            data += response.content
            end, count = 24, 0
            while count < records and end + 4 <= len(data):
                size = abs(struct.unpack(">i", data[end : end + 4])[0])
                if end + 4 + size > len(data):
                    break
                end, count = end + 4 + size, count + 1
            if count == records or len(response.content) < chunk:
                data = data[:end]
                break
    directory.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)
    return data


def _wrap(angle: float) -> float:
    """Wraps an elevation angle (degrees) to [-180, 180)."""
    return angle - 360 if angle >= 180 else angle


def describe_level_2(data: bytes) -> dict:
    """
    Describes a Level II volume: the radar's position, ground (m above mean sea
    level), and feedhorn height (m), its volume coverage pattern (VCP) and the
    commanded elevation angles of its cuts, and, for a whole volume, the lowest
    and highest scanned elevation angles and the reflectivity range (m) of its
    lowest sweep.
    """
    # pylint: disable-next=import-outside-toplevel
    from metpy.io import Level2File

    logging.getLogger("metpy.io.nexrad").setLevel(logging.ERROR)
    volume = Level2File(io.BytesIO(data))
    if not volume.sweeps or not volume.sweeps[0] or volume.vcp_info is None:
        raise ValueError("the volume's records hold no radials")
    site = volume.sweeps[0][0][1]
    description = {
        "latitude": float(site.lat),
        "longitude": float(site.lon),
        "ground": float(site.site_amsl),
        "height": float(site.feedhorn_agl),
        "vcp": int(volume.vcp_info.num),
        "commanded": [round(_wrap(cut.el_angle), 2) for cut in volume.vcp_info.els],
    }
    scanned = [
        float(np.median([_wrap(radial[0].el_angle) for radial in sweep]))
        for sweep in volume.sweeps
        if sweep
    ]
    if scanned:
        description["lowest"] = min(scanned)
        description["highest"] = max(scanned)
        for sweep in volume.sweeps:
            moments = sweep[0][4]
            if b"REF" in moments:
                header = moments[b"REF"][0]
                description["range"] = float(
                    (header.first_gate + header.gate_width * header.num_gates) * 1e3
                )
                break
    return description


def read_ncei_nexrad(directory: pathlib.Path) -> pd.DataFrame:
    """
    Reads NCEI's table of NEXRAD stations (cached in a directory): its WSR-88D
    radars (not the Terminal Doppler Weather Radars it also lists), indexed by
    ICAO identifier, with the antenna's elevation in meters above mean sea level
    (`antenna`).
    """
    path = directory / "nexrad-stations.txt"
    if not path.exists():
        directory.mkdir(parents=True, exist_ok=True)
        path.write_bytes(get_with_retries(NCEI_NEXRAD).content)
    lines = [
        line for line in path.read_text(encoding="utf-8").splitlines() if line.strip()
    ]
    spans = [(m.start(), m.end()) for m in re.finditer(r"-+", lines[1])]
    table = pd.read_fwf(io.StringIO("\n".join(lines)), colspecs=spans, skiprows=[1])
    # the WSR-88D radars, without the FAA's Terminal Doppler Weather Radars
    table = table[table["STNTYPE"] == "NEXRAD"].copy()
    table["antenna"] = table["ELEV"] * FEET
    return table.set_index("ICAO")


def horizontal_errors(table: pd.DataFrame, reference: pd.DataFrame) -> np.ndarray:
    """
    Gets the geodesic distances (m) between the positions of a table's stations
    and those of a reference (both with `latitude` and `longitude`), by index.
    """
    common = table.index.intersection(reference.index)
    _, _, distance = GEOD.inv(
        table.loc[common, "longitude"].to_numpy(),
        table.loc[common, "latitude"].to_numpy(),
        reference.loc[common, "longitude"].to_numpy(),
        reference.loc[common, "latitude"].to_numpy(),
    )
    return np.asarray(distance)


def membership(table_ids: set[str], reporting: set[str]) -> tuple[float, float]:
    """
    Compares a table's stations with the stations reporting data: the recall
    (percent of reporting stations in the table) and precision (percent of the
    table's stations reporting).
    """
    both = len(table_ids & reporting)
    return 100 * both / len(reporting), 100 * both / len(table_ids)

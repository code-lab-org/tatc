"""
Access to reference data: HTTP requests with retries, files read with range
requests, the public archives of the JPSS satellites, a NASA Earthdata login, the
EUMETSAT Data Store, and the element sets of the platforms.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import dataclasses
import functools
import io
import pathlib
import re
import time as clock
from datetime import datetime, timedelta

import h5py
import requests
import yaml

from tatc.catalog import get_entry_ids, get_entry

ELEMENT_SET_FOLDER = pathlib.Path(__file__).resolve().parents[1] / "elements"
ELEMENT_SET_REGISTRY = ELEMENT_SET_FOLDER / "element_sets.yml"

_SESSION = requests.Session()


def login_earthdata():
    """
    Logs in to NASA Earthdata with a token in `EARTHDATA_TOKEN` (or a username and
    password in `EARTHDATA_USERNAME` and `EARTHDATA_PASSWORD`), or else with the
    credentials of urs.earthdata.nasa.gov in `~/.netrc`.
    """
    import earthaccess  # pylint: disable=import-outside-toplevel

    for strategy in ("environment", "netrc"):
        try:
            auth = earthaccess.login(strategy=strategy)
        except Exception:  # pylint: disable=broad-exception-caught
            continue
        if auth is not None and auth.authenticated:
            return auth
    raise RuntimeError(
        "Earthdata login failed: set EARTHDATA_TOKEN, or add urs.earthdata.nasa.gov "
        "to ~/.netrc"
    )


def download_eumetsat(
    collection: str,
    start: datetime,
    end: datetime,
    directory: pathlib.Path,
    *,
    match: str = "",
    entry: str | None = None,
) -> list[pathlib.Path]:
    """
    Downloads the products of a collection of the EUMETSAT Data Store sensed in a
    period (whose names contain `match`), or one file of each (whose name ends with
    `entry`), unless they are cached in a directory. Downloads need a login to the
    Data Store's version 2 APIs (`eumdac auth --interactive`, whose tokens `eumdac`
    refreshes).

    Args:
        collection (str): The identifier of the collection.
        start (datetime.datetime): The start of the period.
        end (datetime.datetime): The end of the period.
        directory (pathlib.Path): Where the files are cached.
        match (str): Text that the names of the products contain.
        entry (str | None): The end of the name of the file of each product to
            download, or `None` for the whole product.

    Returns:
        list[pathlib.Path]: The files, in the order of the products' names.
    """
    # imported here, as only the scripts of EUMETSAT's products need it
    import eumdac  # pylint: disable=import-outside-toplevel
    from eumdac.token import OAuth2Token  # pylint: disable=import-outside-toplevel

    directory.mkdir(parents=True, exist_ok=True)
    store = eumdac.DataStore(OAuth2Token())
    products = sorted(
        (
            product
            for product in store.get_collection(collection).search(
                dtstart=start.replace(tzinfo=None), dtend=end.replace(tzinfo=None)
            )
            if match in str(product)
        ),
        key=str,
    )
    paths = []
    for product in products:
        name = (
            str(product) + ".zip"
            if entry is None
            else next(e for e in product.entries if e.endswith(entry))
        )
        path = directory / pathlib.PurePosixPath(name).name
        if not path.exists():
            for attempt in range(4):
                try:
                    with product.open(entry=None if entry is None else name) as source:
                        temporary = path.with_suffix(path.suffix + ".part")
                        with open(temporary, "wb") as target:
                            while chunk := source.read(1 << 20):
                                target.write(chunk)
                    temporary.replace(path)
                    break
                except Exception:  # pylint: disable=broad-exception-caught
                    if attempt == 3:
                        raise
                    clock.sleep(10 * (attempt + 1))
        paths.append(path)
    return paths


def get_with_retries(url: str, session=None, **kwargs) -> requests.Response:
    """Sends a GET request, retrying transient failures with a backoff."""
    for attempt in range(6):
        try:
            response = (session or _SESSION).get(url, timeout=180, **kwargs)
            response.raise_for_status()
            return response
        except requests.RequestException:
            if attempt == 5:
                raise
            clock.sleep(2**attempt)
    raise RuntimeError("unreachable")


def post_with_retries(url: str, body: dict) -> requests.Response:
    """Posts a JSON body, retrying transient failures with a backoff."""
    for attempt in range(6):
        try:
            response = _SESSION.post(url, json=body, timeout=300)
            response.raise_for_status()
            return response
        except requests.RequestException:
            if attempt == 5:
                raise
            clock.sleep(2**attempt)
    raise RuntimeError("unreachable")


class RemoteFile(io.RawIOBase):
    """A read-only file of a URL, downloaded on demand with range requests."""

    def __init__(self, url: str, session=None, block_size: int = 2**16):
        super().__init__()
        self.session, self.block_size, self.blocks, self.position = (
            session,
            block_size,
            {},
            0,
        )
        response = get_with_retries(url, session, headers={"Range": "bytes=0-0"})
        # later requests go directly to where any redirect led
        self.url = response.url
        self.size = int(response.headers["Content-Range"].split("/")[-1])

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.position

    def seek(self, offset, whence=io.SEEK_SET):
        self.position = {
            io.SEEK_SET: 0,
            io.SEEK_CUR: self.position,
            io.SEEK_END: self.size,
        }[whence] + offset
        return self.position

    def fetch(self, start: int, end: int) -> None:
        """Downloads block-aligned bytes into the cache."""
        response = get_with_retries(
            self.url, self.session, headers={"Range": f"bytes={start}-{end - 1}"}
        )
        for i, block in enumerate(
            range(start // self.block_size, (end - 1) // self.block_size + 1)
        ):
            self.blocks[block] = response.content[
                i * self.block_size : (i + 1) * self.block_size
            ]

    def readinto(self, buffer):
        count = min(len(buffer), self.size - self.position)
        if count <= 0:
            return 0
        first = self.position // self.block_size
        last = (self.position + count - 1) // self.block_size
        missing = [b for b in range(first, last + 1) if b not in self.blocks]
        if missing:
            self.fetch(
                missing[0] * self.block_size,
                min((missing[-1] + 1) * self.block_size, self.size),
            )
        data = b"".join(self.blocks[b] for b in range(first, last + 1))
        start = self.position - first * self.block_size
        buffer[:count] = data[start : start + count]
        self.position += count
        return count


@dataclasses.dataclass(frozen=True)
class Archive:
    """
    The archive of a JPSS satellite on NOAA Open Data Dissemination, from which
    objects and element sets are read, with the prefix that distinguishes its
    cached files.
    """

    url: str
    prefix: str

    def list_keys(self, prefix: str) -> list[str]:
        """Lists the keys of all objects in the archive with a prefix."""
        keys, token = [], None
        while True:
            params = {"list-type": "2", "prefix": prefix}
            if token:
                params["continuation-token"] = token
            listing = get_with_retries(self.url, params=params).text
            keys += re.findall(r"<Key>(.*?)</Key>", listing)
            match = re.search(
                r"<NextContinuationToken>(.*?)</NextContinuationToken>", listing
            )
            if match is None:
                return keys
            token = match.group(1)

    def fetch(self, key: str) -> bytes:
        """Downloads an object of the archive."""
        return get_with_retries(self.url + key).content

    def load_element_set(
        self, day: datetime, directory: pathlib.Path
    ) -> list[str] | None:
        """
        Gets the operational element set with an epoch at the start of a day, cached
        in a directory, or `None` if the archive has none.
        """
        path = directory / f"{self.prefix}_tle_{day:%Y%m%d}.txt"
        if not path.exists():
            key = next(
                (
                    k
                    for k in self.list_keys(f"TLE-AUX/{day:%Y/%m/%d}/")
                    if f"Z_{day:%Y%m%d}000000Z_" in k
                ),
                None,
            )
            if key is None:
                return None
            with h5py.File(io.BytesIO(self.fetch(key))) as f:
                text = bytes(f["All_Data/TLE-AUX_All/Dataset_Array_0"][:]).decode()
            path.write_text(text.strip() + "\n", encoding="utf-8")
        return path.read_text(encoding="utf-8").splitlines()


# the archives of the JPSS satellites, by platform; NOAA's archive has Suomi NPP's
# operational element sets only from 2026-08-19, and earlier days read them from
# the files of element sets
ARCHIVES = {
    "noaa-20": Archive("https://noaa-nesdis-n20-pds.s3.amazonaws.com/", "n20"),
    "noaa-21": Archive("https://noaa-nesdis-n21-pds.s3.amazonaws.com/", "n21"),
    "suomi-npp": Archive("https://noaa-nesdis-snpp-pds.s3.amazonaws.com/", "npp"),
}


def _parse_element_sets(path: pathlib.Path, norad_id: int) -> list[tuple]:
    """Gets the epochs and two-line element sets of a satellite in a file."""
    lines = [
        line.rstrip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith(("1 ", "2 "))
    ]
    sets = []
    for first, second in zip(lines[::2], lines[1::2], strict=True):
        if int(first[2:7]) != norad_id:
            continue
        year, days = int(first[18:20]), float(first[20:32])
        epoch = datetime(2000 + year, 1, 1) + timedelta(days=days - 1)
        sets.append((epoch, [first, second]))
    return sets


def read_element_set_file(
    path: pathlib.Path,
    norad_id: int,
    day: datetime,
    within: timedelta | None = None,
) -> list[str]:
    """
    Gets the two-line element set of a satellite whose epoch is nearest the middle of a
    day, from a file of element sets with or without name lines; or, `within` a time
    of the middle of the day, every such set of a distinct epoch, in order of epoch.
    """
    if not path.exists():
        raise FileNotFoundError(f"element sets for {day:%Y-%m-%d} need {path}")
    sets = _parse_element_sets(path, norad_id)
    if not sets:
        raise ValueError(f"{path} has no element sets of {norad_id}")
    return _select_element_sets(sets, day, within, path)


def _select_element_sets(
    sets: list[tuple], day: datetime, within: timedelta | None, source
) -> list[str]:
    """
    Gets the element set whose epoch is nearest the middle of a day, or every set
    of a distinct epoch `within` a time of it, in order of epoch.
    """
    middle = day.replace(tzinfo=None) + timedelta(hours=12)
    if within is not None:
        # the last set of each epoch
        near = {epoch: tle for epoch, tle in sets if abs(epoch - middle) <= within}
        if not near:
            raise ValueError(
                f"{source} has no element sets within {within} of {middle}"
            )
        return [line for epoch in sorted(near) for line in near[epoch]]
    return min(sets, key=lambda item: abs((item[0] - middle).total_seconds()))[1]


@functools.cache
def _get_norad_ids() -> dict[str, int]:
    """Gets the NORAD catalog number of each platform of the catalog."""
    return {
        platform.id: platform.norad_id
        for entry_id in get_entry_ids()
        for platform in get_entry(entry_id).platforms
    }


@functools.cache
def _read_registry() -> dict[str, dict[str, list[str]]]:
    """Reads the element sets the validations used, by platform and day."""
    with open(ELEMENT_SET_REGISTRY, "r", encoding="utf-8") as f:
        registry = yaml.safe_load(f)
    return {
        platform: {str(day): list(lines) for day, lines in days.items()}
        for platform, days in registry.items()
    }


def get_element_set(
    platform: str,
    day: datetime,
    *,
    file: pathlib.Path | None = None,
    within: timedelta | None = None,
    directory: pathlib.Path | None = None,
) -> list[str]:
    """
    Gets the element set of a platform for a day: from a `file` of element sets, if
    given; or else the one registered for the day in `elements/element_sets.yml`;
    or else NOAA's operational element set of a JPSS satellite, cached in a
    `directory`; or else from the files in `elements/` of the platform's element
    sets. From files, it gets the set whose epoch is nearest the middle of the day,
    or every set `within` a time of it.

    Args:
        platform (str): The identifier of the platform, as in the catalog.
        day (datetime.datetime): The day (UTC).
        file (pathlib.Path | None): A file of element sets, if any.
        within (datetime.timedelta | None): The time from the middle of the day
            within which every set of a file is used, or `None` for the nearest.
        directory (pathlib.Path | None): Where to cache element sets downloaded
            from an archive.

    Returns:
        list[str]: The lines of the element sets.

    Raises:
        FileNotFoundError: If no source has an element set of the platform.
    """
    norad_id = _get_norad_ids()[platform]
    if file is not None:
        return read_element_set_file(file, norad_id, day, within)
    registered = _read_registry().get(platform, {}).get(f"{day:%Y-%m-%d}")
    if registered is not None:
        return list(registered)
    if platform in ARCHIVES and directory is not None:
        tle = ARCHIVES[platform].load_element_set(day, directory)
        if tle is not None:
            return tle
    # the sets of every file that holds the platform's, as several files may
    sets = [
        item
        for path in sorted(ELEMENT_SET_FOLDER.glob("*.3le"))
        for item in _parse_element_sets(path, norad_id)
    ]
    if sets:
        return _select_element_sets(sets, day, within, ELEMENT_SET_FOLDER)
    raise FileNotFoundError(
        f"no element set of {platform} for {day:%Y-%m-%d}: register one in "
        f"{ELEMENT_SET_REGISTRY.name}, or give a file of element sets (--elements)"
    )

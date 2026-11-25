"""
Helpers shared by the scripts that validate catalog entries, by topic:

* `cli`: the shared command line, with the days and platforms of a run from the
  catalog;
* `archives`: reading reference data (HTTP, range requests, the JPSS archives,
  NASA Earthdata, the EUMETSAT Data Store) and the platforms' element sets;
* `geometry`: coordinates, and the lattices of samples that scanners image;
* `scoring`: comparisons through public TAT-C methods;
* `records`: the validation records and reduced test fixtures.
* `footprints`: the footprints and spacing of beam instruments' samples.
* `grids`: the pixels of geostationary imagers' fixed grids.

Scripts import these names from the package.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from .archives import (
    ARCHIVES,
    ELEMENT_SET_FOLDER,
    Archive,
    RemoteFile,
    download_eumetsat,
    get_element_set,
    get_with_retries,
    login_earthdata,
    post_with_retries,
    read_element_set_file,
)
from .cli import Run, parse_day, parse_run
from .footprints import Footprints, measure_along_spacing, validate_footprints
from .grids import FixedGrid, validate_fixed_grid
from .geometry import (
    EARTH_ROTATION,
    TO_ECEF,
    TO_GEODETIC,
    ScanLattice,
    ecef,
    fit_cone,
    iet_to_utc,
    in_lattice_cells,
    inertial_velocity,
    normalize,
    transform_geometry,
    up,
)
from .records import (
    FIXTURES,
    RECORDS,
    REPOSITORY,
    add_earth_orientation,
    source_of,
    update_fixture,
    write_record,
    write_results,
    write_transfer,
)
from .scoring import (
    GRID_STEP,
    PUBLISHED_EPOCH,
    SIZE_TOLERANCE,
    build_grid,
    compute_published_errors,
    compute_spacing_ratios,
    match_passes,
    pick_passes,
    predict_passes,
    recompute_passes,
    score_passes,
    score_published,
    score_swath_area,
    to_observations,
)

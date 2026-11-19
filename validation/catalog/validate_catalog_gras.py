"""
Validates the GRAS catalog entry (the radio occultation receivers of Metop-B and
Metop-C) against UCAR's processing of a day of Metop occultations, in the AWS Registry
of Open Data (read with `awsgnssroutils`), using the exact values of the entry through
public TAT-C methods (see `occultation.py`). The database lists every occultation
recorded on the day, with its receiver, transmitter, and direction, and each
occultation's calibrated phase (Level 1A) file records the receiver's and
transmitter's Earth-fixed positions.

* Receivers and transmitters are identified by the satellite whose element sets
  reproduce their recorded positions, in a sample of files that holds every receiver
  and transmitter of the day, keeping each satellite's element sets within 20 km of
  them.
* Each recorded occultation's time at the sample elevation is computed from these
  orbits, and the occultation is predicted if the entry's occultations (with every
  identified GPS satellite, over the day) include one of the same receiver,
  transmitter, and direction within 60 s (`occultation_recall`). Predicted
  occultations that were not recorded are not errors, as receivers track transmitters
  selectively; their share is reported but not gated.
* In a sample of occultations, the straight line between the recorded positions gives
  the time and tangent point at the sample elevation, which are compared with the
  predicted occultation's time and position (`epoch_p95`, `tangent_position_median`).

The sample elevation and maximum yaw come from `--fit`, on a day of data: the median
elevation of the straight line when the refracted ray's tangent point is at the
reference point of each of a sample of UCAR's refractivity profiles (`PROFILE_POINT`,
see `occultation.py`), rounded to the nearest kilometer, and the largest yaw, from the
receiver's velocity or anti-velocity, of the recorded occultations at that elevation. Running on a day at least 30 days from the fit tests it on data held out
from the fit. The days are those on which the database holds both receivers'
occultations.

Reads anonymous public data, and needs the `validation` extra and `awsgnssroutils`:

    python validation/catalog/validate_catalog_gras.py <data directory> 2026-05-13 --fit
    python validation/catalog/validate_catalog_gras.py <data directory> 2026-06-12

Without a day, it takes that of the entry's record (or, with `--fit`, of the fit). The
element sets of the receivers and transmitters of each day are in the files of
`ELEMENT_FILES`. With `--write`, it writes the entry's validation record and the
reduced fixture from which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import numpy as np

from occultation import validate_occultations
from validation_io import ELEMENT_SET_FOLDER, REPOSITORY, parse_run

NOTEBOOK_ELEMENTS = REPOSITORY / "docs" / "validation" / "limb_and_space" / "elements"
# the files of element sets of each day's receivers and transmitters
ELEMENT_FILES = {
    "2026-05-13": [
        NOTEBOOK_ELEMENTS / "gps_2026-05-14.3le",
        ELEMENT_SET_FOLDER / "metop_2026-05-12.3le",
    ],
    "2026-06-12": [
        ELEMENT_SET_FOLDER / "gps_2026-06-11.3le",
        ELEMENT_SET_FOLDER / "metop_2026-06-11.3le",
    ],
}

# the point of UCAR's refractivity profiles at which the sample elevation is fitted:
# their reference point, which geolocates them
PROFILE_POINT = "reference"


def main() -> None:
    """Validates the entry on a day of data, or fits its maximum yaw."""
    validate_occultations(
        parse_run("gras", __doc__),
        mission="metop",
        receiver_prefix="metop",
        element_files=ELEMENT_FILES,
        reference="UCAR Metop occultations (AWS GNSS RO database, v1.1) and "
        "calibrated phase files (receiver and transmitter positions), and "
        "refractivity profiles",
        profile_point=PROFILE_POINT,
    )


if __name__ == "__main__":
    np.seterr(invalid="ignore")
    main()

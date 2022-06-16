# Tradespace Analysis Toolkit for Constellations (TAT-C)

The Tradespace Analysis Toolkit for Constellations (TAT-C) provides low-level
data structures and functions for systems engineering analysis and design of
Earth-observing space missions.

## Installation

TAT-C uses conda to manage dependencies.

Configure conda package:
```shell
conda-build tatc
conda install --use-local tatc
```

Now, the library `tatc` can be imported from Python scripts running in the
environment. To verify, try:
```shell
python
import tatc
```

## Unit Tests

Run unit tests:
```shell
python -m unittest
```

## Acknowledgements

This project was supported in part by the National Aeronautics and Space
Administration (NASA) Earth Science Division (ESD) Earth Science Technology
Office (ESTO) Advanced Information Systems Technology (AIST) program under
grant numbers: NNX17AE06G, 80NSSC17K0586, 80NSSC20K1118, and 80NSSC21K1515.

### Current Project Team
 * PI: Paul T. Grogan <pgrogan@stevens.edu>
 * I. Josue Tapia-Tamayo <itapiata@stevens.edu>

### Project Alumni
 * Isaac Feldman
 * Hayden Daly
 * Lindsay Portelli
 * Matthew Sabatini
 * Evan Abel
 * Sigfried Hache

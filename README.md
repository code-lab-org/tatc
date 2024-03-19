# Tradespace Analysis Toolkit for Constellations (TAT-C)

The Tradespace Analysis Toolkit for Constellations (TAT-C) provides low-level
data structures and functions for systems engineering analysis and design of
Earth-observing space missions.

Documentation: https://tatc.readthedocs.io

## Installation

TAT-C uses conda and the conda-forge channel for distribution because some of
the underlying libraries are platform dependent.

The simplest way to use TAT-C is to install it with the command:
```shell
conda install tatc -c conda-forge
```
Then, TAT-C is available for use in any Python script by importing:
```python
import tatc
```

Alternatively, to run TAT-C from a local editable source, (e.g., for modifying
TAT-C functionality), clone this repository and create a new conda environment:
```shell
conda env create -f environment.yml
```
(Note that if you need the `dev` or `examples` options below, use `conda env create -f environment-dev.yml` instead).

Then, activate the tatc_env environment:
```shell
conda activate tatc_env
```
And finally install the tatc library in "editable" mode:
```shell
pip install -e .
```

Note: the following optional dependencies are available with bracket notation: 
 * `pip install -e .[dev]`: for development (unit testing, coverage, and linting)
 * `pip install -e .[docs]`: for generating documentation
 * `pip install -e .[examples]`: for running optional examples
 * `pip install -e .[app]`: for running the web application
Multiple optional dependencies can be installed with a comma-separated list (e.g., `pip install -e .[dev,examples]`)

### Faster Installation

For faster dependency solving during installation, consider installing the
mamba package:
```shell
conda install mamba -c conda-forge
```
and replace `conda` with `mamba` in the installation instructions above.

## Development Tools

### Unit Tests

Run unit tests with:
```shell
python -m unittest
```

Optionally, run a test coverage report:
```shell
coverage run -m unittest
```
including html output:
```shell
coverage html
```

### Documentation

Generate documentation from the `docs` directory using the command:
```shell
make html
```

### Code Style

This project uses the black code style, applied from the project root:
```shell
black .
```

## Contact

Paul T. Grogan <paul.grogan@asu.edu>

## Acknowledgements

This project was supported in part by the National Aeronautics and Space
Administration (NASA) Earth Science Division (ESD) Earth Science Technology
Office (ESTO) Advanced Information Systems Technology (AIST) program under
grant numbers: NNX17AE06G, 80NSSC17K0586, 80NSSC20K1118, and 80NSSC21K1515.

Current Project Team
 * PI: Paul T. Grogan <paul.grogan@asu.edu>
 * I. Josue Tapia-Tamayo <josue.tapia@asu.edu>

Project Alumni
 * Isaac Feldman
 * Hayden Daly
 * Lindsay Portelli
 * Matthew Sabatini
 * Evan Abel
 * Sigfried Hache

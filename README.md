# Tradespace Analysis Toolkit for Constellations (TAT-C)

The Tradespace Analysis Toolkit for Constellations (TAT-C) provides low-level
data structures and functions for systems engineering analysis and design of
Earth-observing space missions.

Documentation: https://tatc.readthedocs.io

## Installation

TAT-C uses the pip build system to manage dependencies. Install the tatc library in "editable" mode:
```shell
pip install -e .
```

Note: the following optional dependencies are available with bracket notation: 
 * `pip install -e .[dev]`: for development (unit testing, coverage, and linting)
 * `pip install -e .[docs]`: for generating documentation
 * `pip install -e .[examples]`: for running optional examples
 * `pip install -e .[app]`: for running the web application
Multiple optional dependencies can be installed with a comma-separated list (e.g., `pip install -e .[dev,examples]`)

## Development Tools

### Docker and Docker-Compose

TAT-C includes a Dockerfile (`Dockerfile`) that specifies build targets for a base runtime (the TAT-C library), an application server, and a distributed task worker.

To build images, run:
```shell
docker build -t tatc .
```

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

# Tradespace Analysis Toolkit for Constellations (TAT-C)

[![PyPI](https://img.shields.io/pypi/v/tatc.svg)](https://pypi.org/project/tatc/)
[![Python Versions](https://img.shields.io/pypi/pyversions/tatc.svg)](https://pypi.org/project/tatc/)
[![Unit Tests](https://github.com/code-lab-org/tatc/actions/workflows/unit-test.yml/badge.svg)](https://github.com/code-lab-org/tatc/actions/workflows/unit-test.yml)
[![Documentation](https://readthedocs.org/projects/tatc/badge/?version=latest)](https://tatc.readthedocs.io)
[![License](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17363628.svg)](https://doi.org/10.5281/zenodo.17363628)

The Tradespace Analysis Toolkit for Constellations (TAT-C) provides low-level
data structures and functions for systems engineering analysis and design of
Earth-observing space missions suitable for pre-Phase A concept studies.

Documentation: [https://tatc.readthedocs.io](https://tatc.readthedocs.io)

Repository: [https://github.com/code-lab-org/tatc](https://github.com/code-lab-org/tatc)

## Installation

TAT-C requires Python 3.10&ndash;3.14. Install the latest release from PyPI:
```shell
pip install tatc
```

### Development Installation

To work with the source code, clone the repository and install the tatc library in "editable" mode:
```shell
pip install -e .
```

Note: the following optional dependencies are available with bracket notation: 
 * `pip install -e ".[dev]"`: for development (unit testing, coverage, and linting)
 * `pip install -e ".[docs]"`: for generating documentation
 * `pip install -e ".[examples]"`: for running optional examples
 * `pip install -e ".[osse]"`: for running optional observing system simulation experiment (OSSE) examples
 * `pip install -e ".[preprocess]"`: for preprocessing functions (e.g., terrain masks)
 * `pip install -e ".[validation]"`: for running validation notebooks (includes `examples` and `preprocess`)

Multiple optional dependencies can be installed with a comma-separated list (e.g., `pip install -e ".[dev,examples]"`)

## Development Tools

Development tools are applicable when working with the source code.

### Unit Tests

Run unit tests with:
```shell
python -m unittest
```

Unit tests run numba-compiled functions (see below) compiled, which coverage cannot trace. To include them in a test coverage report, run the tests with the environment variable `NUMBA_DISABLE_JIT=1` (e.g., in bash):
```shell
NUMBA_DISABLE_JIT=1 coverage run -m unittest
```
including html output:
```shell
coverage html
```

### Compiled Functions

Numerical functions in `tatc.utils` are compiled with numba (`njit(cache=True)`), which caches the compiled machine code in `__pycache__` directories next to the source files (or in the directory set by the `NUMBA_CACHE_DIR` environment variable). Numba recompiles a function only if its own source file changes, so after editing values it uses from other modules (e.g., in `tatc/constants.py`) or other compiled functions it calls (e.g., in `tatc/utils/orbital.py`), clear the cache from the project root with:
```shell
python -c "import pathlib; [p.unlink() for p in pathlib.Path('src').rglob('*.nb[ci]')]"
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

Pull requests are also linted with pylint, which must score at least 9.0:
```shell
pylint --rcfile=.pylintrc src/tatc
```

## License

This project is licensed under the BSD 3-Clause License &mdash; see
[LICENSE](LICENSE) for details.

## Citation

If you use TAT-C in your research, please cite it using the metadata in
[CITATION.cff](CITATION.cff), or the DOI badge above for the latest release.

## Contact

Paul T. Grogan <paul.grogan@asu.edu>

## Acknowledgements

This project was supported in part by the National Aeronautics and Space
Administration (NASA) Earth Science Division (ESD) Earth Science Technology
Office (ESTO) Advanced Information Systems Technology (AIST) program. 
Financial support is acknowledged under NASA grant numbers: NNX17AE06G, 
80NSSC17K0586, 80NSSC20K1118, 80NSSC21K1515, 80NSSC22K1705, 80NSSC24K0575, 
80NSSC24K0921; NASA Jet Propulsion Laboratory subcontracts: 1689594, 1686623, 
1704657, 1705655; Texas A \& M University subaward M2403907.

Current Project Team
 * PI: Paul T. Grogan <paul.grogan@asu.edu>

Project Alumni
 * I. Josue Tapia-Tamayo
 * Suvan Kumar
 * Isaac Feldman
 * Hayden Daly
 * Lindsay Portelli
 * Matthew Sabatini
 * Evan Abel
 * Sigfried Hache
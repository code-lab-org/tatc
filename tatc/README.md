# Tradespace Analysis Tool for Constellations (TAT-C) Core Library

This is the core library that supports TAT-C. It provides low-level data
structures and analysis functions.

## Installation

Some of TAT-C's dependencies are only available via the conda package manager
and, specifically, using the conda-forge channel. Therefore, it is easiest to
install dependencies in a new environment created specifically for TAT-C use.

Configure conda environment:
```shell
conda env create -f environment.yml
conda activate tatc_env
```

Within the TAT-C environment, you can install `tatc` library in editable
mode (`-e`) from current directory (`.`):
```shell
python -m pip install -e .
```

Now, the library `tatc` can be import from Python scripts running in the
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

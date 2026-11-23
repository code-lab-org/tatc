.. role:: console(code)
  :language: console

.. role:: python(code)
  :language: python

============
Installation
============

PyPI Installation
------------------

TAT-C uses pip and the PyPI software repository to install the main package along with underlying dependencies and libraries.

TAT-C can be installed from PyPI via the terminal/shell command :console:`pip install tatc`.

Then, TAT-C is available for use in any Python script by importing :python:`import tatc`.

Conda Installation
------------------

TAT-C uses conda and the conda-forge channel for distribution because some of the underlying libraries are platform dependent.

The simplest way to install TAT-C is via the terminal/shell command :console:`conda install conda-forge::tatc`.

Then, TAT-C is available for use in any Python script by importing :python:`import tatc`.

Source Installation
-------------------
Alternatively, TAT-C can be installed from a local editable source, (e.g., for modifying TAT-C functionality). 
This method is required to access run the web-based application.

Clone the project repository using :console:`git clone https://github.com/code-lab-org/tatc.git`.

Install the tatc library in "editable" mode :console:`pip install -e .`

Note: the following optional dependencies are available with bracket notation: 
 * :console:`pip install -e ".[dev]"`: for development functions (unit testing, coverage, and linting)
 * :console:`pip install -e ".[docs]"`: for generating documentation in :console:`docs/`
 * :console:`pip install -e ".[examples]"`: for running optional examples in :console:`docs/examples`
 * :console:`pip install -e ".[osse]"`: for running optional observing system simulation experiment (OSSE) examples
 * :console:`pip install -e ".[preprocess]"`: for preprocessing functions (e.g., terrain masks)
 * :console:`pip install -e ".[airborne]"`: for planning flights and modeling sensors with `HyPlan <https://github.com/ryanpavlick/hyplan>`_, whose outputs the airborne preprocessing functions read (these functions do not require HyPlan)
 * :console:`pip install -e ".[tradespace]"`: for parallel trade studies and searches with optimizers (joblib, Optuna, and pymoo)
 * :console:`pip install -e ".[validation]"`: for running validation notebooks in :console:`docs/validation` and catalog validation scripts in :console:`validation/catalog` (includes :console:`examples` and :console:`preprocess`)

Multiple optional dependencies can be installed with a comma-separated list (e.g., :console:`pip install -e ".[dev,examples]"`)

Faster Package Resolution
-------------------------
Dependency resolution in conda can be very slow.
For faster installation, consider installing the mamba package :console:`conda install mamba -c conda-forge` and replace :console:`conda` with :console:`mamba` in the installation instructions above.

Unit Tests
----------

Run unit tests from the project root directory with :console:`python -m pytest -n auto`, which runs them in parallel on every available core with pytest-xdist (both installed by the :console:`dev` optional dependencies), as continuous integration does.
The tests are written with :console:`unittest`, so :console:`python -m unittest` also runs them (one at a time) without pytest.
For a test coverage report, run :console:`coverage run -m pytest`.

Documentation
-------------

Build the documentation from the :console:`docs/` directory with :console:`make html`.

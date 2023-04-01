.. role:: console(code)
  :language: console

.. role:: python(code)
  :language: python

============
Installation
============

Conda Installation
------------------

TAT-C uses conda and the conda-forge channel for distribution because some of the underlying libraries are platform dependent.

The simplest way to install TAT-C is via the terminal/shell command :console:`conda install tatc -c conda-forge`.

Then, TAT-C is available for use in any Python script by importing :python:`import tatc`.

Source Installation
-------------------
Alternatively, TAT-C can be installed from a local editable source, (e.g., for modifying TAT-C functionality). 
This method is required to access run the web-based application.

Clone the project repository and, from the project root directory, create a new conda environment with the required dependencies :console:`conda env create -f environment.yml`. 
(Note: if needing the :console:`dev`, :console:`docs`, or :console:`examples` options below, use :console:`conda env create -f environment-dev.yml` instead).

Then, activate the tatc_env environment :console:`conda activate tatc_env`.

Finally, install the tatc library within the environment in "editable" mode :console:`pip install -e .`

Note: the following optional dependencies are available with bracket notation: 
 * :console:`pip install -e .[dev]`: for development functions (unit testing, coverage, and linting)
 * :console:`pip install -e .[docs]`: for generating documentation in :console:`docs/`
 * :console:`pip install -e .[examples]`: for running optional examples in :console:`docs/examples`
 * :console:`pip install -e .[app]`: for running the web application in :console:`app/`

Multiple optional dependencies can be installed with a comma-separated list (e.g., :console:`pip install -e .[dev,examples]`)

Faster Package Resolution
-------------------------
Dependency resolution in conda can be very slow.
For faster installation, consider installing the mamba package :console:`conda install mamba -c conda-forge` and replace :console:`conda` with :console:`mamba` in the installation instructions above.

Unit Tests
----------

Run unit tests from the project root directory with :console:`python -m unittest` or :console:`coverage run -m unittest`.

Documentation
-------------

Build the documentation from the :console:`docs/` directory with :console:`make html`.

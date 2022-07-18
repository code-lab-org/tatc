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

The simplest way to install TAT-C is via the terminal/shell command :console:`conda install -c conda-forge tatc`.

Then, TAT-C is available for use in any Python script by importing :python:`import tatc`.

Source Installation
-------------------
Alternatively, TAT-C can be installed from a local editable source, (e.g., for modifying TAT-C functionality).

Clone the project repository and create a new conda environment with the required dependencies :console:`conda env create -f environment.yml`.

Then, activate the tatc_env environment :console:`conda activate tatc_env`.

Finally, install the tatc library in "editable" mode :console:`pip install -e .`

Faster Package Resolution
-------------------------
Dependency resolution in conda can be very slow.
For faster installation, consider installing the mamba package :console:`conda install mamba -c conda-forge` and replace :console:`conda` with :console:`mamba` in the installation instructions above.

Additional Dependencies
-----------------------

Included examples (`docs/examples`) require additional dependencies which can be installed via :console:`conda install -c conda-forge geoplot contextily`.

Recommended development dependencies can be installed with conda via :console:`conda install -c conda-forge coverage sphinx nbsphinx black` and with pip via :console:`pip install sphinx-rtd-theme`.

Unit Tests
----------

Run unit tests with :console:`python -m unittest` or :console:`coverage run -m unittest`.

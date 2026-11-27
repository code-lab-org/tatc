========
Examples
========

The following provide examples of TAT-C use in a Jupyter notebook environment, grouped by topic. Start with the getting started examples, which introduce the analysis functions with one satellite.

.. toctree::
  :maxdepth: 2

  getting_started/index
  coverage/index
  limb_and_space/index
  ground_systems/index
  airborne/index
  osse/index
  response/index
  tradespace/index

Note that running the above examples in a local Python environment requires additional dependencies which can be installed via::

  pip install tatc[examples]

using the pip package manager. The airborne examples require::

  pip install tatc[airborne]
  
the observing system simulation experiments require::

  pip install tatc[osse]

and the trade space tutorial requires::

  pip install tatc[tradespace]

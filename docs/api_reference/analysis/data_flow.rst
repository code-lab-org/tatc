Data Flow
=========

The flow of data from satellites to the ground through ground stations and
relays: the data each satellite generates and stores, the links scheduled
through antennas at rates from link budgets, and the latency (or loss) of each
observation.

Compute Data Flow
-----------------
.. autofunction:: tatc.analysis.compute_data_flow

Analyze Data Flow
-----------------
.. autofunction:: tatc.analysis.analyze_data_flow

.. autoclass:: tatc.analysis.DataFlowResult
  :members: find, satellites

Find Data Flow Downlinks
------------------------
.. autofunction:: tatc.analysis.find_data_flow_downlinks

Compute Generation
------------------
.. autofunction:: tatc.analysis.compute_generation

.. autofunction:: tatc.analysis.get_generation_rate

.. autofunction:: tatc.analysis.is_data_flow_modeled

Sample Link Geometry
--------------------
.. autofunction:: tatc.analysis.sample_link_geometry

Compute Link Rates
------------------
.. autofunction:: tatc.analysis.compute_link_rates

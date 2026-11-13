Instrument Catalog
==================

.. automodule:: tatc.catalog

The entries and their results are listed in :doc:`/validation/catalog`.

Entries
-------

.. autofunction:: tatc.catalog.get_entry

.. autofunction:: tatc.catalog.get_entry_ids

.. autofunction:: tatc.catalog.list_entries

.. autopydantic_model:: tatc.catalog.InstrumentEntry

.. autoclass:: tatc.catalog.InstrumentPayload

Attitudes by Day
----------------

.. autopydantic_model:: tatc.catalog.TiltSchedule

.. autofunction:: tatc.catalog.get_subsolar_latitude

Validation Status
-----------------

See :doc:`/validation/standard`.

.. autoclass:: tatc.catalog.EntryValidation

.. autoclass:: tatc.catalog.ConfigurationValidation

.. autoclass:: tatc.catalog.GateResult

.. autoclass:: tatc.catalog.Tier
   :members:

.. autoclass:: tatc.catalog.Grade
   :members:

.. autoclass:: tatc.catalog.InstrumentClass
   :members:

.. autofunction:: tatc.catalog.get_gates

.. autofunction:: tatc.catalog.get_required_grades

Metrics
-------

.. automodule:: tatc.catalog.metrics
   :members:

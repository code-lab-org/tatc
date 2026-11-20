Instrument Catalog
==================

.. automodule:: tatc.catalog

The entries and their results are listed in :doc:`/validation/catalog`.

Entries
-------

.. autofunction:: tatc.catalog.get_entry

.. autofunction:: tatc.catalog.get_entry_ids

.. autofunction:: tatc.catalog.list_entries

.. autofunction:: tatc.catalog.build_payload

.. autopydantic_model:: tatc.catalog.InstrumentEntry

.. autoclass:: tatc.catalog.InstrumentPayload
   :members:

Validation Status
-----------------

See :doc:`/validation/standard`.

.. autoclass:: tatc.catalog.EntryValidation
   :members:

.. autoclass:: tatc.catalog.ConfigurationValidation
   :members:

.. autoclass:: tatc.catalog.GateResult
   :members:

.. autoclass:: tatc.catalog.Tier
   :members:

.. autoclass:: tatc.catalog.Grade
   :members:

.. autoclass:: tatc.catalog.InstrumentClass
   :members:

Schemas
-------

The schemas of the files of entries.

.. autopydantic_model:: tatc.catalog.schemas.Configuration

.. autopydantic_model:: tatc.catalog.schemas.InstrumentSpecification

.. autopydantic_model:: tatc.catalog.schemas.Provenance

.. autopydantic_model:: tatc.catalog.schemas.TiltSchedule

.. autopydantic_model:: tatc.catalog.schemas.LimbScan

.. autopydantic_model:: tatc.catalog.schemas.ScanStart

.. autopydantic_model:: tatc.catalog.schemas.Occultation

Records
-------

The schemas of the validation records of entries.

.. autopydantic_model:: tatc.catalog.records.ConfigurationRecord

.. autopydantic_model:: tatc.catalog.records.Transfer

.. autofunction:: tatc.catalog.records.get_period_days

Standard
--------

.. automodule:: tatc.catalog.standard
   :members: STANDARD_VERSION, COMPATIBLE_VERSIONS, MIN_HOLDOUT_GAP_DAYS, Gate, get_gates, get_required_grades, is_compatible, serves, evaluate_entry

Metrics
-------

.. automodule:: tatc.catalog.metrics
   :members:

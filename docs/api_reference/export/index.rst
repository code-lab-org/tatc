================
Export Functions
================

Export functions write the observations of an observing system simulation
experiment (OSSE) to the files of other tools, such as data assimilation
systems. They require the optional `osse` dependencies
(``pip install tatc[osse]``).

CF-netCDF Point Files
---------------------

.. autofunction:: tatc.export.to_cf_netcdf

IODA Files of JEDI
------------------

.. autofunction:: tatc.export.to_ioda

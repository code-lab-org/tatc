"""
Functions to export the results of TAT-C to the files of other tools: the
observations of an observing system simulation experiment (OSSE) as
CF-netCDF point files (`to_cf_netcdf`, requires xarray) and as the IODA files
of the Joint Effort for Data assimilation Integration (JEDI, `to_ioda`,
requires netCDF4), both installed with `pip install tatc[osse]`. Each module
is imported when first used, so that it requires only its own dependencies.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from typing import Any

_MODULES = {
    "to_cf_netcdf": "cf",
    "to_ioda": "ioda",
}


def __getattr__(name: str) -> Any:
    """
    Gets a function from its module, imported when first used.
    """
    # pylint: disable=import-outside-toplevel
    module = _MODULES.get(name)
    if module == "cf":
        from . import cf as source
    elif module == "ioda":
        from . import ioda as source
    else:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    return getattr(source, name)


def __dir__() -> list[str]:
    return sorted(list(globals()) + list(_MODULES))


__all__ = sorted(_MODULES)

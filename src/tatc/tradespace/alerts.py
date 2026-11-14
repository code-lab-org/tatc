"""
Warnings of trade studies, which point at the call into the trade space
module that leads to them (rather than at the module's own code), however
deep within the module they are given.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import sys
import warnings

_PACKAGE = __name__.rpartition(".")[0]
"""The name of the package whose frames a warning skips."""


def _warn(message: str, category: type[Warning] = UserWarning) -> None:
    """
    Gives a warning at the first frame of the stack outside the package
    (the user's call into it), as `warnings.warn` would with the stack level
    of that frame.

    Args:
        message (str): The message.
        category (type[Warning]): The category of the warning.
    """
    level, frame = 2, sys._getframe(1)  # pylint: disable=protected-access
    while frame is not None and _is_package(frame.f_globals.get("__name__")):
        frame, level = frame.f_back, level + 1
    warnings.warn(message, category, stacklevel=level)


def _is_package(module: object) -> bool:
    """Whether a module name is of the package (or the package itself)."""
    return isinstance(module, str) and (
        module == _PACKAGE or module.startswith(_PACKAGE + ".")
    )

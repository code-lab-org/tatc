"""
Warnings that refer to the code that called TAT-C.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import os
import sys
import warnings

_PACKAGE = os.path.normcase(os.path.dirname(os.path.abspath(__file__))) + os.sep
"""The directory of the TAT-C package, whose frames a warning skips."""


def _get_external_stacklevel() -> int:
    """
    Gets the stack level, relative to the caller of this function, of the
    first frame outside the TAT-C package: the code that called TAT-C,
    however deeply the warning is issued within it.

    Returns:
        int: the stack level, for `warnings.warn` called by the caller
    """
    frame = sys._getframe(1)  # pylint: disable=protected-access
    level = 1
    while frame is not None and os.path.normcase(
        os.path.abspath(frame.f_code.co_filename)
    ).startswith(_PACKAGE):
        frame = frame.f_back
        level += 1
    return level


def warn(message: str, category: type[Warning] = UserWarning) -> None:
    """
    Issues a warning that refers to the code that called TAT-C (the first
    frame outside the package, see `_get_external_stacklevel`), rather than
    to a line within it, so that it shows where the user called the
    analysis, and warnings issued at different depths within TAT-C for the
    same call refer to the same line (and so are shown once by the default
    filter).

    Args:
        message (str): The message.
        category (type[Warning]): The category of warning.
    """
    warnings.warn(message, category, stacklevel=_get_external_stacklevel())

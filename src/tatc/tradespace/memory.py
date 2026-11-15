"""
Objects that keep results in memory (such as caches of the metrics of the
phases of missions, the targets of replicates, and the caches of their
observations), whose copies sent to the parallel workers of a study start
with empty memories.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations


class Forgetful:
    """
    A mixin of objects that keep results in memory, in the attributes named
    by `_forgotten` (each a container, such as a dictionary), which their
    copies (such as those pickled for the parallel workers of a study) start
    without: empty containers of the same types.
    """

    # pylint: disable=too-few-public-methods

    _forgotten: tuple[str, ...] = ()

    def __getstate__(self) -> dict:
        state = self.__dict__.copy()
        for name in self._forgotten:
            state[name] = type(state[name])()
        return state

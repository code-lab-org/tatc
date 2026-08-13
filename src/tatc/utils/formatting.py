"""
Formatting utility functions.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""


def zero_pad(object_name: str, max_number: int, current_number: int) -> str:
    """
    Appends a zero-padded number to an object name to create a unique,
    sortable label for a member of a numbered collection (e.g. a satellite
    in a constellation). The padding width is derived from `max_number` so
    every generated label consumes the same number of characters, which
    keeps lexicographic (string) sort order consistent with numeric order.

    Note:
        If `current_number` requires more digits than `max_number`, it is
        left at its natural width rather than truncated, so the
        fixed-width guarantee no longer holds for that entry.

    Args:
        object_name (str): Base name to prefix the generated label (e.g. a constellation name).
        max_number (int): The largest number that will be padded; determines the padding width.
        current_number (int): The number to zero-pad and append to `object_name`.

    Returns:
        str: `object_name`, a single space, and `current_number` zero-padded to the width of `max_number`.
    """
    max_length = len(str(max_number))
    return object_name + " " + str(current_number).zfill(max_length)

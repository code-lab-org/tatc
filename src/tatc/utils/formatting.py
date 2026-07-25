"""
Formatting utility functions.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

def zero_pad(object_name: str, max_number: int, current_number: int) -> str:
    """
    Uses length of max_number to zero pad allowing for alphanumeric sorting.

    Args:
        object_name (str): Object name, to be concatenated with zero padded number.
        max_number (int): Maximum number, utilized as reference for zero padding.
        current_number (int): Index number to be zero padded.

    Returns:
        str: The object name with zero padded number appended.
    """
    max_length = len(str(max_number))
    return object_name + " " + str(current_number).zfill(max_length)
"""
Checks that the dependencies of an installed TAT-C are at their minimum
versions.

The unit-test workflow installs TAT-C with the constraints in
`tests/minimum_versions.txt` to test the minimum versions of its dependencies.
This script checks that each dependency with a minimum version (`>=`) in
`pyproject.toml` was installed at that version, so that a minimum missing from
the constraints, or out of step with them, fails rather than going untested.
Run it after installing:

    pip install --prefer-binary -c tests/minimum_versions.txt -e ".[preprocess]"
    python tests/check_minimum_versions.py

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import re
import sys
from importlib.metadata import requires, version

# a required (not optional) dependency: its name, any extras, and specifiers
REQUIREMENT = re.compile(r"^([A-Za-z0-9._-]+)(?:\[[^\]]*\])?\s*([^;]*)$")
MINIMUM = re.compile(r">=\s*(\d+(?:\.\d+)*)")


def get_release(text: str) -> tuple[int, ...]:
    """
    Gets the release numbers of a version, without trailing zeros, so that
    versions such as "1.0" and "1.0.0" are equal.

    Args:
        text (str): The version.

    Returns:
        tuple[int, ...]: The release numbers.
    """
    numbers = [int(number) for number in re.match(r"\d+(?:\.\d+)*", text)[0].split(".")]
    while len(numbers) > 1 and numbers[-1] == 0:
        numbers.pop()
    return tuple(numbers)


def find_mismatches() -> list[str]:
    """
    Finds the dependencies of TAT-C not installed at their minimum versions.

    Returns:
        list[str]: A description of each mismatch.
    """
    mismatches = []
    for requirement in requires("tatc") or []:
        match = REQUIREMENT.match(requirement.strip())
        if match is None:
            continue  # an optional dependency, with a marker
        name, specifiers = match.groups()
        minimum = MINIMUM.search(specifiers)
        if minimum is None:
            continue
        installed = version(name)
        if get_release(installed) != get_release(minimum[1]):
            mismatches.append(f"{name}: minimum {minimum[1]}, installed {installed}")
    return mismatches


def main() -> int:
    """
    Checks the installed dependencies, reporting any not at their minimums.

    Returns:
        int: The exit code (0: all at their minimum versions).
    """
    mismatches = find_mismatches()
    for mismatch in mismatches:
        print(mismatch, file=sys.stderr)
    if mismatches:
        print(
            "pin each minimum version of pyproject.toml in tests/minimum_versions.txt",
            file=sys.stderr,
        )
        return 1
    print("all dependencies at their minimum versions")
    return 0


if __name__ == "__main__":
    sys.exit(main())

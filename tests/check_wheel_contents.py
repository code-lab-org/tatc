"""
Checks that a built wheel contains every packaged resource.

Unit tests run from an editable install, which finds resources in the source
tree, so they cannot detect a file left out of the package data in
`pyproject.toml`. This script compares a wheel with the source tree. It is not
named as a unit test because it needs a built wheel; run it after building:

    python -m build --wheel
    python tests/check_wheel_contents.py dist/*.whl

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import pathlib
import sys
import zipfile

RESOURCES = pathlib.Path(__file__).parents[1] / "src" / "tatc" / "resources"


def find_missing_resources(wheel_path: str) -> list[str]:
    """
    Finds the resource files in the source tree that are not in a wheel.

    Args:
        wheel_path (str): The path to the wheel.

    Returns:
        list[str]: The missing paths, relative to the installed package root.
    """
    expected = {
        pathlib.PurePath("tatc", "resources", path.relative_to(RESOURCES)).as_posix()
        for path in RESOURCES.rglob("*")
        if path.is_file()
        and "__pycache__" not in path.parts
        and path.suffix not in (".py", ".pyc")
    }
    with zipfile.ZipFile(wheel_path) as wheel:
        return sorted(expected.difference(wheel.namelist()))


def main(arguments: list[str]) -> int:
    """
    Checks one wheel, reporting any missing resources.

    Args:
        arguments (list[str]): The command line arguments: the wheel's path.

    Returns:
        int: The exit code (0: all resources found).
    """
    if len(arguments) != 1:
        print("usage: check_wheel_contents.py <wheel>", file=sys.stderr)
        return 2
    missing = find_missing_resources(arguments[0])
    for path in missing:
        print(f"missing from wheel: {path}", file=sys.stderr)
    if missing:
        print("add these files to package data in pyproject.toml", file=sys.stderr)
        return 1
    print(f"all packaged resources found in {arguments[0]}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

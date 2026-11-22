"""
Loads the scripts that write the catalog's documentation pages
(`docs/validation/catalog_page.py` and `systems_page.py`) by their paths, without
putting `docs/validation` on the import path, whose folders (such as `coverage`)
would shadow installed packages of the same names.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import importlib.util
import pathlib
import sys
from types import ModuleType

PAGES = pathlib.Path(__file__).parents[2] / "docs" / "validation"


def load_page_script(name: str) -> ModuleType:
    """
    Loads a page script (and the scripts it imports, `catalog_page` for
    `systems_page`), registered by name.

    Args:
        name (str): The name of the script (such as `systems_page`).

    Returns:
        types.ModuleType: The script.
    """
    if name in sys.modules:
        return sys.modules[name]
    if name == "systems_page":
        load_page_script("catalog_page")
    spec = importlib.util.spec_from_file_location(name, PAGES / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module

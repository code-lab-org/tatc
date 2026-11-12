"""
Unit tests for the structure of the trade space module: the general modules
of trade studies do not depend on those of constellations and missions, and
the modules import each other without cycles.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import ast
import unittest
from pathlib import Path

import tatc.tradespace

PACKAGE = Path(tatc.tradespace.__file__).parent

GENERAL = {
    "attribution",
    "confirmation",
    "criteria",
    "fingerprint",
    "nodes",
    "optimization",
    "pareto",
    "plotting",
    "refinement",
    "scenarios",
    "screening",
    "sensitivity",
    "space",
    "store",
    "study",
    "utility",
    "values",
}
"""The modules of trade studies of any designs."""


def _imports(path: Path) -> set[str]:
    """
    Gets the modules of the package that a module imports when it runs (at
    the top or within functions, not only for type checking).
    """
    found: set[str] = set()

    def visit(node: ast.AST) -> None:
        if (
            isinstance(node, ast.If)
            and getattr(node.test, "id", None) == "TYPE_CHECKING"
        ):
            return
        if isinstance(node, ast.ImportFrom) and node.level == 1 and node.module:
            found.add(node.module.split(".")[0])
        elif isinstance(node, ast.ImportFrom) and node.level == 1:
            found.update(alias.name for alias in node.names)
        for child in ast.iter_child_nodes(node):
            visit(child)

    visit(ast.parse(path.read_text(encoding="utf-8")))
    return found


class TestLayering(unittest.TestCase):
    """
    Unit tests for the dependencies between the modules of trade studies.
    """

    def setUp(self):
        self.modules = {
            path.stem: _imports(path)
            for path in PACKAGE.glob("*.py")
            if path.stem != "__init__"
        }

    def test_general_modules(self):
        """General modules do not import those of constellations and missions."""
        self.assertLessEqual(GENERAL, set(self.modules))
        for module in sorted(GENERAL):
            with self.subTest(module=module):
                self.assertLessEqual(self.modules[module], GENERAL)

    def test_no_cycles(self):
        """Modules import each other without cycles."""
        state: dict[str, str] = {}

        def visit(module: str, path: list[str]) -> None:
            if state.get(module) == "done":
                return
            if state.get(module) == "visiting":
                self.fail(f"import cycle: {' -> '.join(path + [module])}")
            state[module] = "visiting"
            for imported in sorted(self.modules.get(module, ())):
                visit(imported, path + [module])
            state[module] = "done"

        for module in sorted(self.modules):
            visit(module, [])

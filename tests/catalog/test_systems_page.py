"""
Tests of the systems' documentation page, which `docs/conf.py` writes from the
systems and records with `docs/validation/systems_page.py`. CI does not build the
documentation, so these check the page here.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import importlib.util
import pathlib
import unittest

from tatc.catalog import get_system, get_system_ids

SCRIPT = pathlib.Path(__file__).parents[2] / "docs" / "validation" / "systems_page.py"


def _load_script():
    spec = importlib.util.spec_from_file_location("systems_page", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestSystemsPage(unittest.TestCase):
    """
    Test the systems page.
    """

    @classmethod
    def setUpClass(cls):
        cls.script = _load_script()
        cls.page = cls.script.render()

    def test_every_system_listed_with_its_tier(self):
        """
        Test that each system has a section that states its computed tier.
        """
        for system_id in get_system_ids():
            system = get_system(system_id)
            self.assertIn(f".. _system-{system_id}:", self.page)
            self.assertIn(
                f'``get_system("{system_id}")``. **Class:** '
                f"``{system.system_class.value}``. **Tier:** "
                f"`{system.validation.tier.value}`.",
                self.page,
            )

    def test_every_failed_gate_shown(self):
        """
        Test that the page shows as many failed gates as the systems have.
        """
        failed = sum(
            result.passed is False
            for system_id in get_system_ids()
            for configuration in get_system(
                system_id
            ).validation.configurations.values()
            for result in configuration.gates
        )
        self.assertGreater(failed, 0)
        self.assertEqual(self.page.count("     - **fail**"), failed)

    def test_limitations_shown(self):
        """
        Test that every limitation of every configuration is on the page.
        """
        for system_id in get_system_ids():
            system = get_system(system_id)
            texts = [*system.limitations]
            for configuration in system.configurations.values():
                texts += configuration.limitations
            for text in texts:
                self.assertIn(self.script.escape(text), self.page, system_id)


if __name__ == "__main__":
    unittest.main()

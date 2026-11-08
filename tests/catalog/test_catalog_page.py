"""
Tests of the catalog's documentation page, which `docs/conf.py` writes from the
entries and records with `docs/validation/catalog_page.py`. CI does not build the
documentation, so these check the page here.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import importlib.util
import pathlib
import unittest

from tatc.catalog import Gate, get_entry, get_entry_ids

SCRIPT = pathlib.Path(__file__).parents[2] / "docs" / "validation" / "catalog_page.py"


def _load_script():
    spec = importlib.util.spec_from_file_location("catalog_page", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestCatalogPage(unittest.TestCase):
    """
    Test the catalog page.
    """

    @classmethod
    def setUpClass(cls):
        cls.script = _load_script()
        cls.page = cls.script.render()

    def test_every_entry_listed_with_its_tier(self):
        """
        Test that each entry has a section that states its computed tier.
        """
        for entry_id in get_entry_ids():
            entry = get_entry(entry_id)
            self.assertIn(f".. _catalog-{entry_id}:", self.page)
            self.assertIn(
                f'``get_entry("{entry_id}")``. **Class:** '
                f"``{entry.instrument_class.value}``. **Tier:** "
                f"`{entry.validation.tier.value}`.",
                self.page,
            )

    def test_every_failed_gate_shown(self):
        """
        Test that the page shows as many failed gates as the entries have.
        """
        failed = sum(
            result.passed is False
            for entry_id in get_entry_ids()
            for configuration in get_entry(entry_id).validation.configurations.values()
            for result in configuration.gates
        )
        self.assertGreater(failed, 0)
        self.assertEqual(self.page.count("     - **fail**"), failed)

    def test_limitations_shown(self):
        """
        Test that every limitation of every configuration is on the page.
        """
        for entry_id in get_entry_ids():
            entry = get_entry(entry_id)
            texts = [*entry.limitations]
            for configuration in entry.configurations.values():
                texts += configuration.limitations
            for text in texts:
                self.assertIn(self.script.escape(text), self.page, entry_id)

    def test_rounding_keeps_outcome(self):
        """
        Test that a rounded value never appears to meet a gate it fails.
        """
        gate = Gate("footprint_iou_median", 0.95, False, "")
        self.assertEqual(self.script.format_value(0.9486, gate), "0.9486")
        self.assertEqual(self.script.format_value(0.949996, gate), "0.949996")
        gate = Gate("pass_recall", 99, False, "%")
        self.assertEqual(self.script.format_value(98.9996, gate), "98.9996%")
        gate = Gate("epoch_p95", 1, True, "s")
        self.assertEqual(self.script.format_value(1.001, gate), "1.001 s")
        self.assertEqual(self.script.format_value(0.43, gate), "0.43 s")

    def test_markup_escaped(self):
        """
        Test that markup in plain text is shown as written.
        """
        self.assertEqual(
            self.script.escape("a *b* `c` d_ e_f"), r"a \*b\* \`c\` d\_ e_f"
        )


if __name__ == "__main__":
    unittest.main()

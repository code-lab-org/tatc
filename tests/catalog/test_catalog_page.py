"""
Tests of the catalog's documentation page, which `docs/conf.py` writes from the
entries and records with `docs/validation/catalog_page.py`. CI does not build the
documentation, so these check the page here.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import pathlib
import sys
import unittest

from tatc.catalog import get_entry, get_entry_ids

sys.path.insert(0, str(pathlib.Path(__file__).parents[2] / "docs" / "validation"))
import catalog_page  # pylint: disable=wrong-import-position,import-error


class TestCatalogPage(unittest.TestCase):
    """
    Test the catalog page.
    """

    @classmethod
    def setUpClass(cls):
        cls.script = catalog_page
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

    def test_markup_escaped(self):
        """
        Test that markup in plain text is shown as written.
        """
        self.assertEqual(
            self.script.escape("a *b* `c` d_ e_f"), r"a \*b\* \`c\` d\_ e_f"
        )


if __name__ == "__main__":
    unittest.main()

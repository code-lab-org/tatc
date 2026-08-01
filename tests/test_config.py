"""
Unit tests for the tatc.config module.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""
import shutil
import tempfile
import textwrap
import unittest
from pathlib import Path

import yaml
from pydantic import ValidationError

from tatc import config as tatc_config
from tatc.config import ConfigError, RuntimeConfiguration, load_yaml_config


class TestRuntimeConfiguration(unittest.TestCase):
    """
    Unit tests for the tatc.config.RuntimeConfiguration schema.
    """

    def test_defaults(self):
        """
        Test that constructing with no arguments produces the hard-coded
        default value documented for every field.
        """
        rc = RuntimeConfiguration()
        self.assertEqual(rc.footprint_points_elliptical, 32)
        self.assertEqual(rc.footprint_points_rectangular_side, 8)
        self.assertEqual(rc.repeat_cycle_delta_position_m, 10000)
        self.assertEqual(rc.repeat_cycle_delta_velocity_m_per_s, 10)
        self.assertEqual(rc.repeat_cycle_search_elevation_deg, 88)
        self.assertEqual(rc.repeat_cycle_search_duration_days, 30)
        self.assertTrue(rc.repeat_cycle_lazy_load)
        self.assertTrue(rc.repeat_cycle_for_orbit_track)
        self.assertTrue(rc.repeat_cycle_for_observation_events)
        self.assertTrue(rc.gp_orbit_lazy_load)

    def test_footprint_points_elliptical_accepts_minimum(self):
        """
        Test that exactly 4 elliptical footprint points is accepted (the
        boundary of the ge=4 constraint).
        """
        rc = RuntimeConfiguration(footprint_points_elliptical=4)
        self.assertEqual(rc.footprint_points_elliptical, 4)

    def test_footprint_points_elliptical_rejects_below_minimum(self):
        """
        Test that fewer than 4 elliptical footprint points is rejected,
        since a closed elliptical footprint polygon needs at least 4
        vertices.
        """
        with self.assertRaises(ValidationError):
            RuntimeConfiguration(footprint_points_elliptical=3)

    def test_footprint_points_rectangular_side_accepts_minimum(self):
        """
        Test that exactly 1 rectangular footprint side point is accepted
        (the boundary of the ge=1 constraint).
        """
        rc = RuntimeConfiguration(footprint_points_rectangular_side=1)
        self.assertEqual(rc.footprint_points_rectangular_side, 1)

    def test_footprint_points_rectangular_side_rejects_below_minimum(self):
        """
        Test that fewer than 1 rectangular footprint side point is
        rejected.
        """
        with self.assertRaises(ValidationError):
            RuntimeConfiguration(footprint_points_rectangular_side=0)

    def test_repeat_cycle_delta_position_m_rejects_non_positive(self):
        """
        Test that a zero or negative position delta tolerance is rejected,
        since a non-positive distance tolerance cannot discriminate
        repeated positions.
        """
        with self.assertRaises(ValidationError):
            RuntimeConfiguration(repeat_cycle_delta_position_m=0)
        with self.assertRaises(ValidationError):
            RuntimeConfiguration(repeat_cycle_delta_position_m=-1)

    def test_repeat_cycle_delta_velocity_m_per_s_rejects_non_positive(self):
        """
        Test that a zero or negative velocity delta tolerance is rejected.
        """
        with self.assertRaises(ValidationError):
            RuntimeConfiguration(repeat_cycle_delta_velocity_m_per_s=0)
        with self.assertRaises(ValidationError):
            RuntimeConfiguration(repeat_cycle_delta_velocity_m_per_s=-1)

    def test_repeat_cycle_search_elevation_deg_rejects_non_positive(self):
        """
        Test that a zero or negative minimum search elevation angle is
        rejected.
        """
        with self.assertRaises(ValidationError):
            RuntimeConfiguration(repeat_cycle_search_elevation_deg=0)
        with self.assertRaises(ValidationError):
            RuntimeConfiguration(repeat_cycle_search_elevation_deg=-1)


class TestLoadYamlConfig(unittest.TestCase):
    """
    Unit tests for the tatc.config.load_yaml_config function.
    """

    def setUp(self):
        self.tmp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.tmp_dir, ignore_errors=True)

    def _write(self, name: str, content: str) -> Path:
        path = self.tmp_dir / name
        path.write_text(textwrap.dedent(content), encoding="utf-8")
        return path

    def test_missing_file_raises_config_error(self):
        """
        Test that a nonexistent path raises a ConfigError rather than an
        OS-level FileNotFoundError.
        """
        with self.assertRaises(ConfigError):
            load_yaml_config(self.tmp_dir / "does_not_exist.yml")

    def test_valid_yaml_overrides_all_fields(self):
        """
        Test that a fully-specified YAML file is loaded into a matching
        RuntimeConfiguration.
        """
        path = self._write(
            "config.yml",
            """\
            footprint_points_elliptical: 16
            footprint_points_rectangular_side: 4
            repeat_cycle_delta_position_m: 500
            repeat_cycle_delta_velocity_m_per_s: 1
            repeat_cycle_search_elevation_deg: 45
            repeat_cycle_search_duration_days: 7
            repeat_cycle_lazy_load: false
            repeat_cycle_for_orbit_track: false
            repeat_cycle_for_observation_events: false
            gp_orbit_lazy_load: false
            """,
        )
        rc = load_yaml_config(path)
        self.assertIsInstance(rc, RuntimeConfiguration)
        self.assertEqual(rc.footprint_points_elliptical, 16)
        self.assertEqual(rc.footprint_points_rectangular_side, 4)
        self.assertEqual(rc.repeat_cycle_delta_position_m, 500)
        self.assertEqual(rc.repeat_cycle_delta_velocity_m_per_s, 1)
        self.assertEqual(rc.repeat_cycle_search_elevation_deg, 45)
        self.assertEqual(rc.repeat_cycle_search_duration_days, 7)
        self.assertFalse(rc.repeat_cycle_lazy_load)
        self.assertFalse(rc.repeat_cycle_for_orbit_track)
        self.assertFalse(rc.repeat_cycle_for_observation_events)
        self.assertFalse(rc.gp_orbit_lazy_load)

    def test_partial_yaml_falls_back_to_defaults_for_missing_fields(self):
        """
        Test that fields omitted from the YAML file retain their
        RuntimeConfiguration default rather than raising a validation
        error.
        """
        path = self._write("config.yml", "footprint_points_elliptical: 16\n")
        rc = load_yaml_config(path)
        self.assertEqual(rc.footprint_points_elliptical, 16)
        self.assertEqual(rc.footprint_points_rectangular_side, 8)

    def test_unrecognized_key_is_silently_ignored(self):
        """
        Test that an unrecognized top-level key in the YAML file does not
        raise an error. RuntimeConfiguration relies on pydantic's default
        "ignore" behavior for unknown fields, so a typo'd setting name is
        currently dropped without any indication to the user.
        """
        path = self._write("config.yml", "not_a_real_setting: 123\n")
        rc = load_yaml_config(path)
        self.assertEqual(rc, RuntimeConfiguration())

    def test_yaml_value_violating_field_constraint_raises_config_error(self):
        """
        Test that a value violating a field's constraint (here, fewer than
        4 elliptical footprint points) surfaces as a ConfigError rather
        than a raw pydantic ValidationError.
        """
        path = self._write("config.yml", "footprint_points_elliptical: 2\n")
        with self.assertRaises(ConfigError):
            load_yaml_config(path)

    def test_malformed_yaml_parser_error_raises_config_error(self):
        """
        Test that YAML content triggering PyYAML's ParserError (here, an
        unclosed flow sequence) is caught and re-raised as a ConfigError.
        """
        path = self._write("config.yml", "footprint_points_elliptical: [1, 2\n")
        with self.assertRaises(ConfigError):
            load_yaml_config(path)

    def test_malformed_yaml_scanner_error_is_not_converted_to_config_error(self):
        """
        Known gap (flagged during review, fix deferred): load_yaml_config
        only catches yaml.parser.ParserError, not the broader
        yaml.YAMLError hierarchy. Malformed YAML that instead raises a
        ScannerError (e.g. inconsistent indentation) currently propagates
        uncaught rather than becoming a ConfigError. This test pins down
        that actual current behavior so a future fix is a deliberate,
        visible change rather than a silent one.
        """
        path = self._write("config.yml", "a: 1\n  b: 2\n")
        with self.assertRaises(yaml.YAMLError) as ctx:
            load_yaml_config(path)
        self.assertNotIsInstance(ctx.exception, ConfigError)

    def test_empty_yaml_file_raises_type_error(self):
        """
        Known gap (flagged during review, fix deferred): an empty YAML
        file parses to None via yaml.safe_load, and
        RuntimeConfiguration(**None) raises a TypeError rather than a
        pydantic ValidationError, so load_yaml_config does not convert it
        to a ConfigError. This test pins down that actual current
        behavior so a future fix is a deliberate, visible change rather
        than a silent one.
        """
        path = self._write("config.yml", "")
        with self.assertRaises(TypeError):
            load_yaml_config(path)


class TestPackagedDefaults(unittest.TestCase):
    """
    Unit tests for the packaged resources/defaults.yml file and the
    lazily-cached singleton get_rc() populates from it.
    """

    def test_defaults_yaml_ships_with_the_package(self):
        """
        Test that resources/defaults.yml is present on disk relative to
        the installed tatc.config module. Regression test for a packaging
        bug where pyproject.toml's package-data only listed
        "resources/*.bsp", so defaults.yml was silently omitted from
        built wheels/sdists and every real (non-editable) install fell
        back to hard-coded defaults regardless of this file's contents.
        """
        resources_path = Path(tatc_config.__file__).parent / "resources" / "defaults.yml"
        self.assertTrue(
            resources_path.exists(),
            "resources/defaults.yml must ship with the package (see "
            "pyproject.toml [tool.setuptools.package-data])",
        )

    def test_defaults_yaml_matches_expected_values(self):
        """
        Test that resources/defaults.yml contains the exact values
        documented in this repository, guarding against silent drift
        between the shipped defaults file and RuntimeConfiguration's own
        hard-coded defaults.
        """
        resources_path = Path(tatc_config.__file__).parent / "resources" / "defaults.yml"
        with open(resources_path, "r", encoding="utf-8") as f:
            raw = yaml.safe_load(f)
        self.assertEqual(
            raw,
            {
                "footprint_points_elliptical": 32,
                "footprint_points_rectangular_side": 8,
                "repeat_cycle_delta_position_m": 10000,
                "repeat_cycle_delta_velocity_m_per_s": 10,
                "repeat_cycle_search_elevation_deg": 88,
                "repeat_cycle_search_duration_days": 30,
                "repeat_cycle_lazy_load": True,
                "repeat_cycle_for_orbit_track": True,
                "repeat_cycle_for_observation_events": True,
                "gp_orbit_lazy_load": True,
            },
        )

    def test_module_level_rc_is_a_runtime_configuration(self):
        """
        Test that get_rc() (consumed elsewhere in the codebase as
        tatc.config.get_rc(), and via the tatc.config.rc backward-compat
        accessor) returns a valid RuntimeConfiguration instance matching
        the packaged defaults.
        """
        self.assertIsInstance(tatc_config.get_rc(), RuntimeConfiguration)
        self.assertEqual(tatc_config.get_rc(), RuntimeConfiguration())
        self.assertIsInstance(tatc_config.rc, RuntimeConfiguration)


class TestImportTimeFallback(unittest.TestCase):
    """
    Unit tests for the fallback behavior of get_rc() when the packaged
    resources/defaults.yml cannot be loaded.
    """

    def test_fallback_logs_warning_and_still_produces_valid_rc(self):
        """
        Test that if the packaged resources/defaults.yml is missing (e.g.
        the packaging regression this module was patched to fix),
        get_rc() still succeeds by falling back to hard-coded defaults,
        and now logs a warning rather than failing silently.
        """
        resources_path = Path(tatc_config.__file__).parent / "resources" / "defaults.yml"
        backup_dir = Path(tempfile.mkdtemp())
        backup_path = backup_dir / "defaults.yml"
        shutil.copy2(resources_path, backup_path)
        # get_rc() is cached, so force it to reload from disk within this
        # test, and again afterward to pick the real file back up. This
        # only busts the cache, it does not touch class/module identity,
        # so it's safe regardless of test execution order.
        tatc_config.reset_rc()
        try:
            resources_path.unlink()
            with self.assertLogs("tatc.config", level="WARNING") as ctx:
                rc = tatc_config.get_rc()
            self.assertTrue(
                any("hard-coded" in message for message in ctx.output)
            )
            self.assertEqual(rc, RuntimeConfiguration())
        finally:
            shutil.copy2(backup_path, resources_path)
            shutil.rmtree(backup_dir, ignore_errors=True)
            tatc_config.reset_rc()


if __name__ == "__main__":
    unittest.main()

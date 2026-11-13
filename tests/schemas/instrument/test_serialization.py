"""
Unit tests for serializing and loading instruments on the systems that carry them.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
import warnings

from pydantic import ValidationError
from pydantic.json_schema import PydanticJsonSchemaWarning

from tatc.schemas import (
    Attitude,
    CircularOrbit,
    ConicalInstrument,
    Instrument,
    PointedInstrument,
    Satellite,
    SurfacePlatform,
    WalkerConstellation,
)
from tatc.utils import VelocityFrame, ViewGeometry


class TestInstrumentSerialization(unittest.TestCase):
    """
    Unit tests for serializing and loading instruments.
    """

    def setUp(self):
        self.orbit = CircularOrbit(altitude=700e3, inclination=98)
        self.base = Instrument(name="Base", field_of_regard=50)
        self.pointed = PointedInstrument(
            name="Pointed",
            field_of_regard=120,
            cross_track_field_of_view=20,
            along_track_field_of_view=1,
            roll_angle=-33.8,
            roll_angle_profile=[(0, -33), (180, -34)],
            is_rectangular=True,
            view_geometry=ViewGeometry.SCAN,
            cross_track_pixels=96,
            scan_period=8 / 3,
        )
        self.conical = ConicalInstrument(
            name="Conical",
            cone_angle=48.5,
            along_track_field_of_view=1,
            scan_half_width=70,
            beam_width=1.2,
        )
        self.instruments = [self.base, self.pointed, self.conical]

    def test_type_values(self):
        """
        Test that each class has its own type, in its dumped data.
        """
        for instrument, expected in zip(
            self.instruments, ("instrument", "pointed", "conical")
        ):
            self.assertEqual(instrument.type, expected)
            self.assertEqual(instrument.model_dump()["type"], expected)

    def test_satellite_dump_keeps_instrument_fields(self):
        """
        Test that dumping a satellite keeps the fields of each instrument class,
        without serializer warnings.
        """
        satellite = Satellite(name="S", orbit=self.orbit, instruments=self.instruments)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            data = satellite.model_dump(mode="json")
        self.assertEqual(data["instruments"][0], self.base.model_dump(mode="json"))
        self.assertEqual(data["instruments"][1], self.pointed.model_dump(mode="json"))
        self.assertEqual(data["instruments"][2], self.conical.model_dump(mode="json"))

    def test_satellite_round_trip(self):
        """
        Test that a satellite with each instrument class equals itself after a
        round trip through JSON.
        """
        satellite = Satellite(name="S", orbit=self.orbit, instruments=self.instruments)
        copy = Satellite.model_validate_json(satellite.model_dump_json())
        self.assertEqual(copy, satellite)
        for instrument, expected in zip(copy.instruments, self.instruments):
            self.assertIs(type(instrument), type(expected))

    def test_satellite_round_trip_with_attitude(self):
        """
        Test that a satellite with an attitude, which sets the velocity frame of
        the body, equals itself after a round trip through JSON and through
        Python data.
        """
        satellite = Satellite(
            name="S",
            orbit=self.orbit,
            instruments=self.instruments,
            attitude=Attitude(velocity_frame=VelocityFrame.INERTIAL),
        )
        from_json = Satellite.model_validate_json(satellite.model_dump_json())
        from_python = Satellite.model_validate(satellite.model_dump())
        for copy in (from_json, from_python):
            self.assertEqual(copy, satellite)
            self.assertEqual(copy.attitude.velocity_frame, VelocityFrame.INERTIAL)

    def test_satellite_round_trip_with_instrument_velocity_frame(self):
        """
        Test that an instrument's velocity frame survives a round trip if the
        satellite has no attitude.
        """
        pointed = self.pointed.model_copy(
            update={"velocity_frame": VelocityFrame.INERTIAL}
        )
        satellite = Satellite(name="S", orbit=self.orbit, instruments=[pointed])
        copy = Satellite.model_validate_json(satellite.model_dump_json())
        self.assertEqual(copy.instruments[0].velocity_frame, VelocityFrame.INERTIAL)

    def test_explicit_velocity_frame_with_attitude_still_rejected(self):
        """
        Test that an explicit velocity frame on an instrument of a satellite
        with an attitude is still an error.
        """
        data = self.pointed.model_dump(mode="json")
        self.assertIn("velocity_frame", data)
        with self.assertRaises(ValidationError):
            Satellite(
                name="S",
                orbit=self.orbit,
                instruments=[data],
                attitude=Attitude(),
            )

    def test_constellation_round_trip(self):
        """
        Test that a constellation keeps its instrument classes in a round trip,
        and passes them to its members.
        """
        constellation = WalkerConstellation(
            name="W",
            orbit=self.orbit,
            instruments=[self.pointed, self.conical],
            number_satellites=2,
            number_planes=2,
            attitude=Attitude(velocity_frame=VelocityFrame.INERTIAL),
        )
        copy = WalkerConstellation.model_validate_json(constellation.model_dump_json())
        self.assertEqual(copy, constellation)
        self.assertIsInstance(
            copy.generate_members()[0].instruments[1], ConicalInstrument
        )

    def test_surface_platform_round_trip(self):
        """
        Test that a surface platform keeps its pointed instrument in a round trip.
        """
        pointed = PointedInstrument(
            name="Pointed", cross_track_field_of_view=20, along_track_field_of_view=10
        )
        platform = SurfacePlatform(
            name="P", latitude=33, longitude=-112, height=30, instruments=[pointed]
        )
        copy = SurfacePlatform.model_validate_json(platform.model_dump_json())
        self.assertEqual(copy, platform)
        self.assertIsInstance(copy.instruments[0], PointedInstrument)

    def test_instrument_objects_kept(self):
        """
        Test that validating an instrument object keeps the same object.
        """
        satellite = Satellite(name="S", orbit=self.orbit, instruments=self.instruments)
        for instrument, expected in zip(satellite.instruments, self.instruments):
            self.assertIs(instrument, expected)

    def test_type_inferred_from_fields(self):
        """
        Test that data without a type, as written before the field existed,
        loads as the class its fields imply.
        """
        satellite = Satellite(
            name="S",
            orbit=self.orbit,
            instruments=[
                {"name": "A"},
                {"name": "B", "field_of_regard": 50},
                {"cross_track_field_of_view": 20, "along_track_field_of_view": 1},
                {"cone_angle": 48.5, "along_track_field_of_view": 1},
            ],
        )
        self.assertEqual(
            [type(i) for i in satellite.instruments],
            [Instrument, Instrument, PointedInstrument, ConicalInstrument],
        )

    def test_explicit_type_takes_precedence(self):
        """
        Test that an explicit type selects the class, so data that does not
        match it is an error.
        """
        with self.assertRaises(ValidationError):
            Satellite(
                name="S",
                orbit=self.orbit,
                instruments=[{"type": "conical", "cross_track_field_of_view": 20}],
            )

    def test_unknown_type_invalid(self):
        """
        Test that an unknown type is an error.
        """
        for value in ("radar", 1):
            with self.assertRaises(ValidationError):
                Satellite(
                    name="S",
                    orbit=self.orbit,
                    instruments=[{"type": value, "field_of_regard": 50}],
                )

    def test_misspelled_fields_invalid(self):
        """
        Test that misspelled fields are errors, rather than silently giving a
        default field of regard of 180 degrees.
        """
        for data in (
            {"name": "A", "field_of_regardd": 50},
            {"cross_track_field_of_vieww": 20, "along_track_field_of_view": 1},
            {"cone_anglee": 48.5, "along_track_field_of_view": 1},
        ):
            with self.assertRaises(ValidationError):
                Satellite(name="S", orbit=self.orbit, instruments=[data])

    def test_direct_extra_fields_invalid(self):
        """
        Test that each instrument class rejects unknown fields when constructed.
        """
        with self.assertRaises(ValidationError):
            Instrument(field_of_regardd=50)
        with self.assertRaises(ValidationError):
            PointedInstrument(
                cross_track_field_of_view=20, along_track_field_of_view=1, roll=5
            )
        with self.assertRaises(ValidationError):
            ConicalInstrument(cone_angle=48.5, along_track_field_of_view=1, roll=5)

    def test_json_schema(self):
        """
        Test that a JSON schema can be generated for each mode, and that the
        validation schema describes each instrument class.
        """
        with warnings.catch_warnings():
            # defaults of other schemas that are not JSON serializable
            warnings.simplefilter("ignore", PydanticJsonSchemaWarning)
            validation = Satellite.model_json_schema(mode="validation")
            serialization = Satellite.model_json_schema(mode="serialization")
        self.assertIn("PointedInstrument", validation["$defs"])
        self.assertIn("ConicalInstrument", validation["$defs"])
        self.assertEqual(serialization["properties"]["instruments"]["type"], "array")


if __name__ == "__main__":
    unittest.main()

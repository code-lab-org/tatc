"""
Shared fixtures for tatc.analysis unit tests.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""
import unittest

from tatc.schemas import GeneralPerturbationsOrbit, Instrument, Satellite, WalkerConstellation

ISS_TLE = [
    "1 25544U 98067A   22171.11255782  .00008307  00000+0  15444-3 0  9992",
    "2 25544  51.6448 322.0970 0003980 282.3738 231.6559 15.49798078345636",
]


class IssConstellationTestCase(unittest.TestCase):
    """
    Base test case providing a satellite and constellation in an ISS-like
    orbit with a single wide-field instrument.
    """
    def setUp(self):
        self.instrument = Instrument(name="Test", field_of_regard=180.0)
        self.orbit = GeneralPerturbationsOrbit.from_tle(ISS_TLE)
        self.satellite = Satellite(
            name="Test", orbit=self.orbit, instruments=[self.instrument]
        )
        self.constellation = WalkerConstellation(
            name="Test",
            orbit=self.orbit,
            instruments=[self.instrument],
            number_satellites=4,
            number_planes=2,
        )

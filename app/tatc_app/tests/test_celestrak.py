import warnings
from tatc.schemas import TwoLineElements

from .base import TatcTestCase


# suppress `RuntimeError: Event loop is closed` known to exist on Windows
from functools import wraps
from asyncio.proactor_events import _ProactorBasePipeTransport


def silence_event_loop_closed(func):
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        try:
            return func(self, *args, **kwargs)
        except RuntimeError as e:
            if str(e) != "Event loop is closed":
                raise

    return wrapper


_ProactorBasePipeTransport.__del__ = silence_event_loop_closed(
    _ProactorBasePipeTransport.__del__
)
# /suppress `RuntimeError: Event loop is closed` known to exist on Windows


class CelestrakTestCase(TatcTestCase):
    def test_get_tle_valid_user(self):
        # ignore warning for unclosed socket
        warnings.filterwarnings(
            action="ignore", message="unclosed", category=ResourceWarning
        )
        response = self.client.get("/celestrak/tle?name=zarya")
        self.assertEqual(response.status_code, 200)
        self.assertEqual(len(response.text.splitlines()), 3)
        TwoLineElements(tle=response.text.splitlines()[1:3])

    def test_get_tle_empty_response(self):
        # ignore warning for unclosed socket
        warnings.filterwarnings(
            action="ignore", message="unclosed", category=ResourceWarning
        )
        response = self.client.get("/celestrak/tle?name=doesnotexist")
        self.assertEqual(response.status_code, 200)
        self.assertTrue(len(response.text.splitlines()) < 2)

    def test_get_tle_invalid_user(self):
        self.client.app.dependency_overrides = {}
        response = self.client.get("/celestrak/tle?name=zarya")
        self.assertEqual(response.status_code, 401)

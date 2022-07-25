from .base import TatcTestCase

import os


class CesiumTestCase(TatcTestCase):
    def test_get_token_valid_user(self):
        response = self.client.get("/cesium/token")
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.text, os.getenv("TATC_CESIUM_TOKEN"))

    def test_get_token_invalid_user(self):
        self.client.app.dependency_overrides = {}
        response = self.client.get("/cesium/token")
        self.assertEqual(response.status_code, 401)

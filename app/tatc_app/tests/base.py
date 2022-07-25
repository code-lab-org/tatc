import unittest

from uuid import UUID
from ..main import app
from ..utils.users import current_active_user
from ..utils.schemas import UserRead
from fastapi.testclient import TestClient


class TatcTestCase(unittest.TestCase):
    def setUp(self):
        app.dependency_overrides[current_active_user] = lambda: UserRead(
            id=UUID("12345678123456781234567812345678"),
            email="user@example.com",
            is_active=True,
            is_verified=False,
            is_superuser=False,
        )
        self.client = TestClient(app)

    def tearDown(self):
        app.dependency_overrides = {}

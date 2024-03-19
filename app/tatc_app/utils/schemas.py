# -*- coding: utf-8 -*-
"""
Schema specifications for general utilities.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import uuid

from fastapi_users import schemas


class UserRead(schemas.BaseUser[uuid.UUID]):
    """
    User model for read access.
    """


class UserCreate(schemas.BaseUserCreate):
    """
    User model for creating a new user.
    """


class UserUpdate(schemas.BaseUserUpdate):
    """
    User model for updating an existing user.
    """

import uuid

from fastapi_users import schemas


class UserRead(schemas.BaseUser[uuid.UUID]):
    """
    User model for read access.
    """

    pass


class UserCreate(schemas.BaseUserCreate):
    """
    User model for creating a new user.
    """

    pass


class UserUpdate(schemas.BaseUserUpdate):
    """
    User model for updating an existing user.
    """

    pass

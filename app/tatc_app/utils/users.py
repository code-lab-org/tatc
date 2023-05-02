import os
import uuid
from typing import Optional
import contextlib


from fastapi import Depends, Request
from fastapi_users import BaseUserManager, FastAPIUsers, UUIDIDMixin
from fastapi_users.authentication import (
    AuthenticationBackend,
    BearerTransport,
    CookieTransport,
    JWTStrategy,
)
from fastapi_users.db import SQLAlchemyUserDatabase
from fastapi_users.exceptions import UserAlreadyExists

from .db import User, get_async_session, get_user_db
from .schemas import UserCreate


SECRET = os.getenv("TATC_SECRET", "change me")
LOGIN_LIFETIME = int(os.getenv("TATC_LOGIN_LIFETIME_SECONDS", 7200))


# define the user manager
class UserManager(UUIDIDMixin, BaseUserManager[User, uuid.UUID]):
    """
    Manager to bind actions to user management actions.
    """

    async def on_after_register(self, user: User, request: Optional[Request] = None):
        print(f"User {user.id} has registered.")


# define a callback to get the user manager
async def get_user_manager(user_db: SQLAlchemyUserDatabase = Depends(get_user_db)):
    """
    Gets the user manager.

    Args:
        user_db (SQLAlchemyUserDatabase): The user database.

    Returns:
        UserManager: The user manager.
    """
    yield UserManager(user_db)


# define the bearer and cookie transports for login
bearer_transport = BearerTransport(tokenUrl="auth/login")
cookie_transport = CookieTransport(cookie_max_age=LOGIN_LIFETIME)


# define a callback to retrieve the JWT strategy
def get_jwt_strategy() -> JWTStrategy:
    """
    Get the JSON web token (JWT) strategy.

    Returns:
        JWTStrategy: The JWT strategy.
    """
    return JWTStrategy(secret=SECRET, lifetime_seconds=LOGIN_LIFETIME)


# define the JWT and cookie backends
jwt_backend = AuthenticationBackend(
    name="bearer",
    transport=bearer_transport,
    get_strategy=get_jwt_strategy,
)
cookie_backend = AuthenticationBackend(
    name="cookie",
    transport=cookie_transport,
    get_strategy=get_jwt_strategy,
)
# instantiate the FastAPI users module
fastapi_users = FastAPIUsers[User, uuid.UUID](
    get_user_manager,
    [jwt_backend, cookie_backend],
)

# expose context functions
get_async_session_context = contextlib.asynccontextmanager(get_async_session)
get_user_db_context = contextlib.asynccontextmanager(get_user_db)
get_user_manager_context = contextlib.asynccontextmanager(get_user_manager)


# backend function to create a new user
async def create_user(email: str, password: str, is_superuser: bool = False) -> None:
    """
    Create a new user.

    Args:
        email (str): User email.
        password (str): User password.
        is_superuser (bool): `True`, if the user is a superuser, `False` otherwise.
    """
    try:
        async with get_async_session_context() as session:
            async with get_user_db_context(session) as user_db:
                async with get_user_manager_context(user_db) as user_manager:
                    user = await user_manager.create(
                        UserCreate(
                            email=email, password=password, is_superuser=is_superuser
                        )
                    )
                    print(
                        f"Created {'superuser' if is_superuser else 'user'} {email}/{password}."
                    )
    except UserAlreadyExists:
        print(
            f"{'Superuser' if is_superuser else 'User'} {email} already exists, skipping."
        )


# expose active user and superuser dependencies
current_active_user = fastapi_users.current_user(active=True)
current_superuser = fastapi_users.current_user(active=True, superuser=True)

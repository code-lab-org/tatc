import os
from typing import Optional
import contextlib


from fastapi import Depends, Request
from fastapi_users import BaseUserManager, FastAPIUsers
from fastapi_users.authentication import (
    AuthenticationBackend,
    BearerTransport,
    CookieTransport,
    JWTStrategy,
)
from fastapi_users.db import SQLAlchemyUserDatabase
from fastapi_users.manager import UserAlreadyExists

from .db import get_async_session, get_user_db
from .models import User, UserCreate, UserDB, UserUpdate


SECRET = os.getenv("TATC_SECRET", "change me")
LOGIN_LIFETIME = int(os.getenv("TATC_LOGIN_LIFETIME_SECONDS", 7200))


class UserManager(BaseUserManager[UserCreate, UserDB]):
    user_db_model = UserDB

    async def on_after_register(self, user: UserDB, request: Optional[Request] = None):
        print(f"User {user.id} has registered.")


async def get_user_manager(user_db: SQLAlchemyUserDatabase = Depends(get_user_db)):
    yield UserManager(user_db)


bearer_transport = BearerTransport(tokenUrl="auth/login")
cookie_transport = CookieTransport(
    cookie_max_age=LOGIN_LIFETIME
)

def get_jwt_strategy() -> JWTStrategy:
    return JWTStrategy(secret=SECRET, lifetime_seconds=LOGIN_LIFETIME)

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
fastapi_users = FastAPIUsers(
    get_user_manager,
    [jwt_backend, cookie_backend],
    User,
    UserCreate,
    UserUpdate,
    UserDB,
)

get_async_session_context = contextlib.asynccontextmanager(get_async_session)
get_user_db_context = contextlib.asynccontextmanager(get_user_db)
get_user_manager_context = contextlib.asynccontextmanager(get_user_manager)

async def create_user(email: str, password: str, is_superuser: bool = False):
    try:
        async with get_async_session_context() as session:
            async with get_user_db_context(session) as user_db:
                async with get_user_manager_context(user_db) as user_manager:
                    user = await user_manager.create(
                        UserCreate(email=email, password=password, is_superuser=is_superuser)
                    )
                    print(f"Created {'superuser' if is_superuser else 'user'} {email}/{password}.")
    except UserAlreadyExists:
        print(f"{'Superuser' if is_superuser else 'User'} {email} already exists, skipping.")

current_active_user = fastapi_users.current_user(active=True)
current_superuser = fastapi_users.current_user(active=True, superuser=True)

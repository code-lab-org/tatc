from fastapi_users import FastAPIUsers, models
from fastapi_users.authentication import CookieAuthentication, JWTAuthentication
from fastapi_users.db import SQLAlchemyBaseUserTable, SQLAlchemyUserDatabase
import os

from .schemas import User, UserCreate, UserUpdate, UserDB
from .models import UserTable
from .database import Base, database

# load the secret for cookie authentication
SECRET = os.getenv("TATC_SECRET", "change me")
LOGIN_LIFETIME = int(os.getenv("TATC_LOGIN_LIFETIME_SECONDS", 7200))

# configure cookie-based authentication
cookie_authentication = CookieAuthentication(
    secret=SECRET, lifetime_seconds=LOGIN_LIFETIME
)

# configure json web token-based authentication
jwt_authentication = JWTAuthentication(
    secret=SECRET, lifetime_seconds=LOGIN_LIFETIME, tokenUrl="/auth/login"
)

# configure the SQL alchemy database for FastAPI-User
users = UserTable.__table__
user_db = SQLAlchemyUserDatabase(UserDB, database, users)

# configure the FastAPI-User package
fastapi_users = FastAPIUsers(
    user_db,
    [cookie_authentication, jwt_authentication],
    User,
    UserCreate,
    UserUpdate,
    UserDB,
)

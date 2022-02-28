from fastapi_users.db import SQLAlchemyBaseUserTable
from sqlalchemy import Column, String

from .database import Base


class UserTable(Base, SQLAlchemyBaseUserTable):
    name = Column(String)

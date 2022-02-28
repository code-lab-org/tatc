import databases
import os
import sqlalchemy
from sqlalchemy.ext.declarative import DeclarativeMeta, declarative_base
from sqlalchemy.orm import sessionmaker

# define database URL
DATABASE_URL = os.getenv("TATC_DATABASE_URL", "sqlite:///data.db")

# create the database
database = databases.Database(DATABASE_URL)

# define the declarative base
Base: DeclarativeMeta = declarative_base()

# create the database engine
engine = sqlalchemy.create_engine(
    DATABASE_URL, connect_args={"check_same_thread": False}
)

# create a session for route dependencies
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

# function to instantiate a session for route dependencies
def get_db():
    try:
        db = SessionLocal()
        yield db
    finally:
        db.close()

# -*- coding: utf-8 -*-
"""
TAT-C application configuration.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from dotenv import load_dotenv
from fastapi import FastAPI, Depends
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.staticfiles import StaticFiles
from msgpack_asgi import MessagePackMiddleware
import os
import tatc

from .utils.db import create_db_and_tables
from .utils.schemas import UserCreate, UserRead, UserUpdate
from .utils.users import (
    cookie_backend,
    jwt_backend,
    current_active_user,
    current_superuser,
    create_user,
    fastapi_users,
)
from .celery.router import router as celery_router
from .celestrak.router import router as celestrak_router
from .cesium.router import get_cesium_router
from .generation.router import router as generation_router
from .overflight.router import router as overflight_router
from .tracking.router import router as tracking_router
from .coverage.router import router as coverage_router
from .latency.router import router as latency_router

# Load environment variables from the .env file
load_dotenv()
ADMIN_EMAIL = os.getenv("TATC_ADMIN_EMAIL", "admin@example.com")
ADMIN_PASSWORD = os.getenv("TATC_ADMIN_PASSWORD", "admin")
CESIUM_TOKEN = os.getenv("TATC_CESIUM_TOKEN", "")

# Create the FastAPI app
app = FastAPI(
    title="Tradespace Analysis Tool for Constellations (TAT-C)",
    description="Modeling tool for pre-Phase A architecture analysis for Earth science space missions.",
    version=tatc.__version__,
)
# Add MessagePack middleware to allow application/msgpack responses
app.add_middleware(MessagePackMiddleware)
# Add GZip middleware to allow compressed responses
app.add_middleware(GZipMiddleware)
# include the routers for authentication
app.include_router(
    fastapi_users.get_auth_router(cookie_backend), prefix="", tags=["auth"]
)
app.include_router(
    fastapi_users.get_auth_router(jwt_backend), prefix="/auth", tags=["auth"]
)
# include the router for user management
app.include_router(
    fastapi_users.get_users_router(UserRead, UserUpdate),
    prefix="/users",
    tags=["users"],
)
# include the router for user registration
app.include_router(
    fastapi_users.get_register_router(UserRead, UserCreate),
    prefix="/auth",
    tags=["auth"],
    dependencies=[Depends(current_superuser)],
)
# Include the router for the /tasks route
app.include_router(
    celery_router,
    prefix="/tasks",
    dependencies=[Depends(current_active_user)],
)
# Include the router for the /celestrak route
app.include_router(
    celestrak_router,
    prefix="/celestrak",
    dependencies=[Depends(current_active_user)],
)
# Include the router for the /cesium route
app.include_router(
    get_cesium_router(CESIUM_TOKEN),
    prefix="/cesium",
    dependencies=[Depends(current_active_user)],
)
# Include the router for the /generate route
app.include_router(
    generation_router,
    prefix="/generate",
    dependencies=[Depends(current_active_user)],
)
# Include the router for the /analyze route
app.include_router(
    overflight_router,
    prefix="/analyze",
    dependencies=[Depends(current_active_user)],
)
# Include the router for the /analyze route
app.include_router(
    tracking_router,
    prefix="/analyze",
    dependencies=[Depends(current_active_user)],
)
# Include the router for the /analyze route
app.include_router(
    coverage_router,
    prefix="/analyze",
    dependencies=[Depends(current_active_user)],
)
# Include the router for the /analyze route
app.include_router(
    latency_router,
    prefix="/analyze",
    dependencies=[Depends(current_active_user)],
)
# Mount a static directory to the root (/) route for any other requests
app.mount("/", StaticFiles(directory="static", html=True), name="frontend")


# connect to the database and try to create admin user on startup
@app.on_event("startup")
async def startup():
    await create_db_and_tables()
    await create_user(ADMIN_EMAIL, ADMIN_PASSWORD, True)

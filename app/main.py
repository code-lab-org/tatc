from dotenv import load_dotenv
from fastapi import FastAPI, Depends
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi_users.user import UserAlreadyExists
from msgpack_asgi import MessagePackMiddleware
import os

from .utils.database import Base, database, engine
from .utils.dependencies import (
    cookie_authentication,
    jwt_authentication,
    fastapi_users,
    user_db,
)
from .utils.schemas import UserCreate
from .celestrak.router import router as celestrak_router
from .cesium.router import get_cesium_router
from .generation.router import router as generation_router
from .observation.router import router as observation_router
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
    version="3.0.0",
)
# Add MessagePack middleware to allow application/msgpack responses
app.add_middleware(MessagePackMiddleware)
# Add GZip middleware to allow compressed responses
app.add_middleware(GZipMiddleware)
# include the router for cookie authentication
app.include_router(
    fastapi_users.get_auth_router(cookie_authentication),
    prefix="",
    tags=["auth"],
)
# include the router for jwt authentication
app.include_router(
    fastapi_users.get_auth_router(jwt_authentication),
    prefix="/auth",
    tags=["auth"],
)
# include the router for user management
app.include_router(
    fastapi_users.get_users_router(),
    prefix="/users",
    tags=["users"],
)
# include the router for user registration
app.include_router(
    fastapi_users.get_register_router(),
    prefix="/auth",
    tags=["auth"],
    dependencies=[Depends(fastapi_users.current_user(active=True, superuser=True))],
)
# Include the router for the /celestrak route
app.include_router(
    celestrak_router,
    prefix="/celestrak",
    dependencies=[Depends(fastapi_users.current_user(active=True))],
)
# Include the router for the /cesium route
app.include_router(
    get_cesium_router(CESIUM_TOKEN),
    prefix="/cesium",
    dependencies=[Depends(fastapi_users.current_user(active=True))],
)
# Include the router for the /generate route
app.include_router(
    generation_router,
    prefix="/generate",
    dependencies=[Depends(fastapi_users.current_user(active=True))],
)
# Include the router for the /collect route
app.include_router(
    observation_router,
    prefix="/collect",
    dependencies=[Depends(fastapi_users.current_user(active=True))],
)
# Include the router for the /compute route
app.include_router(
    tracking_router,
    prefix="/compute",
    dependencies=[Depends(fastapi_users.current_user(active=True))],
)
# Include the router for the /analyze route
app.include_router(
    coverage_router,
    prefix="/analyze",
    dependencies=[Depends(fastapi_users.current_user(active=True))],
)
# Include the router for the /calculate route
app.include_router(
    latency_router,
    prefix="/compute",
    dependencies=[Depends(fastapi_users.current_user(active=True))],
)
# Mount a static directory to the root (/) route for any other requests
app.mount("/", StaticFiles(directory="static", html=True), name="frontend")
# connect to the database on startup
@app.on_event("startup")
async def startup():
    Base.metadata.create_all(engine)
    await database.connect()
    try:
        await fastapi_users.create_user(
            UserCreate(
                email=ADMIN_EMAIL,
                password=ADMIN_PASSWORD,
                is_superuser=True,
            )
        )
        print(f"Created admin account {ADMIN_EMAIL}/{ADMIN_PASSWORD}.")
    except UserAlreadyExists:
        print(f"Admin account {ADMIN_EMAIL} already exists, skipping.")


# disconnect from the database on shutdown
@app.on_event("shutdown")
async def shutdown():
    await database.disconnect()

from celery import Celery
from skyfield.api import load
import os
from dotenv import load_dotenv

# Load environment variables from the .env file
load_dotenv()

# Create the Celery app
app = Celery(
    "tatc_app",
    broker=os.getenv("TATC_BROKER", "amqp://localhost:5672//"),
    backend=os.getenv("TATC_BACKEND", "redis://localhost:6379/"),
    include=[
        "tatc_app.coverage.tasks",
        "tatc_app.generation.tasks",
        "tatc_app.overflight.tasks",
        "tatc_app.tracking.tasks",
        "tatc_app.utils.tasks",
        "tatc_app.latency.tasks",
    ],
)

app.config_from_object("tatc_app.celeryconfig")

# Pre-load planetary ephemerides on start-up (~15 MB download)
load("de421.bsp")

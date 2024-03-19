# -*- coding: utf-8 -*-
"""
TAT-C worker configuration.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import os
import ssl

from celery import Celery
from skyfield.api import load
from dotenv import load_dotenv

# Load environment variables from the .env file
load_dotenv()

# parse the broker connection
broker_string = os.getenv("TATC_BROKER", "amqp://localhost:5672//")
if "amqps://" in broker_string:
    broker_ssl_option = os.getenv("TATC_BROKER_SSL_CERT_REQS", "NONE")
    broker_ssl_config = {
        "keyfile": os.getenv("TATC_BROKER_SSL_KEYFILE", None),
        "certfile": os.getenv("TATC_BROKER_SSL_CERTFILE", None),
        "ca_certs": os.getenv("TATC_BROKER_SSL_CA_CERTS", None),
        "cert_reqs": (
            ssl.CERT_REQUIRED
            if broker_ssl_option == "REQUIRED"
            else ssl.CERT_OPTIONAL if broker_ssl_option == "OPTIONAL" else ssl.CERT_NONE
        ),
    }
else:
    broker_ssl_config = False

# parse the backend connection
backend_string = os.getenv("TATC_BACKEND", "redis://localhost:6379/")
if "rediss://" in backend_string:
    backend_ssl_option = os.getenv("TATC_BACKEND_SSL_CERT_REQS", "NONE")
    backend_ssl_config = {
        "ssl_keyfile": os.getenv("TATC_BACKEND_SSL_KEYFILE", None),
        "ssl_certfile": os.getenv("TATC_BACKEND_SSL_CERTFILE", None),
        "ssl_ca_certs": os.getenv("TATC_BACKEND_SSL_CA_CERTS", None),
        "ssl_cert_reqs": (
            ssl.CERT_REQUIRED
            if backend_ssl_option == "REQUIRED"
            else (
                ssl.CERT_OPTIONAL if backend_ssl_option == "OPTIONAL" else ssl.CERT_NONE
            )
        ),
    }
else:
    backend_ssl_config = False

# Create the Celery app
app = Celery(
    "tatc_app",
    broker=broker_string,
    broker_use_ssl=broker_ssl_config,
    backend=backend_string,
    redis_backend_use_ssl=backend_ssl_config,
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

# Tradespace Analysis Tool for Constellations (TAT-C) Application

The TAT-C application contains several components:
 * Server: hosts the TAT-C application programming interface (API)
 * Worker: executes low-level tasks
 * Broker (external): coordinates task list between server and worker(s)
 * Backend (external): stores completed task outputs

TAT-C is configured to support the RabbitMQ Broker and Redis Backend.

Documentation: https://tatc.readthedocs.io

## Installation and Usage (Docker)

Due to the installation complexity and platform dependencies for the parallel
task-queuing system (Celery), it is strongly recommended to use Docker
containers to install and run TAT-C.

To enable the Cesium client-side geospatial visualization, your will require
an access token from https://cesium.com/ion/tokens . Paste your access token in
a new file named `.env` in the root of this project in the following format:
```
TATC_CESIUM_TOKEN=your alphanumeric access token here
```
The `.env` file can be read by TAT-C to customize to your account information.

### Docker Compose

Docker Compose can build and operate all four TAT-C applications (server,
worker, broker, and backend). It is suitable for individual workstation or
private use because it does not rely on any networked access between hosts.

Build and run the Docker compose file which provides the server, worker,
broker, and backend components:
```shell
docker-compose up
```

### Docker Images

If you have the broker and backend applications running separately, you can
run the TAT-C server and worker from individual Docker images.

Build the Docker image for the server application.
```shell
docker build --target tatc_server --tag tatc_server .
```
Run a Docker container for the server application.
```shell
docker run -it -p 8000:8000 --env-file=.env tatc_server
```
The client-side graphical user interface is available at http://localhost:8000/

Build the Docker image for the worker application.
```shell
docker build --target tatc_worker --tag tatc_worker .
```
Run a Docker contianer for the worker application.
```shell
docker run -it tatc_worker
```

## Installation and Usage (Direct Python)

This section provides guidance to run the TAT-C server and worker directly in
Python. The conda Python package management system is strongly recommended
because several dependent libraries have binary components.

NOTE: certain libraries have binary dependencies which are incompatible with
Python version 3.8.11. In order to avoid this problem run:
```shell
conda install -c conda-forge python==3.8.10
```

To enable the Cesium client-side geospatial visualization, your will require
an access token from https://cesium.com/ion/tokens . After you create your
account, you *must* add the Asset "Blue Marble Next Generation July, 2004"
from the Asset Depot (ID 3845) to your assets. Paste your access token in
a new file named `.env` in the root of this project in the following format:
```
TATC_CESIUM_TOKEN=your alphanumeric access token here
```
The `.env` file can be read by TAT-C to customize to your account information.

First, create a new conda environment with the TAT-C dependencies and install
the TAT-C core library. Navigate to the `tatc\` directory and run:
```shell
conda env create -f environment.yml
```
Followed by:
```shell
python -m pip install -e .
```

Next, install dependencies for the TAT-C web application. Navigate to the
project root directory and run:
```shell
python -m pip install -r requirements.txt
```

To start the TAT-C server application, run the command
```shell
uvicorn app.main:app --reload --reload-dir app
```

To start the TAT-C worker application, run the command:
```shell
celery -A app.worker worker --loglevel=INFO
```
Note: Celery does not support concurrency on Windows. You can try using `eventlet`:
```shell
pip install eventlet
celery -A app.worker worker --loglevel=INFO --pool=eventlet
```
or the `solo` pool (which does not perform parallel processing):
```shell
celery -A app.worker worker --loglevel=INFO --pool=solo
```

## Environment Variable Configuration

The TAT-C application defines the following environment variables, often
loaded from a `.env` file placed in the project root.

 * TATC_BROKER: broker connection string (default: `amqps://guest:guest@broker:5671//`)
 * TATC_BACKEND: backend connection string (default: `redis://backend:6379/`)
 * TATC_LOGIN_LIFETIME_SECONDS: time until reauthentication required (default: 7200 seconds)
 * TATC_SECRET: cryptographic security secret phrase (default: `change me`)
 * TATC_ADMIN_EMAIL: administrator account email (default: `admin@example.com`)
 * TATC_ADMIN_PASSWORD: administrator account password (default: `admin`)
 * TATC_CESIUM_TOKEN: Cesium access token


## Contact

Principal Investigator: Paul T. Grogan <pgrogan@stevens.edu>

## Acknowledgements

This project was supported in part by the National Aeronautics and Space Administration (NASA) Earth Science Division (ESD) Earth Science Technology Office (ESTO) Advanced Information Systems Technology (AIST) program under grant numbers: NNX17AE06G, 80NSSC17K0586, 80NSSC20K1118, and 80NSSC21K1515.

Current project team:
 * PI: Paul T. Grogan <pgrogan@stevens.edu>
 * I. Josue Tapia-Tamayo
 * Isaac Feldman

Project alumni:
 * Hayden Daly
 * Lindsay Portelli
 * Matthew Sabatini
 * Evan Abel
 * Sigfried Hache


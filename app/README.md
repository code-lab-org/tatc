# TAT-C Application

The TAT-C application contains several components:
 * Server: hosts the TAT-C application programming interface (API)
 * Worker: executes low-level tasks
 * Broker: coordinates task list between server and worker(s)
 * Backend: stores completed task outputs

TAT-C is configured to support the RabbitMQ Broker and Redis Backend.

## Configuration

To enable the Cesium client-side geospatial visualization, your will require an access token from <https://cesium.com/ion/tokens>.
After you create your account, you *must* add the Asset "Blue Marble Next Generation July, 2004" from the Asset Depot (ID 3845) to your assets.
Paste your access token in a new file named `.env` (in this directory) in the following format:
```
TATC_CESIUM_TOKEN=your alphanumeric access token here
```
The `.env` file is read by TAT-C to customize to your account information.

Unless installing via the "Docker Compose" method below, identify the protocol, credentials (if required), URL, and port number for the broker and backend components in the `.env` file.
The default values (below) assume the broker and backend run locally.
```
TATC_BROKER = amqp://guest:guest@broker:5672//
TATC_BACKEND = redis://backend:6379/
```

### Environment Variables

The TAT-C application references the following environment variables, typically loaded from a `.env` file.

 * TATC_BROKER: broker connection string (default: `amqps://guest:guest@broker:5671//`)
 * TATC_BACKEND: backend connection string (default: `redis://backend:6379/`)
 * TATC_LOGIN_LIFETIME_SECONDS: time until reauthentication required (default: `7200` seconds)
 * TATC_SECRET: cryptographic security secret phrase (default: `change me`)
 * TATC_ADMIN_EMAIL: administrator account email (default: `admin@example.com`)
 * TATC_ADMIN_PASSWORD: administrator account password (default: `admin`)
 * TATC_CESIUM_TOKEN: Cesium access token

## Installation and Usage

Due to the installation complexity and platform dependencies for the parallel task-queuing system (Celery), it is strongly recommended to use Docker containers to install and run TAT-C.

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

If you have the broker and backend applications running separately, you can run the TAT-C server and worker from individual Docker images.

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
Run a Docker container for the worker application.
```shell
docker run -it tatc_worker
```

### Direct Python

First, create a new conda environment with the required dependencies
```shell
conda env create -f environment.yml
```
and activate it using
```shell
conda activate tatc_app_env
```

To start the TAT-C server application, run the command
```shell
uvicorn tatc_app.main:app --reload --reload-dir tatc_app
```

To start the TAT-C worker application, run the command:
```shell
celery -A tatc_app.worker worker --loglevel=INFO
```
Note: Celery does not support concurrency on Windows. You can try using `eventlet`:
```shell
conda install -c conda-forge eventlet
celery -A tatc_app.worker worker --loglevel=INFO --pool=eventlet
```
or the `solo` pool (which does not perform parallel processing):
```shell
celery -A tatc_app.worker worker --loglevel=INFO --pool=solo
```

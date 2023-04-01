# Web Application

Starting in version 3.1.3, TAT-C contains an optional web application to expose
analysis functions as RESTful web services, distribute computational load over
multiple worker machines using a parallel task-worker architecture (Celery),
and provide a browser-based graphical user interface (GUI) with geospatial
visualization (Cesium) that allows users to run and visualize results for a
subset of TAT-C analysis functions without requiring Python knowledge.

Running the web application is completely optional and is not enabled by default.

The TAT-C web application contains four major components:
 * Server: hosts the TAT-C application programming interface (API) and graphical user interface (GUI)
 * Worker: executes low-level tasks
 * Broker: coordinates the list of tasks between server and worker(s)
 * Backend: stores completed task outputs

The web application is configured to support the RabbitMQ Broker and Redis Backend.

## Configuration

### Cesium Access Token and Assets

The client-side Cesium geospatial visualization requires an access token
from <https://cesium.com/ion/tokens>. After creating an account, you *must*
add the Asset "Blue Marble Next Generation July, 2004" from the Asset Depot
(ID 3845) to the account assets to enable visualization.

Paste the access token in a new file named `.env` (in the `app/` directory)
with the following format:
```
TATC_CESIUM_TOKEN=your alphanumeric access token here
```
The `.env` file is read by the web application to connect to your account.

### Broker and Backend Connection Strings

Unless installing via the "Docker Compose" method below, identify the protocol,
credentials (if required), domain name/URL, and port number to connect the broker
and backend components in the `.env` file.

The default values (below) assume the broker and backend run locally.
```
TATC_BROKER = amqp://guest:guest@broker:5672//
TATC_BACKEND = redis://backend:6379/
```

### Other Environment Variables

The web application references the following environment variables, typically loaded from a `.env` file.

 * `TATC_BROKER`: broker connection string (default: `amqps://guest:guest@broker:5671//`)
 * `TATC_BACKEND`: backend connection string (default: `redis://backend:6379/`)
 * `TATC_LOGIN_LIFETIME_SECONDS`: time until reauthentication required (default: `7200` seconds)
 * `TATC_SECRET`: cryptographic security secret phrase (default: `change me`)
 * `TATC_ADMIN_EMAIL`: administrator account email (default: `admin@example.com`)
 * `TATC_ADMIN_PASSWORD`: administrator account password (default: `admin`)
 * `TATC_CESIUM_TOKEN`: Cesium access token

## Installation and Usage

Due to the architectural complexity and dependencies for the parallel
task-queuing system (Celery), we strongly recommended using Docker containers
to install and run TAT-C.

### Method 1: Docker Compose (Recommended)

Docker Compose can build and operate all four application components (server,
worker, broker, and backend) locally from one command. It is suitable for
individual workstation use because it does not rely on network access between
components.

From the project root directory (`/`), build and run the Docker compose file
which provides the server, worker, broker, and backend components:
```shell
docker-compose up
```
By default, the client-side GUI is available at <http://localhost:8000/>.


### Method 2: Docker Containers

Docker containers allow the four components (server, worker, broker, and backend)
to run across multiple hosts. Individual containers can be combined with other
services such as reverse proxies to secure connections over the Internet.

The broker and backend services can use the third-party Docker images
[rabbitmq](https://hub.docker.com/_/rabbitmq) and
[redis](https://hub.docker.com/_/redis), respectively. Specify the connection
strings for both in the `.env` file as described above.

For the server application, first build the `tatc_server` Docker image from the
project root (`/`).
```shell
docker build --target tatc_server --tag tatc_server .
```
Next, run a Docker container for the server application (this example maps host
port 8000 to container port 8000).
```shell
docker run -it -p 8000:8000 --env-file=.env tatc_server
```
The client-side GUI is available at <http://localhost:8000/>. Note that the
server needs at least one worker to process tasks.

For the worker application, first build the `tatc_worker` Docker image from the
project root (`/`).
```shell
docker build --target tatc_worker --tag tatc_worker .
```
Next, run a Docker container for the worker application.
```shell
docker run -it --env-file=.env tatc_worker
```
The worker will connect to the broker and wait for new tasks to process.

### Method 3: Python

To use Pythond directly for the server and worker, the broker and backend
components must be running already running (e.g., see Docker Container).
Specify the connection strings for both in the `.env` file as described above.

Follow the TAT-C installation instructions to create and activate a conda
environment `tatc_env` with the TAT-C library installed.

Next, install the additional web application dependencies by installing tatc
with the optional `app` dependency flag:
```shell
pip install -e .[app]
```

To start the TAT-C server application, run the command:
```shell
uvicorn tatc_app.main:app --reload --reload-dir tatc_app
```
The client-side GUI is available at <http://localhost:8000/>.

To start the TAT-C worker application, run the command:
```shell
celery -A tatc_app.worker worker --loglevel=INFO
```
The worker will connect to the broker and wait for new tasks to process. Note
that Celery does not currently support concurrency on Windows. Try using the
`solo` pool option (which does not perform parallel processing):
```shell
celery -A tatc_app.worker worker --loglevel=INFO --pool=solo
```

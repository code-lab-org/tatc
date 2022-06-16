.. _install:

=======================
INSTALLATION AND SET-UP
=======================
The TAT-C application contains several components:
 * Server: hosts the TAT-C application programming interface (API)
 * Worker: executes low-level tasks
 * Broker (external): coordinates task list between server and worker(s)
 * Backend (external): stores completed task outputs

TAT-C is configured to support the RabbitMQ Broker and Redis Backend.


Download
========
The TAT-C files can be downloaded from: https://github.com/code-lab-org/tatc

.. note::
  Certain libraries have binary dependencies which are incompatible with
  Python versions of python newer than 3.8.10. In order to avoid issues please run:

  :code:`conda install -c conda-forge python==3.8.10`


Environment Variable Configuration
==================================
A :code:`.env` file can be created in the root directory with the following information
to be read by TAT-C to customize account information:

* TATC_BROKER: broker connection string (default: :code:`amqps://guest:guest@broker:5671//`)
* TATC_BACKEND: backend connection string (default: :code:`redis://backend:6379/`)
* TATC_LOGIN_LIFETIME_SECONDS: time until reauthentication required (default: 7200 seconds)
* TATC_SECRET: cryptographic security secret phrase (default: :code:`change me`)
* TATC_ADMIN_EMAIL: administrator account email (default: :code:`admin@example.com`)
* TATC_ADMIN_PASSWORD: administrator account password (default: :code:`admin`)

.. note::
  Enabling Cesium client-side geospatial visualization requires
  an access token from https://cesium.com/ion/tokens.  After you create an
  account, you **must** add the Asset "Blue Marble Next Generation July, 2004"
  from the Asset Depot (ID 3845) to your assets. Paste your access token in
  the :code:`.env` file in the root of this project in the following format:

  :code:`TATC_CESIUM_TOKEN=your alphanumeric access token here`



Installation and Usage (Docker)
==================================
Due to the installation complexity and platform dependencies for the parallel
task-queuing system (Celery), it is strongly recommended to use Docker
containers to install and run TAT-C.


Docker Compose
--------------
Docker Compose can build and operate all four TAT-C applications (server,
worker, broker, and backend). It is suitable for individual workstation or
private use because it does not rely on any networked access between hosts.

**Build** and run the Docker compose file which provides the server, worker,
broker, and backend components.

:code:`docker-compose up`


Docker Images
-------------

If you have the broker and backend applications running separately, you can
run the TAT-C server and worker from individual Docker images.

**Build** the Docker image for the server application.

:code:`docker build --target tatc_server --tag tatc_server .`

**Run** a Docker container for the server application.

:code:`docker run -it -p 8000:8000 --env-file=.env tatc_server`

The client-side graphical user interface is available at http://localhost:8000/

**Build** the Docker image for the worker application.

:code:`docker build --target tatc_worker --tag tatc_worker .`

**Run** a Docker contianer for the worker application.

:code:`docker run -it tatc_worker`


Installation and Usage (Direct Python)
======================================
This section provides guidance to run the TAT-C server and worker directly in
Python. The conda Python package management system is strongly recommended
because several dependent libraries have binary components.

NOTE: certain libraries have binary dependencies which are incompatible with
Python version 3.8.11. In order to avoid this problem run:

:code:`conda install -c conda-forge python==3.8.10`

To enable the Cesium client-side geospatial visualization, your will require
an access token from https://cesium.com/ion/tokens . After you create your
account, you *must* add the Asset "Blue Marble Next Generation July, 2004"
from the Asset Depot (ID 3845) to your assets. Paste your access token in
a new file named :code:`.env` in the root of this project in the following format:

:code:`TATC_CESIUM_TOKEN=your alphanumeric access token here`

The :code:`.env` file can be read by TAT-C to customize to your account information.

First, create a new conda environment with the TAT-C dependencies and install
the TAT-C core library. Navigate to the :code:`tatc\` directory and run:

:code:`conda env create -f environment.yml`

Followed by:

:code:`python -m pip install -e .`

Next, install dependencies for the TAT-C web application. Navigate to the
project root directory and run:

:code:`python -m pip install -r requirements.txt``

To start the TAT-C server application, run the command

:code:`uvicorn app.main:app --reload --reload-dir app`

To start the TAT-C worker application, run the command:

:code:`celery -A app.worker worker --loglevel=INFO`

Note: Celery does not support concurrency on Windows. You can try using :code:`eventlet`:

:code:`pip install eventlet`

:code:`celery -A app.worker worker --loglevel=INFO --pool=eventlet`

or the `solo` pool (which does not perform parallel processing):

:code:`celery -A app.worker worker --loglevel=INFO --pool=solo`

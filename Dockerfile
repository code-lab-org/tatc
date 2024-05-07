# This block defines the TAT-C runtime container using the appropriate
# base Python environment.

FROM python:3.10 AS tatc_runtime

WORKDIR /var/tatc
COPY pyproject.toml ./
COPY src src
RUN python -m pip install . --no-cache-dir

# This block defines the TAT-C server container. Using the TAT-C runtime
# container, it installs and starts the server application.

FROM tatc_runtime AS tatc_server

WORKDIR /var/tatc
RUN python -m pip install .[app] --no-cache-dir

COPY app/tatc_app tatc_app
COPY app/static static
COPY app/resources resources

ENV TATC_BROKER=amqp://guest:guest@broker:5672//
ENV TATC_BACKEND=redis://backend:6379/

CMD ["uvicorn", "tatc_app.main:app", "--host", "0.0.0.0", "--port", "8000"]

# This block defines the TAT-C worker container. Using the TAT-C runtime
# container, it installs and starts the worker application.

FROM tatc_runtime AS tatc_worker

WORKDIR /var/tatc
RUN python -m pip install .[app] --no-cache-dir

COPY app/tatc_app tatc_app
COPY app/resources resources

ENV TATC_BROKER=amqp://guest:guest@broker:5672//
ENV TATC_BACKEND=redis://backend:6379/

CMD ["celery", "-A", "tatc_app.worker", "worker", "--uid=nobody", "--gid=nogroup"]

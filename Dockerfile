# This block defines the TAT-C build environment. It assembles the
# requisite Python environment.

FROM condaforge/mambaforge AS tatc_build

COPY tatc/environment.yml .
RUN mamba env create -f environment.yml
RUN mamba install -c conda-forge conda-pack
RUN conda-pack -n tatc_env -o /tmp/env.tar && \
  mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
  rm /tmp/env.tar
RUN /venv/bin/conda-unpack

# This block defines the TAT-C runtime container. Using the built TAT-C
# environment, it installs the TAT-C library.

FROM python:3 AS tatc_runtime
COPY --from=tatc_build /venv /venv
ENV PATH="/venv/bin:$PATH"

WORKDIR /var/tatc
COPY tatc/pyproject.toml tatc/setup.cfg tatc/setup.py ./
COPY tatc/tatc tatc/
RUN python -m pip install . --no-cache-dir --use-feature=in-tree-build

# This block defines the TAT-C server container. Using the TAT-C runtime
# container, it installs and starts the server application.

FROM tatc_runtime AS tatc_server

WORKDIR /var/app

COPY requirements.txt .
RUN python -m pip install --no-cache-dir -r requirements.txt

COPY app app
COPY static static
COPY resources resources

ENV TATC_BROKER=amqp://guest:guest@broker:5672//
ENV TATC_BACKEND=redis://backend:6379/

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]

# This block defines the TAT-C worker container. Using the TAT-C runtime
# container, it installs and starts the worker application.

FROM tatc_runtime AS tatc_worker

WORKDIR /var/app

COPY requirements.txt .
RUN python -m pip install --no-cache-dir -r requirements.txt

COPY app app
COPY resources resources

ENV TATC_BROKER=amqp://guest:guest@broker:5672//
ENV TATC_BACKEND=redis://backend:6379/

CMD ["celery", "-A", "app.worker", "worker", "--uid=nobody", "--gid=nogroup"]

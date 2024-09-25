ARG BASE_IMAGE
ARG BASE_TAG
FROM --platform=linux/amd64 ${BASE_IMAGE}:${BASE_TAG}
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG DEPENDENCY_BLOCKS=test

WORKDIR /app
COPY requirements_data.txt requirements_data.txt
RUN python -m pip install -r requirements_data.txt
COPY src src
COPY pyproject.toml pyproject.toml
RUN --mount=source=.git,target=.git,type=bind \
    --mount=type=secret,id=INDEX_URL \
    PIP_KEYRING_PROVIDER=subprocess python -m pip install \
    --extra-index-url "$(cat /run/secrets/INDEX_URL)" \
    --no-cache-dir '.['"${DEPENDENCY_BLOCKS}"']'

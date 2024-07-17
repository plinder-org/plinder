ARG BASE_IMAGE
ARG BASE_TAG
FROM --platform=linux/amd64 ${BASE_IMAGE}:${BASE_TAG}
ARG MAMBA_DOCKERFILE_ACTIVATE=1
COPY src src
COPY pyproject.toml pyproject.toml
RUN --mount=source=.git,target=.git,type=bind \
    python -m pip install --no-cache-dir .

FROM --platform=linux/amd64 mambaorg/micromamba:git-c160e88-jammy

USER root
ENV DEBIAN_FRONTEND=noninteractive

# Install OS dependencies
RUN apt-get --allow-releaseinfo-change update \
    && apt-get -y install --no-install-recommends \
    apt-utils \
    dialog 2>&1 \
    && apt-get install -y --no-install-recommends \
    git \
    gpg \
    wget \
    man-db \
    procps \
    tree \
    unzip \
    gcc \
    build-essential \
    lsb-release \
    curl \
    vim \
    exuberant-ctags \
    apt-transport-https \
    ca-certificates \
    gnupg \
    sudo \
    libgl1-mesa-glx \
    libxrender1 \
    rsync \
    libtiff-dev \
    && apt-get autoclean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# Install rust toolchain
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y
ENV PATH=/root/.cargo/bin:${PATH}

# Install conda env
COPY env.yml /tmp/env.yml
RUN micromamba install -y -f /tmp/env.yml && \
    micromamba clean --all --yes && \
    rm /tmp/env.yml

# Install gsutil without full gcloud sdk
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | \
    tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
    && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | \
    apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - \
    && apt-get update -y \
    && apt-get install google-cloud-cli -y \
    && apt-get autoclean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

ENV BASH_ENV=/usr/local/bin/_activate_current_env.sh

COPY requirements.txt /tmp/requirements.txt
RUN --mount=type=secret,id=INDEX_URL \
    PIP_KEYRING_PROVIDER=subprocess python -m pip install \
    --extra-index-url "$(cat /run/secrets/INDEX_URL)" \
    --no-cache-dir -r /tmp/requirements.txt && \
    rm /tmp/requirements.txt

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]

CMD ["/bin/bash"]

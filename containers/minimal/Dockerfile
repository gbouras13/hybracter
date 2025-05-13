######################################################
# Minimal contains the bare minimum to run hybracter #
# and requires to run hybracter install before use   #
######################################################
FROM ubuntu:25.04 AS minimal

ENV TERM=xterm-256color

# Install miniforge and hybracter
ARG MINIFORGE_VERSION=25.3.0-1
ARG MINIFORGE_PATH=/opt/miniforge
ENV PATH=${MINIFORGE_PATH}/bin:${PATH}
RUN set -eux; \
    apt-get update \
    && apt-get -y install --no-install-recommends \
        wget=1.24.5-2ubuntu1 \
        ca-certificates=20241223 \
    # Installing miniforge
    # partly based on https://github.com/conda-forge/miniforge-images/blob/8c784163a44af85c5918ac4bede66d3de2074d82/ubuntu/Dockerfilee
    && wget -q -O /tmp/miniforge.sh https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/Miniforge3-${MINIFORGE_VERSION}-Linux-x86_64.sh \
    && bash /tmp/miniforge.sh -b -p ${MINIFORGE_PATH}/ \
    # Installing hybracter
    && mamba install -y \
        bioconda::hybracter==0.11.2 \
    # Cleanup
    && rm -rf /tmp/miniforge.sh \
    && find ${MINIFORGE_PATH} -follow -type f -name '*.a' -delete \
    && find ${MINIFORGE_PATH} -follow -type f -name '*.pyc' -delete \
    && mamba clean -afy \
    && apt-get -y --purge autoremove \
        wget \
        ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

################################################
# Full is built on top of minimal and contains #
# all databases and dependencies               #
################################################
FROM minimal AS full
RUN hybracter install --medaka \
    && hybracter test-hybrid --use-conda --conda-create-envs-only \
    && hybracter test-long --use-conda --conda-create-envs-only \
    && rm -rf hybracter_out

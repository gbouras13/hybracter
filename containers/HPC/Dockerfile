#
# hybracter 
#

FROM --platform=linux/amd64 staphb/unicycler:0.5.1

ARG DEBIAN_FRONTEND="noninteractive"
ARG HYBRACTER_VERSION=0.12.0
ARG MINIFORGE_VERSION=24.11.3-0

# Install required packages and dependencies
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    autoconf \
    automake \
    wget \
    doxygen \
    gnupg \
    git \
    vim-nox \
    gfortran \
    libtool \
    python3-venv \
    ninja-build \
    python3-pip \
    libnuma-dev \
    python3-dev \
    && rm -rf /var/lib/apt/lists/*


#
# Install miniforge
#
RUN set -eux ; \
    wget -q https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/Miniforge3-${MINIFORGE_VERSION}-Linux-x86_64.sh ; \
    bash ./Miniforge3-* -b -p /opt/miniforge3 -s ; \
    rm -rf ./Miniforge3-*
ENV PATH /opt/miniforge3/bin:$PATH

#
# Install conda environment
#
RUN set -eux ; \
    mamba install -y -c conda-forge -c bioconda -c defaults \
    hybracter=${HYBRACTER_VERSION} \
    && mamba clean -af -y

RUN hybracter install --medaka || true && \
    mkdir -p /opt/miniforge3/lib/python3.12/site-packages/hybracter/databases && \
    touch /opt/miniforge3/lib/python3.12/site-packages/hybracter/databases/medaka.flag
RUN ls /opt/miniforge3/lib/python3.12/site-packages/hybracter/databases
RUN hybracter install --medaka
RUN hybracter test-hybrid  --use-conda -t 8 -o hybracter_out_hybrid 
RUN hybracter test-long --use-conda -t 8 -o hybracter_out_long
RUN rm -rf hybracter_out_long 
RUN rm -rf hybracter_out_hybrid


#
# hybracter 
#

FROM --platform=linux/amd64 ubuntu:20.04
FROM staphb/unicycler:0.5.0

ARG DEBIAN_FRONTEND="noninteractive"
ARG LIBFABRIC_VERSION=1.18.1
ARG HYBRACTER_VERSION=0.7.1

# Install required packages and dependencies
RUN apt -y update \
    && apt -y install build-essential \
    wget \
    doxygen \
    gnupg \
    gnupg2 \
    apt-transport-https \
    software-properties-common  \
    git \
    vim \
    gfortran \
    libtool \
    python3-venv \
    ninja-build \
    python3-pip \
    libnuma-dev \
    python3-dev \
    && apt -y remove --purge --auto-remove cmake \
    && wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null\
    | gpg --dearmor - | tee /etc/apt/trusted.gpg.d/kitware.gpg >/dev/null \
    && apt-add-repository -y "deb https://apt.kitware.com/ubuntu/ jammy-rc main" \
    && apt -y update \
    # Build and install libfabric
    && rm -rf /tmp/build \
    && mkdir -p /tmp/build \
    && cd /tmp/build \
    && wget https://github.com/ofiwg/libfabric/archive/refs/tags/v${LIBFABRIC_VERSION}.tar.gz \
    && tar xf v${LIBFABRIC_VERSION}.tar.gz \
    && cd libfabric-${LIBFABRIC_VERSION} \
    && ./autogen.sh \
    && ./configure \
    && make -j 16 \
    && make install \
    && rm -rf /var/lib/apt/lists/*

#
# Install miniforge
#
RUN set -eux ; \
    wget https://github.com/conda-forge/miniforge/releases/${MINIFORGE_VERSION}/download/Miniforge3-Linux-x86_64.sh ; \
    bash ./Miniforge3-* -b -p /opt/miniforge3 -s ; \
    rm -rf ./Miniforge3-*
ENV PATH /opt/miniforge3/bin:$PATH

#
# Install conda environment
# 
RUN set -eux ; \
    mamba install -y -c conda-forge -c bioconda -c defaults \
    hybracter=${HYBRACTER_VERSION}=pyhdfd78af_0 \
    && mamba clean -af -y
ENV PATH /opt/miniforge3/bin:$PATH

RUN hybracter install --medaka \
    & hybracter test-hybrid --threads 8 \
    & hybracter test-long --threads 16  --conda-create-envs-only \
    & rm -rf hybracter_out

#
# hybracter 
#

FROM staphb/unicycler:0.5.0

ARG DEBIAN_FRONTEND="noninteractive"
ARG LIBFABRIC_VERSION=1.18.1
ARG HYBRACTER_VERSION=0.7.1
ARG MINIFORGE_VERSION=24.11.3-0

# Install required packages and dependencies
RUN apt -y update \
    && apt -y install build-essential \
    wget=1.20.3-1ubuntu2.1 \
    doxygen=1.8.17-0ubuntu2  \
    gnupg=2.2.19-3ubuntu2.2  \
    gnupg2=2.2.19-3ubuntu2.2  \
    apt-transport-https=2.0.10  \
    software-properties-common=0.99.9.12  \
    git=1:2.25.1-1ubuntu3.13  \
    vim=2:8.1.2269-1ubuntu5.31  \
    gfortran=4:9.3.0-1ubuntu2  \
    libtool=2.4.6-14 \
    python3-venv=3.8.2-0ubuntu2 \
    ninja-build=1.10.0-1build1 \
    python3-pip=20.0.2-5ubuntu1.11  \
    libnuma-dev=2.0.12-1  \
    python3-dev=3.8.2-0ubuntu2 \
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
    wget https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/Miniforge3-${MINIFORGE_VERSION}-Linux-x86_64.sh ; \
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
    && hybracter test-hybrid --threads 8 \
    && hybracter test-long --threads 16  --conda-create-envs-only \
    && rm -rf hybracter_out

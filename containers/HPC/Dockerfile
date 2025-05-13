#
# hybracter 
#

FROM staphb/unicycler:0.5.1

ARG DEBIAN_FRONTEND="noninteractive"
ARG LIBFABRIC_VERSION=2.0.0
ARG HYBRACTER_VERSION=0.11.2
ARG MINIFORGE_VERSION=24.11.3-0

# Install required packages and dependencies
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
RUN apt-get -y update \
    && apt-get -y install --no-install-recommends \
    build-essential=12.9ubuntu3 \
    autoconf=2.71-2 \
    automake=1:1.16.5-1.3  \
    wget=1.21.2-2ubuntu1.1 \
    doxygen=1.9.1-2ubuntu2 \
    gnupg=2.2.27-3ubuntu2.1 \
    gnupg2=2.2.27-3ubuntu2.1 \
    apt-transport-https=2.4.13 \
    software-properties-common=0.99.22.9 \
    git=1:2.34.1-1ubuntu1.12 \
    vim=2:8.2.3995-1ubuntu2.23 \
    gfortran=4:11.2.0-1ubuntu1 \
    libtool=2.4.6-15build2 \
    python3-venv=3.10.6-1~22.04.1 \
    ninja-build=1.10.1-1 \
    python3-pip=22.0.2+dfsg-1ubuntu0.5 \
    libnuma-dev=2.0.14-3ubuntu2 \
    python3-dev=3.10.6-1~22.04.1 \
    && apt-get -y remove --purge --auto-remove cmake \
    && wget -q -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null\
    | gpg --dearmor - | tee /etc/apt/trusted.gpg.d/kitware.gpg >/dev/null \
    && apt-add-repository -y "deb https://apt.kitware.com/ubuntu/ jammy-rc main" \
    && apt-get -y update \
    && rm -rf /var/lib/apt/lists/* \
    && wget -q https://github.com/ofiwg/libfabric/archive/refs/tags/v${LIBFABRIC_VERSION}.tar.gz \
    && tar xf v${LIBFABRIC_VERSION}.tar.gz

# Build and install libfabric
WORKDIR /data/libfabric-${LIBFABRIC_VERSION}
RUN ./autogen.sh \
    && ./configure \
    && make -j 16 \
    && make install

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

RUN hybracter install --medaka \
    && hybracter test-hybrid --use-conda --conda-create-envs-only \
    && hybracter test-long --use-conda --conda-create-envs-only \
    && rm -rf hybracter_out

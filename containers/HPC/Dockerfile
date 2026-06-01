#
# hybracter
#
# StaPH-B unicycler base gives us:
#   - Ubuntu foundation (apt-get works)
#   - unicycler 0.5.1 installed to /usr/local/bin (system PATH, survives
#     conda env activation so plassembler's sub-env can still find it)

FROM staphb/unicycler:0.5.1

ARG DEBIAN_FRONTEND="noninteractive"
ARG MINIFORGE_VERSION=24.11.3-0

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
ENV PATH="/opt/miniforge3/bin:${PATH}"

#
# Install bioinformatics tools via conda.
# unicycler comes from the base image (/usr/local/bin) — no need to conda-install it.
# hybracter itself is installed from source below (not yet on bioconda for dev builds).
#
RUN set -eux ; \
    conda install -y -c conda-forge -c bioconda -c defaults \
    minimap2 snakemake \
    && conda clean -af -y

#
# Install hybracter from the repo source committed in this build
#
COPY . /opt/hybracter_src
RUN pip install /opt/hybracter_src

RUN HYBRACTER_PKG=$(python3 -c "import hybracter, os; print(os.path.dirname(hybracter.__file__))") && \
    hybracter install --medaka || true && \
    mkdir -p ${HYBRACTER_PKG}/databases && \
    touch ${HYBRACTER_PKG}/databases/medaka.flag
RUN python3 -c "import hybracter, os; print(os.listdir(os.path.join(os.path.dirname(hybracter.__file__), 'databases')))"
RUN hybracter install --medaka
RUN hybracter test-hybrid --use-conda -t 8 -o hybracter_out_hybrid
RUN hybracter test-long --use-conda -t 8 -o hybracter_out_long
RUN rm -rf hybracter_out_long hybracter_out_hybrid

# Minimal contains the bare minimum to run hybracter
# and requires to run hybracter install before use
FROM condaforge/miniforge3:24.11.3-2 AS minimal

ENV DEBIAN_FRONTEND=noninteractive
ENV TERM=xterm-256color

RUN mamba install -y \
        bioconda::hybracter==0.11.2 \
    mamba clean -afy


# Full is built on top of minimal and contains
# all databases and dependencies
FROM minimal AS full
RUN hybracter install --medaka \
    && hybracter test-hybrid --use-conda --conda-create-envs-only \
    && hybracter test-long --use-conda --conda-create-envs-only \
    && rm -rf hybracter_out

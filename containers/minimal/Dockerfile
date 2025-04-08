FROM condaforge/miniforge3:24.11.3-2
RUN mamba install -y \
        bioconda::hybracter==0.11.2 \
        bioconda::unicycler==0.5.1 && \
    mamba clean -afy

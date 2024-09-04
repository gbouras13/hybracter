# Installation

# Google Colab Notebooks

If you don't want to install `plassembler` locally, you can run it without any code using the colab notebook [https://colab.research.google.com/github/gbouras13/plassembler/blob/main/run_plassembler.ipynb](https://colab.research.google.com/github/gbouras13/plassembler/blob/main/run_plassembler.ipynb)

This is only recommend if you have one or a few samples to assemble (it takes a while per sample due to the limited nature of Google Colab resources - probably an hour or two a sample). If you have more than this, a local install is recommended.

# Local Installation

You will need conda and **highly recommended** mamba to run `hybracter`, because it is required for the installation of each compartmentalised environment (e.g. Flye will have its own environment). See the end of this page for steps on how to install mamba.

## Conda

`hybracter` is available to install with `conda`. To install `hybracter` into a conda enviornment called `hybracterENV`:

```bash
mamba create -n hybracterENV  -c bioconda -c conda-forge  hybracter
conda activate hybracterENV
hybracter --help
hybracter install
```

## Pip

`hybracter` is available to install with `pip` . 

You will also need conda or mamba available so `hybracter` can install all the required dependencies. Therefore, it is recommended to install `hybracter` into a conda environment as follows.


```bash
mamba create -n hybracterENV pip
conda activate hybracterENV
pip install hybracter
hybracter --help
hybracter install
```

Mamba is **highly highly** recommend. Please see the [documentation](https://hybracter.readthedocs.io/en/latest/install/) for more details on how to install mamba.

## Source

Alternatively, the development version of `hybracter` (which may include new, untested features) can be installed manually via github. 

```bash
git clone https://github.com/gbouras13/hybracter.git
cd hybracter
pip install -e .
hybracter --help
```

## Docker/Singularity

A Docker/Singularity Linux container image is available for Hybracter (starting from v0.7.1) [here](https://quay.io/repository/gbouras13/hybracter). This will likely be useful for running Hybracter in HPC environments.

* **Note** the container image comes with the database and all environments installed - there is no need to run `hybracter install` or `hybracter test-hybrid`/`hybracter test-long` or to specify a database directory with `-d`.

To install and run v0.8.0 with singularity

```bash

IMAGE_DIR="<the directory you want the .sif file to be in >"
singularity pull --dir $IMAGE_DIR docker://quay.io/gbouras13/hybracter:0.8.0

containerImage="$IMAGE_DIR/hybracter_0.8.0.sif"

# example command with test fastqs
 singularity exec $containerImage    hybracter hybrid-single -l test_data/Fastqs/test_long_reads.fastq.gz \
 -1 test_data/Fastqs/test_short_reads_R1.fastq.gz  -2 test_data/Fastqs/test_short_reads_R2.fastq.gz \
 -o output_test_singularity -t 4 -c 50000
```

# Database Installation

**Note: users (and CI) have reported errors where this does not work due to an MD5 file error. This is due to Zenodo having issues (where the database lives). Waiting a few minutes and trying again usuall works.**

To install the hybracter databases (consisting of a plassembler database) use:

```bash
hybracter install
```

Alternatively, if you would like to specify a different database directory, that can be achieved with `-d` or `--databases`:

```bash
hybracter install -d <database directory>
```

If you would like to download all the Medaka models (e.g. if you are going to run on a HPC node without internet), you can use `-m` or `--medaka`.

```bash
hybracter install -m
```

# Installing Hybracter Environments

One benefit of writing `hybracter` in Snakemake is that it is easy to avoid dependency hell by containerising each environment. 

The containerised environments will automatically be created as required when you run hybracter for the first time. This means the first time you run hybracter, it will take a while as all required enviornments are installed.

For users who will be running `hybracter` offline e.g. on a cluster, this means you should probably run a small test datasets like `hybracter test-hybrid` and `hybracter test-long` where internet is available online (such as on the head node) to ensure all required environments are installed. I am one of these people :)

Finally, but default `hybracter` will use mamba to install your environments. If for some reason you must use conda not mamba, use `--conda-frontend conda` with your hybracter command to force `hybracter` to use conda. 

## Errors with Installing Dependencies & Environments

During installation, `hybracter` will install each separate rule's dependency in a separate conda environment - this is one way to get around the 'dependency hell' problem when you need multiple different tools to work like in `hybracter`. For instance, Flye will have its own environment, as will Dnaapler, Plassembler, Polypolish etc.

However, of course this being bioinformatics, there is a chance `hybracter` may not install some of these environments perfectly on your machine.

If you encounter issues with speficic conda environments (see e.g. [this](https://github.com/gbouras13/hybracter/issues/44) and [this](https://github.com/gbouras13/hybracter/issues/45) which had issues with Plassembler), the way you can troubleshoot is:

* For each rule, there will be a path to `conda-env:` listed in the `hybracter.log` file.

e.g.

```bash
Error in rule plassembler_long:
    jobid: 56
    input: hybracter_out/processing/qc/Sample1_filt_trim.fastq.gz
    output: hybracter_out/processing/plassembler/Sample1/plassembler_plasmids.fasta, hybracter_out/processing/plassembler/Sample1/plassembler_summary.tsv, hybracter_out/versions/Sample1/plassembler.version
    log: hybracter_out/stderr/plassembler_long/Sample1.log (check log file(s) for error details)
    conda-env: <path to conda env>
```

* Activate this environment:

`conda activate <path to conda env>`

* Troubleshoot as desired


# `hybracter` testing

Once you have installed `hybracter` and run `hybracter install`, it is recommended to run a test to make sure hybracter installs all required dependencies. It should take 5-10 minutes.

For example you can run

```bash
hybracter test-hybrid --threads 8
```

and for long (same for Linux and MacOS)

```bash
hybracter test-long --threads 1
```

# Beginner Conda and Mamba Installation

If you are new to using the command-line, please install mamba or conda using either option detailed below.

I would prefer and recommend the miniforge option as it comes with mamba!

##  Install miniforge

##### Install [miniforge](https://github.com/conda-forge/miniforge). This will automatically install conda.

Assuming you are using a Linux x86_64 machine (for other architectures, please replace the URL with the appropriate [here](https://github.com/conda-forge/miniforge)

`curl -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh`

For Mac (Intel, will also work with M1):

`curl -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh`

Install miniforge and follow the prompts.

`sh Miniforge3-Linux-x86_64.sh`


##### After installation is complete, you should add the following channels to your conda configuration:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

##### Finally, I would recommend installing hybracter into a fresh environment. For example to create an environment called hybracterENV with hybracter installed:

```bash
mamba create -n hybracterENV hybracter
conda activate hybracterENV
hybracter -h
```

## Install Miniconda

##### Install [Anaconda](https://www.anaconda.com/products/distribution). I would recommend [miniconda](https://docs.conda.io/en/latest/miniconda.html).

Assuming you are using a Linux x86_64 machine (for other architectures, please replace the URL with the appropriate one on the [miniconda](https://docs.conda.io/en/latest/miniconda.html) website).

`curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`

For Mac (Intel, will also work with M1):

`curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh`

Install miniconda and follow the prompts.

`sh Miniconda3-latest-Linux-x86_64.sh`


##### After installation is complete, you should add the following channels to your conda configuration:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

##### After this, conda should be installed (you may need to restart your terminal). It is recommended that mamba is also installed, as it will solve the enviroment quicker than conda:

```bash
conda install mamba
```

##### Finally, I would recommend installing hybracter into a fresh environment. For example to create an environment called hybracterENV with hybracter installed:

```bash
mamba create -n hybracterENV hybracter
conda activate hybracterENV
hybracter -h
```

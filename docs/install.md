# Installation

You will need conda and **highly recommended** mamba to run `hybracter`, because it is required for the installation of each compartmentalised environment (e.g. Flye will have its own environment). See the end of this page for steps on how to install mamba.

## Pip

`hybracter` is available to install with `pip` (it will be on conda soon). 

You will also need conda or mamba available so `hybracter` can install all the required dependencies. Therefore, it is recommended to install `hybracter` into a conda environment as follows.


```
mamba create -n hybracterENV pip
conda activate hybracterENV
pip install hybracter
hybracter --help
hybracter install
```

Mamba is **highly highly** recommend. Please see the [documentation](https://hybracter.readthedocs.io/en/latest/install/) for more details on how to install mamba.

## Source

Alternatively, the development version of `hybracter` (which may include new, untested features) can be installed manually via github. 

```
git clone https://github.com/gbouras13/hybracter.git
cd hybracter
pip install -e .
hybracter --help
```

# Database Installation

To install the hybracter databases (consisting of a plassembler database) use:

```
hybracter install
```

Alternatively, if you would like to specify a different database directory, that can be achieved with `-d` or `--databases`:

```
hybracter install -d <database directory>
```

# Installing Hybracter Environments

One benefit of writing `hybracter` in Snakemake is that it is easy to avoid dependency hell by containerising each environment. 

The containerised environments will automatically be created as required when you run hybracter for the first time. This means the first time you run hybracter, it will take a while as all required enviornments are installed.

For users who will be running `hybracter` offline e.g. on a cluster, this means you should probably run a small test datasets like `hybracter hybrid-test` and `hybracter long-test` where internet is available online (such as on the head node) to ensure all required environments are installed. I am one of these people :)

Finally, but default `hybracter` will use mamba to install your environments. If for some reason you must use conda not mamba, use `--conda-frontend conda` with your hybracter command to force `hybracter` to use conda. 

# `hybracter` testing

Once you have installed `hybracter` and run `hybracter install`, it is recommended to run a test to make sure hybracter installs all required dependencies. It should take 5-10 minutes.

For example you can run

```
hybracter hybrid-test --threads 8
```

and for long (same for Linux and MacOS)

```
hybracter long-test --threads 1
```

# Beginner Conda and Mamba Installation

If you are new to using the command-line, please install mamba or conda using either option detailed below.

I would prefer and recommend the miniforge option as it comes with mamba!

##  Install miniforge

1. Install [miniforge](https://github.com/conda-forge/miniforge). This will automatically install conda.

Assuming you are using a Linux x86_64 machine (for other architectures, please replace the URL with the appropriate [here](https://github.com/conda-forge/miniforge)

`curl -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh`

For Mac (Intel, will also work with M1):

`curl -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh`

Install miniforge and follow the prompts.

`sh Miniforge3-Linux-x86_64.sh`


2. After installation is complete, you should add the following channels to your conda configuration:

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

3. Finally, I would recommend installing hybracter into a fresh environment. For example to create an environment called hybracterENV with hybracter installed:

```
mamba create -n hybracterENV hybracter
conda activate hybracterENV
hybracter -h
```

## Install Miniconda

1. Install [Anaconda](https://www.anaconda.com/products/distribution). I would recommend [miniconda](https://docs.conda.io/en/latest/miniconda.html).

Assuming you are using a Linux x86_64 machine (for other architectures, please replace the URL with the appropriate one on the [miniconda](https://docs.conda.io/en/latest/miniconda.html) website).

`curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`

For Mac (Intel, will also work with M1):

`curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh`

Install miniconda and follow the prompts.

`sh Miniconda3-latest-Linux-x86_64.sh`


2. After installation is complete, you should add the following channels to your conda configuration:

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

3. After this, conda should be installed (you may need to restart your terminal). It is recommended that mamba is also installed, as it will solve the enviroment quicker than conda:

```
conda install mamba
```

4. Finally, I would recommend installing hybracter into a fresh environment. For example to create an environment called hybracterENV with hybracter installed:

```
mamba create -n hybracterENV hybracter
conda activate hybracterENV
hybracter -h
```

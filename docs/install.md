# Installation

You will need conda and ideally mamba to run `hybracter`, because it is required for the installation of each compartmentalised environment (e.g. Flye will have its own environment). See the end of this page for steps on how to install conda and mamba.

## Installing `hybracter` with Conda

The easiest way to install `hybracter` is via conda. For inexperienced command line users, this method is highly recommended.

```
conda install -c bioconda hybracter
```

This will install all the dependencies along with `hybracter`. The dependencies are listed in environment.yml.

I would recommend using mamba :

```
conda install mamba
mamba install -c bioconda hybracter
```

## Pip

You can also install `hybracter` with pip.

```
pip install hybracter
```

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

Finally, but default `hybracter` will use mamba to install your environments. If you must use conda not mamba, use `--conda-frontend conda` to force `hybracter` to use conda. 

# `hybracter` testing

Once you have installed `hybracter` and run `hybracter install`, it is recommended to run a test to make sure hybracter installs all required compartmentalised environments.

For example you can run

```
hybracter hybrid-test --threads 1
```

or 

```
hybracter long-test --threads 1
```

# Beginner Conda Installation

If you are new to using the command-line, please install conda using the following instructions.

1. Install [Anaconda](https://www.anaconda.com/products/distribution). I would recommend [miniconda](https://docs.conda.io/en/latest/miniconda.html).
2. Assuming you are using a Linux x86_64 machine (for other architectures, please replace the URL with the appropriate one on the [miniconda](https://docs.conda.io/en/latest/miniconda.html) website).

`curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`

For Mac (Intel, will also work with M1):

`curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh`

3. Install miniconda and follow the prompts.

`sh Miniconda3-latest-Linux-x86_64.sh`

4. After installation is complete, you should add the following channels to your conda configuration:

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

5. After this, conda should be installed (you may need to restart your terminal). It is recommended that mamba is also installed, as it will solve the enviroment quicker than conda:

`conda install mamba`

6. Finally, I would recommend installing hybracter into a fresh environment. For example to create an environment called hybracterENV with hybracter installed:

```
mamba create -n hybracterENV hybracter
conda activate hybracterENV
hybracter -h
```


# `hybracter` Advanced Usage

## Snakemake Profiles

It is  highly recommended to run `hybracter` using a Snakemake profile if you are using a HPC and have multiple isolates to assemble. 

An example slurm profile is included in the `profile` [directory](https://github.com/gbouras13/hybracter/tree/main/profiles/slurm), and check out this [link](https://github.com/Snakemake-Profiles) for more detail on other HPC job scheduler profiles.

You can run `hybracter` with a profile using `--profile` e.g.:
```
hybracter hybrid --input <input.csv> --output <output_dir> --threads <threads> --profile profiles/hybracter
```

The following guide below has been copied and modified from the [hecatomb documentation](https://hecatomb.readthedocs.io/en/latest/profiles/)

## About Snakemake profiles

Snakemake profiles are a must-have for running Snakemake pipelines on HPC clusters.
While they can be a pain to set up, you only need to do this once and then life is easy.

For more information, check the [Snakemake documentation on profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles), 
or this recent [blog post on Snakemake profiles](https://fame.flinders.edu.au/blog/2021/08/02/snakemake-profiles-updated).

## Profiles for Hybracter

The example _profile_ `config.yaml` file contains all the Snakemake options for jobs that are submitted to the scheduler.
`hybracter` expects the following in the cluster commands:

 - `resources.time` for time in minutes
 - `resources.mem_mb` for requested memory in Mb
 - `threads` for requested CPUs
 
We have tried to use what we believe is the most common nomenclature for these variables in Snakemake pipelines in the hopes that `hybracter` is compatible with existing Snakemake profiles and the available 
[Cookiecutter profiles for Snakemake](https://github.com/Snakemake-Profiles).

We recommend redirecting STDERR and STDOUT messages to log files using the Snakemake variables `{rule}` and `{jobid}`, 
for instance like this `--output=logs/{rule}/{jobid}.out`.
You should also prepend the scheduler command with a command to make the log directories in case they don't exists
(it can cause errors for some schedulers), in this example like so: `mkdir -p logs/{rule}/ && sbatch ...`. 
This will make troubleshooting easier for jobs that fail due to scheduler issues.

The example profile includes a 'watcher' script.
Snakemake won't always pick up when a scheduler prematurely terminates a job, which is why we need a watcher script.
This line in the config file tells Snakemake how to check on the status of a job:
`cluster-status: ~/.config/snakemake/slurm/slurm-status.py` (be sure to check and update your file path).
The `slurm-status.py` script will query the scheduler with the jobid and report back to Snakemake on the job's status.

## Profile installation examples

We'll walk through two ways to set up a profile for the Slurm workload manager.
If your HPC uses a different workload manager, the process of installing a profile will be similar but different.

## Copy an example profile

We have provided an example profile for the Slurm workload manager that should work for most HPCs using Slurm.
Snakemake will look for profiles in your home directory at:

```bash
~/.config/snakemake/
```

First create a directory for your new profile, we'll call it 'slurm':

```bash
mkdir -p ~/.config/snakemake/slurm
```

Now copy the files for the example slurm profile 
(you can view them [here on GitHub](https://github.com/gbouras13/hybracter/tree/main/profiles/hybracter)):

```bash
# go to your new profile directory
cd ~/.config/snakemake/slurm/
# copy the files (either from GitHub or from where you installed hybracter)
wget https://raw.githubusercontent.com/gbouras13/hybracter/main/profiles/slurm/config.yaml
wget https://raw.githubusercontent.com/gbouras13/hybracter/main/profiles/slurm/slurm_status.py
```

This example includes the necessary `config.yaml` file for profiles and a watcher script called `slurm-status.py`.
Make the watcher script executable:

```bash
chmod +x ~/.config/snakemake/slurm/slurm-status.py
```

Done!
You can now use this profile with hybracter:

```bash
hybracter test-hybrid --profile slurm
```

## Create a profile with cookiecutter

Cookiecutter is a nifty tool for creating projects using a template.
For Snakemake profiles, Cookiecutter takes away a lot of the manual configuration steps involved with setting up a profile for a specific scheduler.
There are currently [Cookiecutter templates for Slurm, SGE, PBS, and several other workload managers](https://github.com/Snakemake-Profiles),
and `hybracter` is intended to be compatible with these profiles.

We will walk through installing the [Slurm profile using Cookiecutter](https://github.com/Snakemake-Profiles/slurm).
To begin, create a new directory for your profile:

```bash
mkdir -p ~/.config/snakemake/slurm
```

Move to this directory and run the Cookiecutter command for this profile:

```bash
cd ~/.config/snakemake/slurm/
cookiecutter https://github.com/Snakemake-Profiles/slurm.git
```

Follow the prompts and you're done!
For our system, we did not need to specify anything, but you may need to specify account information for billing etc.
There is more detail on the [GitHub page for this profile](https://github.com/Snakemake-Profiles/slurm).

Use the profile with Hecatomb:

```bash
hybracter test-hybrid --profile slurm
```


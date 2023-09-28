[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![](https://img.shields.io/static/v1?label=CLI&message=Snaketool&color=blueviolet)](https://github.com/beardymcjohnface/Snaketool)
![GitHub last commit (branch)](https://img.shields.io/github/last-commit/gbouras13/hybracter/dev?color=8a35da)
[![DOI](https://zenodo.org/badge/574521745.svg)](https://zenodo.org/badge/latestdoi/574521745)

# `hybracter`

`hybracter` is an automated long-read first bacterial genome assembly pipeline implemented in Snakemake using [Snaketool](https://github.com/beardymcjohnface/Snaketool).

## Table of Contents

- [`hybracter`](#hybracter)
  - [Table of Contents](#table-of-contents)
  - [Quick Start](#quick-start)
  - [Description](#description)
  - [Documentation](#documentation)
  - [Why Would You Run Hybracter?](#why-would-you-run-hybracter)
  - [Other Options](#other-options)
      - [Trycycler](#trycycler)
      - [Dragonflye](#dragonflye)
  - [Pipeline](#pipeline)
  - [Main Commands](#main-commands)
  - [Input csv](#input-csv)
      - [`hybracter hybrid`](#hybracter-hybrid)
      - [`hybracter long`](#hybracter-long)
  - [Usage](#usage)
      - [`hybracter install`](#hybracter-install)
      - [Installing Dependencies](#installing-dependencies)
      - [`hybracter hybrid`](#hybracter-hybrid-1)
      - [`hybracter hybrid-single`](#hybracter-hybrid-single)
      - [`hybracter long`](#hybracter-long-1)
      - [`hybracter long-single`](#hybracter-long-single)
  - [Outputs](#outputs)
    - [Main Output Files](#main-output-files)
  - [Snakemake Profiles](#snakemake-profiles)
  - [Version Log](#version-log)
  - [System](#system)
  - [Bugs and Suggestions](#bugs-and-suggestions)
- [Citation](#citation)


## Quick Start

`hybracter` is available to install with `pip`. It has no non-python dependencies. 

You will need conda or mamba available so `hybracter` can install all the required dependencies. Therefore, it is recommended to install `hybracter` into a conda environment as follows.

```
mamba create -n hybracterENV pip
conda activate hybracterENV
pip install hybracter
hybracter --help
hybracter install
```

Mamba is **highly highly** recommend. Please see the [documentation](https://hybracter.readthedocs.io/en/latest/install/) for more details on how to install mamba.

When you run `hybracter` for the first time, all the required dependencies will be installed as required, so it will take longer than usual (usually a few minutes). Every time you run it afterwards, it will be a lot faster as the dependenices will be installed.

If you intend to run hybracter offline (e.g. on HPC nodes with no access to the internet), I highly recommend running `hybracter hybrid-test` and/or `hybracter long-test` on a node with internet access so hybracter can download the required dependencies. It should take 5-10 minutes. If your computer/node has internet access, please skip this step.

```
# linux
hybracter hybrid-test --threads 8
# macOS
hybracter hybrid-test --threads 8 --no_polca

# same for both linux and mac
hybracter long-test --threads 8
```

## Description

`hybracter` is designed for assembling bacterial isolate genomes using a long read first assembly approach. 
It scales massively using the embarassingly parallel power of HPC and Snakemake profiles. It is designed for applications where you have isolates with Oxford Nanopore Technologies (ONT) long reads and optionally matched paired-end short reads for polishing.

`hybracter` is desined to straddle the fine line between being as fully feature-rich as possible with as much information as you need to decide upon the best assembly, while also being a one-line automated program. In other words, as awesome as Unicycler, but updated for 2023. Perfect for lazy people like myself.

`hybracter` is largely based off Ryan Wick's [magnificent tutorial](https://github.com/rrwick/Perfect-bacterial-genome-tutorial) and associated [paper](https://doi.org/10.1371/journal.pcbi.1010905). `hybracter` differs in that it adds some additional steps regarding targeted plasmid assembly with [plassembler](https://github.com/gbouras13/plassembler), contig reorientation with [dnaapler](https://github.com/gbouras13/dnaapler) and extra polishing and statistical summaries.

Note: if you have Pacbio reads, as of 2023, you probably can just run [Flye](https://github.com/fenderglass/Flye) or [Dragonflye](https://github.com/rpetit3/dragonflye) (or of course [Trycyler](https://github.com/rrwick/Trycycler) ) and reorient the contigs with [dnaapler](https://github.com/gbouras13/dnaapler) without polishing. See Ryan Wick's [blogpost](https://doi.org/10.5281/zenodo.7703461) for more details. Also, you probably still will get good results with hybracter, but the pre-polished genome will be the highest quality! If you really want this feature to be added, please reach out.

## Documentation

Documentation for `hybracter` is available [here](https://hybracter.readthedocs.io/en/latest/).

## Why Would You Run Hybracter?

* If you want the best possible _automated_ long read only or hybrid bacterial isolate genome assembly.
* In other words, if you love Unicycler like I do, but want something faster and more accurate.
* If you need to assemble many (e.g. 10+) bacterial isolates as efficiently as possible.
* If you want all information about from assembly pipeline such as whether your polishing probably improved the genome, whether your assembly was likely complete, and how many plasmids you probably assembled.

## Other Options

#### Trycycler

If you are looking for the best possible (manual) bacterial assembly for a single isolate, please definitely use [Trycyler](https://github.com/rrwick/Trycycler). 

  * `hybracter` will almost certainly not give you better assemblies than Trycycler. Trycycler is the gold standard for a reason.
  * `hybracter` is automated, scalable, faster and requires less bioinformatics/microbial genomics expertise to run. 
  * If you use Trycycler, I would also highly recommend using (disclaimer: my own program) [plassembler](https://github.com/gbouras13/plassembler) (which is built into hybracter) along side Trycycler to assemble small plasmids if you are especially interested in those, because long read only assemblies often [miss small plasmids](https://doi.org/10.1099/mgen.0.001024).

#### Dragonflye

[Dragonflye](https://github.com/rpetit3/dragonflye) by the awesome @[rpetit3](https://github.com/rpetit3) is a good alternative for automated assembly if `hybracter` doesn't fit your needs, particuarly if you are familiar with [Shovill](https://github.com/tseemann/shovill). Some pros and cons between `hybracter` and `dragonflye` are listed below.

  * `dragonflye` allows for more options with regards to assemblers (it supports [Miniasm](https://github.com/lh3/miniasm) or [Raven](https://github.com/lbcb-sci/raven) as well as Flye).
  * On a single isolate, `dragonflye` should be faster (benchmarking coming but the Plassembler, assessment and extra polishing steps of hybracter should make it slower).
  * `hybracter` should be more accurate, due to the extra round of polishing following reorientation, and integration of Plassembler (benchmarking coming).
  * `hybracter` has the advantage of scalability across multiple samples due to its Snakemake and Snaketool implementation. 
  * So if you have access to a cluster, `hybracter` is for you and likely faster.
  * `hybracter` gives more accurate plasmid assemblies because it uses [plassembler](https://github.com/gbouras13/plassembler)
  * `hybracter` will suggest automatically whether an assembly is 'complete' or 'incomplete'
  * `hybracter` will assess each polishing step and choose the genome most likely to be the best quality.

## Pipeline

<p align="center">
  <img src="img/hybracter.png" alt="Hybracter" height=600>
</p>

1. QC with [Filtlong](https://github.com/rrwick/Filtlong), [Porechop](https://github.com/rrwick/Porechop), [fastp](https://github.com/OpenGene/fastp) and optionally contaminant removal using modules from [trimnami](https://github.com/beardymcjohnface/Trimnami).
2. Long-read assembly with [Flye](https://github.com/fenderglass/Flye). 
3. Determine whether chromosome(s) were assembled (marked as 'complete') or not (marked as 'incomplete') based on the given minimum chromosome length.
4. For complete isolates, plasmid recovery with [Plassembler](https://github.com/gbouras13/plassembler).
5. For all isolates, long read polishing with [Medaka](https://github.com/nanoporetech/medaka).
6. For complete isolates, the chromosome is reorientated to begin with the dnaA gene with [dnaapler](https://github.com/gbouras13/dnaapler).
7. For all isolates, if short reads are provided, short read polishing with [Polypolish](https://github.com/rrwick/Polypolish) and [Polca](https://github.com/alekseyzimin/masurca).
   * **Note: POLCA is not available on MacOS. You must use `--no_polca`**
8. For all isolates, assessment of all assemblies with [ALE](https://github.com/sc932/ALE) for `hybracter hybrid` or [Pyrodigal](https://github.com/althonos/pyrodigal) for `hybracter long`.
9. Selection of the best assembly and output final assembly statistics.

## Main Commands

* `hybracter hybrid`: Assemble multiple genomes from isolates that have long-reads and paired-end short reads.
* `hybracter hybrid-single`: Assembles a single genome from an isolate with long-reads and paired-end short reads. It takes similar parameters to [Unicycler](https://github.com/rrwick/Unicycler).
* `hybracter long`: Assemble multiple genomes from isolates that have long-reads only.
* `hybracter long-single`: Assembles a single genome from an isolate with long-reads only.
* `hybracter install`: Downloads and installs the required `plassembler` database.

```
 _           _                    _            
| |__  _   _| |__  _ __ __ _  ___| |_ ___ _ __ 
| '_ \| | | | '_ \| '__/ _` |/ __| __/ _ \ '__|
| | | | |_| | |_) | | | (_| | (__| ||  __/ |   
|_| |_|\__, |_.__/|_|  \__,_|\___|\__\___|_|   
       |___/


Usage: hybracter [OPTIONS] COMMAND [ARGS]...

  For more options, run: hybracter command --help

Options:
  -h, --help  Show this message and exit.

Commands:
  install        Downloads and installs the plassembler database
  hybrid         Run hybracter with hybrid long and paired end short reads
  hybrid-single  Run hybracter hybrid on 1 isolate
  long           Run hybracter with only long reads
  long-single    Run hybracter long on 1 isolate
  test-hybrid    Test hybracter hybrid
  test-long      Test hybracter long
  config         Copy the system default config file
  citation       Print the citation(s) for hybracter
  version        Print the version for hybracter
```

## Input csv

`hybracter hybrid` and `hybracter long` require an input csv file to be specified with `--input`. No other inputs are required.

* This file requires no headers.
* Other than the reads, `hybracter` requires a value for a lower bound the minimum chromosome length for each isolate in base pairs. It must be an integer.
* `hybracter` will denote contigs about this value as chromosome(s) and if it can recover a chromosome, it will denote the isolate as complete.
* In practice, I suggest choosing 90% of the estimated chromosome size for this value.
* e.g. for _S. aureus_, I'd choose 2500000, _E. coli_, 4000000, _P. aeruginosa_ 5500000.

#### `hybracter hybrid`

* `hybracter hybrid` requires an input csv file with 5 columns. 
* Each row is a sample.
* Column 1 is the sample name you want for this isolate. 
* Column 2 is the long read fastq file.
* Column 3 is the minimum chromosome length for that sample.
* Column 4 is the R1 short read fastq file
* Column 5 is the R2 short read fastq file.

e.g.

```
s_aureus_sample1,sample1_long_read.fastq.gz,2500000,sample1_SR_R1.fastq.gz,sample1_SR_R2.fastq.gz
p_aeruginosa_sample2,sample2_long_read.fastq.gz,5500000,sample2_SR_R1.fastq.gz,sample2_SR_R2.fastq.gz
```

#### `hybracter long`

`hybracter long` also requires an input csv with no headers, but only 3 columns.

* `hybracter long` requires an input csv file with 3 columns. 
* Each row is a sample.
* Column 1 is the sample name you want for this isolate. 
* Column 2 is the long read fastq file.
* Column 3 is the minimum chromosome length for that sample.

e.g.

```
s_aureus_sample1,sample1_long_read.fastq.gz,2500000
p_aeruginosa_sample2,sample2_long_read.fastq.gz,5500000
```

## Usage

#### `hybracter install`

You will first need to install the `hybracter` databases.

```
hybracter install
```

Alternatively, can also specify a particular directory to store them - you will need to specify this with `-d <databases directory>` when you run `hybracter`.

```
hybracter install -d  <databases directory>
```

#### Installing Dependencies

**If you have internet access on the machine or node where you are running hybracter, you can skip this step.**

When you run `hybracter` for the first time, all the required dependencies will be installed as required, so it will take longer than usual (usually a few minutes). Every time you run it afterwards, it will be a lot faster as the dependenices will be installed.

If you intend to run hybracter offline (e.g. on HPC nodes with no access to the internet), I highly recommend running `hybracter hybrid-test` and/or `hybracter long-test` on a node with internet access so hybracter can download the required dependencies. It should take 5-10 minutes.

```
hybracter hybrid-test
hybracter long-test
hybracter --help
```

Once that is done, run `hybracter hybrid` or `hybracter long` as follows.

#### `hybracter hybrid`

```
hybracter hybrid -i <input.csv> -o <output_dir> -t <threads> 
```

* `hybracter hybrid` requires only a CSV file specified with `-i` or `--input`
* `--no_polca` will turn off POLCA polishing. **Must be used on MacOS.**
* Use `--min_length` to specify the minimum long-read length for Filtlong.
* Use `--min_quality` to specify the minimum long-read quality for Filtlong.
* You can specify a FASTA file containing contaminants with `--contaminants`. All long reads that map to contaminants will be filtered out.
  * You can specify Escherichia phage lambda (a common contaminant in Nanopore library preparation) using `--contaminants lambda`.
* `--skip_qc` will skip all read QC steps.
* You can change the `--medakaModel` (all available options are listed in `hybracter hybrid -h`)
* You can change the `--flyeModel` (all available options are listed in `hybracter hybrid -h`)


#### `hybracter hybrid-single`

```
hybracter hybrid-single -l <longread FASTQ> -1 <R1 short reads FASTQ> -2 <R2 short reads FASTQ> -s <sample name> -c <chromosome size> -o <output_dir> -t <threads>  [other arguments]
```

#### `hybracter long`

```
hybracter long -i <input.csv> -o <output_dir> -t <threads> [other arguments]
```

* `hybracter long` requires only a CSV file specified with `-i` or `--input`
* Use `--min_length` to specify the minimum long-read length for Filtlong.
* Use `--min_quality` to specify the minimum long-read quality for Filtlong.
* You can specify a FASTA file containing contaminants with `--contaminants`. All long reads that map to contaminants will be filtered.
  * You can specify Escherichia phage lambda (a common contaminant in Nanopore library preparation) using `--contaminants lambda`.
* `--skip_qc` will skip all read QC steps.
* You can change the `--medakaModel` (all available options are listed in `hybracter long -h`)
* You can change the `--flyeModel` (all available options are listed in `hybracter long -h`)


#### `hybracter long-single`

```
hybracter long-single -l <longread FASTQ> -s <sample name> -c <chromosome size>  -o <output_dir> -t <threads>  [other arguments]
```

## Outputs 

`hybracter` creates a number of output files in different formats. 

For more information about all possible file outputs, please see the documentation here.

### Main Output Files

The main outputs are in the `FINAL_OUTPUT` directory.

This directory will include:

1. `hybracter_summary.tsv` file. This gives the summary statistics for your assemblies with the following columns:


|Sample |Complete (True or False) | Total_assembly_length |	Number_of_contigs |	Most_accurate_polishing_round |	Longest_contig_length |	Number_circular_plasmids |
|--------|-----------------------|-------------------------|-------------------|--------|--|--|


2. `complete` and `incomplete` directories.

All samples that are denoted by hybracter to be complete will have 4 outputs:

   * `sample`_summary.tsv containing the summary statistics for that sample.
   * `sample`_final.fasta containing the final assembly for that sample.
   * `sample`_chromosome.fasta containing only the final chromosome(s) assembly for that sample.
   * `sample`_plasmid.fasta containing only the final plasmid(s) assembly for that sample. Note this may be empty. If this is empty, then that sample had no plasmids.

All samples that are denoted by hybracter to be incomplete will have 2 outputs:

   * `sample`.tsv containing the summary statistics for that sample.
   * `sample`_final.fasta containing the final assembly for that sample.


## Snakemake Profiles

I would highly highly recommend running hybracter using a Snakemake profile. Please see this blog [post](https://fame.flinders.edu.au/blog/2021/08/02/snakemake-profiles-updated) for more details. I have included an example slurm profile in the profile directory, but check out this [link](https://github.com/Snakemake-Profiles) for more detail on other HPC job scheduler profiles.


```
hybracter hybrid --input <input.csv> --output <output_dir> --threads <threads> --profile profiles/hybracter
```

## Version Log

A brief description of what is new in each update of `hybracter` can be found in the HISTORY.md file.

## System

`hybracter` is tested on Linux, and on MacOS (with `--no_polca`). 

## Bugs and Suggestions

If you come across bugs with `hybracter`, or would like to make any suggestions to improve the program, please open an issue or email george.bouras@adelaide.edu.au.

# Citation

Please consider also citing these dependencies (especially my own tools Plassembler and Dnaapler :) ):

Plassembler:
* Bouras G., Sheppard A.E., Mallawaarachchi V., Vreugde S., Plassembler: an automated bacterial plasmid assembly tool, Bioinformatics, Volume 39, Issue 7, July 2023, btad409, https://doi.org/10.1093/bioinformatics/btad409. 

Dnaapler:
* Bouras, G., Grigson., S., Papudeshi., B., Mallawaarachchi V., Roach, M. J. (2023) Dnaapler: A tool to reorient circular microbial genomes https://github.com/gbouras13/dnaapler.

Snaketool:
* Roach MJ, Pierce-Ward NT, Suchecki R, Mallawaarachchi V, Papudeshi B, Handley SA, et al. (2022) Ten simple rules and a template for creating workflows-as-applications. PLoS Comput Biol 18(12): e1010705. https://doi.org/10.1371/journal.pcbi.1010705

Ryan Wick et al's Assembling the perfect bacterial genome paper, which provided the intellectual framework for hybracter:
* Wick RR, Judd LM, Holt KE (2023) Assembling the perfect bacterial genome using Oxford Nanopore and Illumina sequencing. PLoS Comput Biol 19(3): e1010905. https://doi.org/10.1371/journal.pcbi.1010905

Trimnami:
* Roach MJ. (2023) Trimnami. https://github.com/beardymcjohnface/Trimnami.

Filtlong:
* Wick RR (2018) Filtlong. https://github.com/rrwick/Filtlong.

Porechop and Porechop_abi:
* Quentin Bonenfant, Laurent Noé, Hélène Touzet, Porechop_ABI: discovering unknown adapters in Oxford Nanopore Technology sequencing reads for downstream trimming, Bioinformatics Advances, Volume 3, Issue 1, 2023, vbac085, https://doi.org/10.1093/bioadv/vbac085
* Wick RR (2017) https://github.com/rrwick/Porechop.

fastp:
* Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu, fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560. 

Flye:
* Kolmogorov, M., Yuan, J., Lin, Y. et al. Assembly of long, error-prone reads using repeat graphs. Nat Biotechnol 37, 540–546 (2019). https://doi.org/10.1038/s41587-019-0072-8

ALE:
* Scott C. Clark, Rob Egan, Peter I. Frazier, Zhong Wang, ALE: a generic assembly likelihood evaluation framework for assessing the accuracy of genome and metagenome assemblies, Bioinformatics, Volume 29, Issue 4, February 2013, Pages 435–443, https://doi.org/10.1093/bioinformatics/bts723

Medaka:
* Oxford Nanopore Technologies, Medaka. https://github.com/nanoporetech/medaka.

Pyrodigal:
* Larralde, M., (2022). Pyrodigal: Python bindings and interface to Prodigal, an efficient method for gene prediction in prokaryotes. Journal of Open Source Software, 7(72), 4296, https://doi.org/10.21105/joss.04296.

Polypolish:
* Wick RR, Holt KE (2022) Polypolish: Short-read polishing of long-read bacterial genome assemblies. PLoS Comput Biol 18(1): e1009802. https://doi.org/10.1371/journal.pcbi.1009802.

POLCA:
* Aleksey V. Zimin, Guillaume Marçais, Daniela Puiu, Michael Roberts, Steven L. Salzberg, James A. Yorke, The MaSuRCA genome assembler, Bioinformatics, Volume 29, Issue 21, November 2013, Pages 2669–2677, https://doi.org/10.1093/bioinformatics/btt476.

Snakemake:
* Mölder F, Jablonski KP, Letcher B et al. Sustainable data analysis with Snakemake [version 1; peer review: 1 approved, 1 approved with reservations]. F1000Research 2021, 10:33 (https://doi.org/10.12688/f1000research.29032.1).


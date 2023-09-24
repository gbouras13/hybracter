[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![](https://img.shields.io/static/v1?label=CLI&message=Snaketool&color=blueviolet)](https://github.com/beardymcjohnface/Snaketool)

# `hybracter`

Hybrid (and long-only) bacterial assembly pipeline for many isolates using Snakemake and [Snaketool](https://github.com/beardymcjohnface/Snaketool).

## Quick Start

`hybracter` is available to install from source only for now with `pip`.

```
git clone "https://github.com/gbouras13/hybracter.git"
cd hybracter/
pip install -e .
hybracter --help
hybracter run --help
```

## Description

`hybracter` is designed for assembling many bacterial isolate genomes using the embarassingly parallel power of HPC and Snakemake profiles. It is designed for applications where you have a number of isolates with ONT long reads and optionally matched paired end short reads for polishing.

It is largely based off Ryan Wick's[magnificent tutorial](https://github.com/rrwick/Perfect-bacterial-genome-tutorial) and associated [paper](https://doi.org/10.1371/journal.pcbi.1010905). I have added some additional steps regarding targetted plasmid assembly and contig reorientation.

## Why Would You Run Hybracter?

* If you want the best possible _automated_ long read only or hybrid bacterial isolate genome assembly.
* If you need to assembly many (10+) bacterial isolates as efficiently as possible.

## Other Options

#### Trycycler

If you are looking for the best possible (manual) bacterial assembly on a single isolate, please definitely use [Trycyler](https://github.com/rrwick/Trycycler). 

  * `hybracter` will not give you better assemblies than Trycycler. Trycycler is the gold standard for a reason.
  * Instead, `hybracter` is automated, scalable, faster and requires less bioinformatics/microbial genomics expertise to run. 
  * I would also highly recommend using (disclaimer: my own program) [plassembler](https://github.com/gbouras13/plassembler) (which is built into hybracter) along side Trycycler to assemble small plasmids if you are especially interested in those, because long read only assemblies often [miss small plasmids](https://doi.org/10.1099/mgen.0.001024).

#### Dragonflye

[Dragonflye](https://github.com/rpetit3/dragonflye) is a good alternative to `hybracter` for automated assembly, particuarly if you are familiar with [Shovill](https://github.com/tseemann/shovill). Some pros and cons between `hybracter` and `dragonflye` are listed below.

  * `dragonflye` allows for more options with regards to assemblers (it supports [Miniasm](https://github.com/lh3/miniasm) or [Raven](https://github.com/lbcb-sci/raven) as well as Flye).
  * On a single isolate, `dragonflye` should be faster (benchmarking coming but the plassembler, assessment and extra polishing steps of hybracter should make it slower).
  * `hybracter` has the advantage of scalability across multiple samples due to its Snakemake and Snaketool implementation. So if you have access to a cluster, `hybracter` is for you and likely faster.
  * `hybracter` gives more accurate plasmid assemblies because it uses [plassembler](https://github.com/gbouras13/plassembler)
  * `hybracter` will suggest automatically whether an assembly is 'complete' or 'incomplete'


## Pipeline

<p align="center">
  <img src="img/hybracter.png" alt="Hybracter" height=600>
</p>

1. QC with [Filtlong](https://github.com/rrwick/Filtlong), [Porechop](https://github.com/rrwick/Porechop) and [fastp](https://github.com/OpenGene/fastp)
2. Long-read assembly with [Flye](https://github.com/fenderglass/Flye). 
3. Determine whether chromosome(s) were assembled (marked as 'complete') or not (marked as 'incomplete') based on the given minimum chromosome length.
4. For complete isolates, plasmid recovery with [Plassembler](https://github.com/gbouras13/plassembler).
5. Long read polishing with [Medaka](https://github.com/nanoporetech/medaka).
6. For complete isolates, chromosome reorientated to begin with the dnaA gene with [dnaapler](https://github.com/gbouras13/dnaapler).
7. If short reads are provided, short read polishing with [Polypolish](https://github.com/rrwick/Polypolish) and [Polca](https://github.com/alekseyzimin/masurca).
8. Assessment of all assemblies with [ALE](https://github.com/sc932/ALE).
9. Selection of the best assembly and output final assembly statistics.

## Commands

* `hybracter hybrid`: Assembles genomes from isolates that have long-reads and paired-end short reads.
* `hybracter long`: Assembles genomes from isolates that have long-reads only.
* `hybracter download`: Downloads the required `plassembler` database.

## Input

* hybracter requires an input csv file to be specified with `--input`. No other inputs are required.
* This file requires no headers.
* Other than the reads, `hybracter` requires a value for a lower bound the minimum chromosome length for each isolate in base pairs.
* It must be an integer.
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

#### `hybracter download`

You will first need to download the Plassembler databases.

```
hybracter download -d  <Plassembler DB directory>
```

Once that is done, run `hybracter hybrid` or `hybracter long` as follows.

#### `hybracter hybrid`

```
hybracter hybrid -i <input.csv> -o <output_dir> -t <threads> -d <Plassembler DB directory>
```

* `hybracter hybrid` requires only a CSV file specified with `-i` or `--input`
* `--no_polca` will turn off POLCA polishing. Must be used on MacOS.
* `--min_length` will let the minimum long-read length for Filtlong.
* `--min_quality` will let the minimum long-read quality for Filtlong.
* `--skip_qc` will skip all read QC.
* You can change the `--medakaModel` (all available options are listed)
* You can change the `--flyeModel` (all available options are listed)

```
hybracter hybrid --help
...
Usage: hybracter hybrid [OPTIONS] [SNAKE_ARGS]...

  Run hybracter with hybrid long and paired end short reads

Options:
  -i, --input TEXT                Input csv  [required]
  --no_polca                      Do not use Polca to polish assemblies with
                                  short reads
  --min_length INTEGER            min read length for long reads
  --min_quality INTEGER           min read quality for long reads
  --skip_qc                       Do not run porechop, filtlong and fastp to
                                  QC the reads
  -d, --databases PATH            Plassembler Databases directory.
  --medakaModel [r1041_e82_400bps_hac_v4.2.0|r1041_e82_400bps_sup_v4.2.0|r941_sup_plant_g610|r941_min_fast_g507|r941_prom_fast_g507|r941_min_fast_g303|r941_min_high_g303|r941_min_high_g330|r941_prom_fast_g303|r941_prom_high_g303|r941_prom_high_g330|r941_min_high_g344|r941_min_high_g351|r941_min_high_g360|r941_prom_high_g344|r941_prom_high_g360|r941_prom_high_g4011|r10_min_high_g303|r10_min_high_g340|r103_min_high_g345|r103_min_high_g360|r103_prom_high_g360|r103_fast_g507|r103_hac_g507|r103_sup_g507|r104_e81_fast_g5015|r104_e81_sup_g5015|r104_e81_hac_g5015|r104_e81_sup_g610|r1041_e82_400bps_hac_g615|r1041_e82_400bps_fast_g615|r1041_e82_400bps_fast_g632|r1041_e82_260bps_fast_g632|r1041_e82_400bps_hac_g632|r1041_e82_400bps_sup_g615|r1041_e82_260bps_hac_g632|r1041_e82_260bps_sup_g632|r1041_e82_400bps_hac_v4.0.0|r1041_e82_400bps_sup_v4.0.0|r1041_e82_260bps_hac_v4.0.0|r1041_e82_260bps_sup_v4.0.0|r1041_e82_260bps_hac_v4.1.0|r1041_e82_260bps_sup_v4.1.0|r1041_e82_400bps_hac_v4.1.0|r1041_e82_400bps_sup_v4.1.0|r941_min_high_g340_rle|r941_min_hac_g507|r941_min_sup_g507|r941_prom_hac_g507|r941_prom_sup_g507|r941_e81_fast_g514|r941_e81_hac_g514|r941_e81_sup_g514]
                                  Medaka Model.  [default:
                                  r1041_e82_400bps_sup_v4.2.0]
  --flyeModel [--nano-hq|--nano-corr|--nano-raw|--pacbio-raw|--pacbio-corr|--pacbio-hifi]
                                  Flye Assembly Parameter  [default: --nano-
                                  hq]
  -o, --output PATH               Output directory  [default: hybracter.out]
  --configfile TEXT               Custom config file [default:
                                  (outputDir)/config.yaml]
  -t, --threads INTEGER           Number of threads to use  [default: 1]
  --use-conda / --no-use-conda    Use conda for Snakemake rules  [default:
                                  use-conda]
  --conda-prefix PATH             Custom conda env directory
  --snake-default TEXT            Customise Snakemake runtime args  [default:
                                  --rerun-incomplete, --printshellcmds,
                                  --nolock, --show-failed-logs, --conda-
                                  frontend mamba]
  -h, --help                      Show this message and exit.
```


## Snakemake Profiles


I would highly highly recommend running hybracter using a Snakemake profile. Please see this blog [post](https://fame.flinders.edu.au/blog/2021/08/02/snakemake-profiles-updated) for more details. I have included an example slurm profile in the profile directory, but check out this [link](https://github.com/Snakemake-Profiles) for more detail on other HPC job scheduler profiles.


```
hybracter run --input <input.csv> --output <output_dir> --threads <threads> --polca --profile profiles/hybracter
```













Assessing Quality
================

https://github.com/rrwick/Perfect-bacterial-genome-tutorial/wiki/Assessing-quality-without-a-reference

Assessing with ALE

If you have short reads, then running ALE to get an ALE score is the best option. Larger ALE scores are better, and since ALE scores are negative, 'larger' means scores with a smaller magnitude and are closer to zero.

To make this easier, we have created the ale_score.sh script (which you can find in the scripts directory of this repo). Just give it the assembly, the two Illumina read files and a thread count like this:

ale_score.sh assembly.fasta reads/illumina_1.fastq.gz reads/illumina_2.fastq.gz 16
And it will produce a file (the assembly filename with .ale appended) which contains the ALE score.

Assessing with Prodigal

If you don't have short reads, then you can run Prodigal and calculate the mean length of predicted proteins. Larger is better because assembly errors (especially indels) often create stop codons.

To make this easier, we have created the mean_prodigal_length.sh script (which you can find in the scripts directory of this repo). Just give it the assembly like this:

mean_prodigal_length.sh assembly.fasta
And it will produce a file (the assembly filename with .prod appended) which contains the mean Prodigal length
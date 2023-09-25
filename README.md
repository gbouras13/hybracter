[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![](https://img.shields.io/static/v1?label=CLI&message=Snaketool&color=blueviolet)](https://github.com/beardymcjohnface/Snaketool)

# `hybracter`

A modern hybrid (and long-only) bacterial assembly pipeline for many isolates using Snakemake and [Snaketool](https://github.com/beardymcjohnface/Snaketool).

## Quick Start

`hybracter` is available to install from source only for now with `pip`.

```
git clone "https://github.com/gbouras13/hybracter.git"
cd hybracter/
pip install -e .
hybracter install
hybracter hybrid-test
hybracter long-test
hybracter --help
```

## Description

`hybracter` is designed for assembling many bacterial isolate genomes using the embarassingly parallel power of HPC and Snakemake profiles. It is designed for applications where you have a number of isolates with Oxford Nanopore Technologies (ONT) long reads and optionally matched paired-end short reads for polishing.

`hybracter` is desined to straddle the fine line between being as fully feature-rich as possible with as much information as you need to decide upon the best assembly, while also being a one-line automated program. Perfect for lazy people like myself :)

`hybracter` is largely based off Ryan Wick's [magnificent tutorial](https://github.com/rrwick/Perfect-bacterial-genome-tutorial) and associated [paper](https://doi.org/10.1371/journal.pcbi.1010905). `hybracter` differs in that it adds some additional steps regarding targetted plasmid assembly with [plassembler](https://github.com/gbouras13/plassembler), contig reorientation with [dnaapler](https://github.com/gbouras13/dnaapler) and extra polishing and statistical summaries.

Note: if you have Pacbio reads, as of 2023, you probably can just run [Flye](https://github.com/fenderglass/Flye) or [Dragonflye](https://github.com/rpetit3/dragonflye) (or of course [Trycyler](https://github.com/rrwick/Trycycler) ) and reorient the contigs with [dnaapler](https://github.com/gbouras13/dnaapler) without polishing. See Ryan Wick's [blogpost](https://doi.org/10.5281/zenodo.7703461) for more details.


## Why Would You Run Hybracter?

* If you want the best possible _automated_ long read only or hybrid bacterial isolate genome assembly.
* If you need to assemble many (e.g. 10+) bacterial isolates as efficiently as possible.
* If you want all  information about from assembly pipeline such as whether your polishing probably improved the genome, whether your assembly was likely complete, and how many plasmids you probably assembled.

## Other Options

#### Trycycler

If you are looking for the best possible (manual) bacterial assembly on a single isolate, please definitely use [Trycyler](https://github.com/rrwick/Trycycler). 

  * `hybracter` will almost certainly not give you better assemblies than Trycycler. Trycycler is the gold standard for a reason.
  * `hybracter` is automated, scalable, faster and requires less bioinformatics/microbial genomics expertise to run. 
  * If you use Trycler, I would also highly recommend using (disclaimer: my own program) [plassembler](https://github.com/gbouras13/plassembler) (which is built into hybracter) along side Trycycler to assemble small plasmids if you are especially interested in those, because long read only assemblies often [miss small plasmids](https://doi.org/10.1099/mgen.0.001024).

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
5. Long read polishing with [Medaka](https://github.com/nanoporetech/medaka).
6. For complete isolates, chromosome reorientated to begin with the dnaA gene with [dnaapler](https://github.com/gbouras13/dnaapler).
7. If short reads are provided, short read polishing with [Polypolish](https://github.com/rrwick/Polypolish) and [Polca](https://github.com/alekseyzimin/masurca).
   * **Note: POLCA is not available on MacOS. You must use `--no_polca`**
8. Assessment of all assemblies with [ALE](https://github.com/sc932/ALE) for `hybracter hybrid` or [Pyrodigal](https://github.com/althonos/pyrodigal) for `hybracter long`.
9. Selection of the best assembly and output final assembly statistics.

## Commands

* `hybracter hybrid`: Assembles genomes from isolates that have long-reads and paired-end short reads.
* `hybracter hybrid-single`: Assembles a single genome from an isolate with long-reads and paired-end short reads. It takes similar parameters to [Unicycler](https://github.com/rrwick/Unicycler).
* `hybracter long`: Assembles genomes from isolates that have long-reads only.
* `hybracter long-single`: Assembles a single genome from an isolate with long-reads only.
* `hybracter install`: Downloads and installs the required `plassembler` database.

## Input csv

* `hybracter hybrid` and `hybracter long` require an input csv file to be specified with `--input`. No other inputs are required.
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

* It also required an input csv with no headers, but only 3 columns.

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

##### Other Arguments 

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

## Snakemake Profiles


I would highly highly recommend running hybracter using a Snakemake profile. Please see this blog [post](https://fame.flinders.edu.au/blog/2021/08/02/snakemake-profiles-updated) for more details. I have included an example slurm profile in the profile directory, but check out this [link](https://github.com/Snakemake-Profiles) for more detail on other HPC job scheduler profiles.


```
hybracter run --input <input.csv> --output <output_dir> --threads <threads> --profile profiles/hybracter
```

## Version Log

A brief description of what is new in each update of `hybracter` can be found in the HISTORY.md file.

## System

`hybracter` is tested on Linux and MacOS with `--no_polca`. 

## Bugs and Suggestions

If you come across bugs with `hybracter`, or would like to make any suggestions to improve the program, please open an issue or email george.bouras@adelaide.edu.au.

# Citation

Bouras, G. (2023). Hybracter: a modern hybrid and long-only bacterial assembly pipeline for many isolates https://github.com/gbouras13/hybracter.  

Please consider also citing the following dependencies, especially my own tools [plassembler](https://github.com/gbouras13/plassembler) and [dnaapler](https://github.com/gbouras13/dnaapler) :)

Plassembler:
https://doi.org/10.1093/bioinformatics/btad409

Dnaapler:
https://github.com/gbouras13/dnaapler

Snaketool:
https://doi.org/10.31219/osf.io/8w5j3

[Wick](https://github.com/rrwick) et al.'s Assembling the perfect bacterial genome paper (provided the intellectual framework for hybracter):
https://doi.org/10.1371/journal.pcbi.1010905

Trimnami:
https://github.com/beardymcjohnface/Trimnami

Filtlong:
https://github.com/rrwick/Filtlong

Porechop and Porechop_abi:
https://doi.org/10.1093/bioadv/vbac085
https://github.com/rrwick/Porechop

fastp:
https://doi.org/10.1093/bioinformatics/bty560

Flye:
https://doi.org/10.1038/s41587-019-0072-8

ALE:
https://doi.org/10.1093/bioinformatics/bts723

Medaka:
https://github.com/nanoporetech/medaka

Pyrodigal:
https://doi.org/10.21105/joss.04296

Polypolish:
https://doi.org/10.1371/journal.pcbi.1009802

POLCA:
https://doi.org/10.1093/bioinformatics/btt476

Snakemake:
https://doi.org/10.12688/f1000research.29032.1

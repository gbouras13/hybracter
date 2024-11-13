# History

## v0.10.1 Updates (13 November 2024)

* Adds retry functionality for dnaapler with 1 thread if there is an error - for some genomes on some systems, [it was observed](https://github.com/gbouras13/hybracter/issues/54) that using the default resources (8 threads, 16GB RAM) will lead to an error in dnaapler. 
* Thanks @[richardstoeckl](https://github.com/richardstoeckl) for implementing this


## v0.10.0 Updates (18 October 2024)

* Updates Medaka to v2.0.1, implementing the `--bacteria` option by default.
* This is based on the recommendations of Ryan Wick [here](https://rrwick.github.io/2024/10/17/medaka-v2.html) who found it improved assemblies due to (likely) enhanced methylation error correction.
* If you still want to specify a Medaka model, the flag `--medaka_override` has been added. You need to include this along with your model via `--medakaModel`. This is most likely useful for older R9 data.
* Adds `--extra_params_flye` parameter if you want to specify extra commands for the Flye assembly step thanks @pdobbler.

## v0.9.1 Updates (8 October 2024)

* Small change to the `plassembler.yaml` config preventing installation bugs - Unicycler v0.5.1 to be installed in a much simpler fashion via Bioconda

### v0.9.0 Updates (18 September 2024)

**`--auto` for automatic estimation of chromosome size**

* Thanks to an [issue](https://github.com/gbouras13/hybracter/issues/90) and code from @[richardstoeckl](https://github.com/richardstoeckl), Hybracter can now estimate the estimated chromosome size for each sample by passing `--auto`. 
* The implementation uses [kmc](https://github.com/refresh-bio/KMC). Specifically, Hybracter uses kmc to count the number of unique 21mers that appear at least 10 times in your long-read FASTQ file. This is because, for a given assembly of length L,  and a k-mer size of k, the total number of unique possible k-mers  will be given by ( L – k ) + 1, and if L >> k, then it suffices as an estimate of total assembly size
* The estimated chromosome size used by Hybracter will actually be 80% of the number of 21-mers found at least 10 times, as it needs to account for plasmids
* If you aren't sure whether you have enough data for assembly (i.e. coverage lower than 20x), be careful using `--auto`, because the actual assembly size will tend to be larger than the number of unique 21mers found at least 10 times. Therefore, the estimated chromosome size will almost certainly be an underestimate and may lead to Hybracter considering your assembly "complete" when in fact it isn't.

* If you use `--auto`, you do not need to specify the chromosome length in the input. This means you don't need to `-c` with `long-single` or `hybrid-single` and in the input csv sample sheet, you do not need a column with chromosome length.
  
e.g. for `hybracter long` you only need 2 columns with sample name and long-read FASTQ file path:

```bash
s_aureus_sample1,sample1_long_read.fastq.gz
p_aeruginosa_sample2,sample2_long_read.fastq.gz
```

and for `hybracter hybrid` you only need 4 columns with sample name, long-read FASTQ, and R1 and R2 short-read FASTQ file paths:

```bash
s_aureus_sample1,sample1_long_read.fastq.gz,sample1_SR_R1.fastq.gz,sample1_SR_R2.fastq.gz
p_aeruginosa_sample2,sample2_long_read.fastq.gz,sample2_SR_R1.fastq.gz,sample2_SR_R2.fastq.gz
```

**Other changes**

* Hybracter v0.9.0 will automatically support the reorientation of archaeal chromosomes (thanks @[richardstoeckl](https://github.com/richardstoeckl)) to begin with the cog1474 Orc1/cdc6 gene.
* `--datadir` can now also accept 2 paths separated by a comma, if you have long reads and short reads in separate directories e.g. `--datadir "long_read_dir,short_read_dir"` (https://github.com/gbouras13/hybracter/issues/76).
* `--min_depth` parameter added. Hybracter will error out if your QC'd long reads have a coverage lower than `min_depth` for a sample (https://github.com/gbouras13/hybracter/issues/89).

## v0.8.0 (30 August 2024)

* Add `--datadir` that removes the need to add full paths in sample sheet (thanks @oschwengers)
* Update medaka to v1.12.1 to support the newest models (#84)
    * New default medaka model is `r1041_e82_400bps_sup_v5.0.0`
* Adds `--mac` flag if you are running Hybracter on MacOS - it is now recommended from to run Hybracter on Linux if you want the latest Medaka models. 
    * This is because ONT do not support bioconda install anymore and the latest version (v1.12.1) from pip doesn't work on Mac
    * `--mac` will install and run Medaka v1.8.0 as in previous versions and use `r1041_e82_400bps_sup_v4.2.0` as default

## v0.7.3 (4 April 2024)

* Enforce spades>=v3.15.2 in the `plassembler.yaml` environment
* For some reason, the environment on Linux environments was being solved for v3.14.1, which was causing an error with Unicycler within Plassembler for some samples described (https://github.com/rrwick/Unicycler/issues/318) 

## v0.7.2 (2 April 2024)

* Adds 'circualr=True' to chromosome contig headers where Flye has marked these as such. This bug was introduced in v0.7.0.
* Thanks Nicole Lerminiaux for spotting this
 
## v0.7.1 (13 March 2024)

* Fixes bug where `hybracter install -d db_dir` would not work as the `-f` parameter was not being passed to Plassembler. Thanks @npbhavya

## v0.7.0 (04 March 2024)

**Bug fixes**

* Fixes bug where `--configfile` wasn't being passed to Hybracter.
* Fixes bug where `hybracter` would crash if the input long reads were not gzipped #51 thanks @wanyuac.

**Changes to short read polishing.**

    * Logic added to run `polypolish` v0.6.0 with `--careful` and skip pypolca if the SR coverage estimate is below 5x (note: FASTA files for pypolca will be generated in the processing directory to play nice with Snakemake, but these will be identical to the polypolish output). 
    * For 5-25x coverage, `polypolish --careful` and `pypolca` with `--careful` will be run. 
    * For >25x coverage, `polypolish` default and `pypolca` with `--careful` will be run. 
    * A preprint justifying these changes will be available soon.

**`--logic` changes**
* By default, `--logic` defaults to `last` for `hybracter hybrid`, as there we have found that the polishing strategy implemented above never makes the assembly worse. We suggest never using `--logic best` with `hybracter hybrid`.

**Changes for chromosome contigs and circularity.**
* If hybracter assembles a contig that is greater than the minimum chromosome length but not marked as circular by Flye, this will now be denoted as a chromosome, but not circular. The genome will be marked as complete also. 
    * These will usually be assemblies with some issue (e.g. prophages, circularisation issues, heterogeneity) and probably require some more attention.
    * For example, with the _Vibrio cholerae_ larger chromosome described [here](https://rrwick.github.io/2024/02/15/misassemblies.html), the genome will be marked as 'complete' but the contig will not be marked as 'circular' in the `hybracter` output.
    * Such contigs will be polished and be in the final `_chromosome.fasta` output, but they will not be rotated by `dnaapler`. 
    * These were previously being excluded, which was missing assemblies with structural heterogeneity (causing the chromosome not to completely circularise) or even bacteria with linear chromosomes like [_Borrelia_](https://www.nature.com/articles/37551). 

**Adds `--depth_filter`** 
* This is passed to [Plassembler](https://github.com/gbouras13/plassembler) and will filter out all putative plasmid contigs that are lower than this depth fraction compared to the chromosome.
* Defaults to 0.25 like Unicycler's implementation.


## v0.6.0 (18 January 2024)

* Fixes bug with Polypolish v0.6.0 breaking the CLI #48
* Adds `-m` option to download all Medaka models with `hybracter install` - useful for offline use #49
* Adds quick SR coverage estimates (in `processing/qc/coverage`) and other QC stats (using [seqkit](https://bioinf.shenwei.me/seqkit/) ) in `processing/qc/seqkit`. This is calculated as (Total bases / estimated chromosome size) for each sample
* Logic added to run Polypolish and pypolca with `--careful` if the SR coverage estimate is below 25x.


## v0.5.0 (08 January 2024)

Ryan Wick recently ran `hybracter long` on the latest Dorado v0.5.0 Nanopore reads (his [blog post](https://rrwick.github.io/2023/12/18/ont-only-accuracy-update.html)). You can read a write-up of the results [here](https://hybracter.readthedocs.io/en/latest/dorado_ryan_louise_0_5_0/). 

* Adds subsampling using `--subsample_depth` using Filtlong, based on some benchmarking of Dorado v0.5.0. Defaults to 100x of the estimated chromosome size `-c`.
* Adds stricter criteria for complete assemblies (aka ensures that identified chromosomes must be circularised according to Flye). Thanks to Matthew Croxen for pointing this out.

## v0.4.1 (28 November 2023)

* Updates code to work with updated version of Plassembler v1.5.0
* Thanks @[npbhavya](https://github.com/npbhavya) for finding this.

## v0.4.0 (14 November 2023)

* Adds `--logic` parameter. You have 2 choices: `--logic best` (the default) or `--logic last`.
* `--logic best` will run `hybracter` as normal and the best assembly (by ALE or pyrodigal mean length) will be selected as the final assembly.
* `--logic last` will force hybracter to pick the last polished round as the final assembly even if it is not the best as per ALE/pyrodigal. So for `hybracter hybrid` this will default to the pypolca polished round. You may wish to use this if you want all your isolates to be consistently assembled.
* Adds reorientation of pre polished chromosome in case it is selected as the best assembly
* Adds fixes to the chromosome comparisons - now it is much easier to interpret any changes between conditions.

## v0.3.0 (8 November 2023)

* Fixes bug relating to polishing. Prior to v0.3.0, hybracter would only polish the chromosome with the entire readset. Benchmarking revealed that if there was significantly similarity between chromosome and plasmids, polishing would introduce errors. Now the entire assembly (chromosome from Flye + plasmids from Plassembler) is polished in every polishing step with improved results. Upgrading and re-running hybracter is recommended.
* Fixes small bugs in output .tsvs.

## v0.2.1 (1 November 2023)

* Fixes small bugs relating to isolates with multiple chromosomes like ATCC17802
* Add temp files for intermediate FASTQ files and cleans up BAM, SAM and hdf files from Medaka and polypolish to reduce output size

## v0.2.0 (26 October 2023)

* Replacement of POLCA by [pypolca](https://github.com/gbouras13/pypolca) as it will be easier to integrate, install and maintain going forward, and allows for POLCA use on MacOS.
* Adds `--no_medaka` to skip long read polishing with Medaka.
* Adds per contig stats output tsv file.
* Adds `--dnaapler_custom_db`.
* @simone-pignotti fixed some errors with fastp redirection and threads.
* Adds various small improvements (mem for PBS use, various other params and paths, some of the conda envs).

## v0.1.2 (6 October 2023)

* Fixes [bugs](https://github.com/gbouras13/hybracter/issues/13) with bwa index creation and typos in some output files.
* Thanks  @simone-pignotti for detecting and fixing [it](https://github.com/gbouras13/hybracter/pull/14)


## v0.1.1 (5 October 2023)

* Fixes a small [bug](https://github.com/gbouras13/hybracter/issues/9) with samples.smk
* Thanks @npbhavya and @simone-pignotti for detecting it.

## v0.1.0 (28 September 2023)

* Initial release.

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

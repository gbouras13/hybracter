# History

## v0.5.0 (08 January 2024)

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
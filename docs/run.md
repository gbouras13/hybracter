
# `hybracter` Usage

## `hybracter install`

You will first need to install the `hybracter` databases.

```bash
hybracter install
```

Alternatively, can also specify a particular directory to store them - you will need to specify this with `-d <databases directory>` when you run `hybracter`.

```bash
hybracter install -d  <databases directory>
```

## `hybracter hybrid`

You only need to specify a input CSV to run `hybracter hybrid`. It is recommended that you also specify an output directory with `-o` and a thread count with `-t`.

```bash
hybracter hybrid -i <input.csv> -o <output_dir> -t <threads>  [other arguments]
```

### Arguments

* `hybracter hybrid` requires only a CSV file specified with `-i` or `--input`
* `--no_polca` will turn off POLCA polishing with [pypolca](https://github.com/gbouras13/pypolca). 
* Use `--min_length` to specify the minimum long-read length for Filtlong.
* Use `--min_quality` to specify the minimum long-read quality for Filtlong.
* You can specify a FASTA file containing contaminants with `--contaminants`. All long reads that map to contaminants will be filtered out.
  * You can specify Escherichia phage lambda (a common contaminant in Nanopore library preparation) using `--contaminants lambda`.
* `--skip_qc` will skip all read QC steps.
* You can change the `--medakaModel` (all available options are listed in `hybracter hybrid -h`)
* You can change the `--flyeModel` (all available options are listed in `hybracter hybrid -h`)
* You can turn off Medaka polishing using `--no_medaka` - recommended for Q20+ modern Nanopore reads
* You can turn off pypolca polishing using `--no_pypolca` - I wouldn't though!
* You can change the `--depth_filter` from 0.25x chromosome coverage. This will filter out all Plassembler contigs below this depth.
* By default, `hybracter hybrid` takes the last polishing round as the final assembly  (`--logic last`). We would not recommend changing this to `--logic best`, as picking the best polishing round according to ALE with  `--logic best` is not guaranteed to give the most accurate assembly (See our [paper](https://doi.org/10.1099/mgen.0.001244)).
* You can estimate the chromosome size with kmc by using `--auto`
* You can set a minimum long-read depth with `--min_depth`. Hybracter will error out if your estimated long-reads coverage is lower than this.
* If you are running hybracter on a Mac, please use `--mac` (or find a Linux machine). This will make sure Medaka v1.8.0 is installed, as newer versions don't work on Macs.
* From v0.10.0, Hybracter will implement the `--bacteria` flag designed specifically for bacterial genomes. See Ryan Wick's [blogpost](https://rrwick.github.io/2024/10/17/medaka-v2.html) for some more explanation and benchmarking. If you do not want to use `--bactera`, please use `--medaka_override` to make sure hybracter uses your `--medakaModel`. This is likely most useful for R9 data.
* If you have all your FASTQs in a certain directory, you can use `--datadir` to specify these (and omit the directory path in the sample sheet `--input`). You can either specify 1 directory (if long and short FASTQs in the same directory) or 2 (long and short FASTQs in separate directories). If you specify 2, they must be separated by a comma e.g. `--datadir "dirlong,dirshort"`.

```bash
hybracter version 0.9.0


 _           _                    _            
| |__  _   _| |__  _ __ __ _  ___| |_ ___ _ __ 
| '_ \| | | | '_ \| '__/ _` |/ __| __/ _ \ '__|
| | | | |_| | |_) | | | (_| | (__| ||  __/ |   
|_| |_|\__, |_.__/|_|  \__,_|\___|\__\___|_|   
       |___/


Usage: hybracter hybrid [OPTIONS] [SNAKE_ARGS]...

  Run hybracter with hybrid long and paired end short reads

Options:
  -i, --input TEXT                Input csv  [required]
  --datadir TEXT                  Directory/ies where FASTQs are. Can specify
                                  1 directory (long and short FASTQs in the
                                  same directory) or 2 (long and short FASTQs
                                  in separate directories). If you specify 2,
                                  they must be separated by a comma e.g.
                                  dirlong,dirshort. Will be added to the
                                  filenames in the input csv.
  --no_pypolca                    Do not use pypolca to polish assemblies with
                                  short reads
  --logic [best|last]             Hybracter logic to select best assembly. Use
                                  --last to pick the last polishing round. Use
                                  --best to pick best assembly based on ALE
                                  (hybrid).   [default: last]
  -o, --output PATH               Output directory  [default: hybracter_out]
  --configfile TEXT               Custom config file [default: config.yaml]
  -t, --threads INTEGER           Number of threads to use  [default: 1]
  --min_length INTEGER            min read length for long reads  [default:
                                  1000]
  --min_quality INTEGER           min read quality score for long reads in bp.
                                  [default: 9]
  --skip_qc                       Do not run porechop_abi, filtlong and fastp
                                  to QC the reads
  -d, --databases PATH            Plassembler Databases directory.
  --subsample_depth INTEGER       subsampled long read depth to subsample with
                                  Filtlong. By default is 100x.  [default:
                                  100]
  --min_depth INTEGER             minimum long read depth to continue the run.
                                  By default is 0x. Hybracter will error and
                                  exit if a sample has less than
                                  min_depth*chromosome_size bases of long-
                                  reads left AFTER filtlong and porechop-ABI
                                  steps are run.  [default: 0]
  --medakaModel [r1041_e82_400bps_sup_v5.0.0|r1041_e82_400bps_hac_v5.0.0|r1041_e82_400bps_hac_v4.3.0|r1041_e82_400bps_sup_v4.3.0|r1041_e82_400bps_hac_v4.2.0|r1041_e82_400bps_sup_v4.2.0|r941_sup_plant_g610|r941_min_fast_g507|r941_prom_fast_g507|r941_min_fast_g303|r941_min_high_g303|r941_min_high_g330|r941_prom_fast_g303|r941_prom_high_g303|r941_prom_high_g330|r941_min_high_g344|r941_min_high_g351|r941_min_high_g360|r941_prom_high_g344|r941_prom_high_g360|r941_prom_high_g4011|r10_min_high_g303|r10_min_high_g340|r103_min_high_g345|r103_min_high_g360|r103_prom_high_g360|r103_fast_g507|r103_hac_g507|r103_sup_g507|r104_e81_fast_g5015|r104_e81_sup_g5015|r104_e81_hac_g5015|r104_e81_sup_g610|r1041_e82_400bps_hac_g615|r1041_e82_400bps_fast_g615|r1041_e82_400bps_fast_g632|r1041_e82_260bps_fast_g632|r1041_e82_400bps_hac_g632|r1041_e82_400bps_sup_g615|r1041_e82_260bps_hac_g632|r1041_e82_260bps_sup_g632|r1041_e82_400bps_hac_v4.0.0|r1041_e82_400bps_sup_v4.0.0|r1041_e82_260bps_hac_v4.0.0|r1041_e82_260bps_sup_v4.0.0|r1041_e82_260bps_hac_v4.1.0|r1041_e82_260bps_sup_v4.1.0|r1041_e82_400bps_hac_v4.1.0|r1041_e82_400bps_sup_v4.1.0|r941_min_high_g340_rle|r941_min_hac_g507|r941_min_sup_g507|r941_prom_hac_g507|r941_prom_sup_g507|r941_e81_fast_g514|r941_e81_hac_g514|r941_e81_sup_g514]
                                  Medaka Model.  [default:
                                  r1041_e82_400bps_sup_v5.0.0]
  --flyeModel [--nano-hq|--nano-corr|--nano-raw|--pacbio-raw|--pacbio-corr|--pacbio-hifi]
                                  Flye Assembly Parameter  [default: --nano-
                                  hq]
  --contaminants PATH             Contaminants FASTA file to map long
                                  readsagainst to filter out. Choose
                                  --contaminants lambda to filter out phage
                                  lambda long reads.
  --dnaapler_custom_db PATH       Custom amino acid FASTA file of sequences to
                                  be used as a database with dnaapler custom.
  --no_medaka                     Do not polish the long read assembly with
                                  Medaka.
  --auto                          Automatically estimate the chromosome size
                                  using KMC.
  --depth_filter FLOAT            Depth filter to pass to Plassembler. Filters
                                  out all putative plasmid contigs below this
                                  fraction of the chromosome read depth (needs
                                  to be below in both long and short read sets
                                  for hybrid).
  --mac                           If you are running Hybracter on Mac -
                                  installs v1.8.0 of Medaka as higher versions
                                  break.
  --medaka_override               Use this if you do NOT want to use the
                                  --bacteria option with Medaka. Instead your
                                  specified --medakaModel will be used.
  --use-conda / --no-use-conda    Use conda for Snakemake rules  [default:
                                  use-conda]
  --conda-prefix PATH             Custom conda env directory
  --snake-default TEXT            Customise Snakemake runtime args  [default:
                                  --rerun-incomplete, --printshellcmds,
                                  --nolock, --show-failed-logs, --conda-
                                  frontend mamba]
  -h, --help                      Show this message and exit.
```

## `hybracter hybrid-single`

You can also run a single isolate using the same input arguments as Unicycler by specifying `hybracter hybrid-single`. Instead of specifying a CSV with `--input`, use `-l` to specify the long read FASTQ file, `-1` to specify the short read R1 file, `-2` to specify the short read R2 file, `-s` to specify the sample name, `-c` to specify the chromosome size (`-c` can be omitted with `--auto`).

```bash
hybracter hybrid-single -l <longread FASTQ> -1 <R1 short reads FASTQ> -2 <R2 short reads FASTQ> -s <sample name> -c <chromosome size> -o <output_dir> -t <threads>  [other arguments]
```

The other arguments are the same as `hybracter hybrid`

```bash

Usage: hybracter hybrid-single [OPTIONS] [SNAKE_ARGS]...

  Run hybracter hybrid on 1 isolate

Options:
  -l, --longreads TEXT            FASTQ file of longreads  [required]
  -1, --short_one TEXT            R1 FASTQ file of paired end short reads
                                  [required]
  -2, --short_two TEXT            R2 FASTQ file of paired end short reads
                                  [required]
  -s, --sample TEXT               Sample name.  [default: sample]
  -c, --chromosome INTEGER        Approximate lower-bound chromosome length
                                  (in base pairs).  [default: 1000000]
  --no_pypolca                    Do not use pypolca to polish assemblies with
                                  short reads
  --logic [best|last]             Hybracter logic to select best assembly. Use
                                  --best to pick best assembly based on ALE
                                  (hybrid) or pyrodigal mean length (long).
                                  Use --last to pick the last polishing round
                                  regardless.  [default: last]
  -o, --output PATH               Output directory  [default: hybracter_out]
  --configfile TEXT               Custom config file [default: config.yaml]
  -t, --threads INTEGER           Number of threads to use  [default: 1]
  --min_length INTEGER            min read length for long reads  [default:
                                  1000]
  --min_quality INTEGER           min read quality score for long reads in bp.
                                  [default: 9]
  --skip_qc                       Do not run porechop_abi, filtlong and fastp
                                  to QC the reads
  -d, --databases PATH            Plassembler Databases directory.
  --subsample_depth INTEGER       subsampled long read depth to subsample with
                                  Filtlong. By default is 100x.  [default:
                                  100]
  --min_depth INTEGER             minimum long read depth to continue the run.
                                  By default is 0x. Hybracter will error and
                                  exit if a sample has less than
                                  min_depth*chromosome_size bases of long-
                                  reads left AFTER filtlong and porechop-ABI
                                  steps are run.  [default: 0]
  --medakaModel [r1041_e82_400bps_sup_v5.0.0|r1041_e82_400bps_hac_v5.0.0|r1041_e82_400bps_hac_v4.3.0|r1041_e82_400bps_sup_v4.3.0|r1041_e82_400bps_hac_v4.2.0|r1041_e82_400bps_sup_v4.2.0|r941_sup_plant_g610|r941_min_fast_g507|r941_prom_fast_g507|r941_min_fast_g303|r941_min_high_g303|r941_min_high_g330|r941_prom_fast_g303|r941_prom_high_g303|r941_prom_high_g330|r941_min_high_g344|r941_min_high_g351|r941_min_high_g360|r941_prom_high_g344|r941_prom_high_g360|r941_prom_high_g4011|r10_min_high_g303|r10_min_high_g340|r103_min_high_g345|r103_min_high_g360|r103_prom_high_g360|r103_fast_g507|r103_hac_g507|r103_sup_g507|r104_e81_fast_g5015|r104_e81_sup_g5015|r104_e81_hac_g5015|r104_e81_sup_g610|r1041_e82_400bps_hac_g615|r1041_e82_400bps_fast_g615|r1041_e82_400bps_fast_g632|r1041_e82_260bps_fast_g632|r1041_e82_400bps_hac_g632|r1041_e82_400bps_sup_g615|r1041_e82_260bps_hac_g632|r1041_e82_260bps_sup_g632|r1041_e82_400bps_hac_v4.0.0|r1041_e82_400bps_sup_v4.0.0|r1041_e82_260bps_hac_v4.0.0|r1041_e82_260bps_sup_v4.0.0|r1041_e82_260bps_hac_v4.1.0|r1041_e82_260bps_sup_v4.1.0|r1041_e82_400bps_hac_v4.1.0|r1041_e82_400bps_sup_v4.1.0|r941_min_high_g340_rle|r941_min_hac_g507|r941_min_sup_g507|r941_prom_hac_g507|r941_prom_sup_g507|r941_e81_fast_g514|r941_e81_hac_g514|r941_e81_sup_g514]
                                  Medaka Model.  [default:
                                  r1041_e82_400bps_sup_v5.0.0]
  --flyeModel [--nano-hq|--nano-corr|--nano-raw|--pacbio-raw|--pacbio-corr|--pacbio-hifi]
                                  Flye Assembly Parameter  [default: --nano-
                                  hq]
  --contaminants PATH             Contaminants FASTA file to map long
                                  readsagainst to filter out. Choose
                                  --contaminants lambda to filter out phage
                                  lambda long reads.
  --dnaapler_custom_db PATH       Custom amino acid FASTA file of sequences to
                                  be used as a database with dnaapler custom.
  --no_medaka                     Do not polish the long read assembly with
                                  Medaka.
  --auto                          Automatically estimate the chromosome size
                                  using KMC.
  --depth_filter FLOAT            Depth filter to pass to Plassembler. Filters
                                  out all putative plasmid contigs below this
                                  fraction of the chromosome read depth (needs
                                  to be below in both long and short read sets
                                  for hybrid).
  --mac                           If you are running Hybracter on Mac -
                                  installs v1.8.0 of Medaka as higher versions
                                  break.
  --medaka_override               Use this if you do NOT want to use the
                                  --bacteria option with Medaka. Instead your
                                  specified --medakaModel will be used.
  --use-conda / --no-use-conda    Use conda for Snakemake rules  [default:
                                  use-conda]
  --conda-prefix PATH             Custom conda env directory
  --snake-default TEXT            Customise Snakemake runtime args  [default:
                                  --rerun-incomplete, --printshellcmds,
                                  --nolock, --show-failed-logs, --conda-
                                  frontend mamba]
  -h, --help                      Show this message and exit.

```

## `hybracter long`

You only need to specify a input CSV to run  `hybracter long`. It is recommended that you also specify an output directory with `-o` and a thread count with `-t`.

```bash
hybracter long -i <input.csv> -o <output_dir> -t <threads> [other arguments]
```

### Arguments

* `hybracter long` requires only a CSV file specified with `-i` or `--input`
* Use `--min_length` to specify the minimum long-read length for Filtlong.
* Use `--min_quality` to specify the minimum long-read quality for Filtlong.
* You can specify a FASTA file containing contaminants with `--contaminants`. All long reads that map to contaminants will be filtered.
  * You can specify Escherichia phage lambda (a common contaminant in Nanopore library preparation) using `--contaminants lambda`.
* `--skip_qc` will skip all read QC steps.
* You can change the `--medakaModel` (all available options are listed in `hybracter long -h`)
* You can change the `--flyeModel` (all available options are listed in `hybracter long -h`)
* You can turn off Medaka polishing using `--no_medaka` - recommended for Q20+ modern Nanopore and PacBio reads
* You can change the `--depth_filter` from 0.25x chromosome coverage. This will filter out all Plassembler contigs below this depth.
* You can force `hybracter long` to pick the last polishing round (not the best according to pyrodigal mean CDS length) with `--logic last`. `hybracter long` defaults to picking the best i.e. `--logic best`.
* You can estimate the chromosome size with kmc by using `--auto`
* You can set a minimum long-read depth with `--min_depth`. Hybracter will error out if your estimated long-reads coverage is lower than this.
* If you are running hybracter on a Mac, please use `--mac` (or find a Linux machine). This will make sure Medaka v1.8.0 is installed, as newer versions don't work on Macs.
* From v0.10.0, Hybracter will implement the `--bacteria` flag designed specifically for bacterial genomes. See Ryan Wick's [blogpost](https://rrwick.github.io/2024/10/17/medaka-v2.html) for some more explanation and benchmarking. If you do not want to use `--bactera`, please use `--medaka_override` to make sure hybracter uses your `--medakaModel`. This is likely most useful for R9 data.
* If you have all your FASTQs in a certain directory, you can use `--datadir` to specify these (and omit the directory path in the sample sheet `--input`).

```bash
Usage: hybracter long [OPTIONS] [SNAKE_ARGS]...

  Run hybracter with only long reads
Options:
  -i, --input TEXT                Input csv  [required]
  --datadir TEXT                  Directory where FASTQs are. Will be added 
                                  to the filenames in the input csv.
  -o, --output PATH               Output directory  [default: hybracter_out]
  --configfile TEXT               Custom config file [default: config.yaml]
  -t, --threads INTEGER           Number of threads to use  [default: 1]
  --min_length INTEGER            min read length for long reads  [default:
                                  1000]
  --min_quality INTEGER           min read quality score for long reads in bp.
                                  [default: 9]
  --skip_qc                       Do not run porechop_abi, filtlong and fastp
                                  to QC the reads
  -d, --databases PATH            Plassembler Databases directory.
  --subsample_depth INTEGER       subsampled long read depth to subsample with
                                  Filtlong. By default is 100x.  [default:
                                  100]
  --min_depth INTEGER             minimum long read depth to continue the run.
                                  By default is 0x. Hybracter will error and
                                  exit if a sample has less than
                                  min_depth*chromosome_size bases of long-
                                  reads left AFTER filtlong and porechop-ABI
                                  steps are run.  [default: 0]
  --medakaModel [r1041_e82_400bps_sup_v5.0.0|r1041_e82_400bps_hac_v5.0.0|r1041_e82_400bps_hac_v4.3.0|r1041_e82_400bps_sup_v4.3.0|r1041_e82_400bps_hac_v4.2.0|r1041_e82_400bps_sup_v4.2.0|r941_sup_plant_g610|r941_min_fast_g507|r941_prom_fast_g507|r941_min_fast_g303|r941_min_high_g303|r941_min_high_g330|r941_prom_fast_g303|r941_prom_high_g303|r941_prom_high_g330|r941_min_high_g344|r941_min_high_g351|r941_min_high_g360|r941_prom_high_g344|r941_prom_high_g360|r941_prom_high_g4011|r10_min_high_g303|r10_min_high_g340|r103_min_high_g345|r103_min_high_g360|r103_prom_high_g360|r103_fast_g507|r103_hac_g507|r103_sup_g507|r104_e81_fast_g5015|r104_e81_sup_g5015|r104_e81_hac_g5015|r104_e81_sup_g610|r1041_e82_400bps_hac_g615|r1041_e82_400bps_fast_g615|r1041_e82_400bps_fast_g632|r1041_e82_260bps_fast_g632|r1041_e82_400bps_hac_g632|r1041_e82_400bps_sup_g615|r1041_e82_260bps_hac_g632|r1041_e82_260bps_sup_g632|r1041_e82_400bps_hac_v4.0.0|r1041_e82_400bps_sup_v4.0.0|r1041_e82_260bps_hac_v4.0.0|r1041_e82_260bps_sup_v4.0.0|r1041_e82_260bps_hac_v4.1.0|r1041_e82_260bps_sup_v4.1.0|r1041_e82_400bps_hac_v4.1.0|r1041_e82_400bps_sup_v4.1.0|r941_min_high_g340_rle|r941_min_hac_g507|r941_min_sup_g507|r941_prom_hac_g507|r941_prom_sup_g507|r941_e81_fast_g514|r941_e81_hac_g514|r941_e81_sup_g514]
                                  Medaka Model.  [default:
                                  r1041_e82_400bps_sup_v5.0.0]
  --flyeModel [--nano-hq|--nano-corr|--nano-raw|--pacbio-raw|--pacbio-corr|--pacbio-hifi]
                                  Flye Assembly Parameter  [default: --nano-
                                  hq]
  --contaminants PATH             Contaminants FASTA file to map long
                                  readsagainst to filter out. Choose
                                  --contaminants lambda to filter out phage
                                  lambda long reads.
  --dnaapler_custom_db PATH       Custom amino acid FASTA file of sequences to
                                  be used as a database with dnaapler custom.
  --no_medaka                     Do not polish the long read assembly with
                                  Medaka.
  --auto                          Automatically estimate the chromosome size
                                  using KMC.
  --depth_filter FLOAT            Depth filter to pass to Plassembler. Filters
                                  out all putative plasmid contigs below this
                                  fraction of the chromosome read depth (needs
                                  to be below in both long and short read sets
                                  for hybrid).
  --mac                           If you are running Hybracter on Mac -
                                  installs v1.8.0 of Medaka as higher versions
                                  break.
  --medaka_override               Use this if you do NOT want to use the
                                  --bacteria option with Medaka. Instead your
                                  specified --medakaModel will be used.
  --use-conda / --no-use-conda    Use conda for Snakemake rules  [default:
                                  use-conda]
  --conda-prefix PATH             Custom conda env directory
  --snake-default TEXT            Customise Snakemake runtime args  [default:
                                  --rerun-incomplete, --printshellcmds,
                                  --nolock, --show-failed-logs, --conda-
                                  frontend mamba]
  --logic [best|last]             Hybracter logic to select best assembly. Use
                                  --best to pick best assembly based on ALE
                                  (hybrid) or pyrodigal mean length (long).
                                  Use --last to pick the last polishing round
                                  regardless.  [default: best]
  -h, --help                      Show this message and exit.

```

## `hybracter long-single`


Run `hybracter long` on a single isolate. Instead of specifying a CSV with `--input`, use `-l` to specify the long read FASTQ file, `-s` to specify the sample name, `-c` to specify the chromosome size (unless you are using `--auto`, then you don't need to specify `-c`).

```bash
hybracter long-single -l <longread FASTQ> -s <sample name> -c <chromosome size>  -o <output_dir> -t <threads>  [other arguments]
```

```bash
Usage: hybracter long-single [OPTIONS] [SNAKE_ARGS]...

  Run hybracter long on 1 isolate
Options:
  -l, --longreads TEXT            FASTQ file of longreads  [required]
  -s, --sample TEXT               Sample name.  [default: sample]
  -c, --chromosome INTEGER        FApproximate lower-bound chromosome length
                                  (in base pairs).  [default: 1000000]
  -o, --output PATH               Output directory  [default: hybracter_out]
  --configfile TEXT               Custom config file [default: config.yaml]
  -t, --threads INTEGER           Number of threads to use  [default: 1]
  --min_length INTEGER            min read length for long reads  [default:
                                  1000]
  --min_quality INTEGER           min read quality score for long reads in bp.
                                  [default: 9]
  --skip_qc                       Do not run porechop_abi, filtlong and fastp
                                  to QC the reads
  -d, --databases PATH            Plassembler Databases directory.
  --subsample_depth INTEGER       subsampled long read depth to subsample with
                                  Filtlong. By default is 100x.  [default:
                                  100]
  --min_depth INTEGER             minimum long read depth to continue the run.
                                  By default is 0x. Hybracter will error and
                                  exit if a sample has less than
                                  min_depth*chromosome_size bases of long-
                                  reads left AFTER filtlong and porechop-ABI
                                  steps are run.  [default: 0]
  --medakaModel [r1041_e82_400bps_sup_v5.0.0|r1041_e82_400bps_hac_v5.0.0|r1041_e82_400bps_hac_v4.3.0|r1041_e82_400bps_sup_v4.3.0|r1041_e82_400bps_hac_v4.2.0|r1041_e82_400bps_sup_v4.2.0|r941_sup_plant_g610|r941_min_fast_g507|r941_prom_fast_g507|r941_min_fast_g303|r941_min_high_g303|r941_min_high_g330|r941_prom_fast_g303|r941_prom_high_g303|r941_prom_high_g330|r941_min_high_g344|r941_min_high_g351|r941_min_high_g360|r941_prom_high_g344|r941_prom_high_g360|r941_prom_high_g4011|r10_min_high_g303|r10_min_high_g340|r103_min_high_g345|r103_min_high_g360|r103_prom_high_g360|r103_fast_g507|r103_hac_g507|r103_sup_g507|r104_e81_fast_g5015|r104_e81_sup_g5015|r104_e81_hac_g5015|r104_e81_sup_g610|r1041_e82_400bps_hac_g615|r1041_e82_400bps_fast_g615|r1041_e82_400bps_fast_g632|r1041_e82_260bps_fast_g632|r1041_e82_400bps_hac_g632|r1041_e82_400bps_sup_g615|r1041_e82_260bps_hac_g632|r1041_e82_260bps_sup_g632|r1041_e82_400bps_hac_v4.0.0|r1041_e82_400bps_sup_v4.0.0|r1041_e82_260bps_hac_v4.0.0|r1041_e82_260bps_sup_v4.0.0|r1041_e82_260bps_hac_v4.1.0|r1041_e82_260bps_sup_v4.1.0|r1041_e82_400bps_hac_v4.1.0|r1041_e82_400bps_sup_v4.1.0|r941_min_high_g340_rle|r941_min_hac_g507|r941_min_sup_g507|r941_prom_hac_g507|r941_prom_sup_g507|r941_e81_fast_g514|r941_e81_hac_g514|r941_e81_sup_g514]
                                  Medaka Model.  [default:
                                  r1041_e82_400bps_sup_v5.0.0]
  --flyeModel [--nano-hq|--nano-corr|--nano-raw|--pacbio-raw|--pacbio-corr|--pacbio-hifi]
                                  Flye Assembly Parameter  [default: --nano-
                                  hq]
  --contaminants PATH             Contaminants FASTA file to map long
                                  readsagainst to filter out. Choose
                                  --contaminants lambda to filter out phage
                                  lambda long reads.
  --dnaapler_custom_db PATH       Custom amino acid FASTA file of sequences to
                                  be used as a database with dnaapler custom.
  --no_medaka                     Do not polish the long read assembly with
                                  Medaka.
  --auto                          Automatically estimate the chromosome size
                                  using KMC.
  --depth_filter FLOAT            Depth filter to pass to Plassembler. Filters
                                  out all putative plasmid contigs below this
                                  fraction of the chromosome read depth (needs
                                  to be below in both long and short read sets
                                  for hybrid).
  --mac                           If you are running Hybracter on Mac -
                                  installs v1.8.0 of Medaka as higher versions
                                  break.
  --medaka_override               Use this if you do NOT want to use the
                                  --bacteria option with Medaka. Instead your
                                  specified --medakaModel will be used.
  --use-conda / --no-use-conda    Use conda for Snakemake rules  [default:
                                  use-conda]
  --conda-prefix PATH             Custom conda env directory
  --snake-default TEXT            Customise Snakemake runtime args  [default:
                                  --rerun-incomplete, --printshellcmds,
                                  --nolock, --show-failed-logs, --conda-
                                  frontend mamba]
  --logic [best|last]             Hybracter logic to select best assembly. Use
                                  --best to pick best assembly based on ALE
                                  (hybrid) or pyrodigal mean length (long).
                                  Use --last to pick the last polishing round
                                  regardless.  [default: best]
  -h, --help                      Show this message and exit.
```


# `hybracter` Usage

## `hybracter install`

You will first need to install the `hybracter` databases.

```
hybracter install
```

Alternatively, can also specify a particular directory to store them - you will need to specify this with `-d <databases directory>` when you run `hybracter`.

```
hybracter install -d  <databases directory>
```

## `hybracter hybrid`

You only need to specify a input CSV to run `hybracter hybrid`. It is recommended that you also specify an output directory with `-o` and a thread count with `-t`.

```
hybracter hybrid -i <input.csv> -o <output_dir> -t <threads>  [other arguments]
```

### Arguments

* `hybracter hybrid` requires only a CSV file specified with `-i` or `--input`
* `--no_polca` will turn off POLCA polishing. **Must be used on MacOS.**
* Use `--min_length` to specify the minimum long-read length for Filtlong.
* Use `--min_quality` to specify the minimum long-read quality for Filtlong.
* You can specify a FASTA file containing contaminants with `--contaminants`. All long reads that map to contaminants will be filtered out.
  * You can specify Escherichia phage lambda (a common contaminant in Nanopore library preparation) using `--contaminants lambda`.
* `--skip_qc` will skip all read QC steps.
* You can change the `--medakaModel` (all available options are listed in `hybracter hybrid -h`)
* You can change the `--flyeModel` (all available options are listed in `hybracter hybrid -h`)

```
Usage: hybracter [OPTIONS] COMMAND [ARGS]...

  For more options, run: hybracter command --help

Options:
  -h, --help  Show this message and exit.

Commands:
  install        Downloads and installs the plassembler database
  hybrid         Run hybracter with hybrid long and paired end short reads
  hybrid-single  Run hybracter hybrid on 1 isolate
  long           Run hybracter long
  long-single    Run hybracter long on 1 isolate
  test-hybrid    Test hybracter hybrid
  test-long      Test hybracter long
  config         Copy the system default config file
  citation       Print the citation(s) for hybracter
  version        Print the version for hybracter
(hybracter) a1667917@LY0TWV6HTW2 hybracter % hybracter hybrid -h

hybracter version 0.1.0


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
  --no_polca                      Do not use Polca to polish assemblies with
                                  short reads
  -o, --output PATH               Output directory  [default: hybracter.out]
  --configfile TEXT               Custom config file [default:
                                  (outputDir)/config.yaml]
  -t, --threads INTEGER           Number of threads to use  [default: 1]
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
  --contaminants PATH             Contaminants FASTA file to map long
                                  readsagainst to filter out. Choose
                                  --contaminants lambda to filter out phage
                                  lambda long reads.
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

You can also run a single isolate using the same input arguments as Unicycler by specifying `hybracter hybrid-single`. Instead of specifying a CSV with `--input`, use `-l` to specify the long read FASTQ file, `-1` to specify the short read R1 file, `-2` to specify the short read R2 file, `-s` to specify the sample name, `-c` to specify the chromosome size.

```
hybracter hybrid-single -l <longread FASTQ> -1 <R1 short reads FASTQ> -2 <R2 short reads FASTQ> -s <sample name> -c <chromosome size> -o <output_dir> -t <threads>  [other arguments]
```

The other arguments are the same as `hybracter hybrid`

```
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
  --no_polca                      Do not use Polca to polish assemblies with
                                  short reads
  -o, --output PATH               Output directory  [default: hybracter.out]
  --configfile TEXT               Custom config file [default:
                                  (outputDir)/config.yaml]
  -t, --threads INTEGER           Number of threads to use  [default: 1]
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
  --contaminants PATH             Contaminants FASTA file to map long
                                  readsagainst to filter out. Choose
                                  --contaminants lambda to filter out phage
                                  lambda long reads.
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

```
hybracter long -i <input.csv> -o <output_dir> -t <threads> [other arguments]
```

### Other Arguments

* `hybracter long` requires only a CSV file specified with `-i` or `--input`
* Use `--min_length` to specify the minimum long-read length for Filtlong.
* Use `--min_quality` to specify the minimum long-read quality for Filtlong.
* You can specify a FASTA file containing contaminants with `--contaminants`. All long reads that map to contaminants will be filtered.
  * You can specify Escherichia phage lambda (a common contaminant in Nanopore library preparation) using `--contaminants lambda`.
* `--skip_qc` will skip all read QC steps.
* You can change the `--medakaModel` (all available options are listed in `hybracter long -h`)
* You can change the `--flyeModel` (all available options are listed in `hybracter long -h`)

```
Usage: hybracter long [OPTIONS] [SNAKE_ARGS]...

  Run hybracter long

Options:
  -i, --input TEXT                Input csv  [required]
  -o, --output PATH               Output directory  [default: hybracter.out]
  --configfile TEXT               Custom config file [default:
                                  (outputDir)/config.yaml]
  -t, --threads INTEGER           Number of threads to use  [default: 1]
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
  --contaminants PATH             Contaminants FASTA file to map long
                                  readsagainst to filter out. Choose
                                  --contaminants lambda to filter out phage
                                  lambda long reads.
  --use-conda / --no-use-conda    Use conda for Snakemake rules  [default:
                                  use-conda]
  --conda-prefix PATH             Custom conda env directory
  --snake-default TEXT            Customise Snakemake runtime args  [default:
                                  --rerun-incomplete, --printshellcmds,
                                  --nolock, --show-failed-logs, --conda-
                                  frontend mamba]
  -h, --help                      Show this message and exit.
```

## `hybracter long-single`


Run `hybracter long` on a single isolate. Instead of specifying a CSV with `--input`, use `-l` to specify the long read FASTQ file, `-s` to specify the sample name, `-c` to specify the chromosome size.

```
hybracter long-single -l <longread FASTQ> -s <sample name> -c <chromosome size>  -o <output_dir> -t <threads>  [other arguments]
```

```
Usage: hybracter long-single [OPTIONS] [SNAKE_ARGS]...

  Run hybracter long on 1 isolate

Options:
  -l, --longreads TEXT            FASTQ file of longreads  [required]
  -s, --sample TEXT               Sample name.  [default: sample]
  -c, --chromosome INTEGER        FApproximate lower-bound chromosome length
                                  (in base pairs).  [default: 1000000]
  -o, --output PATH               Output directory  [default: hybracter.out]
  --configfile TEXT               Custom config file [default:
                                  (outputDir)/config.yaml]
  -t, --threads INTEGER           Number of threads to use  [default: 1]
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
  --contaminants PATH             Contaminants FASTA file to map long
                                  readsagainst to filter out. Choose
                                  --contaminants lambda to filter out phage
                                  lambda long reads.
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

I would highly recommend running hybracter using a Snakemake profile. Please see this blog [post](https://fame.flinders.edu.au/blog/2021/08/02/snakemake-profiles-updated) for more details. I have included an example slurm profile in the `profile` [directory](https://github.com/gbouras13/hybracter/tree/main/profiles/hybracter), but check out this [link](https://github.com/Snakemake-Profiles) for more detail on other HPC job scheduler profiles.

You can run `hybracter` with a profile using `--profile` e.g.:
```
hybracter hybrid --input <input.csv> --output <output_dir> --threads <threads> --profile profiles/hybracter
```
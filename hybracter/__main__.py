"""
Entrypoint for hybracter

Check out the wiki for a detailed look at customising this file:
https://github.com/beardymcjohnface/Snaketool/wiki/Customising-your-Snaketool
"""

import os

import click

from .util import (
    OrderedCommands,
    all_medaka_models,
    copy_config,
    default_to_ouput,
    print_citation,
    print_version,
    run_snakemake,
    snake_base,
)


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
        click.option(
            "-o",
            "--output",
            help="Output directory",
            type=click.Path(dir_okay=True, writable=True, readable=True),
            default="hybracter_out",
            show_default=True,
        ),
        click.option(
            "--configfile",
            default="config.yaml",
            type=str,
            show_default=False,
            callback=default_to_ouput,
            help="Custom config file [default: config.yaml]",
        ),
        click.option(
            "-t",
            "--threads",
            help="Number of threads to use",
            default=1,
            show_default=True,
        ),
        click.option(
            "--min_length",
            "min_length",
            help="min read length for long reads",
            type=int,
            default=1000,
            show_default=True,
        ),
        click.option(
            "--min_quality",
            "min_quality",
            help="min read quality score for long reads in bp.",
            type=int,
            default=9,
            show_default=True,
        ),
        click.option(
            "--skip_qc",
            help="Do not run porechop_abi, filtlong and fastp to QC the reads",
            is_flag=True,
            default=False,
        ),
        click.option(
            "-d",
            "--databases",
            help="Plassembler Databases directory.",
            type=click.Path(dir_okay=True, readable=True),
        ),
        click.option(
            "--subsample_depth",
            "subsample_depth",
            help="subsampled long read depth to subsample with Filtlong. By default is 100x.",
            type=int,
            default=100,
            show_default=True,
        ),
        click.option(
            "--min_depth",
            "min_depth",
            help="minimum long read depth to continue the run. By default is 0x.\n Hybracter will error and exit if a sample has less than min_depth*chromosome_size bases of long-reads left AFTER filtlong and porechop-ABI steps are run.",
            type=int,
            default=0,
            show_default=True,
        ),
        click.option(
            "--medakaModel",
            "medakaModel",
            help="Medaka Model.",
            default="r1041_e82_400bps_sup_v5.0.0",
            show_default=True,
            type=click.Choice(all_medaka_models),
        ),
        click.option(
            "--flyeModel",
            "flyeModel",
            help="Flye Assembly Parameter",
            show_default=True,
            default="--nano-hq",
            type=click.Choice(
                [
                    "--nano-hq",
                    "--nano-corr",
                    "--nano-raw",
                    "--pacbio-raw",
                    "--pacbio-corr",
                    "--pacbio-hifi",
                ]
            ),
        ),
        click.option(
            "--contaminants",
            help="Contaminants FASTA file to map long readsagainst to filter out. Choose --contaminants lambda to filter out phage lambda long reads.",
            type=click.Path(),
            default="none",
            required=False,
        ),
        click.option(
            "--dnaapler_custom_db",
            help="Custom amino acid FASTA file of sequences to be used as a database with dnaapler custom.",
            type=click.Path(),
            default="none",
            required=False,
        ),
        click.option(
            "--no_medaka",
            help="Do not polish the long read assembly with Medaka.",
            is_flag=True,
            default=False,
        ),
        click.option(
            "--auto",
            help="Automatically estimate the chromosome size using KMC.",
            is_flag=True,
            default=False,
        ),
        click.option(
            "--depth_filter",
            help="Depth filter to pass to Plassembler. Filters out all putative plasmid contigs below this fraction of the chromosome read depth (needs to be below in both long and short read sets for hybrid).",
            type=float,
            default=0.25,
        ),
        click.option(
            "--mac",
            help="If you are running Hybracter on Mac - installs v1.8.0 of Medaka as higher versions break.",
            is_flag=True,
            default=False,
        ),

        click.option(
            "--use-conda/--no-use-conda",
            default=True,
            help="Use conda for Snakemake rules",
            show_default=True,
        ),
        click.option(
            "--conda-prefix",
            default=snake_base(os.path.join("workflow", "conda")),
            help="Custom conda env directory",
            type=click.Path(),
            show_default=False,
        ),
        click.option(
            "--snake-default",
            multiple=True,
            default=[
                "--rerun-incomplete",
                "--printshellcmds",
                "--nolock",
                "--show-failed-logs",
                "--conda-frontend conda",
            ],
            help="Customise Snakemake runtime args",
            show_default=True,
        ),
        click.option(
            "--log",
            default="hybracter.log",
            callback=default_to_ouput,
            type=str,
            show_default=False,
            hidden=True,
        ),
        click.argument("snake_args", nargs=-1),
    ]
    for option in reversed(options):
        func = option(func)
    return func


@click.group(
    cls=OrderedCommands, context_settings=dict(help_option_names=["-h", "--help"])
)
def cli():
    """For more options, run:
    hybracter command --help"""
    pass


def print_splash():
    click.echo(
        """\b

 _           _                    _            
| |__  _   _| |__  _ __ __ _  ___| |_ ___ _ __ 
| '_ \| | | | '_ \| '__/ _` |/ __| __/ _ \ '__|
| | | | |_| | |_) | | | (_| | (__| ||  __/ |   
|_| |_|\__, |_.__/|_|  \__,_|\___|\__\___|_|   
       |___/

"""
    )


help_hybrid_msg_extra = """
\b
CLUSTER EXECUTION:
hybracter hybrid ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           hybracter hybrid --input [file]
Specify output directory:    hybracter hybrid ... --output [directory]
Specify threads:    hybracter hybrid ... --threads [threads]
Disable conda:      hybracter hybrid ... --no-use-conda 
Change defaults:    hybracter hybrid ... --snake-default="-k --nolock"
Add Snakemake args: hybracter hybrid ... --dry-run --keep-going --touch
Specify targets:    hybracter hybrid ... all print_targets
Available targets:
    all             Run everything (default)
    print_targets   List available targets
"""

help_hybrid_single_msg_extra = """
\b
CLUSTER EXECUTION:
hybracter hybrid-single ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           hybracter hybrid-single -l [FASTQ file of longreads]
Required:           hybracter hybrid-single -1 [R1 FASTQ file of paired end short reads]
Required:           hybracter hybrid-single -2 [R2 FASTQ file of paired end short reads]
Specify output directory:    hybracter hybrid-single  ... --output [directory]
Specify threads:    hybracter hybrid-single  ... --threads [threads]
Disable conda:      hybracter hybrid-single  ... --no-use-conda 
Change defaults:    hybracter hybrid-single  ... --snake-default="-k --nolock"
Add Snakemake args: hybracter hybrid-single  ... --dry-run --keep-going --touch
Specify targets:    hybracter hybrid-single  ... all print_targets
Available targets:
    all             Run everything (default)
    print_targets   List available targets
"""

help_long_msg_extra = """
\b
CLUSTER EXECUTION:
hybracter hybrid ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           hybracter long --input [file]
Specify output directory:    hybracter long ... --output [directory]
Specify threads:    hybracter long ... --threads [threads]
Disable conda:      hybracter long ... --no-use-conda 
Change defaults:    hybracter long ... --snake-default="-k --nolock"
Add Snakemake args: hybracter long ... --dry-run --keep-going --touch
Specify targets:    hybracter long ... all print_targets
Available targets:
    all             Run everything (default)
    print_targets   List available targets
"""

help_long_single_msg_extra = """
\b
CLUSTER EXECUTION:
hybracter long-single ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           hybracter long-single -l [FASTQ file of longreads]
Specify output directory:    hybracter long-single  ... --output [directory]
Specify threads:    hybracter long-single  ... --threads [threads]
Disable conda:      hybracter long-single  ... --no-use-conda 
Change defaults:    hybracter long-single  ... --snake-default="-k --nolock"
Add Snakemake args: hybracter long-single  ... --dry-run --keep-going --touch
Specify targets:    hybracter long-single  ... all print_targets
Available targets:
    all             Run everything (default)
    print_targets   List available targets
"""

help_msg_install = """
\b
installs the plassembler database
hybracter install ... 
\b
RUN EXAMPLES:
Database:           hybracter install -d [directory]
"""


help_msg_test_hybrid = """
\b
hybracter test-hybrid  ...
\b
RUN EXAMPLES:
Specify output directory:    hybracter test-hybrid  ... --output [directory]
Specify threads:    hybracter test-hybrid  ... --threads [threads]
Disable conda:      hybracter test-hybrid  ... --no-use-conda 
Change defaults:    hybracter test-hybrid  ... --snake-default="-k --nolock"
Add Snakemake args: hybracter test-hybrid  ... --dry-run --keep-going --touch
Specify targets:    hybracter test-hybrid  ... all print_targets
Available targets:
    all             Run everything (default)
    print_targets   List available targets
"""

help_msg_test_long = """
\b
hybracter test-long ...
\b
RUN EXAMPLES:
Specify output directory:    hybracter test-long  ... --output [directory]
Specify threads:    hybracter test-long  ... --threads [threads]
Disable conda:      hybracter test-long  ... --no-use-conda 
Change defaults:    hybracter test-long  ... --snake-default="-k --nolock"
Add Snakemake args: hybracter test-long  ... --dry-run --keep-going --touch
Specify targets:    hybracter test-long  ... all print_targets
Available targets:
    all             Run everything (default)
    print_targets   List available targets
"""

"""
hybrid
"""


@click.command(
    epilog=help_hybrid_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option("-i", "--input", "_input", help="Input csv", type=str, required=True)
@click.option("--datadir", "datadir", help="Directory/ies where FASTQs are. Can specify 1 directory (long and short FASTQs in the same directory) or 2 (long and short FASTQs in separate directories).\n If you specify 2, they must be separated by a comma e.g. dirlong,dirshort.\n Will be added to the filenames in the input csv.", type=str, default=None)
@click.option(
    "--no_pypolca",
    help="Do not use pypolca to polish assemblies with short reads",
    is_flag=True,
    default=False,
)
@click.option(
            "--logic",
            "logic",
            help="Hybracter logic to select best assembly. Use --last to pick the last polishing round. Use --best to pick best assembly based on ALE (hybrid). ",
            show_default=True,
            default="last",
            type=click.Choice(
                [
                    "best",
                    "last",
                ]
            ),
        )
@common_options
def hybrid(
    _input,
    datadir,
    no_medaka,
    auto,
    no_pypolca,
    skip_qc,
    medakaModel,
    databases,
    subsample_depth,
    min_depth,
    min_quality,
    flyeModel,
    min_length,
    output,
    contaminants,
    dnaapler_custom_db,
    logic,
    depth_filter,
    mac,
    log,
    configfile,
    **kwargs
):
    """Run hybracter with hybrid long and paired end short reads"""
    # Config to add or update in configfile
    merge_config = {
        "args": {
            "input": _input,
            "datadir": datadir,
            "output": output,
            "log": log,
            "min_length": min_length,
            "databases": databases,
            "subsample_depth": subsample_depth,
            "min_depth":min_depth,
            "min_quality": min_quality,
            "no_pypolca": no_pypolca,
            "skip_qc": skip_qc,
            "medakaModel": medakaModel,
            "flyeModel": flyeModel,
            "contaminants": contaminants,
            "dnaapler_custom_db": dnaapler_custom_db,
            "no_medaka": no_medaka,
            "auto": auto,
            "mac": mac,
            "depth_filter": depth_filter,
            "single": False,
            "logic": logic,
            "configfile": configfile
        }
    }

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "hybrid.smk")),
        configfile=configfile,
        merge_config=merge_config,
        log=log,
        **kwargs
    )


"""
hybrid single
"""


@click.command(
    epilog=help_hybrid_single_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option(
    "-l", "--longreads", help="FASTQ file of longreads", type=str, required=True
)
@click.option(
    "-1",
    "--short_one",
    help="R1 FASTQ file of paired end short reads",
    type=str,
    required=True,
)
@click.option(
    "-2",
    "--short_two",
    help="R2 FASTQ file of paired end short reads",
    type=str,
    required=True,
)
@click.option(
    "-s", "--sample", help="Sample name.", type=str, default="sample", show_default=True
)
@click.option(
    "-c",
    "--chromosome",
    help="Approximate lower-bound chromosome length (in base pairs).",
    type=int,
    default=1000000,
    show_default=True,
)
@click.option(
    "--no_pypolca",
    help="Do not use pypolca to polish assemblies with short reads",
    is_flag=True,
    default=False,
)
@click.option(
            "--logic",
            "logic",
            help="Hybracter logic to select best assembly. Use --best to pick best assembly based on ALE (hybrid) or pyrodigal mean length (long). Use --last to pick the last polishing round regardless.",
            show_default=True,
            default="last",
            type=click.Choice(
                [
                    "best",
                    "last",
                ]
            ),
        )
@common_options
def hybrid_single(
    longreads,
    chromosome,
    sample,
    short_one,
    short_two,
    no_pypolca,
    skip_qc,
    medakaModel,
    databases,
    subsample_depth,
    min_depth,
    min_quality,
    flyeModel,
    min_length,
    output,
    contaminants,
    dnaapler_custom_db,
    no_medaka,
    auto,
    mac,
    logic,
    depth_filter,
    log,
    configfile,
    **kwargs
):
    """Run hybracter hybrid on 1 isolate"""
    # Config to add or update in configfile
    merge_config = {
        "args": {
            "longreads": longreads,
            "short_one": short_one,
            "short_two": short_two,
            "sample": sample,
            "chromosome": chromosome,
            "output": output,
            "log": log,
            "min_length": min_length,
            "databases": databases,
            "subsample_depth": subsample_depth,
            "min_depth":min_depth,
            "min_quality": min_quality,
            "no_pypolca": no_pypolca,
            "skip_qc": skip_qc,
            "medakaModel": medakaModel,
            "flyeModel": flyeModel,
            "contaminants": contaminants,
            "dnaapler_custom_db": dnaapler_custom_db,
            "no_medaka": no_medaka,
            "auto": auto,
            "mac": mac,
            "depth_filter": depth_filter,
            "single": True,
            "logic": logic,
            "configfile": configfile
        }
    }

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "hybrid.smk")),
        configfile=configfile,
        merge_config=merge_config,
        log=log,
        **kwargs
    )


"""
long
"""


@click.command(
    epilog=help_long_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option("-i", "--input", "_input", help="Input csv", type=str, required=True)
@click.option("--datadir", "datadir", help="Directory where FASTQs are. Will be added to the filenames in the input csv.", type=str, default=None)
@common_options
@click.option(
            "--logic",
            "logic",
            help="Hybracter logic to select best assembly. Use --best to pick best assembly based on ALE (hybrid) or pyrodigal mean length (long). Use --last to pick the last polishing round regardless.",
            show_default=True,
            default="best",
            type=click.Choice(
                [
                    "best",
                    "last",
                ]
            ),
        )
def long(
    _input,
    datadir,
    medakaModel,
    databases,
    subsample_depth,
    min_depth,
    skip_qc,
    flyeModel,
    min_length,
    output,
    min_quality,
    contaminants,
    dnaapler_custom_db,
    no_medaka,
    auto,
    mac,
    logic,
    depth_filter,
    log,
    configfile,
    **kwargs
):
    """Run hybracter with only long reads"""
    # Config to add or update in configfile
    merge_config = {
        "args": {
            "input": _input,
            "datadir": datadir,
            "output": output,
            "log": log,
            "min_length": min_length,
            "min_quality": min_quality,
            "skip_qc": skip_qc,
            "databases": databases,
            "subsample_depth": subsample_depth,
            "min_depth": min_depth,
            "medakaModel": medakaModel,
            "flyeModel": flyeModel,
            "contaminants": contaminants,
            "dnaapler_custom_db": dnaapler_custom_db,
            "no_medaka": no_medaka,
            "auto": auto,
            "mac": mac,
            "depth_filter": depth_filter,
            "single": False,
            "logic": logic,
            "configfile": configfile
        }
    }

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "long.smk")),
        merge_config=merge_config,
        configfile=configfile,
        log=log,
        **kwargs
    )


"""
long single
"""


@click.command(
    epilog=help_long_single_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option(
    "-l", "--longreads", help="FASTQ file of longreads", type=str, required=True
)
@click.option(
    "-s", "--sample", help="Sample name.", type=str, default="sample", show_default=True
)
@click.option(
    "-c",
    "--chromosome",
    help="FApproximate lower-bound chromosome length (in base pairs).",
    type=int,
    default=1000000,
    show_default=True,
)
@common_options
@click.option(
            "--logic",
            "logic",
            help="Hybracter logic to select best assembly. Use --best to pick best assembly based on ALE (hybrid) or pyrodigal mean length (long). Use --last to pick the last polishing round regardless.",
            show_default=True,
            default="best",
            type=click.Choice(
                [
                    "best",
                    "last",
                ]
            ),
        )
def long_single(
    longreads,
    sample,
    chromosome,
    medakaModel,
    databases,
    subsample_depth,
    min_depth,
    skip_qc,
    flyeModel,
    min_length,
    output,
    min_quality,
    contaminants,
    dnaapler_custom_db,
    log,
    no_medaka,
    auto,
    mac,
    logic,
    depth_filter,
    configfile,
    **kwargs
):
    """Run hybracter long on 1 isolate"""
    # Config to add or update in configfile
    merge_config = {
        "args": {
            "longreads": longreads,
            "sample": sample,
            "chromosome": chromosome,
            "output": output,
            "log": log,
            "min_length": min_length,
            "min_quality": min_quality,
            "skip_qc": skip_qc,
            "databases": databases,
            "subsample_depth": subsample_depth,
            "min_depth": min_depth,
            "medakaModel": medakaModel,
            "flyeModel": flyeModel,
            "contaminants": contaminants,
            "dnaapler_custom_db": dnaapler_custom_db,
            "no_medaka": no_medaka,
            "auto": auto,
            "mac": mac,
            "depth_filter": depth_filter,
            "single": True,
            "logic": logic,
            "configfile": configfile
        }
    }

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "long.smk")),
        configfile=configfile,
        merge_config=merge_config,
        log=log,
        **kwargs
    )


"""
install
"""


@click.command(
    epilog=help_msg_install,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option(
    "--use-conda/--no-use-conda",
    default=True,
    help="Use conda for Snakemake rules",
    show_default=True,
)
@click.option(
    "--snake-default",
    multiple=True,
    default=[
        "--rerun-incomplete",
        "--printshellcmds",
        "--nolock",
        "--show-failed-logs",
        "--conda-frontend conda",
    ],
    help="Customise Snakemake runtime args",
    show_default=True,
)
@click.option(
    "-d",
    "--databases",
    "databases",
    help="Directory where the Plassembler Database will be installed to (optional).",
    show_default=True,
)
@click.option(
    "-m",
    "--medaka",
    "medaka",
    help="Download medaka models.",
    is_flag=True,
    default=False
)
@click.option(
            "--mac",
            help="If you are running Hybracter on Mac - installs v1.8.0 of Medaka as higher versions break.",
            is_flag=True,
            default=False,
)
@click.option(
    "-o",
    "--output",
    "output",
    help="Temporary directory where intermediate files will be stored for hybracter install. \n This will be deleted.",
    type=click.Path(dir_okay=True, writable=True, readable=True),
    default="hybracter_install_intermediate_files",
    show_default=True,
)
@click.option(
    "--configfile",
    default="config.yaml",
    show_default=False,
    callback=default_to_ouput,
    help="Custom config file [default: (outputDir)/config.yaml]",
)
@click.option(
    "--log",
    "log",
    default="hybracter.log",
    callback=default_to_ouput,
    hidden=True,
)
@click.argument("snake_args", nargs=-1)
def install(databases, output, log, medaka, mac, **kwargs):
    # define both together
    """Downloads and installs the plassembler database"""
    merge_config = {"args": {"databases": databases, "output": output, "log": log, "medaka": medaka, "mac": mac}}
    run_snakemake(
        snakefile_path=snake_base(os.path.join("workflow", "install.smk")),
        merge_config=merge_config,
        log=log,
        **kwargs
    )


"""
test hybrid
"""


# Test command
@click.command(
    epilog=help_msg_test_hybrid,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@common_options
@click.option(
    "--no_pypolca",
    help="Do not use pypolca to polish assemblies with short reads",
    is_flag=True,
    default=False,
)
@click.option(
            "--logic",
            "logic",
            help="Hybracter logic to select best assembly. Use --best to pick best assembly based on ALE (hybrid) or pyrodigal mean length (long). Use --last to pick the last polishing round regardless.",
            show_default=True,
            default="last",
            type=click.Choice(
                [
                    "best",
                    "last",
                ]
            ),
        )
def test_hybrid(
    output,
    log,
    no_medaka,
    auto,
    min_length,
    min_quality,
    skip_qc,
    medakaModel,
    flyeModel,
    databases,
    subsample_depth,
    min_depth,
    no_pypolca,
    mac,
    contaminants,
    dnaapler_custom_db,
    logic,
    depth_filter,
    configfile,
    **kwargs
):
    """Test hybracter hybrid"""
    # Config to add or update in configfile
    merge_config = {
        "args": {
            "output": output,
            "log": log,
            "min_length": min_length,
            "min_quality": min_quality,
            "skip_qc": skip_qc,
            "medakaModel": medakaModel,
            "flyeModel": flyeModel,
            "databases": databases,
            "subsample_depth": subsample_depth,
            "min_depth": min_depth,
            "no_pypolca": no_pypolca,
            "contaminants": contaminants,
            "dnaapler_custom_db": dnaapler_custom_db,
            "no_medaka": no_medaka,
            "auto": auto,
            "mac": mac,
            "logic": logic,
            "depth_filter": depth_filter,
            "configfile": configfile
        }
    }
    run_snakemake(
        snakefile_path=snake_base(os.path.join("workflow", "test_hybrid.smk")),
        merge_config=merge_config,
        configfile=configfile,
        log=log,
        **kwargs
    )


"""
test long
"""


# Test command
@click.command(
    epilog=help_msg_test_long,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option(
            "--logic",
            "logic",
            help="Hybracter logic to select best assembly. Use --best to pick best assembly based on ALE (hybrid) or pyrodigal mean length (long). Use --last to pick the last polishing round regardless.",
            show_default=True,
            default="best",
            type=click.Choice(
                [
                    "best",
                    "last",
                ]
            ),
        )
@common_options
def test_long(
    output,
    log,
    min_length,
    min_quality,
    skip_qc,
    medakaModel,
    databases,
    subsample_depth,
    min_depth,
    flyeModel,
    contaminants,
    dnaapler_custom_db,
    no_medaka,
    auto,
    mac,
    logic,
    depth_filter,
    configfile,
    **kwargs
):
    """Test hybracter long"""
    # Config to add or update in configfile
    merge_config = {
        "args": {
            "output": output,
            "log": log,
            "min_length": min_length,
            "min_quality": min_quality,
            "skip_qc": skip_qc,
            "medakaModel": medakaModel,
            "databases": databases,
            "subsample_depth": subsample_depth,
            "min_depth": min_depth,
            "flyeModel": flyeModel,
            "contaminants": contaminants,
            "dnaapler_custom_db": dnaapler_custom_db,
            "no_medaka": no_medaka,
            "auto": auto,
            "mac": mac,
            "logic": logic,
            "depth_filter": depth_filter,
            "configfile": configfile
        }
    }
    run_snakemake(
        snakefile_path=snake_base(os.path.join("workflow", "test_long.smk")),
        configfile=configfile,
        merge_config=merge_config,
        log=log,
        **kwargs
    )


@click.command()
@common_options
def config(configfile, **kwargs):
    """Copy the system default config file"""
    copy_config(configfile)


@click.command()
def citation(**kwargs):
    """Print the citation(s) for hybracter"""
    print_citation()


@click.command()
def version(**kwargs):
    """Print the version for hybracter"""


cli.add_command(install)
cli.add_command(hybrid)
cli.add_command(hybrid_single)
cli.add_command(long)
cli.add_command(long_single)
cli.add_command(test_hybrid)
cli.add_command(test_long)
cli.add_command(config)
cli.add_command(citation)
cli.add_command(version)


def main():
    print_version()
    print_splash()
    cli()


if __name__ == "__main__":
    main()

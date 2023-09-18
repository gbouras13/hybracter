"""
Entrypoint for hybracter

Check out the wiki for a detailed look at customising this file:
https://github.com/beardymcjohnface/Snaketool/wiki/Customising-your-Snaketool
"""

import os

import click

from .util import (
    OrderedCommands,
    copy_config,
    default_to_ouput,
    print_citation,
    print_version,
    run_snakemake,
    snake_base,
    all_medaka_models,
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
            default="hybracter.out",
            show_default=True,
        ),
        click.option(
            "--configfile",
            default="config.yaml",
            show_default=False,
            callback=default_to_ouput,
            help="Custom config file [default: (outputDir)/config.yaml]",
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
        ),
        click.option(
            "--min_quality",
            "min_quality",
            help="min read quality for long reads",
            type=int,
            default=9,
        ),
        click.option(
            "--skip_qc",
            help="Do not run porechop, filtlong and fastp to QC the reads",
            is_flag=True,
            default=False,
        ),
        click.option(
            "-d",
            "--databases",
            help="Plassembler Databases directory.",
            type=click.Path(dir_okay=True, readable=True),
            default="plassembler_DB",
        ),
        click.option(
            "--medakaModel",
            "medakaModel",
            help="Medaka Model.",
            default="r1041_e82_400bps_sup_v4.2.0",
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
                "--conda-frontend mamba",
            ],
            help="Customise Snakemake runtime args",
            show_default=True,
        ),
        click.option(
            "--log",
            default="hybracter.log",
            callback=default_to_ouput,
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
Specify threads:    hybracter hybrid ... --threads [threads]
Disable conda:      hybracter hybrid ... --no-use-conda 
Change defaults:    hybracter hybrid ... --snake-default="-k --nolock"
Add Snakemake args: hybracter hybrid ... --dry-run --keep-going --touch
Specify targets:    hybracter hybrid ... all print_targets
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
Required:           hybracter hybrid --input [file]
Specify threads:    hybracter hybrid ... --threads [threads]
Disable conda:      hybracter hybrid ... --no-use-conda 
Change defaults:    hybracter hybrid ... --snake-default="-k --nolock"
Add Snakemake args: hybracter hybrid ... --dry-run --keep-going --touch
Specify targets:    hybracter hybrid ... all print_targets
Available targets:
    all             Run everything (default)
    print_targets   List available targets
"""

help_msg_download = """
\b
Downloads the plassembler database
hybracter download ... 
\b
RUN EXAMPLES:
Database:           hybracter install --download [directory]
"""

help_msg_ale = """
\b
Compiles ale 
RUN EXAMPLES:
hybracter ale
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
@click.option(
    "--no_polca",
    help="Do not use Polca to polish assemblies with short reads",
    is_flag=True,
    default=False,
)
@common_options
def hybrid(
    _input,
    no_polca,
    skip_qc,
    medakaModel,
    databases,
    min_quality,
    flyeModel,
    min_length,
    output,
    log,
    **kwargs
):
    """Run hybracter with hybrid long and paired end short reads"""
    # Config to add or update in configfile
    merge_config = {
        "args": {
            "input": _input,
            "output": output,
            "log": log,
            "min_length": min_length,
            "databases": databases,
            "min_quality": min_quality,
            "no_polca": no_polca,
            "skip_qc": skip_qc,
            "medakaModel": medakaModel,
            "flyeModel": flyeModel,
        }
    }

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "hybrid.smk")),
        merge_config=merge_config,
        log=log,
        **kwargs
    )


@click.command(
    epilog=help_long_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option("-i", "--input", "_input", help="Input csv", type=str, required=True)
@common_options
def long(
    _input,
    medakaModel,
    databases,
    skip_qc,
    flyeModel,
    min_length,
    output,
    min_quality,
    log,
    **kwargs
):
    """Run hybracter long"""
    # Config to add or update in configfile
    merge_config = {
        "args": {
            "input": _input,
            "output": output,
            "log": log,
            "min_length": min_length,
            "min_quality": min_quality,
            "skip_qc": skip_qc,
            "databases": databases,
            "medakaModel": medakaModel,
            "flyeModel": flyeModel,
        }
    }

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "long.smk")),
        merge_config=merge_config,
        log=log,
        **kwargs
    )


@click.command(
    epilog=help_msg_download,
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
    help="Directory where the Plassembler Database will be downloaded to.",
    show_default=True,
    default="plassembler_DB",
)
def download(databases, log, **kwargs):
    # Config to add or update in configfile
    merge_config = {"args": {"databases": databases, "log": log}}
    """Downloads the plassembler database"""
    run_snakemake(
        snakefile_path=snake_base(os.path.join("workflow", "download.smk")),
        merge_config=merge_config,
        **kwargs
    )


@click.command()
@common_options
def config(configfile, **kwargs):
    """Copy the system default config file"""
    copy_config(configfile)


@click.command()
def citation(**kwargs):
    """Print the citation(s) for this tool"""
    print_citation()


cli.add_command(download)
cli.add_command(hybrid)
cli.add_command(long)
cli.add_command(config)
cli.add_command(citation)


def main():
    print_version()
    print_splash()
    cli()


if __name__ == "__main__":
    main()

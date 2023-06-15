"""
Entrypoint for hybracter

Check out the wiki for a detailed look at customising this file:
https://github.com/beardymcjohnface/Snaketool/wiki/Customising-your-Snaketool
"""

import os
import click

from .util import (
    snake_base,
    print_version,
    default_to_ouput,
    copy_config,
    run_snakemake,
    OrderedCommands,
    print_citation,
)


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
        click.option(
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
            "--threads", help="Number of threads to use", default=1, show_default=True
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
                "--conda-frontend mamba"
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

help_msg_install = """
\b
Installed the plassembler database
hybracter install ... 
\b
RUN EXAMPLES:
Database:           hybracter install --database [file]
"""


@click.command(
    epilog=help_hybrid_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option("--input", "_input", help="Input csv", type=str, required=True)
@click.option('--polca',  help='Use Polca to polish assemblies with short reads', is_flag=True,  default=False)
@click.option('--min_length',  help='min read length for long reads', type=int,  default=False)
@click.option('--plasmids',  help='whether you want to use Plassembler for plasmid recovery.', is_flag=True,  default=False)
@click.option('--no_polish','no_polish',  help='whether you want to turn off Medaka to polishing for your genome.', is_flag=True,  default=False )
@click.option('--medakaModel','medakaModel',  help='Medaka Model.', default='r941_min_sup_g507', show_default=True, type=click.Choice(['r941_min_sup_g507', 'r941_min_hac_g507', 'r941_min_fast_g507', 'r1041_e82_400bps_sup_g615']) )
@click.option('--flyeModel','flyeModel',  help='Flye Assembly Parameter', show_default=True,  default='--nano-hq',type=click.Choice(['--nano-hq', '--nano-corr', '--nano-raw', "--pacbio-raw", "--pacbio-corr", "--pacbio-hifi"]))
@common_options
def hybrid(_input,  polca, medakaModel, plasmids, no_polish, flyeModel, min_length, output, log, **kwargs):
    """Run hybracter"""
    # Config to add or update in configfile
    merge_config = {
        "input": _input, 
        "output": output, 
        "log": log, 
        "min_length": min_length,
        "plasmids": plasmids,
        "no_polish": no_polish,
        "polca": polca, 
        "medakaModel": medakaModel, 
        "flyeModel": flyeModel }

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
@click.option("--input", "_input", help="Input csv", type=str, required=True)
@click.option('--min_length',  help='min read length for long reads', type=int,  default=False)
@click.option('--plasmids',  help='whether you want to use Plassembler for plasmid recovery. Long only mode. Experimental.', is_flag=True,  default=False)
@click.option('--no_polish','no_polish',  help='whether you want to turn off Medaka to polishing for your genome.', is_flag=True,  default=False )
@click.option('--medakaModel','medakaModel',  help='Medaka Model.', default='r941_min_sup_g507', show_default=True, type=click.Choice(['r941_min_sup_g507', 'r941_min_hac_g507', 'r941_e81_fast_g514', 'r1041_e82_400bps_sup_g615']) )
@click.option('--flyeModel','flyeModel',  help='Flye Assembly Parameter', show_default=True,  default='--nano-hq',type=click.Choice(['--nano-hq', '--nano-corr', '--nano-raw', "--pacbio-raw", "--pacbio-corr", "--pacbio-hifi"]))
@common_options
def long(_input, medakaModel, plasmids, no_polish, flyeModel, min_length, output, log, **kwargs):
    """Run hybracter"""
    # Config to add or update in configfile
    merge_config = {
        "input": _input, 
        "output": output, 
        "log": log, 
        "min_length": min_length,
        "plasmids": plasmids,
        "no_polish": no_polish,
        "medakaModel": medakaModel, 
        "flyeModel": flyeModel }

    # run!
    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join("workflow", "long.smk")),
        merge_config=merge_config,
        log=log,
        **kwargs
    )


@click.command(
    epilog=help_msg_install,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ))
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
                "--conda-frontend conda"
            ],
            help="Customise Snakemake runtime args",
            show_default=True,
        )
@click.option('--database','database',  help='Database Directory', show_default=True,  default='Database')
@common_options
def install( database, log,output,  **kwargs):
    # Config to add or update in configfile
    merge_config = { "database": database, "output": output, "log": log }
    """Install databases"""
    run_snakemake(
        snakefile_path=snake_base(os.path.join('workflow','installDB.smk')),
        merge_config=merge_config,
        **kwargs)



@click.command()
@common_options
def config(configfile, **kwargs):
    """Copy the system default config file"""
    copy_config(configfile)


@click.command()
def citation(**kwargs):
    """Print the citation(s) for this tool"""
    print_citation()

cli.add_command(install)
cli.add_command(hybrid)
cli.add_command(long)
cli.add_command(config)
cli.add_command(citation)


def main():
    print_version()
    cli()


if __name__ == "__main__":
    main()

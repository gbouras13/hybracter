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

def print_splash():
    click.echo("""\b

 _           _                    _            
| |__  _   _| |__  _ __ __ _  ___| |_ ___ _ __ 
| '_ \| | | | '_ \| '__/ _` |/ __| __/ _ \ '__|
| | | | |_| | |_) | | | (_| | (__| ||  __/ |   
|_| |_|\__, |_.__/|_|  \__,_|\___|\__\___|_|   
       |___/

""")


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
Installs ALE from source and downloads the plassembler database
hybracter install ... 
\b
RUN EXAMPLES:
Database:           hybracter install --database [file]
"""

# from medaka models.py 
# https://github.com/nanoporetech/medaka/blob/938167e2ff804899d578d04388ba5c15ab339316/medaka/models.py

all_medaka_models = [
    # r1041 e82 (kit14) consensus
    'r1041_e82_400bps_hac_v4.2.0',
    'r1041_e82_400bps_sup_v4.2.0',
    # r1041 variant calling
    # 'r1041_e82_400bps_hac_variant_v4.2.0',
    # 'r1041_e82_400bps_sup_variant_v4.2.0',
    # r9 consensus
    'r941_sup_plant_g610',
    'r941_min_fast_g507', 'r941_prom_fast_g507',
    'r941_min_fast_g303', 'r941_min_high_g303', 'r941_min_high_g330',
    'r941_prom_fast_g303', 'r941_prom_high_g303', 'r941_prom_high_g330',
    'r941_min_high_g344', 'r941_min_high_g351', 'r941_min_high_g360',
    'r941_prom_high_g344', 'r941_prom_high_g360', 'r941_prom_high_g4011',
    # r10 consensus
    'r10_min_high_g303', 'r10_min_high_g340',
    'r103_min_high_g345', 'r103_min_high_g360', 'r103_prom_high_g360',
    'r103_fast_g507', 'r103_hac_g507', 'r103_sup_g507',
    # r104 e81 consensus
    'r104_e81_fast_g5015', 'r104_e81_sup_g5015', 'r104_e81_hac_g5015',
    'r104_e81_sup_g610',
    # r104 e81 variant calling
    # 'r104_e81_fast_variant_g5015', 'r104_e81_hac_variant_g5015',
    # 'r104_e81_sup_variant_g610',
    # r1041 e82 consensus
    'r1041_e82_400bps_hac_g615',  'r1041_e82_400bps_fast_g615',
    'r1041_e82_400bps_fast_g632', 'r1041_e82_260bps_fast_g632',
    'r1041_e82_400bps_hac_g632', 'r1041_e82_400bps_sup_g615',
    'r1041_e82_260bps_hac_g632', 'r1041_e82_260bps_sup_g632',
    'r1041_e82_400bps_hac_v4.0.0', 'r1041_e82_400bps_sup_v4.0.0',
    'r1041_e82_260bps_hac_v4.0.0', 'r1041_e82_260bps_sup_v4.0.0',
    'r1041_e82_260bps_hac_v4.1.0', 'r1041_e82_260bps_sup_v4.1.0',
    'r1041_e82_400bps_hac_v4.1.0', 'r1041_e82_400bps_sup_v4.1.0',
    # r1041 e82 variant calling
    # 'r1041_e82_400bps_hac_variant_g615',
    # 'r1041_e82_400bps_fast_variant_g615',
    # 'r1041_e82_400bps_fast_variant_g632',
    # 'r1041_e82_260bps_fast_variant_g632',
    # 'r1041_e82_400bps_hac_variant_g632',
    # 'r1041_e82_400bps_sup_variant_g615',
    # 'r1041_e82_260bps_hac_variant_g632',
    # 'r1041_e82_260bps_sup_variant_g632',
    # 'r1041_e82_260bps_hac_variant_v4.1.0',
    # 'r1041_e82_260bps_sup_variant_v4.1.0',
    # 'r1041_e82_400bps_hac_variant_v4.1.0',
    # 'r1041_e82_400bps_sup_variant_v4.1.0',
    # snp and variant - flipflop
    # 'r941_prom_snp_g303', 'r941_prom_variant_g303',
    # 'r941_prom_snp_g322', 'r941_prom_variant_g322',
    # 'r941_prom_snp_g360', 'r941_prom_variant_g360',
    # 'r103_prom_snp_g3210', 'r103_prom_variant_g3210',
    # snp and variant - crf guppy507+
    # 'r941_sup_plant_variant_g610',
    # 'r941_min_fast_snp_g507', 'r941_min_fast_variant_g507',
    # 'r941_min_hac_snp_g507',
    # 'r941_min_sup_snp_g507', 'r941_min_sup_variant_g507',
    # 'r941_prom_fast_snp_g507', 'r941_prom_fast_variant_g507',
    # 'r941_prom_hac_snp_g507',
    # 'r941_prom_sup_snp_g507', 'r941_prom_sup_variant_g507',
    # 'r103_fast_snp_g507', 'r103_fast_variant_g507',
    # 'r103_hac_snp_g507', 'r103_hac_variant_g507',
    # 'r103_sup_snp_g507', 'r103_sup_variant_g507',
    # rle consensus
    'r941_min_high_g340_rle',
    # r9 consensus
    'r941_min_hac_g507', 'r941_min_sup_g507',
    'r941_prom_hac_g507', 'r941_prom_sup_g507',
    # r9 variant calling
    # 'r941_min_hac_variant_g507',
    # 'r941_prom_hac_variant_g507',
    # r941 e81 consensus
    'r941_e81_fast_g514', 'r941_e81_hac_g514', 'r941_e81_sup_g514',
    # r941 e81 variant calling
    # 'r941_e81_fast_variant_g514', 'r941_e81_hac_variant_g514',
    # 'r941_e81_sup_variant_g514',
]


"""
hybrid
"""

@click.command(
    epilog=help_hybrid_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option("--input", "_input", help="Input csv", type=str, required=True)
@click.option('--no_polca',  help='Do not use Polca to polish assemblies with short reads', is_flag=True,  default=False)
@click.option('--min_length', 'min_length',  help='min read length for long reads', type=int,  default=1000)
@click.option('--min_quality', 'min_quality',  help='min read quality for long reads', type=int,  default=9)
@click.option('--databases',  help='Plassembler Databases directory.', type=click.Path(dir_okay=True, readable=True),  default='plassembler_DB')
@click.option('--medakaModel','medakaModel',  help='Medaka Model.', default='r1041_e82_400bps_sup_v4.2.0', show_default=True, type=click.Choice(all_medaka_models) )
@click.option('--flyeModel','flyeModel',  help='Flye Assembly Parameter', show_default=True,  default='--nano-hq', type=click.Choice(['--nano-hq', '--nano-corr', '--nano-raw', "--pacbio-raw", "--pacbio-corr", "--pacbio-hifi"]))
@common_options
def hybrid(_input,  no_polca, medakaModel, databases, min_quality, flyeModel, min_length, output, log, **kwargs):
    """Run hybracter with hybrid long and paired end short reads"""
    # Config to add or update in configfile
    merge_config = {
        'args': {
        "input": _input, 
        "output": output, 
        "log": log, 
        "min_length": min_length,
        "databases": databases,
        "min_quality": min_quality,
        "no_polca": no_polca, 
        "medakaModel": medakaModel, 
        "flyeModel": flyeModel } }

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
         'args': {
        "input": _input, 
        "output": output, 
        "log": log, 
        "min_length": min_length,
        "plasmids": plasmids,
        "no_polish": no_polish,
        "medakaModel": medakaModel, 
        "flyeModel": flyeModel } }

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
@click.option('--databases','databases',  help='Plassembler databases Directory', show_default=True,  default='Database')
@common_options
def download( databases, log,output,  **kwargs):
    # Config to add or update in configfile
    merge_config = { "databases": databases, "output": output, "log": log }
    """Downloads the plassembler database"""
    run_snakemake(
        snakefile_path=snake_base(os.path.join('workflow','download.smk')),
        merge_config=merge_config,
        **kwargs)


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
@common_options
def ale(  log,output,  **kwargs):
    # install ale
    merge_config = { "output": output, "log": log }
    """installs ale"""
    run_snakemake(
        snakefile_path=snake_base(os.path.join('workflow','ale.smk')),
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

cli.add_command(download)
cli.add_command(hybrid)
cli.add_command(ale)
cli.add_command(long)
cli.add_command(config)
cli.add_command(citation)


def main():
    print_version()
    print_splash()
    cli()


if __name__ == "__main__":
    main()

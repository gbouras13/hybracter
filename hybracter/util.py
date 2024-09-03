"""MISC FUNCTIONS
You shouldn't need to tweak these much if at all
"""

import os
import subprocess
import sys
from shutil import copyfile
from time import localtime, strftime

import click
import yaml


class OrderedCommands(click.Group):
    """This class will preserve the order of subcommands, which is useful when printing --help"""

    def list_commands(self, ctx: click.Context):
        return list(self.commands)


def snake_base(rel_path):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path)


def print_version():
    with open(snake_base("hybracter.VERSION"), "r") as f:
        version = f.readline()
    echo_click("\n" + "hybracter version " + version + "\n")


def echo_click(msg, log=None):
    click.echo(msg, nl=False, err=True)
    if log:
        # if log already exists
        if os.path.exists(log):
            with open(log, "a") as l:
                l.write(msg)
        # to create it
        else:
            directory_path = os.path.dirname(log)
            # Create the directories recursively
            os.makedirs(directory_path, exist_ok=True)
            with open(log, "w") as l:
                l.write(msg)




def print_citation():
    with open(snake_base("hybracter.CITATION"), "r") as f:
        for line in f:
            echo_click(line)


def msg(err_message, log=None):
    tstamp = strftime("[%Y:%m:%d %H:%M:%S] ", localtime())
    try:
        echo_click(tstamp + err_message + "\n", log=log)
    except FileNotFoundError:
        print("Cleaning up intermediate files.")


def msg_box(splash, errmsg=None, log=None):
    msg("-" * (len(splash) + 4), log=log)
    msg(f"| {splash} |", log=log)
    msg(("-" * (len(splash) + 4)), log=log)
    if errmsg:
        echo_click("\n" + errmsg + "\n", log=log)


def default_to_ouput(ctx, param, value):
    """Callback for --configfile; place value in output directory unless specified"""
    if param.default == value:
        return os.path.join(ctx.params["output"], value)
    return value


def write_config(_config, file, log=None):
    msg(f"Writing config file to {file}", log=log)
    with open(file, "w") as stream:
        yaml.dump(_config, stream)


def update_config(in_config=None, merge=None, output_config=None, log=None):
    """Update config with new values"""

    if output_config is None:
        output_config = in_config

    # read the config
    config = read_config(in_config)

    # merge additional config values
    msg("Updating config file with new values", log=log)
    config.update(merge)

    write_config(config, output_config, log=log)


def copy_config(
    local_config,
    merge_config=None,
    system_config=snake_base(os.path.join("config", "config.yaml")),
    log=None,
):
    if not os.path.isfile(local_config):
        if len(os.path.dirname(local_config)) > 0:
            os.makedirs(os.path.dirname(local_config), exist_ok=True)
        msg(f"Copying system default config to {local_config}", log=log)

        if merge_config:
            update_config(
                in_config=system_config,
                merge=merge_config,
                output_config=local_config,
                log=log,
            )
        else:
            copyfile(system_config, local_config)
    else:
        msg(
            f"Config file {local_config} already exists. Using existing config file.",
            log=log,
        )


def read_config(file):
    """Read a config file to a dictionary

    Args:
        file (str): Filepath to config YAML file for reading

    Returns (dict): Config read from YAML file
    """

    with open(file, "r") as stream:
        config = yaml.safe_load(stream)
    return config


"""RUN A SNAKEFILE
Hopefully you shouldn't need to tweak this function at all.
- You must provide a Snakefile, all else is optional
- Highly recommend supplying a configfile and the default snakemake args"""


def run_snakemake(
    configfile=None,
    snakefile_path=None,
    merge_config=None,
    threads=1,
    use_conda=False,
    conda_prefix=None,
    snake_default=None,
    snake_args=[],
    log=None,
):
    """Run a Snakefile"""
    snake_command = ["snakemake", "-s", snakefile_path]

    # if using a configfile
    if configfile:

        # copy sys default config if not present
        copy_config(configfile, log=log)

        if merge_config:
            update_config(in_config=configfile, merge=merge_config, log=log)

        snake_command += ["--configfile", configfile]

        # display the runtime configuration
        snake_config = read_config(configfile)
        msg_box(
            "Runtime config",
            errmsg=yaml.dump(snake_config, Dumper=yaml.Dumper),
            log=log,
        )

    # add threads
    if not "--profile" in snake_args:
        snake_command += ["--jobs", threads]

    # add conda args if using conda
    if use_conda:
        snake_command += ["--use-conda"]
        if conda_prefix:
            snake_command += ["--conda-prefix", conda_prefix]

    # add snakemake default args
    if snake_default:
        snake_command += snake_default

    # add any additional snakemake commands
    if snake_args:
        snake_command += list(snake_args)

    # Run Snakemake!!!
    snake_command = " ".join(str(s) for s in snake_command)
    msg_box("Snakemake command", errmsg=snake_command, log=log)
    if not subprocess.run(snake_command, shell=True).returncode == 0:
        msg("ERROR: Snakemake failed", log=log)
        sys.exit(1)
    else:
        msg("Snakemake finished successfully", log=log)
    return 0


"""
list of all available medaka models
"""

# from medaka models.py
# https://github.com/nanoporetech/medaka/blob/aa4da2429a5e32d02369ff0b67fad6edc6c4a94e/medaka/options.py

all_medaka_models = [
    
    # default 1.12.1
    'r1041_e82_400bps_sup_v5.0.0',
    # r1041 e82 (kit14) consensus
    'r1041_e82_400bps_hac_v5.0.0',
    'r1041_e82_400bps_hac_v4.3.0',
    'r1041_e82_400bps_sup_v4.3.0',
    "r1041_e82_400bps_hac_v4.2.0",
    "r1041_e82_400bps_sup_v4.2.0",
    # r1041 variant calling
    # 'r1041_e82_400bps_hac_variant_v4.2.0',
    # 'r1041_e82_400bps_sup_variant_v4.2.0',
    # r9 consensus
    "r941_sup_plant_g610",
    "r941_min_fast_g507",
    "r941_prom_fast_g507",
    "r941_min_fast_g303",
    "r941_min_high_g303",
    "r941_min_high_g330",
    "r941_prom_fast_g303",
    "r941_prom_high_g303",
    "r941_prom_high_g330",
    "r941_min_high_g344",
    "r941_min_high_g351",
    "r941_min_high_g360",
    "r941_prom_high_g344",
    "r941_prom_high_g360",
    "r941_prom_high_g4011",
    # r10 consensus
    "r10_min_high_g303",
    "r10_min_high_g340",
    "r103_min_high_g345",
    "r103_min_high_g360",
    "r103_prom_high_g360",
    "r103_fast_g507",
    "r103_hac_g507",
    "r103_sup_g507",
    # r104 e81 consensus
    "r104_e81_fast_g5015",
    "r104_e81_sup_g5015",
    "r104_e81_hac_g5015",
    "r104_e81_sup_g610",
    # r104 e81 variant calling
    # 'r104_e81_fast_variant_g5015', 'r104_e81_hac_variant_g5015',
    # 'r104_e81_sup_variant_g610',
    # r1041 e82 consensus
    "r1041_e82_400bps_hac_g615",
    "r1041_e82_400bps_fast_g615",
    "r1041_e82_400bps_fast_g632",
    "r1041_e82_260bps_fast_g632",
    "r1041_e82_400bps_hac_g632",
    "r1041_e82_400bps_sup_g615",
    "r1041_e82_260bps_hac_g632",
    "r1041_e82_260bps_sup_g632",
    "r1041_e82_400bps_hac_v4.0.0",
    "r1041_e82_400bps_sup_v4.0.0",
    "r1041_e82_260bps_hac_v4.0.0",
    "r1041_e82_260bps_sup_v4.0.0",
    "r1041_e82_260bps_hac_v4.1.0",
    "r1041_e82_260bps_sup_v4.1.0",
    "r1041_e82_400bps_hac_v4.1.0",
    "r1041_e82_400bps_sup_v4.1.0",
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
    "r941_min_high_g340_rle",
    # r9 consensus
    "r941_min_hac_g507",
    "r941_min_sup_g507",
    "r941_prom_hac_g507",
    "r941_prom_sup_g507",
    # r9 variant calling
    # 'r941_min_hac_variant_g507',
    # 'r941_prom_hac_variant_g507',
    # r941 e81 consensus
    "r941_e81_fast_g514",
    "r941_e81_hac_g514",
    "r941_e81_sup_g514",
    # r941 e81 variant calling
    # 'r941_e81_fast_variant_g514', 'r941_e81_hac_variant_g514',
    # 'r941_e81_sup_variant_g514',
]

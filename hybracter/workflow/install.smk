"""
Snakefile for downloading plassembler database 
"""

import os
import glob
import attrmap as ap
import attrmap.utils as au
from pathlib import Path


# config file
configfile: os.path.join(workflow.basedir, "../", "config", "config.yaml")


config = ap.AttrMap(config)

# flag to download all medaka models
MEDAKA_DOWNLOAD = config.args.medaka
MAC = config.args.mac


# directories
include: os.path.join("rules", "preflight", "directories.smk")
# targets
include: os.path.join("rules", "preflight", "targets_download.smk")
# plassembler
include: os.path.join("rules", "download", "download_plassembler_db.smk")
include: os.path.join("rules", "download", "download_medaka_models.smk")


# define the rule all
rule all:
    input:
        TargetFilesDownload,

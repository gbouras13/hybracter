"""
Snakefile for downloading plassembler database 
"""
import os
import glob
import attrmap as ap
import attrmap.utils as au
from pathlib import Path

# Concatenate Snakemake's own log file with the master log file
def copy_log_file():
    files = glob.glob(os.path.join(".snakemake", "log", "*.snakemake.log"))
    if not files:
        return None
    current_log = max(files, key=os.path.getmtime)
    shell("cat " + current_log + " >> " + config['log'])

onsuccess:
    copy_log_file()

onerror:
    copy_log_file()


# config file 
configfile: os.path.join(workflow.basedir, '../', 'config', 'config.yaml')
config = ap.AttrMap(config)

# directories
include: os.path.join("rules", "preflight", "directories.smk")
# targets
include: os.path.join("rules", "preflight", "targets_download.smk")
# plassembler
include: os.path.join("rules", "download", "download_plassembler_db.smk")


# define the rule all
rule all:
    input:
        TargetFilesDownload


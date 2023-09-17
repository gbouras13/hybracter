
import os
import glob
import attrmap as ap
import attrmap.utils as au
from pathlib import Path

# Concatenate Snakemake's own log file with the master log file
# log defined below
def copy_log_file():
    files = glob.glob(os.path.join(".snakemake", "log", "*.snakemake.log"))
    if not files:
        return None
    current_log = max(files, key=os.path.getmtime)
    shell("cat " + current_log + " >> " + LOG)

onsuccess:
    copy_log_file()

onerror:
    copy_log_file()

# config file 
configfile: os.path.join(workflow.basedir, '../', 'config', 'config.yaml')
config = ap.AttrMap(config)


# directories
include: os.path.join("rules", "preflight", "directories.smk")
# functions
include: os.path.join("rules", "preflight", "functions.smk")
# samples
include: os.path.join("rules", "preflight", "samples.smk")
# targets
include: os.path.join("rules", "preflight", "targets_hybrid.smk")


### from config files
#  input as csv
INPUT = config.args.input
OUTPUT = config.args.output
LOG = os.path.join(OUTPUT, 'hybracter.log')
THREADS = config.args.threads
MIN_LENGTH = config.args.min_length
MIN_QUALITY = config.args.min_quality
MEDAKA_MODEL = config.args.medakaModel
FLYE_MODEL = config.args.flyeModel

# Parse the samples and read files

dictReads = parseSamples(INPUT, False) # long flag false
SAMPLES = list(dictReads.keys())
wildcard_constraints:
    sample= '|'.join([re.escape(x) for x in SAMPLES])

##############################
# Import rules and functions
##############################

# qc
include: os.path.join("rules", "processing", "qc.smk")
# assembly
include: os.path.join("rules", "assembly", "assemble.smk")

# get flye stats and combine across all runs
include: os.path.join("rules", "processing", "assembly_statistics.smk")
# extract chrom
include: os.path.join("rules", "processing", "extract_fastas.smk")

# checkpoint
include: os.path.join("rules", "completeness", "aggregate.smk")

# checkpoint here for completeness
# need long read polish files regardless
include: os.path.join("rules", "polishing", "long_read_polish.smk")
include: os.path.join("rules", "polishing", "long_read_polish_incomplete.smk")


#  polish the assemblies
include: os.path.join("rules", "polishing", "short_read_polish.smk")
include: os.path.join("rules", "polishing", "short_read_polish_incomplete.smk")


# plassembler 
include: os.path.join("rules", "assembly", "plassembler.smk")
include: os.path.join("rules", "processing", "combine_plassembler_info.smk")

if config.args.no_polca is False:
    include: os.path.join("rules", "polishing", "short_read_polca.smk")

# ale
include: os.path.join("rules", "assess", "assess_complete.smk")
include: os.path.join("rules", "assess", "assess_incomplete.smk")

# finalse
include: os.path.join("rules", "finalse", "select_best_assembly.smk")

### rule all
rule all:
    input:
        TargetFilesHybrid

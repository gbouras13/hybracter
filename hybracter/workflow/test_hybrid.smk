import glob
import attrmap as ap
import attrmap.utils as au
from pathlib import Path
import os


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
configfile: os.path.join(workflow.basedir, "../", "config", "config.yaml")


config = ap.AttrMap(config)


# directories
include: os.path.join("rules", "preflight", "directories.smk")
# functions
include: os.path.join("rules", "preflight", "functions.smk")


# set plassemblerdb as test
dir.plassemblerdb = os.path.join(dir.test, "Plassembler_DB_Test")

# check db
# from functions.smk
check_db(dir.plassemblerdb)


# include samples
include: os.path.join("rules", "preflight", "samples.smk")
# targets
include: os.path.join("rules", "preflight", "targets_hybrid.smk")


### from config files
#  input as csv
INPUT = config.args.input
OUTPUT = config.args.output
LOG = os.path.join(OUTPUT, "hybracter.log")
THREADS = config.args.threads
MIN_LENGTH = config.args.min_length
MIN_QUALITY = config.args.min_quality
MEDAKA_MODEL = config.args.medakaModel
FLYE_MODEL = config.args.flyeModel

# Parse the samples and read files
# dictReads = parseSamples(INPUT, True)  # long flag true
# SAMPLES = list(dictReads.keys())
# wildcard_constraints:
#     sample="|".join([re.escape(x) for x in SAMPLES]),

# instead of parsing the samples with INPUT, just make the dictionary from scratch

dictReads = {}
dictReads["Sample1"] = {}
dictReads["Sample1"]["LR"] = os.path.join(dir.test_fastqs, "test_long_reads.fastq.gz")
dictReads["Sample1"]["MinChromLength"] = 50000
dictReads["Sample1"]["R1"] = os.path.join(dir.test_fastqs, "test_short_reads_R1.fastq.gz")
dictReads["Sample1"]["R2"] = os.path.join(dir.test_fastqs, "test_short_reads_R2.fastq.gz")

SAMPLES = ["Sample1"]


##############################
# Import rules and functions
##############################

# qc and host
# depends on whehter --contaminants has been specified and --skip_qc flag activiated
if config.args.contaminants != "none":  # where --contaminants specified
    CONTAM = (
        check_host()
    )  # from functions.smk to make sure the specified file is lambda or a FASTA

    include: os.path.join("rules", "processing", "remove_contaminants_qc.smk")

else:  # where no contaminants to be removed
    if config.args.skip_qc is True:

        include: os.path.join("rules", "processing", "skip_qc.smk")

    else:

        include: os.path.join("rules", "processing", "qc.smk")


# assembly
include: os.path.join("rules", "assembly", "assemble.smk")
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
include: os.path.join("rules", "finalise", "select_best_assembly.smk")


### rule all
rule all:
    input:
        TargetFilesHybrid,

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


# samples
include: os.path.join("rules", "preflight", "samples.smk")
# targets
include: os.path.join("rules", "preflight", "targets_long.smk")


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
LOGIC = config.args.logic
DEPTH_FILTER = config.args.depth_filter
SUBSAMPLE_DEPTH = config.args.subsample_depth
MIN_DEPTH = config.args.min_depth
AUTO = config.args.auto
MAC = config.args.mac

# MAC medaka

new_models = [
    "r1041_e82_400bps_sup_v5.0.0",
    "r1041_e82_400bps_hac_v5.0.0",
    "r1041_e82_400bps_hac_v4.3.0",
    "r1041_e82_400bps_sup_v4.3.0",
]

if MAC:
    if MEDAKA_MODEL in new_models:
        print(
            f"{MEDAKA_MODEL} is not available in medaka v1.8.0 as it is too new. If you want this model, try Hybracter on a Linux machine."
        )
        print(f"Changing the medaka model to r1041_e82_400bps_sup_v4.2.0.")
        MEDAKA_MODEL = "r1041_e82_400bps_sup_v4.2.0"

# Parse the samples and read files
# dictReads = parseSamples(INPUT, True)  # long flag true
# SAMPLES = list(dictReads.keys())
# wildcard_constraints:
#     sample="|".join([re.escape(x) for x in SAMPLES]),

# instead of parsing the samples with INPUT, just make the dictionary from scratch

dictReads = {}
dictReads["Sample1"] = {}
dictReads["Sample1"]["LR"] = os.path.join(dir.test_fastqs, "test_long_reads.fastq.gz")
if not AUTO:
    dictReads["Sample1"]["MinChromLength"] = 50000
dictReads["Sample1"]["TargetBases"] = SUBSAMPLE_DEPTH * 50000
dictReads["Sample1"]["MinBases"] = MIN_DEPTH * 50000

dictReads["Sample2"] = {}
dictReads["Sample2"]["LR"] = os.path.join(dir.test_fastqs, "test_long_reads.fastq.gz")
if not AUTO:
    dictReads["Sample2"]["MinChromLength"] = 100000
dictReads["Sample2"]["TargetBases"] = SUBSAMPLE_DEPTH * 100000
dictReads["Sample2"]["MinBases"] = MIN_DEPTH * 100000

SAMPLES = ["Sample1", "Sample2"]


##############################
# Import rules and functions
##############################

# kmc - needs to be included due to the chckpointing
include: os.path.join("rules", "processing", "estimate_chromosome.smk")

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


# for seqkit
include: os.path.join("rules", "processing", "coverage.smk")
# assembly
include: os.path.join("rules", "assembly", "assemble.smk")
# extract chrom
include: os.path.join("rules", "processing", "extract_fastas.smk")
# checkpoint
# needs its own rules for long
include: os.path.join("rules", "completeness", "aggregate_long.smk")
# plassembler long and info
include: os.path.join("rules", "assembly", "plassembler_long.smk")
include: os.path.join("rules", "processing", "combine_plassembler_info.smk")


# checkpoint here for completeness

## medaka vs no medaka
if config.args.no_medaka is False:  # standard - uses Medaka

    # need long read polish files regardless
    include: os.path.join("rules", "polishing", "long_read_polish.smk")
    include: os.path.join("rules", "polishing", "long_read_polish_incomplete.smk")

    # dnaapler or dnaapler custom
    if config.args.dnaapler_custom_db == "none":  # standard - no custom

        include: os.path.join("rules", "reorientation", "dnaapler.smk")

    else:

        include: os.path.join("rules", "reorientation", "dnaapler_custom.smk")
    # finalse & pyrodigal
    include: os.path.join("rules", "finalise", "select_best_assembly_long.smk")


# --no_medaka is chosen
else:
    # dnaapler or dnaapler custom
    if config.args.dnaapler_custom_db == "none":  # standard - no custom

        include: os.path.join("rules", "reorientation", "dnaapler_no_medaka.smk")

    else:

        include: os.path.join("rules", "reorientation", "dnaapler_custom_no_medaka.smk")
    # finalse & pyrodigal
    include: os.path.join("rules", "finalise", "select_best_assembly_long_no_medaka.smk")


### rule all
rule all:
    input:
        TargetFilesLong,

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


configfile: os.path.join(workflow.basedir, "../", "config", "config.yaml")


config = ap.AttrMap(config)


# directories
CHECKDB = True  # to check db installations inside directories.smk


include: os.path.join("rules", "preflight", "directories.smk")
# functions
include: os.path.join("rules", "preflight", "functions.smk")


# check db
# from functions.smk
check_db(dir.plassemblerdb)


# samples
include: os.path.join("rules", "preflight", "samples.smk")
# targets
include: os.path.join("rules", "preflight", "targets_hybrid.smk")


### from config files
#  input as csv
INPUT = config.args.input
DATADIR = config.args.datadir
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

# for hybracter hybrid
if config.args.single is False:
    dictReads = parseSamples(INPUT, False, SUBSAMPLE_DEPTH, DATADIR, MIN_DEPTH, AUTO)  # long flag false
    SAMPLES = list(dictReads.keys())
# for hybracter hybrid-single
else:
    dictReads = {}
    dictReads[config.args.sample] = {}
    dictReads[config.args.sample]["LR"] = config.args.longreads
    dictReads[config.args.sample]["MinChromLength"] = config.args.chromosome
    dictReads[config.args.sample]["TargetBases"] = SUBSAMPLE_DEPTH * config.args.chromosome
    dictReads[config.args.sample]["MinBases"] = MIN_DEPTH * config.args.chromosome
    dictReads[config.args.sample]["R1"] = config.args.short_one
    dictReads[config.args.sample]["R2"] = config.args.short_two
    SAMPLES = [config.args.sample]


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


# for short read polishing --careful
include: os.path.join("rules", "processing", "coverage.smk")

# kmc - needs to be included due to the checkpointing
include: os.path.join("rules", "processing", "estimate_chromosome.smk")

### medaka vs no medaka
# default - medaka will be run
if config.args.no_medaka is False:

    # checkpoint here for completeness
    # need long read polish files regardless
    include: os.path.join("rules", "polishing", "long_read_polish.smk")
    include: os.path.join("rules", "polishing", "long_read_polish_incomplete.smk")

    # dnaapler or dnaapler custom
    if config.args.dnaapler_custom_db == "none":  # standard - no custom

        include: os.path.join("rules", "reorientation", "dnaapler.smk")

    else:

        include: os.path.join("rules", "reorientation", "dnaapler_custom.smk")
    #  polish the assemblies
    include: os.path.join("rules", "polishing", "short_read_polish.smk")
    include: os.path.join("rules", "polishing", "short_read_polish_incomplete.smk")
    # plassembler
    include: os.path.join("rules", "assembly", "plassembler.smk")
    include: os.path.join("rules", "processing", "combine_plassembler_info.smk")

    if config.args.no_pypolca is False:

        include: os.path.join("rules", "polishing", "short_read_pypolca.smk")
    # ale
    include: os.path.join("rules", "assess", "assess_complete.smk")
    include: os.path.join("rules", "assess", "assess_incomplete.smk")
    # finalse
    include: os.path.join("rules", "finalise", "select_best_assembly.smk")


# --no_medaka
else:
    # dnaapler or dnaapler custom
    if config.args.dnaapler_custom_db == "none":  # standard - no custom

        include: os.path.join("rules", "reorientation", "dnaapler_no_medaka.smk")

    else:

        include: os.path.join("rules", "reorientation", "dnaapler_custom_no_medaka.smk")
    #  polish the assemblies with short reads
    include: os.path.join("rules", "polishing", "short_read_polish_no_medaka.smk")
    include: os.path.join("rules", "polishing", "short_read_polish_incomplete_no_medaka.smk")
    # plassembler
    include: os.path.join("rules", "assembly", "plassembler.smk")
    include: os.path.join("rules", "processing", "combine_plassembler_info.smk")

    if config.args.no_pypolca is False:

        include: os.path.join("rules", "polishing", "short_read_pypolca_no_medaka.smk")
    # ale
    include: os.path.join("rules", "assess", "assess_complete_no_medaka.smk")
    include: os.path.join("rules", "assess", "assess_incomplete_no_medaka.smk")
    # finalse
    include: os.path.join("rules", "finalise", "select_best_assembly_no_medaka.smk")


### rule all
rule all:
    input:
        TargetFilesHybrid,

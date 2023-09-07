
import os
import glob
import attrmap as ap
import attrmap.utils as au

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
configfile: os.path.join(workflow.basedir, '../', 'config', 'config_hybrid.yaml')
config = ap.AttrMap(config)


### from config files
#  input as csv
CSV = config.input
OUTPUT = config.output
THREADS = config.threads


# Config

PLASMIDS = config['plasmids']
MIN_LENGTH = config['min_length']
MIN_QUALITY = config['min_quality']
POLCA_FLAG = config['polca']
NO_POLISH = config['no_polish']
MEDAKA_MODEL = config['medakaModel']
FLYE_MODEL = config['flyeModel']

# define functions
# get long reads
def get_input_lr_fastqs(wildcards):
    return dictReads[wildcards.sample]["LR"]

# get min chrom length (define chrom size)
def getMinChromLength(wildcards):
    return dictReads[wildcards.sample]["MinChromLength"]


### Include Directories
include: "rules/directories.smk"

# Parse the samples and read files
include: "rules/samples.smk"
dictReads = parseSamples(CSV, False)
SAMPLES = list(dictReads.keys())

wildcard_constraints:
    sample= '|'.join([re.escape(x) for x in SAMPLES])

##############################
# Import rules and functions
##############################

# import targets
include: "rules/targets.smk"

# qc
include: "rules/qc.smk"
# assembly
include: "rules/assemble.smk"
# get flye stats and combine across all runs
include: "rules/assembly_statistics.smk"
# extract chrom
include: "rules/extract_fastas.smk"

# checkpoint
include: "rules/aggregate.smk"

# checkpoint here for completeness
# need long read polish files regardless
include: "rules/long_read_polish.smk"
include: "rules/long_read_incomplete.smk"

#  polish the assemblies
include: "rules/short_read_polish.smk"
include: "rules/short_read_polish_incomplete.smk"

if POLCA_FLAG == True:
    include: "rules/short_read_polca.smk"

# plassembler if PLASMIDS true
if PLASMIDS is True:
    include: "rules/plassembler.smk"
    include: "rules/combine_plassembler_info.smk"

### need to get final chromosome and plasmids

rule all:
    input:
        TargetFilesHybrid


### DEFAULT CONFIG FILE
import os
configfile: os.path.join(workflow.basedir, '../', 'config', 'config.yaml')

### from config files
#  input as csv
CSV = config['input']
OUTPUT = config['output']
THREADS = config['threads']

# snakemake params 
BigJobMem = config["BigJobMem"]
BigJobCpu = config["BigJobCpu"]
SmallJobMem = config["SmallJobMem"]
SmallJobCpu = config["SmallJobCpu"]

SmallTime = config["SmallTime"]
BigTime = config["BigTime"]
MediumTime = config["MediumTime"]

# plassembler DB
PlassemblerDatabase = config["plassemblerDatabase"]

# LR Only flag 
PLASMIDS = config['plasmids']
MIN_LENGTH = config['min_length']
MIN_QUALITY = config['min_quality']
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
dictReads = parseSamples(CSV, True)
SAMPLES = list(dictReads.keys())

wildcard_constraints:
    sample= '|'.join([re.escape(x) for x in SAMPLES])

##############################
# Import rules and functions
##############################

# import targets
include: "rules/targets_long.smk"

# qc
include: "rules/qc.smk"
# assembly
include: "rules/assemble.smk"
# get flye stats and combine across all runs
include: "rules/assembly_statistics.smk"
# extract chrom
include: "rules/extract_fastas.smk"

include: "rules/aggregate.smk"

# checkpoint here for completeness

include: "rules/long_read_polish.smk"
include: "rules/long_read_incomplete.smk"


# plassembler if PLASMIDS true
if PLASMIDS is True:
    include: "rules/plassembler.smk"
    include: "rules/combine_plassembler_info.smk"

### need to get final chromosome and plasmids

rule all:
    input:
        TargetFilesLong

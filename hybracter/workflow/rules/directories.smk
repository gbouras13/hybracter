"""
Ensures consistent variable names and file locations for the pipeline.
"""

### OUTPUT DIRs
FLAGS = os.path.join(OUTPUT, 'FLAGS')
PROCESSING = os.path.join(OUTPUT, 'PROCESSING')
RESULTS = os.path.join(OUTPUT, 'RESULTS')

DELETE = os.path.join(OUTPUT, 'DELETE_LOGS')

# qc.smk
QC = os.path.join(OUTPUT, 'QC')

# assemble.smk
ASSEMBLIES = os.path.join(PROCESSING, 'ASSEMBLIES')

# assembly_statistics.smk 
ASSEMBLY_STATISTICS = os.path.join(PROCESSING, 'ASSEMBLY_STATISTICS')
ASSEMBLY_SUMMARY = os.path.join(RESULTS, 'ASSEMBLY_SUMMARY')

#extract_fastas.smk
CHROMOSOME_PRE_POLISH = os.path.join(PROCESSING, 'CHROMOSOME_PRE_POLISH')
FLYE_PLASMIDS = os.path.join(PROCESSING, 'FLYE_PLASMIDS')

COMPLETENESS_FLAG = os.path.join(PROCESSING, 'COMPLETENESS_FLAG')

INCOMPLETE_PRE_POLISH = os.path.join(PROCESSING, 'INCOMPLETE_PRE_POLISH')
AGGREGATED = os.path.join(AGGREGATED, 'INCOMPLETE_PRE_POLISH')



# long_read_polish.smk
MEDAKA = os.path.join(PROCESSING, 'MEDAKA')
MEDAKA_RD_1 = os.path.join(MEDAKA, 'MEDAKA_RD_1')
MEDAKA_RD_2 = os.path.join(MEDAKA, 'MEDAKA_RD_2')
DNAAPLER = os.path.join(PROCESSING, 'DNAAPLER')

# short_read_polish.smk 
FASTP = os.path.join(PROCESSING, 'FASTP')
BWA = os.path.join(PROCESSING, 'BWA')
POLYPOLISH = os.path.join(PROCESSING, 'POLYPOLISH')

# short_read_polca.smk
POLCA = os.path.join(PROCESSING, 'POLCA')


# PLASSEMBLER DIR
PLASSEMBLER = os.path.join(PROCESSING, 'PLASSEMBLER')
PLASSEMBLER_FASTAS = os.path.join(OUTPUT, 'PLASSEMBLER_FASTAS')
PLASSEMBLER_SUMMARIES = os.path.join(OUTPUT, 'PLASSEMBLER_SUMMARIES')




# SUMMARY_OUT = os.path.join(OUTPUT, 'SUMMARY')
# CHROMOSOME_PRE_POLISH = os.path.join(PROCESSING, 'CHROMOSOME_PRE_POLISH')
# CHROMOSOME_POST_POLISHING = os.path.join(OUTPUT, 'CHROMOSOME_POST_POLISHING')
# PLASMIDS = os.path.join(PROCESSING, 'PLASMIDS')
# ASSEMBLY_INFO = os.path.join(OUTPUT, 'ASSEMBLY_INFO')
# PLASMID_COVERAGE = os.path.join(OUTPUT, 'ASSEMBLY_INFO')






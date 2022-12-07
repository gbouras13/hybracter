"""
Ensures consistent variable names and file locations for the pipeline.
"""

### OUTPUT DIRs
FLAGS = os.path.join(OUTPUT, 'FLAGS')
PROCESSING = os.path.join(OUTPUT, 'PROCESSING')
RESULTS = os.path.join(OUTPUT, 'RESULTS')
DELETE = os.path.join(OUTPUT, 'DELETE_LOGS')
COMPLETE = os.path.join(PROCESSING, 'COMPLETE')
INCOMPLETE = os.path.join(PROCESSING, 'INCOMPLETE')

# qc.smk
QC = os.path.join(OUTPUT, 'QC')

# assemble.smk
ASSEMBLIES = os.path.join(PROCESSING, 'ASSEMBLIES')

# assembly_statistics.smk 
ASSEMBLY_STATISTICS = os.path.join(PROCESSING, 'ASSEMBLY_STATISTICS')
ASSEMBLY_SUMMARY = os.path.join(RESULTS, 'ASSEMBLY_SUMMARY')

#extract_fastas.smk
CHROMOSOME_PRE_POLISH = os.path.join(PROCESSING, 'CHROMOSOME_PRE_POLISH')
INCOMPLETE_PRE_POLISH = os.path.join(PROCESSING, 'INCOMPLETE_PRE_POLISH')
COMPLETENESS_FLAG = os.path.join(PROCESSING, 'COMPLETENESS_FLAG')





# aggregation dirs
AGGREGATE_LONG_READ_POLISH = os.path.join(FLAGS, 'AGGREGATE_LONG_READ_POLISH')
AGGREGATE_SHORT_READ_POLISH = os.path.join(FLAGS, 'AGGREGATE_SHORT_READ_POLISH')
AGGREGATE_POLCA_POLISH = os.path.join(FLAGS, 'AGGREGATE_POLCA_POLISH')

# long_read_polish.smk
MEDAKA = os.path.join(COMPLETE, 'MEDAKA')
MEDAKA_RD_1 = os.path.join(MEDAKA, 'MEDAKA_RD_1')
MEDAKA_RD_2 = os.path.join(MEDAKA, 'MEDAKA_RD_2')
DNAAPLER = os.path.join(COMPLETE, 'DNAAPLER')

# long_Read_polish_incomplete.smk
MEDAKA_INCOMPLETE = os.path.join(INCOMPLETE, 'MEDAKA_INCOMPLETE')

# short_read_polish.smk 
FASTP = os.path.join(PROCESSING, 'FASTP')
BWA = os.path.join(COMPLETE, 'BWA')
POLYPOLISH = os.path.join(COMPLETE, 'POLYPOLISH')

# short_read_polish_incomplete.smk 
BWA_INCOMPLETE = os.path.join(INCOMPLETE, 'BWA_INCOMPLETE')
POLYPOLISH_INCOMPLETE = os.path.join(INCOMPLETE, 'POLYPOLISH_INCOMPLETE')

# short_read_polca.smk (has complete and incomplete)
POLCA = os.path.join(COMPLETE, 'POLCA')
POLCA_INCOMPLETE = os.path.join(INCOMPLETE, 'POLCA_INCOMPLETE')


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






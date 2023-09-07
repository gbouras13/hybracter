"""
All target output files are declared here
"""

# need long read polish aggr file regardless even if no polish selected 
# because of dnaapler
long_read_polish_files = os.path.join(FLAGS, "aggr_long_read_polish.txt")


# Polca
if POLCA_FLAG == True:
    polca_files = os.path.join(FLAGS, "aggr_polca.txt")
else:
    polca_files = []

# plassembler

if PLASMIDS is True:
    plassembler_files = [os.path.join(FLAGS, "aggr_plassembler.txt"),
    os.path.join(FLAGS, "aggr_combine_plassembler_info.txt")]
else:
    plassembler_files = []

# Preprocessing files
TargetFilesHybrid = [
    os.path.join(FLAGS, "aggr_qc.txt"),
    os.path.join(FLAGS, "aggr_assemble.txt"),
    os.path.join(FLAGS, "aggr_assembly_statistics.txt"),
    os.path.join(FLAGS, "aggr_short_read_polish.txt"),
    os.path.join(FLAGS, "aggr_long_read_polish.txt"),
    polca_files,
    plassembler_files
]



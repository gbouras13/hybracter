"""
All target output files are declared here
"""

# need long read polish aggr file regardless even if no polish selected 
long_read_polish_files = os.path.join(FLAGS, "aggr_long_read_polish.txt")

# Short reads
if LR_ONLY == False:
    short_read_polish_files = os.path.join(FLAGS, "aggr_short_read_polish.txt")
else:
    short_read_polish_files = []

# Polca
if POLCA_FLAG == True:
    polca_files = os.path.join(FLAGS, "aggr_polca.txt")
else:
    polca_files = []

# plassembler

if LR_ONLY is False and PLASMIDS is True:
    plassembler_files = [os.path.join(FLAGS, "aggr_plassembler.txt"),
    os.path.join(FLAGS, "aggr_combine_plassembler_info.txt")]
else:
    plassembler_files = []




# Preprocessing files
TargetFiles = [
    os.path.join(FLAGS, "aggr_qc.txt"),
    os.path.join(FLAGS, "aggr_assemble.txt"),
    os.path.join(FLAGS, "aggr_assembly_statistics.txt"),
    long_read_polish_files,
    short_read_polish_files,
    polca_files,
    plassembler_files
]

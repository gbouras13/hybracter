"""
All target output files are declared here
"""

# need long read polish aggr file regardless even if no polish selected 
# because of dnaapler

# plassembler

if PLASMIDS is True:
    plassembler_files = [os.path.join(FLAGS, "aggr_plassembler.txt"),
    os.path.join(FLAGS, "aggr_combine_plassembler_info.txt")]
else:
    plassembler_files = []




# Preprocessing files
TargetFilesLong = [
    os.path.join(FLAGS, "aggr_qc.txt"),
    os.path.join(FLAGS, "aggr_assemble.txt"),
    os.path.join(FLAGS, "aggr_assembly_statistics.txt"),
    os.path.join(FLAGS, "aggr_long_read_polish.txt"),
    plassembler_files
]

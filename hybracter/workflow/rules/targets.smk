"""
All target output files are declared here
"""

# SR
if LR_ONLY == False:
    short_read_polish_files = os.path.join(FLAGS, "aggr_short_read_polish.txt")
else:
    short_read_polish_files = []

# Polca
if POLCA_FLAG == False:
    polca_files = os.path.join(FLAGS, "aggr_polca.txt")
else:
    polca_files = []

# plassembler

if LR_ONLY == False:
    plassembler_files = [os.path.join(LOGS, "aggr_plasembler.txt"),
    os.path.join(LOGS, "aggr_combine_plassembler_info.txt")]
else:
    plassembler_files = []



# Preprocessing files
TargetFiles = [
    os.path.join(FLAGS, "aggr_qc.txt"),
    os.path.join(FLAGS, "aggr_assemble.txt"),
    os.path.join(FLAGS, "aggr_assembly_statistics.txt"),
    os.path.join(FLAGS, "aggr_chr_plas.txt"),
    #os.path.join(FLAGS, "aggr_long_read_polish.txt"),
    short_read_polish_files,
    polca_files,
    plassembler_files



]

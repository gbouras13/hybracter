"""
All target output files are declared here
"""

# need long read polish aggr file regardless even if no polish selected
# because of dnaapler
"""
polca depends on --no_pypolca flag
"""

# Pypolca
if config.args.no_pypolca == True:
    pypolca_files = []

else:
    pypolca_files = os.path.join(dir.out.flags, "aggr_pypolca.flag")


"""
hybrid
"""

TargetFilesHybrid = [
    os.path.join(dir.out.flags, "aggr_long_qc.flag"),
    os.path.join(dir.out.flags, "aggr_short_qc.flag"),
    os.path.join(dir.out.flags, "aggr_kmc.flag"),
    os.path.join(dir.out.flags, "aggr_seqkit_short.flag"),
    os.path.join(dir.out.flags, "aggr_seqkit_long.flag"),
    os.path.join(dir.out.flags, "aggr_assemble.flag"),
    os.path.join(dir.out.flags, "aggr_short_read_polish.flag"),
    os.path.join(dir.out.flags, "aggr_long_read_polish.flag"),
    pypolca_files,
    os.path.join(dir.out.flags, "aggr_plassembler.flag"),
    os.path.join(dir.out.flags, "aggr_combine_plassembler_info.flag"),
    os.path.join(dir.out.flags, "aggr_ale.flag"),
    os.path.join(dir.out.flags, "aggr_final.flag"),
]

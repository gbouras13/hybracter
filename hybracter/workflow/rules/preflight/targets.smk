"""
All target output files are declared here
"""

# need long read polish aggr file regardless even if no polish selected 
# because of dnaapler


# Polca
if config.no_polca == True:
    polca_files = []
    
else:
    polca_files = os.path.join(dir.out.flags, "aggr_polca.flag")

# plassembler

if config.plasmids is True:
    plassembler_files = [os.path.join(dir.out.flags, "aggr_plassembler.flag"),
    os.path.join(dir.out.flags, "aggr_combine_plassembler_info.flag")]
else:
    plassembler_files = []



TargetFilesHybrid = [
    os.path.join(dir.out.flags, "aggr_qc.flag"),
    os.path.join(dir.out.flags, "aggr_assemble.flag"),
    os.path.join(dir.out.flags, "aggr_assembly_statistics.flag"),
    os.path.join(dir.out.flags, "aggr_short_read_polish.flag"),
    os.path.join(dir.out.flags, "aggr_long_read_polish.flag"),
    polca_files,
    plassembler_files
]


"""
long 
"""

# need long read polish aggr file regardless even if no polish selected 
# because of dnaapler

TargetFilesLong = [
    os.path.join(dir.out.flags, "aggr_qc.flag"),
    os.path.join(dir.out.flags, "aggr_assemble.flag"),
    os.path.join(dir.out.flags, "aggr_assembly_statistics.flag"),
    os.path.join(dir.out.flags, "aggr_long_read_polish.flag"),
    plassembler_files
]



"""
download 
"""

TargetFilesDownload = [
    os.path.join(dir.plassemblerdb  ,'plsdb.msh'),
    os.path.join(dir.plassemblerdb, 'plsdb.tsv')
]


"""
ale
"""

TargetFilesAle = [os.path.join(dir.env, "ALE-master", "src", "ALE")]
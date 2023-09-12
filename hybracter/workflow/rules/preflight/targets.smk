"""
All target output files are declared here
"""

# need long read polish aggr file regardless even if no polish selected 
# because of dnaapler
"""
polca depends on --no_polca flag
"""

# Polca
if config.args.no_polca == True:
    polca_files = []
    
else:
    polca_files = os.path.join(dir.out.flags, "aggr_polca.flag")


"""
plassembler
"""
# plassembler files
plassembler_files = [
    os.path.join(dir.out.flags, "aggr_plassembler.flag"),
    os.path.join(dir.out.flags, "aggr_combine_plassembler_info.flag")
    ]


"""
hybrid
"""

TargetFilesHybrid = [
    os.path.join(dir.out.flags, "aggr_qc.flag"),
    os.path.join(dir.out.flags, "aggr_assemble.flag"),
    os.path.join(dir.out.flags, "aggr_assembly_statistics.flag"),
    os.path.join(dir.out.flags, "aggr_short_read_polish.flag"),
    os.path.join(dir.out.flags, "aggr_long_read_polish.flag"),
    polca_files,
    plassembler_files,
    os.path.join(dir.out.flags, "aggr_ale.flag") # ale
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
    os.path.join(dir.plassemblerdb,'plsdb.msh'),
    os.path.join(dir.plassemblerdb, 'plsdb.tsv')
]


"""
ale
"""

TargetFilesAle = [os.path.join(dir.env, "ALE-master", "src", "ALE")]
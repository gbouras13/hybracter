"""
All target output files are declared here
"""


"""
long 
"""


# plassembler files
if config.args.no_plasmids == True:
    plassembler_files = []
else:
    plassembler_files = 
    [ 
        os.path.join(dir.out.flags, "aggr_plassembler.flag"),
        os.path.join(dir.out.flags, "aggr_combine_plassembler_info.flag")
    ]

TargetFilesLong = [
    os.path.join(dir.out.flags, "aggr_qc.flag"),
    os.path.join(dir.out.flags, "aggr_assemble.flag"),
    os.path.join(dir.out.flags, "aggr_assembly_statistics.flag"),
    os.path.join(dir.out.flags, "aggr_long_read_polish.flag"),
    plassembler_files,
    os.path.join(dir.out.flags, "aggr_final.flag")
]


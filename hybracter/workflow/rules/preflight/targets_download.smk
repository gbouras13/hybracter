"""
All target output files are declared here
"""

TargetFilesDownload = [
    os.path.join(dir.plassemblerdb, "plsdb_2023_11_03_v2.msh"),
    os.path.join(dir.plassemblerdb, "plsdb_2023_11_03_v2.tsv"),
    os.path.join(dir.plassemblerdb, "cleanup.flag"),
]

if MEDAKA_DOWNLOAD is True:
    TargetFilesDownload.append(os.path.join(dir.plassemblerdb, "medaka.flag"))

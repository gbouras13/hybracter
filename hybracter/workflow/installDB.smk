"""
Snakefile for downloading transcriptomes 
"""
import os

# load default config
configfile: os.path.join(workflow.basedir, '../', 'config', 'config.yaml')

# config
DatabaseDir = config["plassemblerDatabase"]

# snakemake params 
BigJobMem = config["BigJobMem"]
BigJobCpu = config["BigJobCpu"]
SmallJobMem = config["SmallJobMem"]
SmallJobCpu = config["SmallJobCpu"]

SmallTime = config["SmallTime"]
BigTime = config["BigTime"]
MediumTime = config["MediumTime"]


rule all:
    input:
        os.path.join(DatabaseDir,'plsdb.msh'),
        os.path.join(DatabaseDir, 'plsdb.tsv')


rule download_db:
    """Rule to Download plassembler db."""
    params:
        DatabaseDir
    conda:
        os.path.join( 'envs','plassembler.yaml')
    output:
        os.path.join(DatabaseDir,'plsdb.msh'),
        os.path.join(DatabaseDir, 'plsdb.tsv')
    threads:
        1
    resources:
        mem_mb=SmallJobMem,
        time=BigTime
    shell:
        """
        install_database.py -d {params[0]}
        """

# download plassembler db
rule download_db:
    """Rule to Download plassembler db."""
    params:
        config.databases
    conda:
        os.path.join(dir.env,'plassembler.yaml')
    output:
        os.path.join(config.databases,'plsdb.msh'),
        os.path.join(config.databases, 'plsdb.tsv')
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        install_database.py -d {params[0]}
        """

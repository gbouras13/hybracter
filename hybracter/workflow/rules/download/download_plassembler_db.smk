# download plassembler db
rule download_db:
    """Rule to Download plassembler db."""
    params:
        db = dir.plassemblerdb,
        intermediate_outputs = dir.out.base
    conda:
        os.path.join(dir.env, "plassembler.yaml")
    output:
        mash = os.path.join(dir.plassemblerdb, "plsdb.msh"),
        tsv = os.path.join(dir.plassemblerdb, "plsdb.tsv"),
        flag = os.path.join(dir.plassemblerdb, "cleanup.flag")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        plassembler download -d {params.db}
        rm -r {params.intermediate_outputs}
        touch {output.flag}
        """

# download plassembler db
rule download_db:
    """Rule to Download plassembler db."""
    params:
        db=dir.plassemblerdb,
        intermediate_outputs=dir.out.base,
    conda:
        os.path.join(dir.env, "plassembler.yaml")
    output:
        mash=os.path.join(dir.plassemblerdb, "plsdb_2023_11_03_v2.msh"),
        tsv=os.path.join(dir.plassemblerdb, "plsdb_2023_11_03_v2.tsv"),
        flag=os.path.join(dir.plassemblerdb, "cleanup.flag"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.sml.cpu
    shell:
        """
        plassembler download -d {params.db} -f
        rm -r {params.intermediate_outputs}
        touch {output.flag}
        """

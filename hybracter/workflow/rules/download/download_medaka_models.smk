# download medaka models
rule download_medaka_models:
    """Rule to Download plassembler db."""
    conda:
        (
            os.path.join(dir.env, "medaka_mac.yaml")
            if MAC
            else os.path.join(dir.env, "medaka.yaml")
        )
    output:
        flag=os.path.join(dir.plassemblerdb, "medaka.flag"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.sml.cpu
    shell:
        """
        medaka tools download_models
        touch {output.flag}
        """

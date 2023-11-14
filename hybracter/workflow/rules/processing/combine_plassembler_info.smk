
rule combine_plassembler_info:
    input:
        os.path.join(dir.out.flags, "aggr_plassembler.flag"),  # need the flag as in the put - to be run after plassembler
    output:
        out=os.path.join(
            dir.out.plassembler_all_summary, "plassembler_assembly_info.tsv"
        ),
    params:
        summary_dir=dir.out.plassembler_individual_summaries,  # all the samples where plassembler ran (complete)
    conda:
        os.path.join(dir.env, "scripts.yaml")
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    script:
        os.path.join(dir.scripts, "combine_plassembler_info.py")


rule aggr_combine_plassembler_info:
    """Aggregate."""
    input:
        os.path.join(dir.out.plassembler_all_summary, "plassembler_assembly_info.tsv"),
    output:
        flag=os.path.join(dir.out.flags, "aggr_combine_plassembler_info.flag"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        touch {output.flag}
        """

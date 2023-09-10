
rule combine_plassembler_info:
    input:
        summary_dir = directory(dir.out.plassembler_individual_summaries), # all the samples where plassembler ran (complete)
        flag_dir = directory(dir.out.plassembler_individual_summaries) # all the samples where plassembler was skipped (incomplete)
    output:
        out = os.path.join(dir.out.plassembler_all_summary,"plassembler_assembly_info.txt")
    conda:
        os.path.join(dir.env,'scripts.yaml')
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    script:
        os.path.join(dir.scripts,  'combine_plassembler_info.py')


rule aggr_combine_plassembler_info:
    """Aggregate."""
    input:
        os.path.join(dir.out.plassembler_all_summary,"plassembler_assembly_info.txt")
    output:
        flag = os.path.join(dir.out.flags, "aggr_combine_plassembler_info.flag")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        touch {output.flag}
        """

rule add_sample_plassembler:
    input:
        inp = os.path.join(dir.out.plassembler_individual_summaries , "{sample}.tsv")
    output:
        out = os.path.join(dir.out.plassembler_individual_summaries ,"{sample}_with_sample.tsv")
    conda:
        os.path.join(dir.env,'scripts.yaml')
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    script:
        '../scripts/add_sample_plassembler.py'


rule combine_plassembler_info:
    input:
        summaries = expand(os.path.join(dir.out.plassembler_individual_summaries , "{sample}_with_sample.tsv"), sample = SAMPLES)
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
        '../scripts/combine_plassembler_info.py'


rule aggr_combine_plassembler_info:
    """Aggregate."""
    input:
        os.path.join(dir.out.plassembler_all_summary,"plassembler_assembly_info.txt")
    output:
        os.path.join(dir.out.flags, "aggr_combine_plassembler_info.flag")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        touch {output[0]}
        """
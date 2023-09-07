
rule add_sample_plassembler:
    input:
        inp = os.path.join(PLASSEMBLER_SUMMARIES, "{sample}.tsv")
    output:
        out = os.path.join(PLASSEMBLER_SUMMARIES,"{sample}_with_sample.tsv")
    conda:
        os.path.join('..', 'envs','scripts.yaml')
    resources:
        mem_mb=config.resources.small.mem,
        time=config.resources.small.time
    threads:
        config.resources.small.cpu
    script:
        '../scripts/add_sample_plassembler.py'


rule combine_plassembler_info:
    input:
        summaries = expand(os.path.join(PLASSEMBLER_SUMMARIES, "{sample}_with_sample.tsv"), sample = SAMPLES)
    output:
        out = os.path.join(PLASSEMBLER_SUMMARY_ALL,"plassembler_assembly_info.txt")
    conda:
        os.path.join('..', 'envs','scripts.yaml')
    resources:
        mem_mb=config.resources.small.mem,
        time=config.resources.small.time
    threads:
        config.resources.small.cpu
    script:
        '../scripts/combine_plassembler_info.py'


rule aggr_combine_plassembler_info:
    """Aggregate."""
    input:
        os.path.join(PLASSEMBLER_SUMMARY_ALL,"plassembler_assembly_info.txt")
    output:
        os.path.join(FLAGS, "aggr_combine_plassembler_info.txt")
    resources:
        mem_mb=config.resources.small.mem,
        time=config.resources.small.time
    threads:
        config.resources.small.cpu
    shell:
        """
        touch {output[0]}
        """
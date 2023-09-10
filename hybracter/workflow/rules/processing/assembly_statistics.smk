
rule extract_summary_assembly_stats:
    input:
        os.path.join(dir.out.assemblies,"{sample}", "assembly_info.txt")
    output:
        os.path.join(dir.out.assembly_statistics,"{sample}_clean_assembly_info.csv"),
        os.path.join(dir.out.assembly_statistics,"{sample}_summary.csv")
    conda:
        os.path.join(dir.env,'scripts.yaml')
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    script:
        os.path.join(dir.scripts,  'extract_summary_assembly_stats.py')

rule combine_summary_stats:
    input:
        summaries = expand(os.path.join(dir.out.assembly_statistics,"{sample}_summary.csv"), sample = SAMPLES)
    output:
        out = os.path.join(dir.out.assembly_summary,"total_assembly_summary.csv")
    conda:
        os.path.join(dir.env,'scripts.yaml')
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    script:
        os.path.join(dir.scripts,  'combine_summary_assembly_stats.py')


rule aggr_statistics:
    """Aggregate."""
    input:
        expand(os.path.join(dir.out.assembly_statistics,"{sample}_clean_assembly_info.csv"), sample = SAMPLES),
        expand(os.path.join(dir.out.assembly_statistics,"{sample}_summary.csv"), sample = SAMPLES),
        os.path.join(dir.out.assembly_summary,"total_assembly_summary.csv")
    output:
        flag = os.path.join(dir.out.flags, "aggr_assembly_statistics.flag")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        touch {output.flag}
        """

rule extract_summary_assembly_stats:
    input:
        os.path.join(ASSEMBLIES,"{sample}", "assembly_info.txt")
    output:
        os.path.join(ASSEMBLY_STATISTICS,"{sample}_clean_assembly_info.csv"),
        os.path.join(ASSEMBLY_STATISTICS,"{sample}_summary.csv")
    conda:
        os.path.join('..', 'envs','scripts.yaml')
    resources:
        mem_mb=SmallJobMem,
        time=SmallTime
    threads:
        1
    script:
        '../scripts/extract_summary_assembly_stats.py'

rule combine_summary_stats:
    input:
        summaries = expand(os.path.join(ASSEMBLY_STATISTICS,"{sample}_summary.csv"), sample = SAMPLES)
    output:
        out = os.path.join(ASSEMBLY_SUMMARY,"total_assembly_summary.csv")
    conda:
        os.path.join('..', 'envs','scripts.yaml')
    resources:
        mem_mb=SmallJobMem,
        time=SmallTime
    threads:
        1
    script:
        '../scripts/combine_summary_assembly_stats.py'


rule aggr_statistics:
    """Aggregate."""
    input:
        expand(os.path.join(ASSEMBLY_STATISTICS,"{sample}_clean_assembly_info.csv"), sample = SAMPLES),
        expand(os.path.join(ASSEMBLY_STATISTICS,"{sample}_summary.csv"), sample = SAMPLES),
        os.path.join(ASSEMBLY_SUMMARY,"total_assembly_summary.csv")
    output:
        os.path.join(FLAGS, "aggr_assembly_statistics.txt")
    resources:
        mem_mb=SmallJobMem,
        time=SmallTime
    threads:
        1
    shell:
        """
        touch {output[0]}
        """
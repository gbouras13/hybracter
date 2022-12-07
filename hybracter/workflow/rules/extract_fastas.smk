
def getMinChromLength(wildcards):
    return dictReads[wildcards.sample]["MinChromLength"]


checkpoint check_completeness:
    input:
        os.path.join(ASSEMBLIES,"{sample}", "assembly.fasta")
    output:
        os.path.join(COMPLETENESS_FLAG,"{sample}.txt")
    script:
        '../scripts/check_completeness.py'


rule extract_complete:
    input:
        os.path.join(ASSEMBLIES,"{sample}", "assembly.fasta"),
        os.path.join(COMPLETENESS_FLAG,"{sample}.txt")
    output:
        os.path.join(CHROMOSOME_PRE_POLISH,"{sample}.fasta")
    params:
        getMinChromLength
    conda:
        os.path.join('..', 'envs','scripts.yaml')
    resources:
        mem_mb=SmallJobMem,
        time=SmallTime
    threads:
        1
    script:
        '../scripts/extract_chromosome.py'

rule extract_incomplete:
    input:
        os.path.join(ASSEMBLIES,"{sample}", "assembly.fasta"),
        os.path.join(COMPLETENESS_FLAG,"{sample}.txt")
    output:
        os.path.join(INCOMPLETE_PRE_POLISH,"{sample}.fasta")
    params:
        getMinChromLength
    conda:
        os.path.join('..', 'envs','scripts.yaml')
    resources:
        mem_mb=SmallJobMem,
        time=SmallTime
    threads:
        1
    script:
        '../scripts/extract_incomplete.py'

# input function for the rule aggregate
def aggregate_input(wildcards):
    # decision based on content of output file
    # Important: use the method open() of the returned file!
    # This way, Snakemake is able to automatically download the file if it is generated in
    # a cloud environment without a shared filesystem.
    with checkpoints.check_completeness.get(sample=wildcards.sample).output[0].open() as f:
        if f.read().strip() == "C":
            return os.path.join(CHROMOSOME_PRE_POLISH,"{sample}.fasta")
        else:
            return os.path.join(INCOMPLETE_PRE_POLISH,"{sample}.fasta")

rule aggregate_input:
    input:
        aggregate_input
    output:
        os.path.join(AGGREGATED,"{sample}.txt")
    shell:
        """
        touch {output}
        """



rule aggr_chromosome:
    """Aggregate."""
    input:
        expand(os.path.join(AGGREGATED,"{sample}.txt"), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "aggr_chr_plas.txt")
    resources:
        mem_mb=SmallJobMem,
        time=SmallTime
    threads:
        1
    shell:
        """
        touch {output[0]}
        """







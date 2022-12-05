
def getMinChromLength(wildcards):
    return dictReads[wildcards.sample]["MinChromLength"]

rule extract_chromosome:
    input:
        os.path.join(ASSEMBLIES,"{sample}", "assembly.fasta")
    output:
        os.path.join(CHROMOSOME_PRE_POLISH,"{sample}.fasta"),
        os.path.join(FLYE_PLASMIDS,"{sample}.fasta")
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


rule aggr_chromosome:
    """Aggregate."""
    input:
        expand(os.path.join(CHROMOSOME_PRE_POLISH,"{sample}.fasta"), sample = SAMPLES)
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
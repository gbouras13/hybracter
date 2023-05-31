
checkpoint check_completeness:
    input:
        os.path.join(ASSEMBLIES,"{sample}", "assembly.fasta")
    output:
        os.path.join(COMPLETENESS_FLAG,"{sample}.txt")
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


### no aggr rule - flows into the long_read_polish rules
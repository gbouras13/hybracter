
checkpoint check_completeness:
    input:
        os.path.join(dir.out.assemblies,"{sample}", "assembly.fasta")
    output:
        os.path.join(dir.out.completeness,"{sample}.txt")
    params:
        getMinChromLength
    conda:
        os.path.join(dir.env,'scripts.yaml')
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    script:
        '../scripts/check_completeness.py'

rule extract_complete:
    input:
        os.path.join(dir.out.assemblies,"{sample}", "assembly.fasta"),
        os.path.join(dir.out.completeness,"{sample}.txt")
    output:
        os.path.join(dir.out.chrom_pre_polish,"{sample}.fasta")
    params:
        getMinChromLength
    conda:
        os.path.join(dir.env,'scripts.yaml')
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    script:
        '../scripts/extract_chromosome.py'

rule extract_incomplete:
    input:
        os.path.join(dir.out.assemblies,"{sample}", "assembly.fasta"),
        os.path.join(dir.out.completeness ,"{sample}.txt")
    output:
        os.path.join(dir.out.incomp_pre_polish,"{sample}.fasta")
    params:
        getMinChromLength
    conda:
        os.path.join(dir.env,'scripts.yaml')
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    script:
        '../scripts/extract_incomplete.py'


### no aggr rule - flows into the long_read_polish rules
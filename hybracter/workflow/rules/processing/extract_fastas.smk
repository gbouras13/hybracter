
checkpoint check_completeness:
    """
    adds checkpoint to determine whether the Flye assembly recovered the complete chromosome (or not)
    """
    input:
        fasta=os.path.join(
            dir.out.intermediate_assemblies, "{sample}", "{sample}_flye.fasta"
        ),
    output:
        completeness_check=os.path.join(dir.out.completeness, "{sample}.txt"),
    params:
        min_chrom_length=getMinChromLength,
    conda:
        os.path.join(dir.env, "scripts.yaml")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    script:
        os.path.join(dir.scripts, "check_completeness.py")


rule extract_chromosome_complete:
    """
    extracts the chromosome for complete samples
    """
    input:
        fasta=os.path.join(dir.out.assemblies, "{sample}", "assembly.fasta"),
        completeness_check=os.path.join(dir.out.completeness, "{sample}.txt"),
    output:
        fasta=os.path.join(dir.out.chrom_pre_polish, "{sample}.fasta"),
    params:
        min_chrom_length=getMinChromLength,
    conda:
        os.path.join(dir.env, "scripts.yaml")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    script:
        os.path.join(dir.scripts, "extract_chromosome.py")


rule extract_incomplete:
    """
    extracts the chromosome for complete samples
    """
    input:
        fasta=os.path.join(dir.out.assemblies, "{sample}", "assembly.fasta"),
        completeness_check=os.path.join(dir.out.completeness, "{sample}.txt"),
    output:
        fasta=os.path.join(dir.out.incomp_pre_polish, "{sample}.fasta"),
    params:
        min_chrom_length=getMinChromLength,
    conda:
        os.path.join(dir.env, "scripts.yaml")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    script:
        os.path.join(dir.scripts, "extract_incomplete.py")


### no aggr rule - it flows into the long_read_polish rules

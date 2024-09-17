
checkpoint check_completeness:
    """
    adds checkpoint to determine whether the Flye assembly recovered the complete chromosome (or not)
    """
    input:
        fasta=os.path.join(dir.out.assemblies, "{sample}", "assembly.fasta"),
        info=os.path.join(dir.out.assemblies, "{sample}", "assembly_info.txt"),
        kmc=os.path.join(dir.out.kmc,"{sample}", "{sample}_kmcLOG.txt"),
    output:
        completeness_check=os.path.join(dir.out.completeness, "{sample}.txt"),
    params:
        min_chrom_length=lambda wildcards: str(getMinChromLength(kmc_log_path=os.path.join(dir.out.kmc, f"{wildcards.sample}", f"{wildcards.sample}_kmcLOG.txt"), sample=wildcards.sample,auto=AUTO)),
    conda:
        os.path.join(dir.env, "scripts.yaml")
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    script:
        os.path.join(dir.scripts, "check_completeness.py")


rule extract_chromosome_complete:
    """
    extracts the chromosome for complete samples
    ignore_list contains any contigs > min_chrom_length that are not circular as marked by Flye 
    """
    input:
        fasta=os.path.join(dir.out.assemblies, "{sample}", "assembly.fasta"),
        info=os.path.join(dir.out.assemblies, "{sample}", "assembly_info.txt"),
        completeness_check=os.path.join(dir.out.completeness, "{sample}.txt"),
        kmc=os.path.join(dir.out.kmc,"{sample}", "{sample}_kmcLOG.txt")
    output:
        fasta=os.path.join(dir.out.chrom_pre_polish, "{sample}_chromosome.fasta"),
        ignore_list=os.path.join(dir.out.chrom_pre_polish, "{sample}_ignore_list.txt"),
    params:
        min_chrom_length=lambda wildcards: str(getMinChromLength(kmc_log_path=os.path.join(dir.out.kmc, f"{wildcards.sample}", f"{wildcards.sample}_kmcLOG.txt"), sample=wildcards.sample,auto=AUTO)),
        polypolish_flag=False,
    conda:
        os.path.join(dir.env, "scripts.yaml")
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    script:
        os.path.join(dir.scripts, "extract_chromosome.py")


rule copy_flye_intermediate_chrom_assembly:
    """
    copies the flye chromosome for to the intermediate chromosome assemblies directory
    also copies the 
    """
    input:
        fasta=os.path.join(dir.out.chrom_pre_polish, "{sample}_chromosome.fasta"),
    output:
        fasta=os.path.join(
            dir.out.intermediate_assemblies, "{sample}", "{sample}_flye.fasta"
        ),
    params:
        min_chrom_length=lambda wildcards: str(getMinChromLength(kmc_log_path=os.path.join(dir.out.kmc, f"{wildcards.sample}", f"{wildcards.sample}_kmcLOG.txt"), sample=wildcards.sample,auto=AUTO)),
    conda:
        os.path.join(dir.env, "scripts.yaml")
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        cp {input.fasta} {output.fasta} 
        """


rule concatenate_chrom_plassembler:
    """
    concatenates chrom and plassembler outputs
    """
    input:
        chrom_fasta=os.path.join(dir.out.chrom_pre_polish, "{sample}_chromosome.fasta"),
        plasmid_fasta=os.path.join(
            dir.out.plassembler, "{sample}", "plassembler_plasmids.fasta"
        ),
    output:
        combo_fasta=os.path.join(
            dir.out.chrom_pre_polish, "{sample}_chromosome_plus_plasmids.fasta"
        ),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        cat {input.chrom_fasta} {input.plasmid_fasta} > {output.combo_fasta}
        """


rule extract_incomplete:
    """
    extracts the chromosome for complete samples
    """
    input:
        fasta=os.path.join(dir.out.assemblies, "{sample}", "assembly.fasta"),
        completeness_check=os.path.join(dir.out.completeness, "{sample}.txt"),
        kmc=os.path.join(dir.out.kmc,"{sample}", "{sample}_kmcLOG.txt")
    output:
        fasta=os.path.join(dir.out.incomp_pre_polish, "{sample}.fasta"),
    params:
        min_chrom_length=lambda wildcards: str(getMinChromLength(kmc_log_path=os.path.join(dir.out.kmc, f"{wildcards.sample}", f"{wildcards.sample}_kmcLOG.txt"), sample=wildcards.sample,auto=AUTO)),
    conda:
        os.path.join(dir.env, "scripts.yaml")
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    script:
        os.path.join(dir.scripts, "extract_incomplete.py")


### no aggr rule - it flows into the long_read_polish rules

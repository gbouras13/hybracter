rule medaka_round_1:
    """
    Runs medaka round 1
    cleans up BAM and hdf for disc space
    """
    input:
        fasta=os.path.join(
            dir.out.chrom_pre_polish, "{sample}_chromosome_plus_plasmids.fasta"
        ),
        fastq=os.path.join(dir.out.qc, "{sample}_filt_trim.fastq.gz"),
    output:
        fasta=os.path.join(dir.out.medaka_rd_1, "{sample}", "consensus.fasta"),
        version=os.path.join(dir.out.versions, "{sample}", "medaka_complete.version"),
    conda:
        (
            os.path.join(dir.env, "medaka_mac.yaml")
            if MAC
            else os.path.join(dir.env, "medaka.yaml")
        )
    params:
        model=MEDAKA_MODEL,
        dir=os.path.join(dir.out.medaka_rd_1, "{sample}"),
        bam=os.path.join(dir.out.medaka_rd_1, "{sample}", "calls_to_draft.bam"),
        hdf=os.path.join(dir.out.medaka_rd_1, "{sample}", "consensus_probs.hdf"),
    resources:
        mem_mb=config.resources.big.mem,
        mem=str(config.resources.big.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.big.cpu
    benchmark:
        os.path.join(dir.out.bench, "medaka_round_1", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "medaka_round_1", "{sample}.log"),
    shell:
        """
        medaka_consensus -i {input.fastq} -d {input.fasta} -o {params.dir} -m {params.model}  -t {threads} 2> {log}
        medaka --version > {output.version}
        touch {params.bam}
        rm {params.bam}
        touch {params.hdf}
        rm {params.hdf}
        """


rule medaka_round_1_extract_intermediate_assembly:
    """
    extracts the chromosome intermediate assembly
    """
    input:
        fasta=os.path.join(dir.out.medaka_rd_1, "{sample}", "consensus.fasta"),
        completeness_check=os.path.join(dir.out.completeness, "{sample}.txt"),
        info=os.path.join(dir.out.assemblies, "{sample}", "assembly_info.txt"),
        kmc=os.path.join(dir.out.kmc,"{sample}", "{sample}_kmcLOG.txt")
    output:
        fasta=os.path.join(
            dir.out.intermediate_assemblies, "{sample}", "{sample}_medaka_rd_1.fasta"
        ),
        ignore_list=os.path.join(dir.out.medaka_rd_1, "{sample}_ignore_list.txt"),
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


rule compare_assemblies_medaka_round_1:
    """
    compare chromosome assemblies between medaka and pre-polished chromosome
    take the one that is dnaaplered in finalise folder so we can compare properly
    """
    input:
        reference=os.path.join(
            dir.out.dnaapler, "{sample}_pre_chrom", "{sample}_reoriented.fasta"
        ),
        assembly=os.path.join(dir.out.dnaapler, "{sample}", "{sample}_reoriented.fasta"),
    output:
        diffs=os.path.join(
            dir.out.differences, "{sample}", "medaka_round_1_vs_pre_polish.txt"
        ),
    conda:
        os.path.join(dir.env, "scripts.yaml")
    params:
        reference_polishing_round="pre_polish",
        query_polishing_round="medaka_round_1",
    resources:
        mem_mb=config.resources.med.mem,
        mem=str(config.resources.med.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.sml.cpu
    script:
        os.path.join(dir.scripts, "compare_assemblies.py")


rule concatenate_intermediate_dnaapler_with_plassembler:
    """
    concatenates dnaapler (medaka rd 1) and plassembler outputs
    """
    input:
        chrom_fasta=os.path.join(
            dir.out.dnaapler, "{sample}", "{sample}_reoriented.fasta"
        ),
        plasmid_fasta=os.path.join(
            dir.out.plassembler, "{sample}", "plassembler_plasmids.fasta"
        ),
    output:
        combo_fasta=os.path.join(
            dir.out.dnaapler,
            "{sample}",
            "{sample}_reoriented_chrom_and_plasmids.fasta",
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


rule medaka_round_2:
    """
    runs medaka round 2 on the reoriented genome
    note: can't run compare_assemblies.py as the order of the genome has completely changed.
    cleans up BAM and hdf for disc space
    """
    input:
        fasta=os.path.join(
            dir.out.dnaapler,
            "{sample}",
            "{sample}_reoriented_chrom_and_plasmids.fasta",
        ),
        fastq=os.path.join(dir.out.qc, "{sample}_filt_trim.fastq.gz"),
    output:
        fasta=os.path.join(dir.out.medaka_rd_2, "{sample}", "consensus.fasta"),
    conda:
        (
            os.path.join(dir.env, "medaka_mac.yaml")
            if MAC
            else os.path.join(dir.env, "medaka.yaml")
        )
    params:
        model=MEDAKA_MODEL,
        dir=os.path.join(dir.out.medaka_rd_2, "{sample}"),
        bam=os.path.join(dir.out.medaka_rd_2, "{sample}", "calls_to_draft.bam"),
        hdf=os.path.join(dir.out.medaka_rd_2, "{sample}", "consensus_probs.hdf"),
    resources:
        mem_mb=config.resources.big.mem,
        mem=str(config.resources.big.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.big.cpu
    benchmark:
        os.path.join(dir.out.bench, "medaka_round_2", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "medaka_round_2", "{sample}.log"),
    shell:
        """
        medaka_consensus -i {input.fastq} -d {input.fasta} -o {params.dir} -m {params.model}  -t {threads} 2> {log}
        touch {params.bam}
        rm {params.bam}
        touch {params.hdf}
        rm {params.hdf}
        """


rule medaka_round_2_extract_intermediate_assembly:
    """
    extracts the chromosome intermediate assembly
    Note: polished plasmids are not kept - unpolished plassembler plasmids are only kept
    Maybe something to add later, but not easy to keep plasmid information as medaka changes headers
    """
    input:
        fasta=os.path.join(dir.out.medaka_rd_2, "{sample}", "consensus.fasta"),
        completeness_check=os.path.join(dir.out.completeness, "{sample}.txt"),
        info=os.path.join(dir.out.assemblies, "{sample}", "assembly_info.txt"),
    output:
        fasta=os.path.join(
            dir.out.intermediate_assemblies, "{sample}", "{sample}_medaka_rd_2.fasta"
        ),
        ignore_list=os.path.join(dir.out.medaka_rd_2, "{sample}_ignore_list.txt"),
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


rule compare_assemblies_medaka_round_2:
    """
    compare chromosome assemblies between medaka rd 1 (post dnaapler) and rd 2
    """
    input:
        reference=os.path.join(
            dir.out.dnaapler, "{sample}", "{sample}_reoriented.fasta"
        ),
        assembly=os.path.join(
            dir.out.intermediate_assemblies, "{sample}", "{sample}_medaka_rd_2.fasta"
        ),
        diffs=os.path.join(
            dir.out.differences, "{sample}", "medaka_round_1_vs_pre_polish.txt"
        ),
    output:
        diffs=os.path.join(
            dir.out.differences, "{sample}", "medaka_round_2_vs_medaka_round_1.txt"
        ),
    conda:
        os.path.join(dir.env, "scripts.yaml")
    params:
        reference_polishing_round="medaka_round_1",
        query_polishing_round="medaka_round_2",
    resources:
        mem_mb=config.resources.med.mem,
        mem=str(config.resources.med.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.sml.cpu
    script:
        os.path.join(dir.scripts, "compare_assemblies.py")

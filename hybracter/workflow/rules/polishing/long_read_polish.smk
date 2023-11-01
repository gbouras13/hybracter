rule medaka_round_1:
    """
    Runs medaka round 1
    cleans up BAM and hdf for disc space
    """
    input:
        fasta=os.path.join(dir.out.chrom_pre_polish, "{sample}.fasta"),
        fastq=os.path.join(dir.out.qc, "{sample}_filt_trim.fastq.gz"),
    output:
        fasta=os.path.join(dir.out.medaka_rd_1, "{sample}", "consensus.fasta"),
        version=os.path.join(dir.out.versions, "{sample}", "medaka_complete.version"),
        copy_fasta=os.path.join(
            dir.out.intermediate_assemblies, "{sample}", "{sample}_medaka_rd_1.fasta"
        ),
    conda:
        os.path.join(dir.env, "medaka.yaml")
    params:
        model=MEDAKA_MODEL,
        dir=os.path.join(dir.out.medaka_rd_1, "{sample}"),
        bam=os.path.join(dir.out.medaka_rd_1, "{sample}", "calls_to_draft.bam"),
        hdf=os.path.join(dir.out.medaka_rd_1, "{sample}", "consensus_probs.hdf")
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
        cp {output.fasta} {output.copy_fasta}
        touch {params.bam}
        rm {params.bam}
        touch {params.hdf}
        rm {params.hdf}
        """


rule compare_assemblies_medaka_round_1:
    """
    compare assemblies between medaka and pre-polished chromosome
    """
    input:
        reference=os.path.join(dir.out.chrom_pre_polish, "{sample}.fasta"),
        assembly=os.path.join(dir.out.medaka_rd_1, "{sample}", "consensus.fasta"),
    output:
        diffs=os.path.join(
            dir.out.differences, "{sample}", "medaka_round_1_vs_pre_polish.txt"
        ),
    conda:
        os.path.join(dir.env, "scripts.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        mem=str(config.resources.med.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.sml.cpu
    script:
        os.path.join(dir.scripts, "compare_assemblies.py")


rule medaka_round_2:
    """
    runs medaka round 2 on the reoriented genome
    note: can't run compare_assemblies.py as the order of the genome has completely changed.
    cleans up BAM and hdf for disc space
    """
    input:
        fasta=os.path.join(dir.out.dnaapler, "{sample}", "{sample}_reoriented.fasta"),
        fastq=os.path.join(dir.out.qc, "{sample}_filt_trim.fastq.gz"),
    output:
        fasta=os.path.join(dir.out.medaka_rd_2, "{sample}", "consensus.fasta"),
        copy_fasta=os.path.join(
            dir.out.intermediate_assemblies, "{sample}", "{sample}_medaka_rd_2.fasta"
        ),
    conda:
        os.path.join(dir.env, "medaka.yaml")
    params:
        model=MEDAKA_MODEL,
        dir=os.path.join(dir.out.medaka_rd_2, "{sample}"),
        bam=os.path.join(dir.out.medaka_rd_2, "{sample}", "calls_to_draft.bam"),
        hdf=os.path.join(dir.out.medaka_rd_2, "{sample}", "consensus_probs.hdf")
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
        cp {output.fasta} {output.copy_fasta}
        touch {params.bam}
        rm {params.bam}
        touch {params.hdf}
        rm {params.hdf}
        """

rule medaka_incomplete:
    input:
        fasta=os.path.join(dir.out.incomp_pre_polish, "{sample}.fasta"),
        fastq=os.path.join(dir.out.qc, "{sample}_filt_trim.fastq.gz"),
    output:
        fasta=os.path.join(dir.out.medaka_incomplete, "{sample}", "consensus.fasta"),
        version=os.path.join(dir.out.versions, "{sample}", "medaka_incomplete.version"),
        copy_fasta=os.path.join(
            dir.out.intermediate_assemblies_incomplete,
            "{sample}",
            "{sample}_medaka.fasta",
        ),
    conda:
        (
            os.path.join(dir.env, "medaka_mac.yaml")
            if MAC
            else os.path.join(dir.env, "medaka.yaml")
        )
    params:
        model=MEDAKA_MODEL,
        dir=os.path.join(dir.out.medaka_incomplete, "{sample}"),
        bam=os.path.join(dir.out.medaka_incomplete, "{sample}", "calls_to_draft.bam"),
        hdf=os.path.join(dir.out.medaka_incomplete, "{sample}", "consensus_probs.hdf"),
    resources:
        mem_mb=config.resources.big.mem,
        mem=str(config.resources.big.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.big.cpu
    benchmark:
        os.path.join(dir.out.bench, "medaka_incomplete", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "medaka_incomplete", "{sample}.log"),
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

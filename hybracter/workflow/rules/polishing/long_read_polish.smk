rule medaka_round_1:
    """
    Runs medaka round 1
    """
    input:
        fasta=os.path.join(dir.out.chrom_pre_polish, "{sample}.fasta"),
        fastq=os.path.join(dir.out.qc, "{sample}_filt.fastq.gz"),
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
    resources:
        mem_mb=config.resources.big.mem,
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
        rm {log}
        """


rule dnaapler:
    """
    Runs dnaapler to begin chromosome with dnaa
    """
    input:
        fasta=os.path.join(dir.out.medaka_rd_1, "{sample}", "consensus.fasta"),
    output:
        fasta=os.path.join(dir.out.dnaapler, "{sample}", "{sample}_reoriented.fasta"),
        version=os.path.join(dir.out.versions, "{sample}", "dnaapler.version"),
    conda:
        os.path.join(dir.env, "dnaapler.yaml")
    params:
        dir=os.path.join(dir.out.dnaapler, "{sample}"),
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time,
    threads: config.resources.med.cpu
    benchmark:
        os.path.join(dir.out.bench, "dnaapler", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "dnaapler", "{sample}.log"),
    shell:
        """
        dnaapler chromosome -i {input.fasta} -o {params.dir} -p {wildcards.sample} -t {threads} -f 2> {log}
        dnaapler --version > {output.version}
        rm {log}
        """


rule medaka_round_2:
    input:
        fasta=os.path.join(dir.out.dnaapler, "{sample}", "{sample}_reoriented.fasta"),
        fastq=os.path.join(dir.out.qc, "{sample}_filt.fastq.gz"),
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
    resources:
        mem_mb=config.resources.big.mem,
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
        rm {log}
        """

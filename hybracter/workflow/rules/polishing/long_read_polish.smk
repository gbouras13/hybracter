rule medaka_round_1:
    """
    Runs medaka round 1
    """
    input:
        fasta = os.path.join(dir.out.chrom_pre_polish,"{sample}.fasta"),
        fastq = os.path.join(dir.out.qc,"{sample}_filt.fastq.gz")
    output:
        fasta = os.path.join(dir.out.medaka_rd_1 ,"{sample}", "consensus.fasta"),
        version = os.path.join(dir.out.versions, "{sample}", "medaka_complete.version")
    conda:
        os.path.join(dir.env,'medaka.yaml')
    params:
        model = MEDAKA_MODEL,
        dir = os.path.join(dir.out.medaka_rd_1 ,"{sample}")
    resources:
        mem_mb=config.resources.big.mem,
        time=config.resources.med.time
    threads:
        config.resources.big.cpu
    shell:
        """
        medaka_consensus -i {input.fastq} -d {input.fasta} -o {params.dir} -m {params.model}  -t {threads}
        medaka --version > {output.version}
        """

rule dnaapler:
    """
    Runs dnaapler to begin chromosome with dnaa
    """
    input:
        fasta = os.path.join(dir.out.medaka_rd_1 ,"{sample}", "consensus.fasta")
    output:
        fasta = os.path.join(dir.out.dnaapler , "{sample}", "{sample}_reoriented.fasta"),
        version = os.path.join(dir.out.versions, "{sample}", "dnaapler.version")
    conda:
        os.path.join(dir.env,'dnaapler.yaml')
    params:
        dir = os.path.join(dir.out.dnaapler, "{sample}")
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time
    threads:
        config.resources.med.cpu
    shell:
        """
        dnaapler chromosome -i {input.fasta} -o {params.dir} -p {wildcards.sample} -t {threads} -f
        dnaapler --version > {output.version}
        """

rule medaka_round_2:
    input:
        fasta = os.path.join(dir.out.dnaapler , "{sample}", "{sample}_reoriented.fasta"),
        fastq = os.path.join(dir.out.qc,"{sample}_filt.fastq.gz")
    output:
        fasta = os.path.join(dir.out.medaka_rd_2,"{sample}", "consensus.fasta")
    conda:
        os.path.join(dir.env,'medaka.yaml')
    params:
        model = MEDAKA_MODEL,
        dir = os.path.join(dir.out.medaka_rd_2,"{sample}")
    resources:
        mem_mb=config.resources.big.mem,
        time=config.resources.med.time
    threads:
        config.resources.big.cpu
    shell:
        """
        medaka_consensus -i {input.fastq} -d {input.fasta} -o {params.dir} -m {params.model}  -t {threads}
        """


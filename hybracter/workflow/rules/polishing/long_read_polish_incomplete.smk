rule medaka_incomplete:
    input:
        fasta = os.path.join(dir.out.incomp_pre_polish,"{sample}.fasta"),
        fastq = os.path.join(dir.out.qc,"{sample}_filt.fastq.gz")
    output:
        fasta = os.path.join(dir.out.medaka_incomplete,"{sample}", "consensus.fasta"),
        version = os.path.join(dir.out.versions, "{sample}", "medaka_incomplete.version")
    conda:
        os.path.join(dir.env,'medaka.yaml')
    params:
        model = MEDAKA_MODEL,
        dir = os.path.join(dir.out.medaka_incomplete,"{sample}")
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



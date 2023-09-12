
rule filtlong:
    """
    runs filtlong to filter my quality and length
    """
    input:
        fastq = get_input_lr_fastqs
    output:
        fastq = os.path.join(dir.out.qc,"{sample}_filt.fastq.gz"),
        version = os.path.join(dir.out.versions, "{sample}", "filtlong.version")
    conda:
        os.path.join(dir.env,'qc.yaml')
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time
    threads:
        config.resources.sml.cpu
    params:
        qual = config.args.min_quality,
        length = config.args.min_length
    shell:
        """
        filtlong --min_mean_q {params.qual} --min_length {params.length} {input.fastq} | pigz > {output.fastq}
        filtlong --version > {output.version}
        """

rule porechop:
    """
    runs porechop to trim adapters
    """
    input:
        fastq = os.path.join(dir.out.qc,"{sample}_filt.fastq.gz")
    output:
        fastq = os.path.join(dir.out.qc,"{sample}_filt_trim.fastq.gz"),
        version = os.path.join(dir.out.versions, "{sample}", "porechop.version")
    conda:
        os.path.join(dir.env,'qc.yaml')
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        porechop -i {input.fastq}  -o {output.fastq} -t {threads}
        porechop --version > {output.version}
        """

rule aggr_qc:
    """
    aggregates over all samples
    """
    input:
        expand(os.path.join(dir.out.qc,"{sample}_filt_trim.fastq.gz"), sample = SAMPLES),
        expand(os.path.join(dir.out.versions, "{sample}", "filtlong.version"), sample = SAMPLES),
        expand(os.path.join(dir.out.versions, "{sample}", "porechop.version"), sample = SAMPLES)
    output:
        flag = os.path.join(dir.out.flags, "aggr_qc.flag")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        touch {output.flag}
        """

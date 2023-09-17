
rule filtlong:
    """
    runs filtlong to filter quality and length
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
    benchmark:
        os.path.join(dir.out.bench, "filtlong", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "filtlong", "{sample}.log")    
    shell:
        """
        filtlong --min_mean_q {params.qual} --min_length {params.length} {input.fastq} | pigz > {output.fastq} 2> {log}
        filtlong --version > {output.version}
        rm {log}
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
        config.resources.med.cpu
    benchmark:
        os.path.join(dir.out.bench, "porechop", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "porechop", "{sample}.log")    
    shell:
        """
        porechop -i {input.fastq}  -o {output.fastq} -t {threads} 2> {log}
        porechop --version > {output.version}
        rm {log}
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

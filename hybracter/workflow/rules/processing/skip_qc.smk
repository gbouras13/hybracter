
rule skip_qc_long:
    """
    copies long reads
    """
    input:
        fastq=get_input_lr_fastqs,
    output:
        fastq=os.path.join(dir.out.qc, "{sample}_filt_trim.fastq.gz"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        cp {input.fastq} {output.fastq}
        """


rule skip_qc_short:
    """
    copies short reads
    """
    input:
        r1=get_input_r1,
        r2=get_input_r2,
    output:
        r1=os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"),
        r2=os.path.join(dir.out.fastp, "{sample}_2.fastq.gz"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        cp {input.r1} {output.r1}
        cp {input.r2} {output.r2}
        """


rule aggr_long_skip_qc:
    """
    aggregates over all samples long
    """
    input:
        expand(os.path.join(dir.out.qc, "{sample}_filt_trim.fastq.gz"), sample=SAMPLES),
    output:
        flag=os.path.join(dir.out.flags, "aggr_long_qc.flag"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        touch {output.flag}
        """


rule aggr_short_skip_qc:
    """ 
    aggregates over all samples short
    """
    input:
        expand(os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"), sample=SAMPLES),
        expand(os.path.join(dir.out.fastp, "{sample}_2.fastq.gz"), sample=SAMPLES),
    output:
        flag=os.path.join(dir.out.flags, "aggr_short_qc.flag"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        touch {output.flag}
        """

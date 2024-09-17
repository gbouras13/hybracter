
rule filtlong:
    """
    runs filtlong to filter quality and length
    """
    input:
        fastq=get_input_lr_fastqs,
        kmc=os.path.join(dir.out.kmc,"{sample}", "{sample}_kmcLOG.txt")
    output:
        fastq=temp(os.path.join(dir.out.qc, "{sample}_filt.fastq.gz")),
        version=os.path.join(dir.out.versions, "{sample}", "filtlong.version"),
    conda:
        os.path.join(dir.env, "filtlong.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        mem=str(config.resources.med.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.sml.cpu
    params:
        qual=config.args.min_quality,
        length=config.args.min_length,
        target_bases=lambda wildcards: str(getTargetBases(kmc_log_path=os.path.join(dir.out.kmc,f"{wildcards.sample}", f"{wildcards.sample}_kmcLOG.txt"),sample=wildcards.sample, auto=AUTO, target_depth=SUBSAMPLE_DEPTH)),
    benchmark:
        os.path.join(dir.out.bench, "filtlong", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "filtlong", "{sample}.log"),
    shell:
        """
        filtlong --target_bases {params.target_bases} --min_mean_q {params.qual} --min_length {params.length} {input.fastq} | pigz > {output.fastq} 2> {log}
        filtlong --version > {output.version}
        """


rule porechop_abi:
    """
    runs porechop_abi to trim adapters
    """
    input:
        fastq=os.path.join(dir.out.qc, "{sample}_filt.fastq.gz"),
    output:
        fastq=os.path.join(dir.out.qc, "{sample}_filt_trim.fastq.gz"),
        version=os.path.join(dir.out.versions, "{sample}", "porechop.version"),
    conda:
        os.path.join(dir.env, "porechop_abi.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        mem=str(config.resources.med.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.med.cpu
    benchmark:
        os.path.join(dir.out.bench, "porechop", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "porechop", "{sample}.log"),
    shell:
        """
        porechop_abi -i {input.fastq}  -o {output.fastq} -t {threads} 2> {log}
        porechop_abi --version > {output.version}
        """


"""
fastp rule for both short_read_polish.smk short_read_polish_incomplete.smk
"""


rule fastp:
    """
    runs fastp on the paired end short reads
    """
    input:
        r1=get_input_r1,
        r2=get_input_r2,
    output:
        r1=os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"),
        r2=os.path.join(dir.out.fastp, "{sample}_2.fastq.gz"),
        html=os.path.join(dir.out.fastp, "{sample}.html"),
        json=os.path.join(dir.out.fastp, "{sample}.json"),
        version=os.path.join(dir.out.versions, "{sample}", "fastp.version"),
    conda:
        os.path.join(dir.env, "fastp.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        mem=str(config.resources.med.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.sml.cpu
    benchmark:
        os.path.join(dir.out.bench, "fastp", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "fastp", "{sample}.log"),
    shell:
        """
        fastp --in1 {input.r1} --in2 {input.r2} --out1 {output.r1} --out2 {output.r2} --html {output.html} --json {output.json} --thread {threads} 2> {log}
        fastp --version 2> {output.version}
        """


rule aggr_long_qc:
    """
    aggregates over all samples for long read qc
    """
    input:
        expand(os.path.join(dir.out.qc, "{sample}_filt_trim.fastq.gz"), sample=SAMPLES),
        expand(
            os.path.join(dir.out.versions, "{sample}", "filtlong.version"),
            sample=SAMPLES,
        ),
        expand(
            os.path.join(dir.out.versions, "{sample}", "porechop.version"),
            sample=SAMPLES,
        ),
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


rule aggr_short_qc:
    """
    aggregates over all samples
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

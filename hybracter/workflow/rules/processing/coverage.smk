rule run_seqkit_short:
    input:
        r1=os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"),
        r2=os.path.join(dir.out.fastp, "{sample}_2.fastq.gz"),
    output:
        r1=os.path.join(dir.out.seqkit, "{sample}_r1.txt"),
        r2=os.path.join(dir.out.seqkit, "{sample}_r2.txt"),
    conda:
        os.path.join(dir.env, "seqkit.yaml")
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    benchmark:
        os.path.join(dir.out.bench, "seqkit", "{sample}_short.txt")
    log:
        os.path.join(dir.out.stderr, "seqkit", "{sample}_short.log"),
    shell:
        """
        seqkit stats {input.r1} -T  > {output.r1}
        seqkit stats {input.r2} -T  > {output.r2}
        """


rule run_seqkit_long:
    input:
        long=os.path.join(dir.out.qc, "{sample}_filt_trim.fastq.gz"),
    output:
        long=os.path.join(dir.out.seqkit, "{sample}_long.txt"),
    conda:
        os.path.join(dir.env, "seqkit.yaml")
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    benchmark:
        os.path.join(dir.out.bench, "seqkit", "{sample}_long.txt")
    log:
        os.path.join(dir.out.stderr, "seqkit", "{sample}_long.log"),
    shell:
        """
        seqkit stats {input.long} -T -N 50 -N 90 > {output.long}
        """


rule aggr_seqkit_short:
    """
    aggregates the seqkit stats over all samples
    """
    input:
        expand(os.path.join(dir.out.seqkit, "{sample}_r1.txt"), sample=SAMPLES),
        expand(os.path.join(dir.out.seqkit, "{sample}_r2.txt"), sample=SAMPLES),
    output:
        flag=os.path.join(dir.out.flags, "aggr_seqkit_short.flag"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        touch {output.flag}
        """


rule aggr_seqkit_long:
    """
    aggregates the seqkit stats over all samples
    """
    input:
        expand(os.path.join(dir.out.seqkit, "{sample}_long.txt"), sample=SAMPLES),
    output:
        flag=os.path.join(dir.out.flags, "aggr_seqkit_long.flag"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        touch {output.flag}
        """


rule estimate_sr_coverage:
    """
    quickly estimates short read coverage for careful
    """
    input:
        r1=os.path.join(dir.out.seqkit, "{sample}_r1.txt"),
        r2=os.path.join(dir.out.seqkit, "{sample}_r2.txt"),
    output:
        sr_coverage=os.path.join(dir.out.coverage, "{sample}.txt"),
    params:
        chromlen=getMinChromLength,
    conda:
        os.path.join(dir.env, "scripts.yaml")
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    script:
        os.path.join(dir.scripts, "sr_coverage.py")

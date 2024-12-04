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


rule estimate_lr_coverage:
    """
    quickly estimates long read coverage for min_depth
    """
    input:
        long_bases=os.path.join(dir.out.seqkit, "{sample}_long.txt"),
        lrge=os.path.join(
            dir.out.chrom_size, "{sample}_lrge_estimated_chromosome_size.txt"
        ),
    output:
        lr_coverage=os.path.join(dir.out.coverage, "{sample}_lr.txt"),
    params:
        min_chrom_length=lambda wildcards: str(
            getMinChromLength(
                lrge_path=os.path.join(
                    dir.out.chrom_size,
                    f"{wildcards.sample}_lrge_estimated_chromosome_size.txt",
                ),
                sample=wildcards.sample,
                auto=AUTO,
            )
        ),
        min_bases=lambda wildcards: str(
            getMinBases(
                lrge_path=os.path.join(
                    dir.out.chrom_size,
                    f"{wildcards.sample}_lrge_estimated_chromosome_size.txt",
                ),
                sample=wildcards.sample,
                auto=AUTO,
                min_depth=MIN_DEPTH,
            )
        ),
    conda:
        os.path.join(dir.env, "scripts.yaml")
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    script:
        os.path.join(dir.scripts, "lr_coverage.py")


rule aggr_seqkit_long:
    """
    aggregates the seqkit stats and coverage check over all samples
    """
    input:
        expand(os.path.join(dir.out.seqkit, "{sample}_long.txt"), sample=SAMPLES),
        expand(os.path.join(dir.out.coverage, "{sample}_lr.txt"), sample=SAMPLES),
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
    quickly estimates short read coverage for pypolca and polypolish careful
    """
    input:
        r1=os.path.join(dir.out.seqkit, "{sample}_r1.txt"),
        r2=os.path.join(dir.out.seqkit, "{sample}_r2.txt"),
        lrge=os.path.join(
            dir.out.chrom_size, "{sample}_lrge_estimated_chromosome_size.txt"
        ),
    output:
        sr_coverage=os.path.join(dir.out.coverage, "{sample}_sr.txt"),
    params:
        chromlen=lambda wildcards: str(
            getMinChromLength(
                lrge_path=os.path.join(
                    dir.out.chrom_size,
                    f"{wildcards.sample}_lrge_estimated_chromosome_size.txt",
                ),
                sample=wildcards.sample,
                auto=AUTO,
            )
        ),
    conda:
        os.path.join(dir.env, "scripts.yaml")
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    script:
        os.path.join(dir.scripts, "sr_coverage.py")

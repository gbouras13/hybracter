rule plassembler_hybrid:
    """
    runs plassembler for hybrid
    need the if statement due to unicycler conda installation being a bit broken.
    if unicycler --version works, then all good
    otherwise will install unicycler from pip
    """
    input:
        l=os.path.join(dir.out.qc, "{sample}_filt_trim.fastq.gz"),
        r1=os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"),
        r2=os.path.join(dir.out.fastp, "{sample}_2.fastq.gz"),
        kmc=os.path.join(dir.out.kmc,"{sample}", "{sample}_kmcLOG.txt")
    output:
        fasta=os.path.join(
            dir.out.plassembler, "{sample}", "plassembler_plasmids.fasta"
        ),
        summary=os.path.join(dir.out.plassembler, "{sample}", "plassembler_summary.tsv"),
        version=os.path.join(dir.out.versions, "{sample}", "plassembler.version"),
    params:
        db=dir.plassemblerdb,
        chromlen = lambda wildcards: str(getMinChromLength(kmc_log_path=os.path.join(dir.out.kmc, f"{wildcards.sample}", f"{wildcards.sample}_kmcLOG.txt"), sample=wildcards.sample,auto=AUTO)),
        outdir=os.path.join(dir.out.plassembler, "{sample}"),
        flye_dir=os.path.join(dir.out.assemblies, "{sample}"),
        depth_filter=DEPTH_FILTER,
    conda:
        os.path.join(dir.env, "plassembler.yaml")
    resources:
        mem_mb=config.resources.big.mem,
        mem=str(config.resources.big.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.big.cpu
    benchmark:
        os.path.join(dir.out.bench, "plassembler_hybrid", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "plassembler_hybrid", "{sample}.log"),
    shell:
        """
        if unicycler --version ; then
            plassembler run -l {input.l} -o {params.outdir} -1 {input.r1} -2 {input.r2} -d {params.db} -t {threads} -c {params.chromlen} --skip_qc --flye_directory {params.flye_dir} --depth_filter {params.depth_filter} -f 2> {log}
        else
            pip install git+https://github.com/rrwick/Unicycler.git
            plassembler run -l {input.l} -o {params.outdir} -1 {input.r1} -2 {input.r2} -d {params.db} -t {threads} -c {params.chromlen} --skip_qc --flye_directory {params.flye_dir} --depth_filter {params.depth_filter} -f 2> {log}
        fi
        touch {output.fasta}
        touch {output.summary}
        plassembler --version > {output.version}
        """


rule add_sample_plassembler:
    input:
        inp=os.path.join(dir.out.plassembler, "{sample}", "plassembler_summary.tsv"),
    output:
        out=os.path.join(
            dir.out.plassembler_individual_summaries,
            "{sample}_plassembler_summary.tsv",
        ),
    conda:
        os.path.join(dir.env, "scripts.yaml")
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    script:
        os.path.join(dir.scripts, "add_sample_plassembler.py")


# aggr rule in the aggregate.smk rule - because don't run plassembler on incomplete samples


rule plassembler_incomplete:
    """
    touch the output flag for the plassembler incomplete samples (i.e. don't run plassembler at all)
    """
    output:
        flag=os.path.join(dir.out.plassembler_incomplete, "{sample}.flag"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        touch {output.flag}
        """

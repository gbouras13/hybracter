
rule assemble:
    """
    LR assembly with Flye
    """
    input:
        fastq=os.path.join(dir.out.qc, "{sample}_filt_trim.fastq.gz"),
        lr_coverage=os.path.join(dir.out.coverage, "{sample}_lr.txt"), # make sure coverage is run before assembly in case below min_depth
    output:
        fasta=os.path.join(dir.out.assemblies, "{sample}", "assembly.fasta"),
        info=os.path.join(dir.out.assemblies, "{sample}", "assembly_info.txt"),
        version=os.path.join(dir.out.versions, "{sample}", "flye.version"),
        infocopy=os.path.join(dir.out.assembly_statistics, "{sample}_assembly_info.txt"),
    conda:
        os.path.join(dir.env, "flye.yaml")
    resources:
        mem_mb=config.resources.big.mem,
        mem=str(config.resources.big.mem) + "MB",
        time=config.resources.med.time,
    params:
        model=FLYE_MODEL,
        dir=os.path.join(dir.out.assemblies, "{sample}"),
    threads: config.resources.big.cpu
    benchmark:
        os.path.join(dir.out.bench, "assemble", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "assemble", "{sample}.log"),
    shell:
        """
        flye {params.model} {input.fastq} -t {threads}  --out-dir {params.dir} 2> {log}
        flye --version > {output.version}
        cp {output.info} {output.infocopy}
        """


rule aggr_assemble:
    input:
        expand(
            os.path.join(dir.out.assemblies, "{sample}", "assembly.fasta"),
            sample=SAMPLES,
        ),
        expand(
            os.path.join(dir.out.assemblies, "{sample}", "assembly_info.txt"),
            sample=SAMPLES,
        ),
        expand(
            os.path.join(dir.out.assembly_statistics, "{sample}_assembly_info.txt"),
            sample=SAMPLES,
        ),
        expand(
            os.path.join(dir.out.versions, "{sample}", "flye.version"), sample=SAMPLES
        ),
    output:
        flag=os.path.join(dir.out.flags, "aggr_assemble.flag"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        touch {output.flag}
        """

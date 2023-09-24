
rule assemble:
    """
    LR assembly with Flye
    """
    input:
        fastq=os.path.join(dir.out.qc, "{sample}_filt_trim.fastq.gz"),
    output:
        fasta=os.path.join(dir.out.assemblies, "{sample}", "assembly.fasta"),
        info=os.path.join(dir.out.assemblies, "{sample}", "assembly_info.txt"),
        version=os.path.join(dir.out.versions, "{sample}", "flye.version"),
        copy_fasta=os.path.join(
            dir.out.intermediate_assemblies, "{sample}", "{sample}_flye.fasta"
        ),
    conda:
        os.path.join(dir.env, "flye.yaml")
    resources:
        mem_mb=config.resources.big.mem,
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
        cp {output.fasta} {output.copy_fasta}
        """


rule extract_flye_assembly_info:
    input:
        info=os.path.join(dir.out.assemblies, "{sample}", "assembly_info.txt"),
    output:
        info=os.path.join(dir.out.assembly_statistics, "{sample}_assembly_info.txt"),
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        cp {input.info} {output.info}
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
        mem_mb=config.resources.big.mem,
        time=config.resources.big.time,
    threads: config.resources.sml.cpu
    shell:
        """
        touch {output.flag}
        """

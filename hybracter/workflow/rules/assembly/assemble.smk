
rule assemble:
    """
    LR assembly with Flye
    """
    input:
        fastq = os.path.join(dir.out.qc,"{sample}_filt_trim.fastq.gz")
    output:
        fasta = os.path.join(dir.out.assemblies ,"{sample}", "assembly.fasta"),
        info = os.path.join(dir.out.assemblies ,"{sample}", "assembly_info.txt"),
        version = os.path.join(dir.out.versions, "{sample}", "flye.version")
    conda:
        os.path.join(dir.env,'flye.yaml')
    resources:
        mem_mb=config.resources.med.cpu,
        time=config.resources.big.time
    params:
        model = FLYE_MODEL,
        dir = os.path.join(dir.out.assemblies ,"{sample}")
    threads:
        config.resources.big.cpu
    shell:
        """
        flye {params.model} {input.fastq} -t {threads}  --out-dir {params.dir}
        flye --version > {output.version}
        """

rule aggr_assemble:
    input:
        expand(os.path.join(dir.out.assemblies,"{sample}", "assembly.fasta"), sample = SAMPLES),
        expand(os.path.join(dir.out.assemblies,"{sample}", "assembly_info.txt"), sample = SAMPLES),
        expand(os.path.join(dir.out.versions, "{sample}", "flye.version"), sample = SAMPLES)
    output:
        flag = os.path.join(dir.out.flags, "aggr_assemble.flag")
    resources:
        mem_mb=config.resources.big.mem,
        time=config.resources.big.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        touch {output.flag}
        """

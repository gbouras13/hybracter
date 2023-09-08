
rule assemble:
    input:
        fastq = os.path.join(dir.out.qc,"{sample}_filt_trim.fastq.gz")
    output:
        os.path.join(dir.out.assemblies ,"{sample}", "assembly.fasta"),
        os.path.join(dir.out.assemblies ,"{sample}", "assembly_info.txt")
    conda:
        os.path.join(dir.env,'flye.yaml')
    resources:
        mem_mb=config.resources.med.cpu,
        time=config.resources.big.time
    params:
        FLYE_MODEL,
        os.path.join(dir.out.assemblies ,"{sample}")
    threads:
        config.resources.big.cpu
    shell:
        """
        flye {params[0]} {input.fastq} -t {threads}  --out-dir {params[1]}
        """

rule aggr_assemble:
    input:
        expand(os.path.join(dir.out.assemblies,"{sample}", "assembly.fasta"), sample = SAMPLES),
        expand(os.path.join(dir.out.assemblies,"{sample}", "assembly_info.txt"), sample = SAMPLES)
    output:
        os.path.join(dir.out.flags, "aggr_assemble.flag")
    resources:
        mem_mb=config.resources.big.mem,
        time=config.resources.big.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        touch {output[0]}
        """


rule assemble:
    input:
        os.path.join(QC,"{sample}_filt_trim.fastq.gz")
    output:
        directory(os.path.join(ASSEMBLIES,"{sample}")),
        os.path.join(ASSEMBLIES,"{sample}", "assembly.fasta"),
        os.path.join(ASSEMBLIES,"{sample}", "assembly_info.txt")
    conda:
        os.path.join('..', 'envs','flye.yaml')
    resources:
        mem_mb=BigJobMem,
        time=BigTime
    params:
        FLYE_MODEL
    threads:
        BigJobCpu
    shell:
        """
        flye {params[0]} {input[0]} -t {threads}  --out-dir {output[0]}
        """

rule aggr_assemble:
    input:
        expand(os.path.join(ASSEMBLIES,"{sample}", "assembly.fasta"), sample = SAMPLES),
        expand(os.path.join(ASSEMBLIES,"{sample}", "assembly_info.txt"), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "aggr_assemble.txt")
    resources:
        mem_mb=SmallJobMem,
        time=SmallTime
    threads:
        1
    shell:
        """
        touch {output[0]}
        """

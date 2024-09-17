"""
@richardstoeckl with a nice implementation
https://github.com/gbouras13/hybracter/issues/90
"""
rule sanitiseInputReads:
    input:
        fastq=get_input_lr_fastqs
    output:
        fastq=temp(os.path.join(dir.out.qc, "{sample}_sanitised.fastq.gz"))
    conda:
        os.path.join(dir.env, "seqkit.yaml")
    threads:
        1
    shell:
        """
        seqkit replace -p "\t" -r '_' {input.fastq} -o {output.fastq}
        """

rule kmc:
    """
    need params.prefix so kmc writes into the output dir
    """
    input:
        fastq = rules.sanitiseInputReads.output.fastq,
    output:
        kmcLOG = os.path.join(dir.out.kmc,"{sample}", "{sample}_kmcLOG.txt"),
        dir = directory(os.path.join(dir.out.kmc,"{sample}")),
    conda:
        os.path.join(dir.env, "kmc.yaml")
    threads:
        config.resources.med.cpu
    resources:
        mem_mb=config.resources.med.mem,
        mem=str(config.resources.med.mem) + "MB",
        time=config.resources.med.time,
    params:
        prefix=os.path.join(dir.out.kmc,"{sample}","{sample}")
    shell:
        "kmc -k21 -ci10 -t{threads} -v -fq {input.fastq} {params.prefix} {output.dir} > {output.kmcLOG} 2>&1"

rule write_chrom_size:
    """
    need params.prefix so kmc writes into the output dir
    """
    input:
        kmcLOG = os.path.join(dir.out.kmc,"{sample}", "{sample}_kmcLOG.txt")
    output:
        chrom_size = os.path.join(dir.out.chrom_size,"{sample}", "{sample}_kmc_estimated_chrom_size.txt"),
    conda:
        os.path.join(dir.env, "scripts.yaml")
    threads:
        config.resources.sml.cpu
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    script:
        os.path.join(dir.scripts, "write_chrom_size.py")

rule aggr_kmc:
    """
    aggregates over all samples
    """
    input:
        expand(os.path.join(dir.out.chrom_size,"{sample}", "{sample}_kmc_estimated_chrom_size.txt"), sample=SAMPLES),
        expand(os.path.join(dir.out.kmc,"{sample}", "{sample}_kmcLOG.txt"), sample=SAMPLES)
    output:
        flag=os.path.join(dir.out.flags, "aggr_kmc.flag"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        touch {output.flag}
        """


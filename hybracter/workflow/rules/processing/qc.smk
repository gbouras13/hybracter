
rule filtlong:
    input:
        get_input_lr_fastqs
    output:
        os.path.join(dir.out.qc,"{sample}_filt.fastq.gz")
    conda:
        os.path.join(dir.env,'qc.yaml')
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time
    threads:
        config.resources.sml.cpu
    params:
        config.min_quality,
        config.min_length
    shell:
        """
        filtlong --min_mean_q {params[0]} --min_length {params[1]} {input[0]} | gzip > {output[0]}
        """

rule porechop:
    input:
        os.path.join(dir.out.qc,"{sample}_filt.fastq.gz")
    output:
        os.path.join(dir.out.qc,"{sample}_filt_trim.fastq.gz")
    conda:
        os.path.join(dir.env,'qc.yaml')
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        porechop -i {input[0]}  -o {output[0]} -t {threads}
        """

rule aggr_qc:
    input:
        expand(os.path.join(dir.out.qc,"{sample}_filt_trim.fastq.gz"), sample = SAMPLES)
    output:
        os.path.join(dir.out.flags, "aggr_qc.flag")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        touch {output[0]}
        """

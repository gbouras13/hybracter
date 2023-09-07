
rule filtlong:
    input:
        get_input_lr_fastqs
    output:
        os.path.join(QC,"{sample}_filt.fastq.gz")
    conda:
        os.path.join('..', 'envs','qc.yaml')
    resources:
        mem_mb=BigJobMem,
        time=MediumTime
    params:
        MIN_QUALITY,
        MIN_LENGTH
    threads: 
        1
    shell:
        """
        filtlong --min_mean_q {params[0]} --min_length {params[1]} {input[0]} | gzip > {output[0]}
        """

rule porechop:
    input:
        os.path.join(QC,"{sample}_filt.fastq.gz")
    output:
        os.path.join(QC,"{sample}_filt_trim.fastq.gz")
    conda:
        os.path.join('..', 'envs','qc.yaml')
    resources:
        mem_mb=BigJobMem,
        time=BigTime
    threads:
        BigJobCpu
    shell:
        """
        porechop -i {input[0]}  -o {output[0]} -t {threads}
        """

rule aggr_qc:
    input:
        expand(os.path.join(QC,"{sample}_filt_trim.fastq.gz"), sample = SAMPLES)
    output:
        os.path.join(FLAGS, "aggr_qc.txt")
    resources:
        mem_mb=SmallJobMem,
        time=SmallTime
    threads:
        1
    shell:
        """
        touch {output[0]}
        """

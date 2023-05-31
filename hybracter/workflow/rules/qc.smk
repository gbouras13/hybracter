
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

# for really large sets, can include rasusa - no need I think
# no rasusa 
# according to Ryan Wick, more depth is better (up to 400x at least)
# https://rrwick.github.io/2021/08/10/accuracy-vs-depth.html
# rule rasusa:
#     input:
#         get_input_fastqs
#     output:
#         os.path.join(QC,"{sample}_filt_ras.fastq.gz")
#     conda:
#         os.path.join('..', 'envs','rasusa.yaml')
#     resources:
#         mem_mb=SmallJobMem,
#         time=30
#     params:
#         MIN_CHROM_LENGTH
#     shell:
#         """
#         rasusa -i {input[0]} --coverage 100 --genome-size {params[0]} -o {output[0]}
#         """

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

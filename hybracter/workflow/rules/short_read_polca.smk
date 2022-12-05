
rule polca:
    input:
        os.path.join(POLYPOLISH,"{sample}.fasta")
    output:
        os.path.join(POLCA,"{sample}.fasta.PolcaCorrected.fa")
    threads:
        BigJobCpu
    params:
        os.path.join(FASTP,"{sample}_1.fastq.gz"),
        os.path.join(FASTP,"{sample}_2.fastq.gz"),
        POLCA
    conda:
        os.path.join('..', 'envs','polca.yaml')
    resources:
        mem_mb=BigJobMem,
        time=MediumTime
    shell:
        """
        cd {params[2]}
        polca.sh -a {input[0]}  -r "{params[0]} {params[1]}" -t {resources.th}
        """



rule aggr_polca:
    input:
        expand(os.path.join(POLCA,"{sample}.fasta.PolcaCorrected.fa")), sample = SAMPLES )
    output:
        os.path.join(FLAGS, "aggr_polca.txt")
    threads:
        1
    resources:
        mem_mb=SmallJobMem,
        time=SmallTime
    shell:
        """
        touch {output[0]}
        """
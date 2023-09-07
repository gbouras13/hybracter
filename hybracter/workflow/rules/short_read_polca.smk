
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
        # real struggle running polca honestly
        cd {params[2]}
        polca.sh -a ../../../../{input[0]}  -r "../../../../{params[0]} ../../../../{params[1]}" -t {threads}
        """


rule polca_incomplete:
    input:
        os.path.join(POLYPOLISH_INCOMPLETE,"{sample}.fasta")
    output:
        os.path.join(POLCA_INCOMPLETE,"{sample}.fasta.PolcaCorrected.fa")
    threads:
        BigJobCpu
    params:
        os.path.join(FASTP,"{sample}_1.fastq.gz"),
        os.path.join(FASTP,"{sample}_2.fastq.gz"),
        POLCA_INCOMPLETE
    conda:
        os.path.join('..', 'envs','polca.yaml')
    resources:
        mem_mb=BigJobMem,
        time=MediumTime
    shell:
        """
        # real struggle running polca honestly
        cd {params[2]}
        polca.sh -a ../../../../{input[0]}  -r "../../../../{params[0]} ../../../../{params[1]}" -t {threads}
        """







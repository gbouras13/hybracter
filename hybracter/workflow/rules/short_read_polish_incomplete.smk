def get_input_r1(wildcards):
    return dictReads[wildcards.sample]["R1"]
def get_input_r2(wildcards):
    return dictReads[wildcards.sample]["R2"]
def get_input_fastqs(wildcards):
    return dictReads[wildcards.sample]["LR"]




rule bwa_index_incomplete:
    input:
        os.path.join(MEDAKA_INCOMPLETE,"{sample}", "consensus.fasta")
    output:
        os.path.join(MEDAKA_INCOMPLETE,"{sample}", "consensus.fasta.bwt")
    threads:
        SmallJobCpu
    conda:
        os.path.join('..', 'envs','short_read_polish.yaml')
    resources:
        mem_mb=BigJobMem,
        time=SmallTime
    shell:
        """
        bwa index {input}
        """

rule bwa_mem_incomplete:
    input:
        os.path.join(MEDAKA_INCOMPLETE,"{sample}", "consensus.fasta"),
        os.path.join(FASTP,"{sample}_1.fastq.gz"),
        os.path.join(FASTP,"{sample}_2.fastq.gz"),
        os.path.join(MEDAKA_INCOMPLETE,"{sample}", "consensus.fasta.bwt")
    output:
        os.path.join(BWA,"{sample}_1.sam"),
        os.path.join(BWA,"{sample}_2.sam")
    threads:
        BigJobCpu
    conda:
        os.path.join('..', 'envs','short_read_polish.yaml')
    resources:
        mem_mb=BigJobMem,
        time=BigTime
    shell:
        """
        bwa mem -t {threads} -a {input[0]} {input[1]} > {output[0]}
        bwa mem -t {threads} -a {input[0]} {input[2]} > {output[1]}
        """

rule polypolish_incomplete:
    input:
        os.path.join(MEDAKA_INCOMPLETE,"{sample}", "consensus.fasta"),
        os.path.join(BWA,"{sample}_1.sam"),
        os.path.join(BWA,"{sample}_2.sam")
    output:
        os.path.join(POLYPOLISH_INCOMPLETE,"{sample}.fasta")
    threads:
        BigJobCpu
    conda:
        os.path.join('..', 'envs','polypolish.yaml')
    resources:
        mem_mb=SmallJobMem,
        time=MediumTime
    shell:
        """
        polypolish {input[0]} {input[1]} {input[2]} > {output[0]}
        """

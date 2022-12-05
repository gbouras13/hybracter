def get_input_r1(wildcards):
    return dictReads[wildcards.sample]["R1"]
def get_input_r2(wildcards):
    return dictReads[wildcards.sample]["R2"]
def get_input_fastqs(wildcards):
    return dictReads[wildcards.sample]["LR"]


rule fastp:
    input:
        get_input_r1,
        get_input_r2
    output:
        os.path.join(FASTP,"{sample}_1.fastq.gz"),
        os.path.join(FASTP,"{sample}_2.fastq.gz")
    threads:
        SmallJobCpu
    conda:
        os.path.join('..', 'envs','short_read_polish.yaml')
    resources:
        mem_mb=BigJobMem,
        time=MediumTime
    shell:
        """
        fastp --in1 {input[0]} --in2 {input[1]} --out1 {output[0]} --out2 {output[1]} 
        """

rule bwa_index:
    input:
        os.path.join(MEDAKA_RD_2,"{sample}", "consensus.fasta")
    output:
        os.path.join(MEDAKA_RD_2,"{sample}", "consensus.fasta.bwt")
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

rule bwa_mem:
    input:
        os.path.join(MEDAKA_RD_2,"{sample}", "consensus.fasta"),
        os.path.join(FASTP,"{sample}_1.fastq.gz"),
        os.path.join(FASTP,"{sample}_2.fastq.gz"),
        os.path.join(MEDAKA_RD_2,"{sample}", "consensus.fasta.bwt")
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

rule polypolish:
    input:
        os.path.join(MEDAKA_RD_2,"{sample}", "consensus.fasta"),
        os.path.join(BWA,"{sample}_1.sam"),
        os.path.join(BWA,"{sample}_2.sam")
    output:
        os.path.join(POLYPOLISH,"{sample}.fasta")
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


rule aggr_short_read_polish:
    input:
        expand(os.path.join(POLYPOLISH,"{sample}.fasta"), sample = SAMPLES )
    output:
        os.path.join(FLAGS, "aggr_short_read_polish.txt")
    threads:
        1
    resources:
        mem_mb=SmallJobMem,
        time=SmallTime
    shell:
        """
        touch {output[0]}
        """
"""
fastp rule for both short_read_polish.smk short_read_polish_incomplete.smk
"""

rule fastp:
    input:
        get_input_r1,
        get_input_r2
    output:
        os.path.join(dir.out.fastp,"{sample}_1.fastq.gz"),
        os.path.join(dir.out.fastp,"{sample}_2.fastq.gz")
    conda:
        os.path.join(dir.env,'short_read_polish.yaml')
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time
    threads:
        config.resources.med.cpu
    shell:
        """
        fastp --in1 {input[0]} --in2 {input[1]} --out1 {output[0]} --out2 {output[1]} 
        """

rule bwa_index:
    input:
        os.path.join(dir.out.medaka_rd_2 ,"{sample}", "consensus.fasta")
    output:
        os.path.join(dir.out.medaka_rd_2 ,"{sample}", "consensus.fasta.bwt")
    conda:
        os.path.join(dir.env,'short_read_polish.yaml')
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        bwa index {input}
        """

rule bwa_mem:
    input:
        os.path.join(dir.out.medaka_rd_2 ,"{sample}", "consensus.fasta"),
        os.path.join(dir.out.fastp ,"{sample}_1.fastq.gz"),
        os.path.join(dir.out.fastp ,"{sample}_2.fastq.gz"),
        os.path.join(dir.out.medaka_rd_2 ,"{sample}", "consensus.fasta.bwt")
    output:
        os.path.join(dir.out.bwa ,"{sample}_1.sam"),
        os.path.join(dir.out.bwa ,"{sample}_2.sam")
    conda:
        os.path.join(dir.env,'short_read_polish.yaml')
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time
    threads:
        config.resources.med.cpu
    shell:
        """
        bwa mem -t {threads} -a {input[0]} {input[1]} > {output[0]}
        bwa mem -t {threads} -a {input[0]} {input[2]} > {output[1]}
        """

rule polypolish:
    input:
        os.path.join(dir.out.medaka_rd_2,"{sample}", "consensus.fasta"),
        os.path.join(dir.out.bwa ,"{sample}_1.sam"),
        os.path.join(dir.out.bwa ,"{sample}_2.sam")
    output:
        os.path.join(dir.out.polypolish,"{sample}.fasta")
    conda:
        os.path.join(dir.env,'polypolish.yaml')
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time
    threads:
        config.resources.med.cpu
    shell:
        """
        polypolish {input[0]} {input[1]} {input[2]} > {output[0]}
        """


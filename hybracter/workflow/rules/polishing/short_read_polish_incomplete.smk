rule bwa_index_incomplete:
    input:
        os.path.join(dir.out.medaka_incomplete ,"{sample}", "consensus.fasta")
    output:
        os.path.join(dir.out.medaka_incomplete ,"{sample}", "consensus.fasta.bwt")
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

rule bwa_mem_incomplete:
    input:
        os.path.join(dir.out.medaka_incomplete ,"{sample}", "consensus.fasta"),
        os.path.join(dir.out.fastp ,"{sample}_1.fastq.gz"),
        os.path.join(dir.out.fastp ,"{sample}_2.fastq.gz"),
        os.path.join(dir.out.medaka_incomplete ,"{sample}", "consensus.fasta.bwt")
    output:
        os.path.join(dir.out.bwa_incomplete ,"{sample}_1.sam"),
        os.path.join(dir.out.bwa_incomplete ,"{sample}_2.sam")
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

rule polypolish_incomplete:
    input:
        os.path.join(dir.out.medaka_incomplete ,"{sample}", "consensus.fasta"),
        os.path.join(dir.out.bwa_incomplete ,"{sample}_1.sam"),
        os.path.join(dir.out.bwa_incomplete ,"{sample}_2.sam")
    output:
        os.path.join(dir.out.polypolish_incomplete,"{sample}.fasta")
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

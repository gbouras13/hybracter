rule bwa_index_incomplete:
    input:
        fasta = os.path.join(dir.out.medaka_incomplete ,"{sample}", "consensus.fasta")
    output:
        index = os.path.join(dir.out.medaka_incomplete ,"{sample}", "consensus.fasta.bwt"),
        version = os.path.join(dir.out.versions, "{sample}", "bwa_incomplete.version")
    conda:
        os.path.join(dir.env,'short_read_polish.yaml')
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        bwa index {input.fasta}
        bwa  > {output.version}
        """

rule bwa_mem_incomplete:
    input:
        fasta = os.path.join(dir.out.medaka_incomplete ,"{sample}", "consensus.fasta"),
        r1 = os.path.join(dir.out.fastp ,"{sample}_1.fastq.gz"),
        r2 = os.path.join(dir.out.fastp ,"{sample}_2.fastq.gz"),
        index = os.path.join(dir.out.medaka_incomplete ,"{sample}", "consensus.fasta.bwt")
    output:
        sam1 = os.path.join(dir.out.bwa_incomplete ,"{sample}_1.sam"),
        sam2 = os.path.join(dir.out.bwa_incomplete ,"{sample}_2.sam")
    conda:
        os.path.join(dir.env,'short_read_polish.yaml')
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time
    threads:
        config.resources.med.cpu
    shell:
        """
        bwa mem -t {threads} -a {input.fasta} {input.r1} > {output.sam1}
        bwa mem -t {threads} -a {input.fasta} {input.r2} > {output.sam2}
        """

rule polypolish_incomplete:
    input:
        fasta = os.path.join(dir.out.medaka_incomplete ,"{sample}", "consensus.fasta"),
        sam1 = os.path.join(dir.out.bwa_incomplete ,"{sample}_1.sam"),
        sam2 = os.path.join(dir.out.bwa_incomplete ,"{sample}_2.sam")
    output:
        fasta = os.path.join(dir.out.polypolish_incomplete,"{sample}.fasta"),
        version = os.path.join(dir.out.versions, "{sample}", "polypolish_incomplete.version")
    conda:
        os.path.join(dir.env,'polypolish.yaml')
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time
    threads:
        config.resources.med.cpu
    shell:
        """
        polypolish {input.fasta} {input.sam1} {input,sam2} > {output.fasta}
        polypolish --version > {output.version}
        """

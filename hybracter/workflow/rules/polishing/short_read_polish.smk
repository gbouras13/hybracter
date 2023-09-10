"""
fastp rule for both short_read_polish.smk short_read_polish_incomplete.smk
"""

rule fastp:
    """
    runs fastp on the paired end short reads
    """
    input:
        r1 = get_input_r1,
        r2 = get_input_r2
    output:
        r1 = os.path.join(dir.out.fastp,"{sample}_1.fastq.gz"),
        r2 = os.path.join(dir.out.fastp,"{sample}_2.fastq.gz"),
        version = os.path.join(dir.out.versions, "{sample}", "fastp.version")
    conda:
        os.path.join(dir.env,'short_read_polish.yaml')
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        fastp --in1 {input.r1} --in2 {input.r2} --out1 {output.r1} --out2 {output.r2} 
        fastp --version > {output.version}
        """

rule bwa_index:
    input:
        fasta = os.path.join(dir.out.medaka_rd_2 ,"{sample}", "consensus.fasta")
    output:
        index = os.path.join(dir.out.medaka_rd_2 ,"{sample}", "consensus.fasta.bwt")
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
        """

rule bwa_mem:
    input:
        fasta = os.path.join(dir.out.medaka_rd_2 ,"{sample}", "consensus.fasta"),
        r1 = os.path.join(dir.out.fastp ,"{sample}_1.fastq.gz"),
        r2 = os.path.join(dir.out.fastp ,"{sample}_2.fastq.gz"),
        index = os.path.join(dir.out.medaka_rd_2 ,"{sample}", "consensus.fasta.bwt")
    output:
        sam1 = os.path.join(dir.out.bwa ,"{sample}_1.sam"),
        sam2 = os.path.join(dir.out.bwa ,"{sample}_2.sam")
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

rule polypolish:
    input:
        fasta = os.path.join(dir.out.medaka_rd_2,"{sample}", "consensus.fasta"),
        sam1 = os.path.join(dir.out.bwa ,"{sample}_1.sam"),
        sam2 = os.path.join(dir.out.bwa ,"{sample}_2.sam")
    output:
        fasta = os.path.join(dir.out.polypolish,"{sample}.fasta"),
        version = os.path.join(dir.out.versions, "{sample}", "polypolish_complete.version")
    conda:
        os.path.join(dir.env,'polypolish.yaml')
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time
    threads:
        config.resources.med.cpu
    shell:
        """
        polypolish {input.fasta} {input.sam1} {input.sam2} > {output.fasta}
        polypolish --version > {output.version}
        """


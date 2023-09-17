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
        os.path.join(dir.env,'fastp.yaml')
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time
    threads:
        config.resources.sml.cpu
    benchmark:
        os.path.join(dir.out.bench, "fastp", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "fastp", "{sample}.log")
    shell:
        """
        fastp --in1 {input.r1} --in2 {input.r2} --out1 {output.r1} --out2 {output.r2}  2> {log}
        fastp --version > {output.version}
        rm {log}
        """

rule bwa_index:
    input:
        fasta = os.path.join(dir.out.medaka_rd_2 ,"{sample}", "consensus.fasta")
    output:
        index = os.path.join(dir.out.medaka_rd_2 ,"{sample}", "consensus.fasta.bwt")
    conda:
        os.path.join(dir.env,'polypolish.yaml')
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        bwa index {input.fasta}
        """

rule polypolish:
    input:
        fasta = os.path.join(dir.out.medaka_rd_2,"{sample}", "consensus.fasta"),
        r1 = os.path.join(dir.out.fastp ,"{sample}_1.fastq.gz"),
        r2 = os.path.join(dir.out.fastp ,"{sample}_2.fastq.gz"),
        index = os.path.join(dir.out.medaka_rd_2 ,"{sample}", "consensus.fasta.bwt")
    output:
        fasta = os.path.join(dir.out.polypolish,"{sample}.fasta"),
        version = os.path.join(dir.out.versions, "{sample}", "polypolish_incomplete.version"),
        sam1 = temp(os.path.join(dir.out.bwa ,"{sample}_1.sam")),
        sam2 = temp(os.path.join(dir.out.bwa ,"{sample}_2.sam")).
        copy_fasta = os.path.join(dir.out.intermediate_assemblies , "{sample}_polypolish.fasta")
    conda:
        os.path.join(dir.env,'polypolish.yaml')
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time
    threads:
        config.resources.med.cpu
    benchmark:
        os.path.join(dir.out.bench, "polypolish", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "polypolish", "{sample}.log")
    shell:
        """
        bwa mem -t {threads} -a {input.fasta} {input.r1} > {output.sam1} 2> {log}
        bwa mem -t {threads} -a {input.fasta} {input.r2} > {output.sam2} 2> {log}
        polypolish {input.fasta} {output.sam1} {output.sam2} > {output.fasta} 2> {log}
        cp {output.fasta} {output.copy_fasta}
        polypolish --version > {output.version}
        rm {log}
        """


rule bwa_index_incomplete:
    input:
        fasta=os.path.join(dir.out.medaka_incomplete, "{sample}", "consensus.fasta"),
    output:
        index=os.path.join(dir.out.medaka_incomplete, "{sample}", "consensus.fasta.bwt"),
    conda:
        os.path.join(dir.env, "polypolish.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        bwa index {input.fasta} 
        """


rule bwa_mem_incomplete:
    input:
        fasta=os.path.join(dir.out.medaka_incomplete, "{sample}", "consensus.fasta"),
        r1=os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"),
        r2=os.path.join(dir.out.fastp, "{sample}_2.fastq.gz"),
        index=os.path.join(dir.out.medaka_incomplete, "{sample}", "consensus.fasta.bwt"),
    output:
        sam1=os.path.join(dir.out.bwa_incomplete, "{sample}_1.sam"),
        sam2=os.path.join(dir.out.bwa_incomplete, "{sample}_2.sam"),
    conda:
        os.path.join(dir.env, "short_read_polish.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time,
    threads: config.resources.med.cpu
    benchmark:
        os.path.join(dir.out.bench, "bwa_mem_incomplete", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "bwa_mem_incomplete", "{sample}.log"),
    shell:
        """
        bwa mem -t {threads} -a {input.fasta} {input.r1} > {output.sam1} 2> {log}
        bwa mem -t {threads} -a {input.fasta} {input.r2} > {output.sam2} 2> {log}
        """


rule polypolish_incomplete:
    input:
        fasta=os.path.join(dir.out.medaka_incomplete, "{sample}", "consensus.fasta"),
        r1=os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"),
        r2=os.path.join(dir.out.fastp, "{sample}_2.fastq.gz"),
        index=os.path.join(dir.out.medaka_incomplete, "{sample}", "consensus.fasta.bwt"),
    output:
        fasta=os.path.join(dir.out.polypolish_incomplete, "{sample}.fasta"),
        version=os.path.join(
            dir.out.versions, "{sample}", "polypolish_incomplete.version"
        ),
        sam1=os.path.join(dir.out.bwa_incomplete, "{sample}_1.sam"),
        sam2=os.path.join(dir.out.bwa_incomplete, "{sample}_2.sam"),
        copy_fasta=os.path.join(
            dir.out.intermediate_assemblies, "{sample}", "{sample}_polypolish.fasta"
        ),
    conda:
        os.path.join(dir.env, "polypolish.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time,
    threads: config.resources.med.cpu
    benchmark:
        os.path.join(dir.out.bench, "polypolish_incomplete", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "polypolish_incomplete", "{sample}.log"),
    shell:
        """
        bwa mem -t {threads} -a {input.fasta} {input.r1} > {output.sam1} 2> {log}
        bwa mem -t {threads} -a {input.fasta} {input.r2} > {output.sam2} 2> {log}
        polypolish {input.fasta} {output.sam1} {output.sam2} > {output.fasta} 2> {log}
        polypolish --version > {output.version}
        cp {output.fasta} {output.copy_fasta}
        """

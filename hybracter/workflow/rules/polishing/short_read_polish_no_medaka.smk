rule bwa_index:
    input:
        fasta=os.path.join(dir.out.dnaapler, "{sample}", "{sample}_reoriented.fasta"),
    output:
        index=os.path.join(
            dir.out.dnaapler, "{sample}", "{sample}_reoriented.fasta.bwt"
        ),
    conda:
        os.path.join(dir.env, "bwa.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        mem=str(config.resources.med.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        bwa index {input.fasta}
        """


rule bwa_mem:
    input:
        fasta=os.path.join(dir.out.dnaapler, "{sample}", "{sample}_reoriented.fasta"),
        r1=os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"),
        r2=os.path.join(dir.out.fastp, "{sample}_2.fastq.gz"),
        index=os.path.join(
            dir.out.dnaapler, "{sample}", "{sample}_reoriented.fasta.bwt"
        ),
    output:
        sam1=os.path.join(dir.out.bwa, "{sample}_1.sam"),
        sam2=os.path.join(dir.out.bwa, "{sample}_2.sam"),
    conda:
        os.path.join(dir.env, "bwa.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        mem=str(config.resources.med.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.med.cpu
    benchmark:
        os.path.join(dir.out.bench, "bwa_mem", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "bwa_mem", "{sample}.log"),
    shell:
        """
        bwa mem -t {threads} -a {input.fasta} {input.r1} > {output.sam1} 2> {log}
        bwa mem -t {threads} -a {input.fasta} {input.r2} > {output.sam2} 2> {log}
        """


rule polypolish:
    input:
        fasta=os.path.join(dir.out.dnaapler, "{sample}", "{sample}_reoriented.fasta"),
        sam1=os.path.join(dir.out.bwa, "{sample}_1.sam"),
        sam2=os.path.join(dir.out.bwa, "{sample}_2.sam"),
    output:
        fasta=os.path.join(dir.out.polypolish, "{sample}.fasta"),
        version=os.path.join(dir.out.versions, "{sample}", "polypolish.version"),
        copy_fasta=os.path.join(
            dir.out.intermediate_assemblies, "{sample}", "{sample}_polypolish.fasta"
        ),
    conda:
        os.path.join(dir.env, "polypolish.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        mem=str(config.resources.med.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.med.cpu
    benchmark:
        os.path.join(dir.out.bench, "polypolish", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "polypolish", "{sample}.log"),
    shell:
        """
        polypolish {input.fasta} {input.sam1} {input.sam2} > {output.fasta} 2> {log}
        cp {output.fasta} {output.copy_fasta}
        polypolish --version > {output.version}
        """


rule compare_assemblies_polypolish_vs_prechrom:
    """
    compare assemblies 
    """
    input:
        reference=os.path.join(
            dir.out.dnaapler, "{sample}", "{sample}_reoriented.fasta"
        ),
        assembly=os.path.join(dir.out.polypolish, "{sample}.fasta"),
    output:
        diffs=os.path.join(
            dir.out.differences, "{sample}", "polypolish_vs_pre_polish.txt"
        ),
    conda:
        os.path.join(dir.env, "scripts.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        mem=str(config.resources.med.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.sml.cpu
    script:
        os.path.join(dir.scripts, "compare_assemblies.py")

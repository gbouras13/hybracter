rule bwa_index:
    input:
        fasta=os.path.join(dir.out.medaka_rd_2, "{sample}", "consensus.fasta"),
    output:
        index=os.path.join(dir.out.medaka_rd_2, "{sample}", "consensus.fasta.bwt"),
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
        fasta=os.path.join(dir.out.medaka_rd_2, "{sample}", "consensus.fasta"),
        r1=os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"),
        r2=os.path.join(dir.out.fastp, "{sample}_2.fastq.gz"),
        index=os.path.join(dir.out.medaka_rd_2, "{sample}", "consensus.fasta.bwt"),
    output:
        sam1=temp(os.path.join(dir.out.bwa, "{sample}_1.sam")),
        sam2=temp(os.path.join(dir.out.bwa, "{sample}_2.sam")),
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
        fasta=os.path.join(dir.out.medaka_rd_2, "{sample}", "consensus.fasta"),
        sam1=os.path.join(dir.out.bwa, "{sample}_1.sam"),
        sam2=os.path.join(dir.out.bwa, "{sample}_2.sam"),
    output:
        fasta=os.path.join(dir.out.polypolish, "{sample}.fasta"),
        version=os.path.join(dir.out.versions, "{sample}", "polypolish.version"),
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
        polypolish --version > {output.version}
        """


rule polypolish_extract_intermediate_assembly:
    """
    extracts the chromosome intermediate assembly
    """
    input:
        fasta=os.path.join(dir.out.polypolish, "{sample}.fasta"),
        completeness_check=os.path.join(dir.out.completeness, "{sample}.txt"),
        info=os.path.join(dir.out.assemblies, "{sample}", "assembly_info.txt")
    output:
        fasta=os.path.join(
            dir.out.intermediate_assemblies, "{sample}", "{sample}_polypolish.fasta"
        ),
    params:
        min_chrom_length=getMinChromLength,
        polypolish_flag=True
    conda:
        os.path.join(dir.env, "scripts.yaml")
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    script:
        os.path.join(dir.scripts, "extract_chromosome.py")


rule compare_assemblies_polypolish_vs_medaka_round_2:
    """
    compare chrom assemblies 
    """
    input:
        reference=os.path.join(
            dir.out.intermediate_assemblies, "{sample}", "{sample}_medaka_rd_2.fasta"
        ),
        assembly=os.path.join(
            dir.out.intermediate_assemblies, "{sample}", "{sample}_polypolish.fasta"
        ),
        diffs=os.path.join(
            dir.out.differences, "{sample}", "medaka_round_2_vs_medaka_round_1.txt"
        ),
    output:
        diffs=os.path.join(
            dir.out.differences, "{sample}", "polypolish_vs_medaka_round_2.txt"
        ),
    conda:
        os.path.join(dir.env, "scripts.yaml")
    params:
        reference_polishing_round="medaka_round_2",
        query_polishing_round="polypolish",
    resources:
        mem_mb=config.resources.med.mem,
        mem=str(config.resources.med.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.sml.cpu
    script:
        os.path.join(dir.scripts, "compare_assemblies.py")


rule pypolca:
    input:
        polypolish_fasta=os.path.join(dir.out.polypolish, "{sample}.fasta"),
        r1=os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"),
        r2=os.path.join(dir.out.fastp, "{sample}_2.fastq.gz"),
    output:
        fasta=os.path.join(dir.out.pypolca, "{sample}", "{sample}_corrected.fasta"),
        version=os.path.join(dir.out.versions, "{sample}", "pypolca_complete.version"),
    params:
        pypolca_dir=os.path.join(dir.out.pypolca, "{sample}"),
        version=os.path.join(dir.out.versions, "{sample}", "pypolca_complete.version"),
        copy_fasta=os.path.join(
            dir.out.intermediate_assemblies, "{sample}", "{sample}_pypolca.fasta"
        ),
    conda:
        os.path.join(dir.env, "pypolca.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        mem=str(config.resources.med.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.med.cpu
    benchmark:
        os.path.join(dir.out.bench, "pypolca", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "pypolca", "{sample}.log"),
    shell:
        """
        pypolca run -a {input.polypolish_fasta} -1 {input.r1} -2 {input.r2} -o {params.pypolca_dir} -t {threads} -f -p {wildcards.sample} 2> {log}
        pypolca --version > {params.version}
        cp {output.fasta} {params.copy_fasta} 
        """


rule compare_assemblies_pypolca_vs_polypolish:
    """
    compare assemblies 
    """
    input:
        reference=os.path.join(dir.out.polypolish, "{sample}.fasta"),
        assembly=os.path.join(dir.out.pypolca, "{sample}", "{sample}_corrected.fasta"),
        diffs=os.path.join(
            dir.out.differences, "{sample}", "polypolish_vs_pre_polish.txt"
        ),
    output:
        diffs=os.path.join(dir.out.differences, "{sample}", "pypolca_vs_polypolish.txt"),
    conda:
        os.path.join(dir.env, "scripts.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        mem=str(config.resources.med.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.sml.cpu
    script:
        os.path.join(dir.scripts, "compare_assemblies.py")


rule pypolca_incomplete:
    input:
        polypolish_fasta=os.path.join(dir.out.polypolish_incomplete, "{sample}.fasta"),
        r1=os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"),
        r2=os.path.join(dir.out.fastp, "{sample}_2.fastq.gz"),
    output:
        fasta=os.path.join(
            dir.out.pypolca_incomplete, "{sample}", "{sample}_corrected.fasta"
        ),
        version=os.path.join(dir.out.versions, "{sample}", "pypolca_incomplete.version"),
    params:
        pypolca_dir=os.path.join(dir.out.pypolca_incomplete, "{sample}"),
        version=os.path.join(dir.out.versions, "{sample}", "pypolca_incomplete.version"),
        copy_fasta=os.path.join(
            dir.out.intermediate_assemblies, "{sample}", "{sample}_pypolca.fasta"
        ),
    conda:
        os.path.join(dir.env, "pypolca.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        mem=str(config.resources.med.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.med.cpu
    benchmark:
        os.path.join(dir.out.bench, "pypolca_incomplete", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "pypolca_incomplete", "{sample}.log"),
    shell:
        """
        pypolca run -a {input.polypolish_fasta} -1 {input.r1} -2 {input.r2} -o {params.pypolca_dir} -t {threads} -f -p {wildcards.sample} 2> {log}
        pypolca --version > {params.version}
        cp {output.fasta} {params.copy_fasta} 
        """

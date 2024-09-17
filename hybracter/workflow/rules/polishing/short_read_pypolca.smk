
rule pypolca:
    """
    if coverage > 5, run pypolca --careful
    otherwise don't run pypolca at all (depth < 5) - copies file to make it work in snakemake for later
    """
    input:
        polypolish_fasta=os.path.join(dir.out.polypolish, "{sample}.fasta"),
        r1=os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"),
        r2=os.path.join(dir.out.fastp, "{sample}_2.fastq.gz"),
        coverage=os.path.join(dir.out.coverage, "{sample}_sr.txt"),
    output:
        fasta=os.path.join(dir.out.pypolca, "{sample}", "{sample}_corrected.fasta"),
        version=os.path.join(dir.out.versions, "{sample}", "pypolca_complete.version"),
    params:
        pypolca_dir=os.path.join(dir.out.pypolca, "{sample}"),
        version=os.path.join(dir.out.versions, "{sample}", "pypolca_complete.version"),
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
        coverage=$(head -n 1 {input.coverage})
        if [ "$coverage" -gt 5 ]; then
            pypolca run -a {input.polypolish_fasta} -1 {input.r1} -2 {input.r2} -o {params.pypolca_dir} -t {threads} -f -p {wildcards.sample} --careful 2> {log}
        else
            cp {input.polypolish_fasta} {output.fasta}
        fi

        pypolca --version > {params.version}
        """


rule pypolca_extract_intermediate_assembly:
    """
    extracts the chromosome intermediate assembly from pypolca
    """
    input:
        fasta=os.path.join(dir.out.pypolca, "{sample}", "{sample}_corrected.fasta"),
        completeness_check=os.path.join(dir.out.completeness, "{sample}.txt"),
        info=os.path.join(dir.out.assemblies, "{sample}", "assembly_info.txt"),
        kmc=os.path.join(dir.out.kmc,"{sample}", "{sample}_kmcLOG.txt")
    output:
        fasta=os.path.join(
            dir.out.intermediate_assemblies, "{sample}", "{sample}_pypolca.fasta"
        ),
        ignore_list=os.path.join(dir.out.pypolca, "{sample}_ignore_list.txt"),
    params:
        min_chrom_length=lambda wildcards: str(getMinChromLength(kmc_log_path=os.path.join(dir.out.kmc, f"{wildcards.sample}", f"{wildcards.sample}_kmcLOG.txt"), sample=wildcards.sample,auto=AUTO)),
        polypolish_flag=True,
    conda:
        os.path.join(dir.env, "scripts.yaml")
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    script:
        os.path.join(dir.scripts, "extract_chromosome.py")


rule compare_assemblies_pypolca_vs_polypolish:
    """
    compare chrom assemblies 
    """
    input:
        reference=os.path.join(
            dir.out.intermediate_assemblies, "{sample}", "{sample}_polypolish.fasta"
        ),
        assembly=os.path.join(
            dir.out.intermediate_assemblies, "{sample}", "{sample}_pypolca.fasta"
        ),
        diffs=os.path.join(
            dir.out.differences, "{sample}", "polypolish_vs_medaka_round_2.txt"
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
    params:
        reference_polishing_round="polypolish",
        query_polishing_round="pypolca",
    script:
        os.path.join(dir.scripts, "compare_assemblies.py")


rule pypolca_incomplete:
    """
    if coverage > 5, run pypolca --careful
    otherwise don't run pypolca at all (depth < 5) - copies file to make it work in snakemake for later
    """
    input:
        polypolish_fasta=os.path.join(dir.out.polypolish_incomplete, "{sample}.fasta"),
        r1=os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"),
        r2=os.path.join(dir.out.fastp, "{sample}_2.fastq.gz"),
        coverage=os.path.join(dir.out.coverage, "{sample}_sr.txt"),
    output:
        fasta=os.path.join(
            dir.out.pypolca_incomplete, "{sample}", "{sample}_corrected.fasta"
        ),
        version=os.path.join(dir.out.versions, "{sample}", "pypolca_incomplete.version"),
    params:
        pypolca_dir=os.path.join(dir.out.pypolca_incomplete, "{sample}"),
        version=os.path.join(dir.out.versions, "{sample}", "pypolca_incomplete.version"),
        copy_fasta=os.path.join(
            dir.out.intermediate_assemblies_incomplete,
            "{sample}",
            "{sample}_pypolca.fasta",
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
        coverage=$(head -n 1 {input.coverage})
        if [ "$coverage" -gt 5 ]; then
            pypolca run -a {input.polypolish_fasta} -1 {input.r1} -2 {input.r2} -o {params.pypolca_dir} -t {threads} -f -p {wildcards.sample} --careful 2> {log}
        else
            cp {input.polypolish_fasta} {output.fasta}
        fi

        pypolca --version > {params.version}
        cp {output.fasta} {params.copy_fasta} 
        """

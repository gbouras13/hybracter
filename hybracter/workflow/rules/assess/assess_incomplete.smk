"""
scores need to be sequential so that it can easily be aggregated based on a polca flag or not - otherwise I need another rule
"""


rule assess_incomp_pre_polish:
    """
    Run ALE on incomp_pre_polish
    """
    input:
        r1=os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"),
        fasta=os.path.join(dir.out.incomp_pre_polish, "{sample}.fasta"),
    output:
        ale=temp(
            os.path.join(dir.out.ale_out_files, "{sample}", "incomp_pre_polish.ale")
        ),
        sam1=temp(os.path.join(dir.out.ale_sams, "{sample}_incomp_pre_polish_1.sam")),
        score=os.path.join(
            dir.out.ale_scores_incomplete, "{sample}", "incomp_pre_polish.score"
        ),
    conda:
        os.path.join(dir.env, "ale.yaml")
    resources:
        mem_mb=config.resources.med.cpu,
        mem=str(config.resources.med.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.big.cpu
    benchmark:
        os.path.join(dir.out.bench, "ale", "{sample}_incomp_pre_polish.txt")
    log:
        os.path.join(dir.out.stderr, "ale", "{sample}_incomp_pre_polish.log"),
    shell:
        """

        bwa index {input.fasta}
        bwa mem -t {threads} -a {input.fasta} {input.r1} > {output.sam1} 2> {log}
        ALE {output.sam1} {input.fasta} {output.ale} 2> {log}
        grep "# ALE_score: " {output.ale} | sed 's/# ALE_score: //' > {output.score}
        """


rule assess_medaka_incomplete:
    """
    Run ALE on medaka_incomplete
    """
    input:
        r1=os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"),
        fasta=os.path.join(dir.out.medaka_incomplete, "{sample}", "consensus.fasta"),
        index=os.path.join(dir.out.medaka_incomplete, "{sample}", "consensus.fasta.bwt"),
        score=os.path.join(
            dir.out.ale_scores_incomplete, "{sample}", "incomp_pre_polish.score"
        ),
    output:
        ale=temp(
            os.path.join(dir.out.ale_out_files, "{sample}", "medaka_incomplete.ale")
        ),
        sam1=temp(os.path.join(dir.out.ale_sams, "{sample}_medaka_incomplete_1.sam")),
        score=os.path.join(
            dir.out.ale_scores_incomplete, "{sample}", "medaka_incomplete.score"
        ),
    conda:
        os.path.join(dir.env, "ale.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        mem=str(config.resources.med.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.big.cpu
    benchmark:
        os.path.join(dir.out.bench, "ale", "{sample}_incomp_medaka.txt")
    log:
        os.path.join(dir.out.stderr, "ale", "{sample}_incomp_medaka.log"),
    shell:
        """
        bwa mem -t {threads} -a {input.fasta} {input.r1} > {output.sam1} 2> {log}
        ALE {output.sam1} {input.fasta} {output.ale} 2> {log}
        grep "# ALE_score: " {output.ale} | sed 's/# ALE_score: //' > {output.score}
        """


rule assess_polypolish_incomplete:
    """
    Run ALE on polyplosh incomplete
    """
    input:
        r1=os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"),
        fasta=os.path.join(dir.out.polypolish_incomplete, "{sample}.fasta"),
        score=os.path.join(
            dir.out.ale_scores_incomplete, "{sample}", "medaka_incomplete.score"
        ),
    output:
        ale=temp(
            os.path.join(
                dir.out.ale_out_files, "{sample}", "polypolish_incomplete.ale"
            )
        ),
        sam1=temp(
            os.path.join(dir.out.ale_sams, "{sample}_polypolish_incomplete_1.sam")
        ),
        score=os.path.join(
            dir.out.ale_scores_incomplete, "{sample}", "polypolish_incomplete.score"
        ),
    conda:
        os.path.join(dir.env, "ale.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        mem=str(config.resources.med.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.big.cpu
    benchmark:
        os.path.join(dir.out.bench, "ale", "{sample}_incomp_polypolish.txt")
    log:
        os.path.join(dir.out.stderr, "ale", "{sample}_incomp_polypolish.log"),
    shell:
        """

        bwa index {input.fasta}
        bwa mem -t {threads} -a {input.fasta} {input.r1} > {output.sam1} 2> {log}
        ALE {output.sam1} {input.fasta} {output.ale} 2> {log}
        grep "# ALE_score: " {output.ale} | sed 's/# ALE_score: //' > {output.score}
        """


rule assess_polca_incomplete:
    """
    Run ALE on polca incomplete
    """
    input:
        r1=os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"),
        fasta=os.path.join(
            dir.out.pypolca_incomplete, "{sample}", "{sample}_corrected.fasta"
        ),
        score=os.path.join(
            dir.out.ale_scores_incomplete, "{sample}", "polypolish_incomplete.score"
        ),
    output:
        ale=temp(
            os.path.join(dir.out.ale_out_files, "{sample}", "pypolca_incomplete.ale")
        ),
        sam1=temp(os.path.join(dir.out.ale_sams, "{sample}_pypolca_incomplete_1.sam")),
        score=os.path.join(
            dir.out.ale_scores_incomplete, "{sample}", "pypolca_incomplete.score"
        ),
    conda:
        os.path.join(dir.env, "ale.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        mem=str(config.resources.med.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.big.cpu
    benchmark:
        os.path.join(dir.out.bench, "ale", "{sample}_incomp_pypolca.txt")
    log:
        os.path.join(dir.out.stderr, "ale", "{sample}_incomp_pypolca.log"),
    shell:
        """

        bwa index {input.fasta}
        bwa mem -t {threads} -a {input.fasta} {input.r1} > {output.sam1} 2> {log}
        ALE {output.sam1} {input.fasta} {output.ale} 2> {log}
        grep "# ALE_score: " {output.ale} | sed 's/# ALE_score: //' > {output.score}
        """

"""
scores need to be sequential so that it can easily be aggregated based on a polca flag or not - otherwise I need another rule
"""


rule assess_chrom_pre_polish:
    """
    Run ALE on chrom_pre_polish
    """
    input:
        r1=os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"),
        fasta=os.path.join(dir.out.chrom_pre_polish, "{sample}.fasta"),
    output:
        score=os.path.join(
            dir.out.ale_scores_complete, "{sample}", "chrom_pre_polish.score"
        ),
        ale=temp(
            os.path.join(dir.out.ale_out_files, "{sample}", "chrom_pre_polish.ale")
        ),
        sam1=temp(os.path.join(dir.out.ale_sams, "{sample}_chrom_pre_polish_1.sam")),
    conda:
        os.path.join(dir.env, "ale.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.big.time,
    threads: config.resources.big.cpu
    benchmark:
        os.path.join(dir.out.bench, "ale", "{sample}_chrom_pre_polish.txt")
    log:
        os.path.join(dir.out.stderr, "ale", "{sample}_chrom_pre_polish.log"),
    shell:
        """
        bwa index {input.fasta}
        bwa mem -t {threads} -a {input.fasta} {input.r1} > {output.sam1} 2> {log}
        ALE {output.sam1} {input.fasta} {output.ale} 2> {log}
        grep "# ALE_score: " {output.ale} | sed 's/# ALE_score: //' > {output.score}
        rm {log}
        """


rule assess_medaka_rd_1:
    """
    Run ALE on medaka_rd_1
    """
    input:
        r1=os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"),
        fasta=os.path.join(dir.out.medaka_rd_1, "{sample}", "consensus.fasta"),
        score=os.path.join(
            dir.out.ale_scores_complete, "{sample}", "chrom_pre_polish.score"
        ),
    output:
        score=os.path.join(dir.out.ale_scores_complete, "{sample}", "medaka_rd_1.score"),
        ale=temp(os.path.join(dir.out.ale_out_files, "{sample}", "medaka_rd_1.ale")),
        sam1=temp(os.path.join(dir.out.ale_sams, "{sample}medaka_rd_1.sam")),
    conda:
        os.path.join(dir.env, "ale.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.big.time,
    threads: config.resources.big.cpu
    benchmark:
        os.path.join(dir.out.bench, "ale", "{sample}_medaka_rd_1.txt")
    log:
        os.path.join(dir.out.stderr, "ale", "{sample}_medaka_rd_1.log"),
    shell:
        """
        bwa index {input.fasta}
        bwa mem -t {threads} -a {input.fasta} {input.r1} > {output.sam1} 2> {log}
        ALE {output.sam1} {input.fasta} {output.ale} 2> {log}
        grep "# ALE_score: " {output.ale} | sed 's/# ALE_score: //' > {output.score}
        rm {log}
        """


rule assess_medaka_rd_2:
    """
    Run ALE on medaka_rd_2
    """
    input:
        r1=os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"),
        fasta=os.path.join(dir.out.medaka_rd_2, "{sample}", "consensus.fasta"),
        score=os.path.join(dir.out.ale_scores_complete, "{sample}", "medaka_rd_1.score"),
    output:
        score=os.path.join(dir.out.ale_scores_complete, "{sample}", "medaka_rd_2.score"),
        ale=temp(os.path.join(dir.out.ale_out_files, "{sample}", "medaka_rd_2.ale")),
        sam1=temp(os.path.join(dir.out.ale_sams, "{sample}medaka_rd_2.sam")),
    conda:
        os.path.join(dir.env, "ale.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.big.time,
    threads: config.resources.big.cpu
    benchmark:
        os.path.join(dir.out.bench, "ale", "{sample}_medaka_rd_2.txt")
    log:
        os.path.join(dir.out.stderr, "ale", "{sample}_medaka_rd_2.log"),
    shell:
        """
        bwa index {input.fasta}
        bwa mem -t {threads} -a {input.fasta} {input.r1} > {output.sam1} 2> {log}
        ALE {output.sam1} {input.fasta} {output.ale} 2> {log}
        grep "# ALE_score: " {output.ale} | sed 's/# ALE_score: //' > {output.score}
        rm {log}
        """


rule assess_polypolish:
    """
    Run ALE on polypolish
    """
    input:
        r1=os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"),
        fasta=os.path.join(dir.out.polypolish, "{sample}.fasta"),
        score=os.path.join(dir.out.ale_scores_complete, "{sample}", "medaka_rd_2.score"),
    output:
        score=os.path.join(dir.out.ale_scores_complete, "{sample}", "polypolish.score"),
        ale=temp(os.path.join(dir.out.ale_out_files, "{sample}", "polypolish.ale")),
        sam1=temp(os.path.join(dir.out.ale_sams, "{sample}polypolish.sam")),
    conda:
        os.path.join(dir.env, "ale.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.big.time,
    threads: config.resources.big.cpu
    benchmark:
        os.path.join(dir.out.bench, "ale", "{sample}_polypolish.txt")
    log:
        os.path.join(dir.out.stderr, "ale", "{sample}_polypolish.log"),
    shell:
        """
        bwa index {input.fasta}
        bwa mem -t {threads} -a {input.fasta} {input.r1} > {output.sam1} 2> {log}
        ALE {output.sam1} {input.fasta} {output.ale} 2> {log}
        grep "# ALE_score: " {output.ale} | sed 's/# ALE_score: //' > {output.score}
        rm {log}
        """


rule assess_polca:
    """
    Run ALE on polypolish
    """
    input:
        r1=os.path.join(dir.out.fastp, "{sample}_1.fastq.gz"),
        fasta=os.path.join(
            dir.out.polca, "{sample}", "{sample}.fasta.PolcaCorrected.fa"
        ),
        score=os.path.join(dir.out.ale_scores_complete, "{sample}", "polypolish.score"),
    output:
        score=os.path.join(dir.out.ale_scores_complete, "{sample}", "polca.score"),
        ale=temp(os.path.join(dir.out.ale_out_files, "{sample}", "polca.ale")),
        sam1=temp(os.path.join(dir.out.ale_sams, "{sample}polca.sam")),
    conda:
        os.path.join(dir.env, "ale.yaml")
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.big.time,
    threads: config.resources.big.cpu
    benchmark:
        os.path.join(dir.out.bench, "ale", "{sample}_polca.txt")
    log:
        os.path.join(dir.out.stderr, "ale", "{sample}_polca.log"),
    shell:
        """
        bwa index {input.fasta}
        bwa mem -t {threads} -a {input.fasta} {input.r1} > {output.sam1} 2> {log}
        ALE {output.sam1} {input.fasta} {output.ale} 2> {log}
        grep "# ALE_score: " {output.ale} | sed 's/# ALE_score: //' > {output.score}
        rm {log}
        """

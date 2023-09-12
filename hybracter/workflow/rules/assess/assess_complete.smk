rule assess_chrom_pre_polish:
    """
    Run ALE on chrom_pre_polish
    """
    input:
        ALE = os.path.join(dir.env, "ALE-master", "src", "ALE"),
        r1 = os.path.join(dir.out.fastp ,"{sample}_1.fastq.gz"),
        fasta = os.path.join(dir.out.chrom_pre_polish,"{sample}.fasta")
    output:
        score = os.path.join(dir.out.ale_scores_complete ,"{sample}", "chrom_pre_polish.score"),
        ale = temp(os.path.join(dir.out.ale_out_files ,"{sample}", "chrom_pre_polish.ale")),
        sam1 = temp(os.path.join(dir.out.ale_sams ,"{sample}_chrom_pre_polish_1.sam"))
    conda:
        os.path.join(dir.env,'polca.yaml')
    resources:
        mem_mb=config.resources.med.cpu,
        time=config.resources.big.time
    threads:
        config.resources.big.cpu
    benchmark:
        os.path.join(dir.out.bench, "ale", "{sample}_chrom_pre_polish.txt")
    log:
        os.path.join(dir.out.stderr, "ale", "{sample}_chrom_pre_polish.log")
    shell:
        """
        bwa index {input.fasta}
        bwa mem -t {threads} -a {input.fasta} {input.r1} > {output.sam1} 2> {log}
        {input.ALE} {output.sam1} {input.fasta} {output.ale} 2> {log}
        grep "# ALE_score: " {output.ale} | sed 's/# ALE_score: //' > {output.score}
        """

rule assess_medaka_rd_1:
    """
    Run ALE on medaka_rd_1
    """
    input:
        ALE = os.path.join(dir.env, "ALE-master", "src", "ALE"),
        r1 = os.path.join(dir.out.fastp ,"{sample}_1.fastq.gz"),
        fasta =  os.path.join(dir.out.medaka_rd_1 ,"{sample}", "consensus.fasta"),
        score = os.path.join(dir.out.ale_scores_complete ,"{sample}", "chrom_pre_polish.score")
    output:
        score = os.path.join(dir.out.ale_scores_complete ,"{sample}", "medaka_rd_1.score"),
        ale = temp(os.path.join(dir.out.ale_out_files ,"{sample}", "medaka_rd_1.ale")),
        sam1 = temp(os.path.join(dir.out.ale_sams ,"{sample}medaka_rd_1.sam"))
    conda:
        os.path.join(dir.env,'polca.yaml')
    resources:
        mem_mb=config.resources.med.cpu,
        time=config.resources.big.time
    threads:
        config.resources.big.cpu
    benchmark:
        os.path.join(dir.out.bench, "ale", "{sample}_medaka_rd_1.txt")
    log:
        os.path.join(dir.out.stderr, "ale", "{sample}_medaka_rd_1.log")
    shell:
        """
        bwa index {input.fasta}
        bwa mem -t {threads} -a {input.fasta} {input.r1} > {output.sam1} 2> {log}
        {input.ALE} {output.sam1} {input.fasta} {output.ale} 2> {log}
        grep "# ALE_score: " {output.ale} | sed 's/# ALE_score: //' > {output.score}
        """

rule assess_incomp_pre_polish:
    """
    Run ALE on incomp_pre_polish
    """
    input:
        ALE = os.path.join(dir.env, "ALE-master", "src", "ALE"),
        r1 = os.path.join(dir.out.fastp ,"{sample}_1.fastq.gz"),
        fasta = os.path.join(dir.out.incomp_pre_polish,"{sample}.fasta")
    output:
        ale = temp(os.path.join(dir.out.ale_out_files ,"{sample}", "incomp_pre_polish.ale")),
        sam1 = temp(os.path.join(dir.out.ale_sams ,"{sample}_incomp_pre_polish_1.sam"))
        score = os.path.join(dir.out.ale_scores_incomplete ,"{sample}", "incomp_pre_polish.score")
    conda:
        os.path.join(dir.env,'polca.yaml')
    resources:
        mem_mb=config.resources.med.cpu,
        time=config.resources.big.time
    threads:
        config.resources.big.cpu
    benchmark:
        os.path.join(dir.out.bench, "ale", "{sample}_incomp_pre_polish.txt")
    log:
        os.path.join(dir.out.stderr, "ale", "{sample}_incomp_pre_polish.log")
    shell:
        """

        bwa index {input.fasta}
        bwa mem -t {threads} -a {input.fasta} {input.r1} > {output.sam1} 2> {log}
        {input.ALE} {output.sam1} {input.fasta} {output.ale} 2> {log}
        grep "# ALE_score: " {output.ale} | sed 's/# ALE_score: //' > {output.score}
        """

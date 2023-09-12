rule assess_chrom_pre_polish:
    """
    Run ALE on chrom_pre_polish
    """
    input:
        ALE = os.path.join(dir.env, "ALE-master", "src", "ALE"),
        r1 = os.path.join(dir.out.fastp ,"{sample}_1.fastq.gz"),
        fasta = os.path.join(dir.out.chrom_pre_polish,"{sample}.fasta")
    params:
        ale = os.path.join(dir.out.ale_out_files ,"{sample}", "chrom_pre_polish.ale"),
        sam1 = os.path.join(dir.out.ale_sams ,"{sample}_chrom_pre_polish_1.sam")
    output:
        score = os.path.join(dir.out.ale_scores_complete ,"{sample}", "chrom_pre_polish.score")
    conda:
        os.path.join(dir.env,'polypolish.yaml')
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
        bwa mem -t {threads} -a {input.fasta} {input.r1} > {params.sam1} 2> {log}
        {input.ALE} {params.sam1} {input.fasta} {params.ale} 2> {log}
        grep "# ALE_score: " {params.ale} | sed 's/# ALE_score: //' > {output.score}
        rm {params.ale}
        """

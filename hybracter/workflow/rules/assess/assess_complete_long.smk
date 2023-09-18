"""
scores need to be sequential so that it can easily be aggregated based on a polca flag or not - otherwise I need another rule
"""

rule assess_long_complete_rule:
    """
    Run pyrodigal to assess all chromosomes
    """
    input:
        pre_polish_fasta = os.path.join(dir.out.chrom_pre_polish,"{sample}.fasta"),
        medaka_rd_1_fasta = os.path.join(dir.out.medaka_rd_1 ,"{sample}", "consensus.fasta"),
        medaka_rd_2_fasta = os.path.join(dir.out.medaka_rd_2 ,"{sample}", "consensus.fasta"),
    output:
        score = os.path.join(dir.out.pyrodigal_mean_lengths_complete ,"{sample}", "chrom_pre_polish.score"),
    conda:
        os.path.join(dir.env,'pyrodigal.yaml')
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    benchmark:
        os.path.join(dir.out.bench, "pyrodigal", "{sample}_assess_all_chroms.txt")
    log:
        os.path.join(dir.out.stderr, "pyrodigal", "{sample}_assess_all_chroms.log")
    script:
        os.path.join(dir.scripts,  'assess_long_complete.py')


"""
aggregate the pyrdigal mean length cds and finalise
"""


def aggregate_finalise(wildcards):
    # decision based on content of output file
    # Important: use the method open() of the returned file!
    # This way, Snakemake is able to automatically download the file if it is generated in
    # a cloud environment without a shared filesystem.
    with checkpoints.check_completeness.get(sample=wildcards.sample).output[
        0
    ].open() as f:
        if f.read().strip() == "C":  # complete
            return os.path.join(dir.out.ale_scores_complete, "{sample}", "polca.score")
        else:  # incomplete
            return os.path.join(
                dir.out.ale_scores_incomplete, "{sample}", "polca_incomplete.score"
            )


### from the aggregate_finalise function - so it dynamic
rule aggregate_finalise_complete:
    input:
        chrom_pre_polish_fasta=os.path.join(dir.out.chrom_pre_polish, "{sample}.fasta"),
        medaka_rd_1_fasta=os.path.join(
            dir.out.medaka_rd_1, "{sample}", "consensus.fasta"
        ),
        medaka_rd_2_fasta=os.path.join(
            dir.out.medaka_rd_2, "{sample}", "consensus.fasta"
        ),
        final_plasmid_fasta=os.path.join(
            dir.out.final_contigs_complete, "{sample}_plasmid.fasta"
        ),
    output:
        chromosome_fasta=os.path.join(
            dir.out.final_contigs_complete, "{sample}_chromosome.fasta"
        ),
        total_fasta=os.path.join(dir.out.final_contigs_complete, "{sample}_final.fasta"),
        pyrodigal_summary=os.path.join(
            dir.out.pyrodigal_summary, "complete", "{sample}_summary.tsv"
        ),
        hybracter_summary=os.path.join(dir.out.final_summaries_complete, "{sample}.tsv"),
    params:
        complete_flag=True,
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.med.time,
    conda:
        os.path.join(dir.env, "pyrodigal.yaml")
    threads: config.resources.sml.cpu
    script:
        os.path.join(dir.scripts, "select_best_chromosome_assembly_long_complete.py")


### from the aggregate_finalise function - so it dynamic
rule aggregate_finalise_incomplete:
    input:
        pre_polish_fasta=os.path.join(dir.out.incomp_pre_polish, "{sample}.fasta"),
        medaka_fasta=os.path.join(
            dir.out.medaka_incomplete, "{sample}", "consensus.fasta"
        ),
    output:
        total_fasta=os.path.join(
            dir.out.final_contigs_incomplete, "{sample}_final.fasta"
        ),
        pyrodigal_summary=os.path.join(
            dir.out.pyrodigal_summary, "incomplete", "{sample}_summary.tsv"
        ),
        hybracter_summary=os.path.join(
            dir.out.final_summaries_incomplete, "{sample}.tsv"
        ),
    params:
        pre_polish_fasta=os.path.join(dir.out.incomp_pre_polish, "{sample}.fasta"),
        medaka_fasta=os.path.join(
            dir.out.medaka_incomplete, "{sample}", "consensus.fasta"
        ),
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.med.time,
    conda:
        os.path.join(dir.env, "pyrodigal.yaml")
    threads: config.resources.sml.cpu
    script:
        os.path.join(dir.scripts, "select_best_chromosome_assembly_long_incomplete.py")

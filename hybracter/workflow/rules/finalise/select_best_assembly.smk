"""
aggregate the ALE scores and pick the best one
"""

# to import aggregate_ale_input


# input function for the rule aggregate polca
def aggregate_ale_input_finalise(wildcards):
    # decision based on content of output file
    # Important: use the method open() of the returned file!
    # This way, Snakemake is able to automatically download the file if it is generated in
    # a cloud environment without a shared filesystem.
    with checkpoints.check_completeness.get(sample=wildcards.sample).output[
        0
    ].open() as f:
        if f.read().strip() == "C":  # complete
            if config.args.no_polca is False:  # with polca
                return os.path.join(
                    dir.out.ale_scores_complete, "{sample}", "polca.score"
                )
            else:  # with polca, best is polypolish
                return os.path.join(
                    dir.out.ale_scores_complete, "{sample}", "polypolish.score"
                )
        else:  # incomplete
            if config.args.no_polca is False:  # with polca
                return os.path.join(
                    dir.out.ale_scores_incomplete, "{sample}", "polca_incomplete.score"
                )
            else:
                return os.path.join(
                    dir.out.ale_scores_incomplete,
                    "{sample}",
                    "polypolish_incomplete.score",
                )


### from the aggregate_ale_input function - so it dynamic
# also calculates the summary
rule select_best_chromosome_assembly_complete:
    input:
        ale_input=aggregate_ale_input_finalise,
        aggr_ale_flag=os.path.join(dir.out.aggr_ale, "{sample}.txt"),  # to make sure ale has finished
        plassembler_fasta=os.path.join(
            dir.out.plassembler, "{sample}", "plassembler_plasmids.fasta"
        ),
    output:
        chromosome_fasta=os.path.join(
            dir.out.final_contigs_complete, "{sample}_chromosome.fasta"
        ),
        plasmid_fasta=os.path.join(
            dir.out.final_contigs_complete, "{sample}_plasmid.fasta"
        ),
        total_fasta=os.path.join(dir.out.final_contigs_complete, "{sample}_final.fasta"),
        ale_summary=os.path.join(dir.out.ale_summary, "complete", "{sample}.tsv"),
        hybracter_summary=os.path.join(dir.out.final_summaries_complete, "{sample}.tsv"),
    params:
        ale_dir=os.path.join(dir.out.ale_scores_complete, "{sample}"),
        chrom_pre_polish_fasta=os.path.join(dir.out.chrom_pre_polish, "{sample}.fasta"),
        medaka_rd_1_fasta=os.path.join(
            dir.out.medaka_rd_1, "{sample}", "consensus.fasta"
        ),
        medaka_rd_2_fasta=os.path.join(
            dir.out.medaka_rd_2, "{sample}", "consensus.fasta"
        ),
        polypolish_fasta=os.path.join(dir.out.polypolish, "{sample}.fasta"),
        polca_fasta=os.path.join(dir.out.polca, "{sample}", "{sample}.fasta"),
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.med.time,
    conda:
        os.path.join(dir.env, "scripts.yaml")
    threads: config.resources.sml.cpu
    script:
        os.path.join(dir.scripts, "select_best_chromosome_assembly_complete.py")


### from the aggregate_ale_input function - so it dynamic
#  also calculates the summary
rule select_best_chromosome_assembly_incomplete:
    input:
        aggregate_ale_input,
        os.path.join(dir.out.aggr_ale, "{sample}.txt"),  # to make sure ale has finished
    output:
        fasta=os.path.join(dir.out.final_contigs_incomplete, "{sample}_final.fasta"),
        ale_summary=os.path.join(dir.out.ale_summary, "incomplete", "{sample}.tsv"),
        hybracter_summary=os.path.join(
            dir.out.final_summaries_incomplete, "{sample}.tsv"
        ),
    params:
        ale_dir=os.path.join(dir.out.ale_scores_incomplete, "{sample}"),
        pre_polish_fasta=os.path.join(dir.out.incomp_pre_polish, "{sample}.fasta"),
        medaka_fasta=os.path.join(
            dir.out.medaka_incomplete, "{sample}", "consensus.fasta"
        ),
        polypolish_fasta=os.path.join(dir.out.polypolish_incomplete, "{sample}.fasta"),
        polca_fasta=os.path.join(dir.out.polca_incomplete, "{sample}", "{sample}.fasta"),
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.med.time,
    conda:
        os.path.join(dir.env, "scripts.yaml")
    threads: config.resources.sml.cpu
    script:
        os.path.join(dir.scripts, "select_best_chromosome_assembly_incomplete.py")

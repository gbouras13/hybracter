"""

Stores all the aggregation rules and checkpoints post assembly - due to the split between complete and incomplete

"""

"""
plassembler
"""


# input function for the rule aggregate
def aggregate_plassembler_input(wildcards):
    # decision based on content of output file
    # Important: use the method open() of the returned file!
    # This way, Snakemake is able to automatically download the file if it is generated in
    # a cloud environment without a shared filesystem.
    with checkpoints.check_completeness.get(sample=wildcards.sample).output[
        0
    ].open() as f:
        if f.read().strip() == "C":
            return os.path.join(
                dir.out.plassembler, "{sample}", "plassembler_plasmids.fasta"
            )
        else:  # if incomplete
            return os.path.join(dir.out.plassembler_incomplete, "{sample}.flag")


### from the long_read_polishing


rule aggregate_plassembler_input_rule:
    input:
        aggregate_plassembler_input,
    output:
        os.path.join(dir.out.aggr_plassembler, "{sample}.txt"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        echo {input[0]}
        touch {output}
        """


rule aggr_plassembler_flag:
    """Aggregate."""
    input:
        expand(os.path.join(dir.out.aggr_plassembler, "{sample}.txt"), sample=SAMPLES),
    output:
        flag=os.path.join(dir.out.flags, "aggr_plassembler.flag"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        touch {output.flag}
        """


"""
long read polishing
"""


# input function for the rule aggregate
def aggregate_long_read_polish_input(wildcards):
    # decision based on content of output file
    # Important: use the method open() of the returned file!
    # This way, Snakemake is able to automatically download the file if it is generated in
    # a cloud environment without a shared filesystem.
    with checkpoints.check_completeness.get(sample=wildcards.sample).output[0].open() as f:
        if config.args.no_medaka is False:  # default
            if f.read().strip() == "C":
                return os.path.join(
                    dir.out.differences, "{sample}", "medaka_round_2_vs_medaka_round_1.txt"
                )
            else:  # if incomplete
                return os.path.join(
                    dir.out.medaka_incomplete, "{sample}", "consensus.fasta"
                )
        else:  # no medaka
            if f.read().strip() == "C":
                return os.path.join(
                    dir.out.dnaapler, "{sample}", "{sample}_reoriented.fasta"
                )
            else:  # if incomplete
                return os.path.join(dir.out.incomp_pre_polish, "{sample}.fasta")


### from the long_read_polishing


rule aggregate_long_read_polish_input_rule:
    input:
        aggregate_long_read_polish_input,
    output:
        os.path.join(dir.out.aggr_lr_polish, "{sample}.txt"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        echo {input[0]}
        touch {output}
        """


rule aggr_long_read_polish_flag:
    """Aggregate."""
    input:
        expand(os.path.join(dir.out.aggr_lr_polish, "{sample}.txt"), sample=SAMPLES),
    output:
        flag=os.path.join(dir.out.flags, "aggr_long_read_polish.flag"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        touch {output.flag}
        """


####### pyrodigal ##########

"""
ale

hybrid

"""


# input function for aggresgating pyrodigal plassembler analysis
def aggregate_pyrodigal_plassembler_input(wildcards):
    # decision based on content of output file
    # Important: use the method open() of the returned file!
    # This way, Snakemake is able to automatically download the file if it is generated in
    # a cloud environment without a shared filesystem.
    with checkpoints.check_completeness.get(sample=wildcards.sample).output[
        0
    ].open() as f:
        if f.read().strip() == "C":  # complete
            return os.path.join(
                dir.out.plassembler_individual_summaries, "{sample}_with_sample.tsv"
            )
        else:
            return os.path.join(dir.out.plassembler_incomplete, "{sample}.flag")


### from the aggregate_pyrodigal_plassembler_input
rule aggregate_pyrodigal_plassembler_input_rule:
    input:
        aggregate_pyrodigal_plassembler_input,
    output:
        os.path.join(dir.out.aggr_pyrodigal_plassembler, "{sample}.txt"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        touch {output}
        """


rule aggr_pyrodigal_plassembler_flag:
    """Aggregate."""
    input:
        expand(
            os.path.join(dir.out.aggr_pyrodigal_plassembler, "{sample}.txt"),
            sample=SAMPLES,
        ),
    output:
        flag=os.path.join(dir.out.flags, "aggr_pyrodigal_plassembler.flag"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        touch {output.flag}
        """


"""
finalise

long

"""


# input function for the rule aggregate polca
def aggregate_finalised_assemblies(wildcards):
    # decision based on content of output file
    # Important: use the method open() of the returned file!
    # This way, Snakemake is able to automatically download the file if it is generated in
    # a cloud environment without a shared filesystem.
    with checkpoints.check_completeness.get(sample=wildcards.sample).output[0].open() as f:
        if f.read().strip() == "C":  # complete
            return os.path.join(dir.out.final_contigs_complete, "{sample}_final.fasta")
        else:  # incomplete
            return os.path.join(dir.out.final_contigs_incomplete, "{sample}_final.fasta")


### from the aggregate_ale_input
rule aggregate_finalised_assemblies_rule:
    input:
        aggregate_finalised_assemblies,
    output:
        os.path.join(dir.out.aggr_final, "{sample}.txt"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        touch {output}
        """


"""
to create the very final summaries
"""


rule create_final_summary:
    input:
        expand(os.path.join(dir.out.aggr_final, "{sample}.txt"), sample=SAMPLES),
    output:
        hybracter_summary=os.path.join(dir.out.final_summaries, "hybracter_summary.tsv"),
    params:
        complete_summaries_dir=dir.out.final_summaries_complete,
        incomplete_summaries_dir=dir.out.final_summaries_incomplete,
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    conda:
        os.path.join(dir.env, "scripts.yaml")
    script:
        os.path.join(dir.scripts, "create_final_hybracter_summary.py")


rule aggr_final_flag:
    """Aggregate."""
    input:
        hybracter_summary=os.path.join(dir.out.final_summaries, "hybracter_summary.tsv"),
    output:
        flag=os.path.join(dir.out.flags, "aggr_final.flag"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        touch {output.flag}
        """

"""
aggregate the ALE scores and pick the best one
"""

# to import aggregate_ale_input
include: os.path.join("..", "completeness", "aggregate.smk")


### from the aggregate_ale_input function - so it dynamic whether or not polca is selected
rule select_best_chromosome_assembly_complete:
    input:
        ale_input = aggregate_ale_input,
        aggr_ale_flag = os.path.join(dir.out.aggr_ale,"{sample}.txt"), # to make sure ale has finished
        plassembler_fasta = os.path.join(dir.out.plassembler ,"{sample}", "plassembler_plasmids.fasta")
    output:
        chromosome_fasta = os.path.join(dir.out.final_contigs_complete,"{sample}_chromosome.fasta"),
        plasmid_fasta = os.path.join(dir.out.final_contigs_complete,"{sample}_plasmid.fasta"),
        total_fasta = os.path.join(dir.out.final_contigs_complete,"{sample}_final.fasta"),
        ale_summary = os.path.join(dir.out.ale_summary,"complete","{sample}.csv")
    params:
        ale_dir = dir.out.ale_scores_complete,
        chrom_pre_polish_fasta = os.path.join(dir.out.chrom_pre_polish,"{sample}.fasta"),
        medaka_rd_1_fasta = os.path.join(dir.out.medaka_rd_1 ,"{sample}", "consensus.fasta"),
        medaka_rd_2_fasta = os.path.join(dir.out.medaka_rd_2 ,"{sample}", "consensus.fasta"),
        polypolish_fasta = os.path.join(dir.out.polypolish,"{sample}.fasta"),
        polca_fasta = os.path.join(dir.out.polca, "{sample}", "{sample}.fasta")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.med.time
    threads:
        config.resources.sml.cpu
    script:
        os.path.join(dir.scripts,  'select_best_chromosome_assembly_complete.py')



### from the aggregate_ale_input function - so it dynamic whether or not polca is selected
rule select_best_chromosome_assembly_incomplete:
    input:
        aggregate_ale_input,
        os.path.join(dir.out.aggr_ale,"{sample}.txt") # to make sure ale has finished
    output:
        fasta = os.path.join(dir.out.final_contigs_incomplete,"{sample}_final.fasta"),
        ale_summary = os.path.join(dir.out.ale_summary, "incomplete", "{sample}.csv")
    params:
        ale_dir = dir.out.ale_scores_incomplete,
        pre_polish_fasta = os.path.join(dir.out.incomp_pre_polish,"{sample}.fasta"),
        medaka_fasta = os.path.join(dir.out.medaka_incomplete ,"{sample}", "consensus.fasta"),
        polypolish_fasta = os.path.join(dir.out.polypolish_incomplete,"{sample}.fasta"),
        polca_fasta = os.path.join(dir.out.polca_incomplete, "{sample}", "{sample}.fasta")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.med.time
    threads:
        config.resources.sml.cpu
    script:
        os.path.join(dir.scripts,  'select_best_chromosome_assembly_incomplete.py')



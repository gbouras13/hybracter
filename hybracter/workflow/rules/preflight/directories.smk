"""
Ensures consistent variable names and file locations for the pipeline.

A lot taken and modified from hecatomb
"""

import attrmap as ap

dir = ap.AttrMap()

### DATABASE LOCATION


try:
    assert (ap.utils.to_dict(config.args)["databases"]) is not None
    dir.plassemblerdb = config.args.databases
except (KeyError, AssertionError):
    dir.plassemblerdb = os.path.join(workflow.basedir, "..", "databases")

### OUTPUT LOCATION
try:
    assert (ap.utils.to_dict(config.args)["output"]) is not None
    dir.out.base = config.args.output
except (KeyError, AssertionError):
    dir.out.base = "hybracter_out"


### WORKFLOW DIRs
dir.env = os.path.join(workflow.basedir, "envs")
dir.rules = os.path.join(workflow.basedir, "rules")
dir.scripts = os.path.join(workflow.basedir, "scripts")

dir.scripts_no_medaka = os.path.join(workflow.basedir, "scripts", "no_medaka")

### OUTPUT DIRs
dir.out.results = os.path.join(dir.out.base, "supplementary_results")
dir.out.flags = os.path.join(dir.out.base, "flags")
dir.out.versions = os.path.join(dir.out.base, "versions")
dir.out.processing = os.path.join(dir.out.base, "processing")
dir.out.complete = os.path.join(dir.out.processing, "complete")
dir.out.incomplete = os.path.join(dir.out.processing, "incomplete")
dir.out.qc = os.path.join(dir.out.processing, "qc")
dir.out.kmc = os.path.join(dir.out.processing, "kmc")
# auto only - write chrom size to file
dir.out.chrom_size = os.path.join(dir.out.processing, "chrom_size")

# logs and benchmarks
dir.out.bench = os.path.join(dir.out.base, "benchmarks")
dir.out.stderr = os.path.join(dir.out.base, "stderr")

# assemble.smk and assembly_statistics.smk
dir.out.assemblies = os.path.join(dir.out.processing, "assemblies")
dir.out.assembly_statistics = os.path.join(dir.out.results, "flye_individual_summaries")
dir.out.assembly_summary = os.path.join(dir.out.results, "flye_all_assembly_summary")

# extract_fastas.smk
dir.out.chrom_pre_polish = os.path.join(dir.out.processing, "pre_polish")
dir.out.incomp_pre_polish = os.path.join(dir.out.processing, "incomp_pre_polish")
dir.out.completeness = os.path.join(dir.out.base, "completeness")

# aggregation dirs
dir.out.aggr_lr_polish = os.path.join(dir.out.flags, "aggr_long_read_polish")
dir.out.aggr_sr_polish = os.path.join(dir.out.flags, "aggr_short_read_polish")
# dir.out.aggr_polca_polish = os.path.join(dir.out.flags, "aggr_polca_polish")
dir.out.aggr_pypolca_polish = os.path.join(dir.out.flags, "aggr_pypolca_polish")
dir.out.aggr_plassembler = os.path.join(dir.out.flags, "aggr_plassembler")
dir.out.aggr_ale = os.path.join(dir.out.flags, "aggr_ale")
dir.out.aggr_pyrodigal = os.path.join(dir.out.flags, "aggr_pyrodigal")
dir.out.aggr_pyrodigal_plassembler = os.path.join(
    dir.out.flags, "aggr_pyrodigal_plassembler"
)
dir.out.aggr_final = os.path.join(dir.out.flags, "aggr_final")

# long_read_polish.smk
dir.out.medaka_rd_1 = os.path.join(dir.out.complete, "medaka_rd_1")
dir.out.medaka_rd_2 = os.path.join(dir.out.complete, "medaka_rd_2")

# dnaapler
dir.out.dnaapler = os.path.join(dir.out.complete, "dnaapler")

# long_Read_polish_incomplete.smk
dir.out.medaka_incomplete = os.path.join(dir.out.incomplete, "medaka_incomplete")

# for both complete and incomplete
dir.out.fastp = os.path.join(dir.out.qc, "fastp")

# for coverage
dir.out.seqkit = os.path.join(dir.out.qc, "seqkit")
dir.out.coverage = os.path.join(dir.out.qc, "coverage")

# short_read_polish.smk
dir.out.bwa = os.path.join(dir.out.complete, "bwa")
dir.out.polypolish = os.path.join(dir.out.complete, "polypolish")
# short_read_polish_incomplete.smk
dir.out.bwa_incomplete = os.path.join(dir.out.incomplete, "bwa_incomplete")
dir.out.polypolish_incomplete = os.path.join(
    dir.out.incomplete, "polypolish_incomplete"
)

# short_read_pypolca.smk (has complete and incomplete)
# dir.out.polca = os.path.join(dir.out.complete, "polca")
# dir.out.polca_incomplete = os.path.join(dir.out.incomplete, "polca")
dir.out.pypolca = os.path.join(dir.out.complete, "pypolca")
dir.out.pypolca_incomplete = os.path.join(dir.out.incomplete, "pypolca")


# plassembler dirs
dir.out.plassembler = os.path.join(dir.out.processing, "plassembler")
dir.out.plassembler_incomplete = os.path.join(
    dir.out.processing, "plassembler_incomplete"
)
dir.out.plassembler_fastas = os.path.join(dir.out.processing, "plassembler_fastas")
dir.out.plassembler_individual_summaries = os.path.join(
    dir.out.results, "plassembler_individual_summaries"
)
dir.out.plassembler_all_summary = os.path.join(
    dir.out.results, "plassembler_all_assembly_summary"
)

# ale
dir.out.ale_scores_complete = os.path.join(dir.out.processing, "ale_scores_complete")
dir.out.ale_scores_incomplete = os.path.join(
    dir.out.processing, "ale_scores_incomplete"
)
dir.out.ale_sams = os.path.join(dir.out.processing, "ale_sams")
dir.out.ale_out_files = os.path.join(dir.out.processing, "ale_out_files")
dir.out.ale_summary = os.path.join(dir.out.results, "ale_score_summaries")

# pyrodigal
dir.out.pyrodigal_mean_lengths_complete = os.path.join(
    dir.out.processing, "pyrodigal_mean_lengths_complete"
)
dir.out.pyrodigal_mean_lengths_incomplete = os.path.join(
    dir.out.processing, "pyrodigal_mean_lengths_incomplete"
)
dir.out.pyrodigal_summary = os.path.join(
    dir.out.results, "pyrodigal_mean_length_summaries"
)
dir.out.pyrodigal_summary_plassembler = os.path.join(
    dir.out.results, "pyrodigal_mean_length_summaries_plassembler"
)

# intermediate assemblies
dir.out.intermediate_assemblies = os.path.join(
    dir.out.results, "intermediate_chromosome_assemblies"
)

dir.out.intermediate_assemblies_incomplete = os.path.join(
    dir.out.results, "intermediate_incomplete_assemblies"
)

# final contigs
dir.out.final_contigs = os.path.join(dir.out.base, "FINAL_OUTPUT")
dir.out.final_contigs_complete = os.path.join(dir.out.final_contigs, "complete")
dir.out.final_contigs_incomplete = os.path.join(dir.out.final_contigs, "incomplete")

# final summaries
dir.out.final_summaries = os.path.join(dir.out.base, "FINAL_OUTPUT")
dir.out.final_summaries_complete = os.path.join(dir.out.final_summaries, "complete")
dir.out.final_summaries_incomplete = os.path.join(dir.out.final_summaries, "incomplete")

# Test dirs
dir.test = os.path.join(workflow.basedir, "../", "test_data")
dir.test_fastqs = os.path.join(dir.test, "Fastqs")
dir.env = os.path.join(workflow.basedir, "envs")

# contaminants
dir.contaminant_genomes = os.path.join(workflow.basedir, "../", "contaminant_genomes")
dir.out.contaminant_index = os.path.join(dir.out.processing, "contaminant_index")
dir.out.contaminant_removal = os.path.join(dir.out.processing, "contaminant_removal")

# compare assemblies
dir.out.differences = os.path.join(dir.out.results, "comparisons")

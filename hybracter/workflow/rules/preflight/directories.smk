"""
Ensures consistent variable names and file locations for the pipeline.

A lot taken and modified from hecatomb
"""

import attrmap as ap

dir = ap.AttrMap()

### DATABASE LOCATION
try:
    assert(config.args.databases) is not None
    dir.plassemblerdb = config.args.databases
except (KeyError,AssertionError):
    dir.plassemblerdb = os.path.join(workflow.basedir,"..","databases")


### OUTPUT LOCATION
try:
    assert(ap.utils.to_dict(config.args)["output"]) is not None
    dir.out.base = config.args.output
except (KeyError, AssertionError):
    dir.out.base = "hybracter.out"


### WORKFLOW DIRs
dir.env     = os.path.join(workflow.basedir, "envs")
dir.rules = os.path.join(workflow.basedir, "rules")
dir.scripts = os.path.join(workflow.basedir, "scripts")

### OUTPUT DIRs
dir.out.results  = os.path.join(dir.out.base, "results")
dir.out.flags  = os.path.join(dir.out.base, "flags")
dir.out.processing  = os.path.join(dir.out.base, "processing")
dir.out.complete  = os.path.join(dir.out.results, "complete")
dir.out.incomplete  = os.path.join(dir.out.results, "incomplete")
dir.out.qc  = os.path.join(dir.out.processing, "qc")
# assemble.smk and assembly_statistics.smk
dir.out.assemblies  = os.path.join(dir.out.processing, "assemblies")
dir.out.assembly_statistics  = os.path.join(dir.out.results, "assembly_statistics")
dir.out.assembly_summary  = os.path.join(dir.out.results, "assembly_summary")

#extract_fastas.smk
dir.out.chrom_pre_polish  = os.path.join(dir.out.results, "chrom_pre_polish")
dir.out.incomp_pre_polish  = os.path.join(dir.out.results, "incomp_pre_polish")
dir.out.completeness  = os.path.join(dir.out.base, "completeness")

# aggregation dirs
dir.out.aggr_lr_polish  = os.path.join(dir.out.flags, "aggr_long_read_polish")
dir.out.aggr_sr_polish  = os.path.join(dir.out.flags, "aggr_short_read_polish")
dir.out.aggr_polca_polish  = os.path.join(dir.out.flags, "aggr_polca_polish")


# long_read_polish.smk
dir.out.medaka_rd_1  = os.path.join(dir.out.complete, "medaka_rd_1")
dir.out.medaka_rd_2  = os.path.join(dir.out.complete, "medaka_rd_2")
dir.out.dnaapler  = os.path.join(dir.out.complete, "dnaapler")

# long_Read_polish_incomplete.smk
dir.out.medaka_incomplete  = os.path.join(dir.out.incomplete, "medaka_incomplete")

# for both complete and incomplete
dir.out.fastp  = os.path.join(dir.out.qc, "fastp")
# short_read_polish.smk 
dir.out.bwa  = os.path.join(dir.out.complete, "bwa")
dir.out.polypolish  = os.path.join(dir.out.complete, "polypolish")
# short_read_polish_incomplete.smk 
dir.out.bwa_incomplete  = os.path.join(dir.out.complete, "bwa_incomplete")
dir.out.polypolish_incomplete  = os.path.join(dir.out.complete, "polypolish_incomplete")


# short_read_polca.smk (has complete and incomplete)
dir.out.polca  = os.path.join(dir.out.complete, "polca")
dir.out.polca_incomplete  = os.path.join(dir.out.incomplete, "polca")

# PLASSEMBLER DIR
dir.out.plassembler =  os.path.join(dir.out.processing, 'plassembler')
dir.out.plassembler_fastas =  os.path.join(dir.out.results, 'plassembler_fastas')
dir.out.plassembler_individual_summaries =  os.path.join(dir.out.results, 'plassembler_individual_summaries')
dir.out.plassembler_all_summary =  os.path.join(dir.out.results, 'plassembler_all_summary')









rule plassembler_long:
    """
    runs plassembler for long
    """
    input:
        l=os.path.join(dir.out.qc, "{sample}_filt_trim.fastq.gz"),
    output:
        fasta=os.path.join(
            dir.out.plassembler, "{sample}", "plassembler_plasmids.fasta"
        ),
        summary=os.path.join(dir.out.plassembler, "{sample}", "plassembler_summary.tsv"),
        version=os.path.join(dir.out.versions, "{sample}", "plassembler.version"),
    params:
        db=dir.plassemblerdb,
        chromlen=getMinChromLength,
        outdir=os.path.join(dir.out.plassembler, "{sample}"),
        flye_dir=os.path.join(dir.out.assemblies, "{sample}"),
    conda:
        os.path.join(dir.env, "plassembler.yaml")
    resources:
        mem_mb=config.resources.big.mem,
        time=config.resources.med.time,
    threads: config.resources.big.cpu
    benchmark:
        os.path.join(dir.out.bench, "plassembler_long", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "plassembler_long", "{sample}.log"),
    shell:
        """
        plassembler long -l {input.l} -o {params.outdir} -d {params.db} -t {threads} -c {params.chromlen} --skip_qc --flye_directory {params.flye_dir} -f 2> {log}
        touch {output.fasta}
        touch {output.summary}
        plassembler --version > {output.version}
        """


rule add_sample_plassembler:
    input:
        inp=os.path.join(dir.out.plassembler, "{sample}", "plassembler_summary.tsv"),
    output:
        out=os.path.join(
            dir.out.plassembler_individual_summaries, "{sample}_with_sample.tsv"
        ),
    conda:
        os.path.join(dir.env, "scripts.yaml")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    script:
        os.path.join(dir.scripts, "add_sample_plassembler.py")


rule plassembler_polish_medaka:
    """
    runs 1 round of medaka on plassembler output to polish them up
    """
    input:
        fasta=os.path.join(
            dir.out.plassembler, "{sample}", "plassembler_plasmids.fasta"
        ),
        fastq=os.path.join(dir.out.qc, "{sample}_filt_trim.fastq.gz"),
        individual_summary=os.path.join(
            dir.out.plassembler_individual_summaries, "{sample}_with_sample.tsv"
        ),
        # to make it sequential for the aggregate step
    output:
        fasta=os.path.join(dir.out.plassembler, "{sample}", "medaka", "consensus.fasta"),
    conda:
        os.path.join(dir.env, "medaka.yaml")
    params:
        model=MEDAKA_MODEL,
        medaka_dir=os.path.join(dir.out.plassembler, "{sample}", "medaka"),
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time,
    threads: config.resources.med.cpu
    benchmark:
        os.path.join(dir.out.bench, "medaka_plassembler", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "medaka_plassembler", "{sample}.log"),
    shell:
        """
        plass_out={input.fasta}
        if [ -s "$plass_out" ]
        then
            medaka_consensus -i {input.fastq} -d {input.fasta} -o {params.medaka_dir} -m {params.model}  -t {threads} 2> {log}
        else
            touch {output.fasta} 2> {log}
        fi
        rm {log}
        """


rule plassembler_assess_polish:
    """
    runs pyrodigal to check
    """
    input:
        plassembler_fasta=os.path.join(
            dir.out.plassembler, "{sample}", "plassembler_plasmids.fasta"
        ),
        medaka_fasta=os.path.join(
            dir.out.plassembler, "{sample}", "medaka", "consensus.fasta"
        ),
    output:
        final_plasmid_fasta=os.path.join(
            dir.out.final_contigs_complete, "{sample}_plasmid.fasta"
        ),
        plassembler_prodigal_summary=os.path.join(
            dir.out.pyrodigal_summary_plassembler,
            "complete",
            "{sample}.tsv",
        ),
    conda:
        os.path.join(dir.env, "pyrodigal.yaml")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    benchmark:
        os.path.join(dir.out.bench, "pyrodigal_plassembler", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "pyrodigal_plassembler", "{sample}.log"),
    script:
        os.path.join(dir.scripts, "assess_plassembler_long_complete.py")


# aggr rule in the aggregate_long.smk - because don't run plassembler on incomplete samples


rule plassembler_incomplete:
    """
    touch the output flag for the plassembler incomplete samples (i.e. don't run plassembler at all)
    """
    output:
        flag=os.path.join(dir.out.plassembler_incomplete, "{sample}.flag"),
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        touch {output.flag}
        """

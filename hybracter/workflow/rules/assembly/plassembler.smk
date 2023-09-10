rule plassembler_hybrid:
    input:
        l = os.path.join(dir.out.qc,"{sample}_filt_trim.fastq.gz"),
        r1 = os.path.join(dir.out.fastp ,"{sample}_1.fastq.gz"),
        r2 = os.path.join(dir.out.fastp ,"{sample}_2.fastq.gz")
    output:
        fasta = os.path.join(dir.out.plassembler ,"{sample}", "plassembler_plasmids.fasta"),
        summary = os.path.join(dir.out.plassembler ,"{sample}", "plassembler_summary.tsv")
    params:
        db = dir.plassemblerdb,
        chromlen = getMinChromLength,
        outdir = os.path.join(dir.out.plassembler ,"{sample}")
    conda:
        os.path.join(dir.env,'plassembler.yaml')
    resources:
        mem_mb=config.resources.big.mem,
        time=config.resources.med.time
    threads:
        config.resources.big.cpu
    shell:
        """
        plassembler run -l {input.l} -o {params.outdir} -1 {input.r1} -2 {input.r2} -d {params.db} -t {threads} -c {params.chromlen} --skip_qc -f
        """

rule plassembler_move_fastas:
    input:
        fasta = os.path.join(dir.out.plassembler,"{sample}", "plassembler_plasmids.fasta")
    output:
        fasta = os.path.join(dir.out.plassembler_fastas, "{sample}.fasta")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        cp {input.fasta} {output.fasta}
        """
        
rule plassembler_move_summaries:
    input:
        tsv = os.path.join(dir.out.plassembler,"{sample}", "plassembler_summary.tsv")
    output:
        tsv = os.path.join(dir.out.plassembler_individual_summaries, "{sample}.tsv")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        cp {input.tsv} {output.tsv}
        """

rule aggr_plassembler:
    """Aggregate."""
    input:
        expand(os.path.join(dir.out.plassembler_fastas, "{sample}.fasta"), sample = SAMPLES),
        expand(os.path.join(dir.out.plassembler_individual_summaries, "{sample}.tsv"), sample = SAMPLES)
    output:
        flag = os.path.join(dir.out.flags, "aggr_plassembler.flag")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        touch {output.flag}
        """

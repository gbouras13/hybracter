rule plassembler:
    input:
        os.path.join(dir.out.qc,"{sample}_filt_trim.fastq.gz"),
        os.path.join(dir.out.fastp ,"{sample}_1.fastq.gz"),
        os.path.join(dir.out.fastp ,"{sample}_2.fastq.gz")
    output:
        os.path.join(dir.out.plassembler ,"{sample}", "plassembler_plasmids.fasta"),
        os.path.join(dir.out.plassembler ,"{sample}", "plassembler_copy_number_summary.tsv")
    params:
        db = dir.dbs.plassemblerdb,
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
        plassembler.py -l {input[0]} -o {params.outdir} -1 {input[1]} -2 {input[1]} -d {params.db} -m 500 -t {threads} -c {params.chromlen} --skip_qc -f
        """

rule plassembler_move_fastas:
    input:
        os.path.join(dir.out.plassembler,"{sample}", "plassembler_plasmids.fasta")
    output:
        os.path.join(dir.out.plassembler_fastas, "{sample}.fasta")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        cp {input[0]} {output[0]}
        """
        
rule plassembler_move_summaries:
    input:
        os.path.join(dir.out.plassembler,"{sample}", "plassembler_copy_number_summary.tsv")
    output:
        os.path.join(dir.out.plassembler_individual_summaries, "{sample}.tsv")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        cp {input[0]} {output[0]}
        """

rule aggr_plassembler:
    """Aggregate."""
    input:
        expand(os.path.join(dir.out.plassembler_fastas, "{sample}.fasta"), sample = SAMPLES),
        expand(os.path.join(dir.out.plassembler_individual_summaries, "{sample}.tsv"), sample = SAMPLES)
    output:
        os.path.join(dir.out.flags, "aggr_plassembler.flag")
    resources:
        mem_mb=config.resources.sml.mem,
        time=config.resources.sml.time
    threads:
        config.resources.sml.cpu
    shell:
        """
        touch {output[0]}
        """

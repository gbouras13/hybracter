"""
dnaapler no medaka
"""


rule dnaapler_no_medaka:
    """
    Runs dnaapler to begin chromosome with dnaa when medaka isn't called
    """
    input:
        fasta=os.path.join(dir.out.chrom_pre_polish, "{sample}_chromosome.fasta"),
        ignore_list=os.path.join(dir.out.chrom_pre_polish, "{sample}_ignore_list.txt"),
    output:
        fasta=os.path.join(dir.out.dnaapler, "{sample}", "{sample}_reoriented.fasta"),
        version=os.path.join(dir.out.versions, "{sample}", "dnaapler.version"),
    conda:
        os.path.join(dir.env, "dnaapler.yaml")
    params:
        dir=os.path.join(dir.out.dnaapler, "{sample}"),
    resources:
        mem_mb=config.resources.med.mem,
        mem=str(config.resources.med.mem) + "MB",
        time=config.resources.med.time,
    threads: config.resources.med.cpu
    benchmark:
        os.path.join(dir.out.bench, "dnaapler", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "dnaapler", "{sample}.log"),
    shell:
        """
        dnaapler all -i {input.fasta} -o {params.dir} --ignore {input.ignore_list} -p {wildcards.sample} -t {threads} -a nearest --db dnaa,repa,cog1474 -f 2> {log}
        dnaapler --version > {output.version}
        """


rule dnaapler_combine_assembly_with_plasmids_assembly:
    """
    extracts the reoriented prepolished chrom assembly
    """
    input:
        chromosome_fasta=os.path.join(
            dir.out.dnaapler, "{sample}", "{sample}_reoriented.fasta"
        ),
        plasmid_fasta=os.path.join(
            dir.out.plassembler, "{sample}", "plassembler_plasmids.fasta"
        ),
    output:
        combined_fasta=os.path.join(
            dir.out.dnaapler,
            "{sample}",
            "{sample}_reoriented_chromosome_plasmids.fasta",
        ),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        cat {input.chromosome_fasta} {input.plasmid_fasta} > {output.combined_fasta}
        """

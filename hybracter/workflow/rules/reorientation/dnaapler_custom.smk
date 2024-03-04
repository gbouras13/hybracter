"""
dnaapler custom
"""


rule dnaapler_custom:
    """
    Runs dnaapler to begin  with custom db
    """
    input:
        fasta=os.path.join(
            dir.out.intermediate_assemblies, "{sample}", "{sample}_medaka_rd_1.fasta"
        ),
        ignore_list=os.path.join(dir.out.chrom_pre_polish, "{sample}_ignore_list.txt"),
    output:
        fasta=os.path.join(dir.out.dnaapler, "{sample}", "{sample}_reoriented.fasta"),
        version=os.path.join(dir.out.versions, "{sample}", "dnaapler.version"),
    conda:
        os.path.join(dir.env, "dnaapler.yaml")
    params:
        dir=os.path.join(dir.out.dnaapler, "{sample}"),
        custom_db=config.args.dnaapler_custom_db,
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
        dnaapler all -i {input.fasta} -o {params.dir} --ignore {input.ignore_list} -p {wildcards.sample} -t {threads} -a nearest -c {params.custom_db}  -f 2> {log}
        dnaapler --version > {output.version}
        """

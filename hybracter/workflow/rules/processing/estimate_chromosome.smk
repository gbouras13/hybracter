"""
Michael Hall's new tool lrge is awesome
Thanks to Michael Hall for some code to determine the best way to run it too
"""


rule lrge:
    input:
        fastq=get_input_lr_fastqs,
    output:
        lrge=os.path.join(
            dir.out.chrom_size, "{sample}_lrge_estimated_chromosome_size.txt"
        ),
    conda:
        os.path.join(dir.env, "lrge.yaml")
    threads: config.resources.med.cpu
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    benchmark:
        os.path.join(dir.out.bench, "lrge", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "lrge", "{sample}.log"),
    shell:
        """
        if [[ {input.fastq} =~ \.gz$ ]]; then
        INPUT_CMD="zcat {input.fastq}"
        else
        INPUT_CMD="cat {input.fastq}"
        fi

        # Now run the awk command on the decompressed file
        if $INPUT_CMD | awk '{{if (NR % 4 == 1) count++}} END {{exit count > 10000 ? 0 : 1}}'; then
            # There are more than 10,000 reads, run lrge with defaults
            lrge -t {threads} -s 42 {input.fastq} -o {output.lrge}
        else
            # There are less than 10,000 reads, use the all-vs-all strategy with all available reads
            lrge -t {threads} -s 42 -n 10000 {input.fastq} -o {output.lrge}
        fi

        """


rule aggr_lrge:
    """
    aggregates over all samples
    """
    input:
        expand(
            os.path.join(
                dir.out.chrom_size, "{sample}_lrge_estimated_chromosome_size.txt"
            ),
            sample=SAMPLES,
        ),
    output:
        flag=os.path.join(dir.out.flags, "aggr_lrge.flag"),
    resources:
        mem_mb=config.resources.sml.mem,
        mem=str(config.resources.sml.mem) + "MB",
        time=config.resources.sml.time,
    threads: config.resources.sml.cpu
    shell:
        """
        touch {output.flag}
        """

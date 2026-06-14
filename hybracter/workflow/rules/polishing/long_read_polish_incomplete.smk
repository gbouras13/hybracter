rule medaka_incomplete:
    input:
        fasta=os.path.join(dir.out.incomp_pre_polish, "{sample}.fasta"),
        fastq=os.path.join(dir.out.qc, "{sample}_filt_trim.fastq.gz"),
    output:
        fasta=os.path.join(dir.out.medaka_incomplete, "{sample}", "consensus.fasta"),
        version=os.path.join(dir.out.versions, "{sample}", "medaka_incomplete.version"),
        copy_fasta=os.path.join(
            dir.out.intermediate_assemblies_incomplete,
            "{sample}",
            "{sample}_medaka.fasta",
        ),
    conda:
        os.path.join(dir.env, "medaka.yaml")
    params:
        model=MEDAKA_MODEL,
        bacteria=BACTERIA,
        dir=os.path.join(dir.out.medaka_incomplete, "{sample}"),
        bam=os.path.join(dir.out.medaka_incomplete, "{sample}", "calls_to_draft.bam"),
        hdf=os.path.join(dir.out.medaka_incomplete, "{sample}", "consensus_probs.hdf"),
    resources:
        mem_mb=config.resources.big.mem,
        time=config.resources.med.time,
    threads: config.resources.big.cpu
    benchmark:
        os.path.join(dir.out.bench, "medaka_incomplete", "{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "medaka_incomplete", "{sample}.log"),
    shell:
        """
        export CUDA_VISIBLE_DEVICES=""
        if [ "{params.bacteria}" = "True" ]; then
           # from v 0.10.0, hybracter supports --bacteria
            if medaka tools resolve_model --auto_model consensus_bacteria {input.fastq} > /dev/null 2>> {log}; then
                medaka_consensus -i {input.fastq} -d {input.fasta} -o {params.dir} --bacteria -t {threads} 2>> {log}
            else
                echo "hybracter: medaka --bacteria cannot auto-select a model (reads carry no basecaller model tag); falling back to -m r1041_e82_400bps_bacterial_methylation" | tee -a {log}
                medaka_consensus -i {input.fastq} -d {input.fasta} -o {params.dir} -m r1041_e82_400bps_bacterial_methylation -t {threads} 2>> {log}
            fi
            medaka --version > {output.version}
            touch {params.bam}
            rm {params.bam}
            touch {params.hdf}
            rm {params.hdf}
        else
           
            medaka_consensus -i {input.fastq} -d {input.fasta} -o {params.dir} -m {params.model} -t {threads} 2> {log}
            medaka --version > {output.version}
            touch {params.bam}
            rm {params.bam}
            touch {params.hdf}
            rm {params.hdf}
        fi
        cp {output.fasta} {output.copy_fasta}
        """

# hybracter
Hybrid (and long-only) bacterial assembly pipeline for many isolates using Snakemake and Snaketool

```
git clone "https://github.com/gbouras13/hybracter.git"
cd hybracter/
pip install -e .
hybracter --help
hybracter run --help
```

Designed for ONT long reads with paired end short reads on the same bacterial isolate.

Pipeline
==========

1. Long-read QC: filtlong for quality and porechop for adapters (qc.smk)
2. Long-read only assembly with Flye (assemble.smk)
3. Combine flye assembly stats (assembly_statistics.smk)
4. Extract chromosome if over minChromLength (extract_fastas.smk)
5. Polish chromosome with medaka, reorient with dnaapler so it begins with dnaA gene, then polish again with medaka (long_read_polish.smk)

If you have short reads and step 4 succeeded in getting a full chromosome:

6. Polish with short reads with Polypolish (short_read_polish.smk).

If you choose to polish with polca:

7. Polish with short reads with polca (short_read_polish.smk).

If you have short reads and step 4 succeeded in getting a full chromosome:

8. Use plassembler to assembly plasmids and calculate copy number (plassembler.smk).
9. Combine plassembler output for each sample into one file (combine_plassembler_info.smk).

To be Done
=========

* Option to allow for incomplete assemblies - would require steps 4 -> 7 to proceed without extracting chromosome.

Input
=======

* hybracter requires an input csv file with 5 columns if running hybrid, or 3 columns if running in long-only mode
* Each row is a sample
* Column 1 is the sample name, column 2 is the long read ONT fastq file , column 3 is the minChromLength for that sample
* For hybrid, column 4 is the R1 short read fastq file, column 5 is the R2 short read fastq file.
* minChromLength is the minimum chromosome length (bp) that hybracter will consider as a "complete" chromosome for that sample. It is required to be an integer e.g. 2000000 for 2Mbp. 
* No headers

e.g.

Hybrid

sample1,sample1_long_read.fastq.gz,2000000,sample1_SR_R1.fastq.gz,sample1_SR_R2.fastq.gz
sample2,sample2_long_read.fastq.gz,1500000,sample2_SR_R1.fastq.gz,sample2_SR_R2.fastq.gz

Long Only 

sample1,sample1_long_read.fastq.gz,2000000\n
sample2,sample2_long_read.fastq.gz,1500000


Usage
=======

```
hybracter run --input <input.csv> --output <output_dir> --threads <threads> 
```

hybracter will run in hybrid mode by default. To specify long-read only, specify `--long`

```
hybracter run --input <input.csv> --output <output_dir> --threads <threads> --long
```

hybracter will not use polca as a polisher by default. To use polca specify `--polca`

```
hybracter run --input <input.csv> --output <output_dir> --threads <threads> --polca
```

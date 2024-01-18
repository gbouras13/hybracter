`hybracter` creates a number of output files in different formats.

# Main Output

The main outputs are in the `FINAL_OUTPUT` directory.

This directory will include:

## Summary File

1. `hybracter_summary.tsv` file. This gives the summary statistics for your assemblies with the following columns:

|Sample |Complete (True or False) | Total_assembly_length |	Number_of_contigs |	Most_accurate_polishing_round |	Longest_contig_length | Longest_contig_coverage|	Number_circular_plasmids |
|--------|-----------------------|-------------------------|-------------------|--------|--|--|--|


## Summary Assemblies

2. The `complete` and `incomplete` directories will contain the summary assemblies for all samples.

All samples that are denoted by hybracter to be complete will have 5 outputs in the `complete` directory:

   * `sample`_summary.tsv containing the summary statistics for that sample.
   * `sample`_per_contig_stats.tsv containing the contig names, lengths, GC% and whether the contig is circular.
   * `sample`_final.fasta containing the final assembly for that sample.
   * `sample`_chromosome.fasta containing only the final chromosome(s) assembly for that sample.
   * `sample`_plasmid.fasta containing only the final plasmid(s) assembly for that sample. Note this may be empty. If this is empty, then that sample had no plasmids. 

All samples that are denoted by hybracter to be incomplete will have 3 outputs in the `incomplete` directory:

   * `sample`_summary.tsv containing the summary statistics for that sample.
   * `sample`_per_contig_stats.tsv containing the contig names, lengths, GC% and whether the contig is circular.
   * `sample`_final.fasta containing the final assembly for that sample.

# Other Outputs

## `supplementary_results` directory

The `supplementary_results` directory contains a number of supplementary results that you might find useful:

##### 1. `comparisons` directory

* This directory contains visual representations comparing the effect of each polishing round for each sample using a modified version of Ryan Wick's [compare_assemblies.py script](https://github.com/rrwick/Perfect-bacterial-genome-tutorial/blob/main/scripts/compare_assemblies.py). An example is below 

```
    contig_1 37368-37398: ACCATTTTTGTTTTATTTTTTGTAAAGACAC
    contig_1 37368-37397: ACCATTTTTGTTTTA-TTTTTGTAAAGACAC
                                         *              

    contig_1 43247-43277: CAACGTTGTTTTCCCTGAGCCTAAATAACCA
    contig_1 43246-43276: CAACGTTGTTTTCCCCGAGCCTAAATAACCA
                                         *              

    contig_1 44658-44688: CTTGATCTTTATCTATGATTTCATTAATACT
    contig_1 44657-44687: CTTGATCTTTATCTACGATTTCATTAATACT
                                         *              
```

* If this file is empty, there are no differences between assemblies

##### 2. `intermediate_chromosome_assemblies` directory

* This directory contains intermediate chromosome assemblies for all polishing rounds for each sample.

##### 3. `flye_individual_summaries` directory

* This directory contains individual sample summaries from Flye for all samples.

##### 4. `plassembler_individual_summaries` directory

* This directory contains individual sample summaries from Plassembler for each sample.

##### 5. `plassembler_all_assembly_summary` directory

* This directory contains individual sample summaries from Plassembler for all samples.

##### 6. `pyrodigal_mean_length_summaries` directory

*  For `long`, this directory contains pyrodigal mean CDS length summary files for each polishing round for each sample.

##### 7. `pyrodigal_mean_length_summaries_plassembler` directory

*  For `long`, this directory contains pyrodigal mean CDS length summary files for each polishing round for each sample for the plassembler assembled plasmids.

## `processing` directory

The `processing` directory will contain a number of intermediate directories whose information you might find useful:

##### 1. `flye` directory

* This directory will contain the Flye assembly output and associated intermediate files for each sample

##### 2. `qc` directory

* This directory will contain the filtered, trimmed and contaminant removed FASTQ reads (where applicable) for each sample. From v0.6.0, it will also contain `seqkit` directory with Seqkit outputs for each read set, and `coverage` with calculated quick short read coverage (for `hybracter hybrid`).

##### 3. `plassembler` directory

* This directory will contain the Plassembler assembly output and associated intermediate files for each sample

##### 4. `chrom_pre_polish` directory

  * This directory will contain the pre-polished chromosome assemblies for complete isolates

##### 5. `complete` and `incomplete` directories

  * These directories will contain the medaka, polypolish and pypolca polishing and dnaapler reorientation intermediate files for each sample

##### 6. `ale_out_files` directory

  * For `hybrid`, this directory will intermediate ALE files for each assembly polishing round internal to `hybracter` (so can be ignored).

##### 7. `ale_scores_complete` and `ale_scores_incomplete`  directories

  *  These directories will containin ALE scores for each assembly polishing round.

## `stderr` directory

* This will contain log files for each program in `hybracter`.

## `versions` directory

* This will contain the specific versions used for each program in `hybracter`.

## `flags` directory

* This will contain flag files internal to `hybracter` (so can be ignored).

## `completeness` directory

* This will contain flag files internal to determine completeness internal to `hybracter` (so can be ignored).

## `benchmarks` directory

* This will contain benchmarking time and memory usage statistics for each program in hybracter.


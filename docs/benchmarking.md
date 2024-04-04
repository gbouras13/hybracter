# `hybracter` Benchmarks

## Summary of Main Results

* `hybracter hybrid` assembled the most accurate chromosomes overall, particularly in terms of insertions and deletions (Dragonflye hybrid was second, Unicycler was the worst). For long reads only, `hybracter long` was superior to dragonflye long.
* `hybracter hybrid` was the best tool at recovering plasmids, recovering 65/69 possibly plasmids completely, while Unicycler recovered 60/69. Accuracy for both was similar. Dragonflye hybrid recovered only 44/69. 
* `hybracter long` recovered 60/69 plasmids from only long read data, comparable to the hybrid methods. Dragonflye recovered only 44/69 plasmids. 

## Recommendations

* `Hybracter hybrid` is superior to Unicycler in terms of accuracy, time taken and (slighly) in terms of plasmid recovery. It should be preferred to Unicycler.
* You should use `hybracter long` if you care about plasmids and have only long reads. It performs similarly to hybrid methods and its inclusion of [Plassembler](https://github.com/gbouras13/plassembler) largely solves the [problem of long read assemblers recovering small plasmids](https://doi.org/10.1099/mgen.0.001024).
* `Hybracter` in both modes is inferior to Dragonflye in terms of time and though better in terms of chromosome accuracy. 
* If you want the fastest possible chromosome assemblies for applications like species ID or sequence typing that retain a high level of accuracy, Dragonflye is a good option.
* Dragonflye should not be used if you care about recovering plasmids.

## Sample Selection

For the full details on how `hybracter` was benchmarked, please see the [repository](https://github.com/gbouras13/hybracter_benchmarking). You can find all the benchmarking output and FASTQ files used [here](https://doi.org/10.5281/zenodo.10906937).

`Hybracter` v0.7.0 was benchmarked in both hybrid and long modes (specifically using the `hybrid-single` and `long-single` commands) against [Unicycler](https://github.com/rrwick/Unicycler) v0.5.0 and [Dragonflye](https://github.com/rpetit3/dragonflye) v1.1.2.

For Dragonflye, both 'long' (default assembly with Flye, default polishing with Racon) and 'hybrid' (default assembly with Flye, default polishing with Racon followed by Polypolish ) were tested.

`Hybracter` was benchmarked on a panel of 30 isolates from 5 studies with accompanying reference genomes.

These samples contain varying levels of long-read quality (reflecting improvements in Oxford Nanopore Technologies long read technology), with the median Q score of long read sets ranging from a minimum of 10.6 to maximum of 26.8. The four studies are:

1.	Five ATCC strain isolates (ATCC-10708 _Salmonella enterica_, ATCC-17802 _Vibrio paragaemolyticus_, ATCC-25922 _Escherichia coli_, ATCC-33560 _Campylobacter jejuni_ and ATCC-BAA-679 _Listeria monocytogenes_ ) with R10 chemistry super-accuracy model basecalled simplex long-reads sequenced by Louise Judd from Ryan Wick's [blogpost](https://rrwick.github.io/2023/05/05/ont-only-accuracy-with-r10.4.1.html) and made available as a part of the Hybracter preprint.
2.	The same 5 ATCC isolates with R10 chemistry fast model basecalled long-reads, and R10 chemistry super-accuracy model basecalled duplex long-reads from [Hall et al](https://www.biorxiv.org/content/10.1101/2024.03.15.585313v2).
3.	Twelve diverse carbapenemase-producing Gram negative bacteria from [Lerminiaux et al](https://doi.org/10.1139/cjm-2023-0175).
4.	Staphylococcus aureus JKD6159 with both R9 and R10 chemistry long read sets from [Wick et al](https://doi.org/10.1128/mra.01129-22).
5.	Mycobacterium tuberculosis HR37v from [Chitale et al](https://www.nature.com/articles/s41467-022-34853-x).

## Wallclock Time Results

To compare the performance of Hybracter, we compared wall-clock runtime consumption on a machine with an Intel® Core™ i9-13900K CPU @ 5.60 GHz on a machine running Ubuntu 20.04.6 LTS with a total of 32 available threads (24 cores). We ran all tools with 8 and 16 threads to reflect performance and runtimes that are comparable to commonly available consumer hardware, with a cap at 32GB of RAM requested for each job. Hybracter hybrid and long were run with ‘hybracter hybrid-single’ and ‘hybracter long-single’ for each isolate in order to generate a comparable per sample runtime for comparison with the other tools. 

| **Tool**              | **Type** | **8 Threads**                                                           | **16 Threads**                                                          |
| --------------------- | -------- | ----------------------------------------------------------------------- | ----------------------------------------------------------------------- |
| **Hybracter hybrid**  | Hybrid   | Median = 00h:15m:03s Min = 00h:04m:29s Max = 00h:54m:41s | Median = 00h:13m:44s Min = 00h:03m:27s Max = 00h:44m:36s      |
| **Dragonflye hybrid** | Hybrid   | Median = 00h:04m:34s Min = 00h:01m:32s Max = 00h:07m:27s      | Median = 00h:03m:46s Min = 00h:01m:22s Max = 00h:06m:01s      |
| **Unicycler**         | Hybrid   | Median = 00h:50m:25s Min =  00h:12m:04s Max = 01h:13m:32s     | Median = 00h:34m:20s Min = 00h:08m:36s Max = 00h:48m:23s       |
| **Hybracter long**    | Long     | Median = 00h:11m:46s Min = 00h:03m:26s Max = 00h:36m:09s | Median = 00h:10m:20s Min = 00h:03m:17s Max = 00h:29m:50s |
| **Dragonflye long**   | Long     | Median = 00h:04m:10s Min = 00h:01m:22s Max = 00h:06m:01s      | Median = 00h:04m:34s Min = 00h:01m:32s Max = 00h:07m:27s |

## Chromosome Accuracy Results

The chromosome assemblies for all samples were complete and circularised, and were compared against the reference with [Dnadiff](https://github.com/marbl/MUMmer3/blob/master/docs/dnadiff.README) v1.3. Small (<60bp) insertions and deletions (indels), small nucleotide variants (SNVs) and large (>60bps) were compared. 

The summary results are presented below:

| **Tool**              | **Type** | **Small Indels**                         | **SNVs**                                 | **Small Indels + SNVs**               | **Large Indels**                                           |
| --------------------- | -------- | --------------------------------------------------- | --------------------------------------------------- | ------------------------------------------------ | ---------------------------------------------------------- |
| **Hybracter hybrid**  | Hybrid   | Median = 0 Min = 0 Max = 41    | Median = 0 Min = 0 Max = 26        | Median = 1 Min = 0 Max = 67 | Total = 9 Median = 0 Min = 0 Max = 2 |
| **Dragonflye hybrid** | Hybrid   | Median = 2.5 Min = 0 Max = 112   | Median = 0 Min = 0  Max = 64    | Median = 4.5  Min = 0  Max = 154     | Total = 70 Median = 2 Min = 0  Max = 12      |
| **Unicycler**         | Hybrid   | Median = 11 Min = 0  Max = 125      | Median = 34 Min = 0 Max = 165       | Median = 57.5 Min = 3 Max = 290   | Total = 87 Median = 1  Min = 0  Max = 16     |
| **Hybracter long**    | Long     | Median = 16 Min = 1 Max = 743 | Median = 21.5 Min = 0 Max = 156 | Median = 54 Min = 1 Max = 852   | Total = 11 Median = 1  Min = 0  Max = 3        |
| **Dragonflye long**   | Long     | Median = 125 Min = 2 Max = 4814       | Median = 34.5 Min = 0 Max = 2172 | Median = 170.5 Min = 2 Max = 6332    | Total = 68 Median = 2 Min = 0 Max = 12    |


## Plasmid Recovery Results

All samples were analysed using the 4-step approach outlined below using summary length and GC% statistics for all contigs and the output of Dnadiff v1.3 comparisons generated for each sample and tool combination against the reference genome plasmids. 

1.	The number of circularised plasmid contigs recovered for each isolate was compared to the reference genome. If the tool recovered a circularised contig homologous to that in the reference (based on the genome lengths, GC% and Dnadiff output), it was denoted as completely recovered. For Dragonflye assemblies, some plasmids were multiplicated due to known issues with the long-read first assembly approach for small plasmids. Any circularised contigs that were multiplicated compared to the reference plasmid were therefore denoted as misassembled. 
2.	For additional circularised contigs not found in the reference recovered, these were tested for homology with NCBI nt database using the web version of blastn. If there was a hit to a plasmid, the Plassembler output within Hybracter was checked for Mash hits to the PLSDB. If there was a hit there too, the contig was denoted as an additional recovered plasmid. There were 2 in total, both from Lerminiaux Isolate G.
3.	Plasmids with contigs that were either not circularised but homologous to a reference plasmid, or circularised but incomplete (based on genome length and Dnadiff output) were denoted as partially recovered or misassembled. 
4.	Reference plasmids without any homologous contigs in the assembly were denoted as missed.

The summary results are presented below:

| **Tool**              | **Complete Plasmids Recovered** | **Total Plasmids Partially Recovered or Misassembled** | **Total Plasmids Missed** | **Additional Plasmids Recovered not in Reference** | **Additional Non-Plasmid Contigs Recovered** |
| --------------------- | ------------------------------- | ------------------------------------------------------ | ------------------------- | -------------------------------------------------- | -------------------------------------------- |
| **Hybracter hybrid**  | 65                              | 4                                                      | 0                         | 2                                                  | 10                                            |
| **Unicycler**         | 60                              | 6                                                      | 3                         | 1                                                  | 2                                            |
| **Dragonflye hybrid** | 44                              | 16                                                     | 9                         | 1                                                  | 10                                            |
| **Hybracter long**    | 60                              | 5                                                      | 4                         | 2                                                  | 3                                            |
| **Dragonflye long**   | 44                              | 16                                                     | 9                         | 1                                                  | 10                                            |

## Plasmid Accuracy Results

The plasmids assemblies for all samples were complete and circularised, and were compared against the reference with [Dnadiff](https://github.com/marbl/MUMmer3/blob/master/docs/dnadiff.README) v1.3. Small (<60bp) insertions and deletions (indels), small nucleotide variants (SNVs) and large (>60bps) were compared. SNVs and small indels were normalised per 100kbp of reference sequence.

| Tool              | Type   | Min SNVs per 100kbp | Median SNVs per 100kbp | Max SNVs per 100kbp | Min Small Indels per 100kbp | Median Small Indels per 100kbp | Max Small Indels per 100kbp | Min SNVs + Small Indels per 100kbp | Median SNVs + Small Indels per 100kbp | Max SNVs + Small Indels per 100kbp | Min  Large Indels | Median  Large Indels | Max  Large Indels | Total Large Indels |
| ----------------- | ------ | ----------------------- | ---------------------- | ----------------------- | ------------------------------- | ------------------------------ | ------------------------------- | -------------------------------------- | ------------------------------------- | -------------------------------------- | --------------------- | -------------------- | --------------------- | ------------------ |
| dragonflye_hybrid | hybrid | 0                       | 0                      | 0.7                       | 0                               | 0                              | 31.75                           | 0                                      | 0.255                                     | 31.75                                  | 1                     | 5.5                  | 18                    | 116                |
| dragonflye_long   | long   | 0                       | 1.285                   | 55.29                   | 0                             | 6.255                           | 129.07                           | 0                                    | 7.655                                 | 167.17                                  | 1                     | 5.5                  | 18                    | 116                |
| hybracter_hybrid  | hybrid | 0                       | 0.255                    | 46.19                   | 0                               | 0                              | 16.19                           | 0                                      | 1.62                                  | 47.34                                  | 0                     | 2                    | 9                     | 44                 |
| hybracter_long    | long   | 0                       | 4.47                  | 57.87                   | 0                               | 2.995                           | 189.57                          | 0                                      | 8.735                                | 227.89                                 | 0                     | 1                    | 27                     | 59                 |
| unicycler         | hybrid | 0                       | 0.855                  | 10.99                   | 0                               | 0                              | 3.45                            | 0                                      | 2.02                                  | 11.28                                  | 0                     | 2                    | 16                    | 63                 |
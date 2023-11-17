# `hybracter` Benchmarks

## Summary of Main Results

* `hybracter` assembled the most accurate chromosomes overall, particularly in terms of insertions and deletions (Dragonflye hybrid was second, Unicycler was the worst). For long reads only, `hybracter` was superior to dragonflye long.
* `hybracter hybrid` was the best tool at recovering plasmids, recovering 57/59 possibly plasmids completely, while Unicycler recovered 54/59. Accuracy for both was similar. Dragonflye hybrid recovered only 34/59. 
* `hybracter long` recovered 54/59 plasmids from only long read data, comparable to the hybrid methods. Dragonflye recovered only 35/59 plasmids. 

## Recommendations

* `Hybracter hybrid` is superior to Unicycler in terms of accuracy, time taken and (slighly) in terms of plasmid recovery. It should be preferred to Unicycler.
* You should use `hybracter long` if you care about plasmids and have only long reads. It performs similarly to hybrid methods and its inclusion of [Plassembler](https://github.com/gbouras13/plassembler) largely solves the [problem of long read assemblers recovering small plasmids](https://doi.org/10.1099/mgen.0.001024).
* `Hybracter` in both modes is inferior to Dragonflye in terms of time and though better in terms of chromosome accuracy. 
* If you want the fastest possible chromosome assemblies for applications like species ID or sequence typing that retain a high level of accuracy, Dragonflye is a good option.
* Dragonflye should not be used if you care about recovering plasmids.

## Sample Selection

For the full details on how `hybracter` was benchmarked, please see the [repository](https://github.com/gbouras13/hybracter_benchmarking).

`Hybracter` was benchmarked in both hybrid and long modes (specifically using the `hybrid-single` and `long-single` commands) against [Unicycler](https://github.com/rrwick/Unicycler) v0.5.0 and [Dragonflye](https://github.com/rpetit3/dragonflye) v1.1.2.

For Dragonflye, both 'long' (default assembly with Flye,, default polishing with Racon) and 'hybrid' (default assembly with Flye, default polishing with Racon followed by Polypolish ) were tested.

`Hybracter` was benchmarked on a panel of 20 isolates form 4 studies with accompanying reference genomes.

These samples contain varying levels of long-read quality (reflecting improvements in Oxford Nanopore Technologies long read technology), with the median Q score of long read sets ranging from a minimum of 12.3 to maximum of 18.3. The four studies are:

1.	Five ATCC strain isolates (ATCC-10708 Salmonella enterica, ATCC-17802 Vibrio paragaemolyticus, ATCC-25922 Escherichia coli, ATCC-33560 Campylobacter jejuni and ATCC-BAA-679 Listeria monocytogenes ) sequenced by Louise Judd from Ryan Wick's [blogpost](https://rrwick.github.io/2023/05/05/ont-only-accuracy-with-r10.4.1.html).
2.	Twelve diverse carbapenemase-producing Gram negative bacteria from [Lerminiaux et al](https://www.biorxiv.org/content/10.1101/2023.09.25.559359v1 ).
3.	Staphylococcus aureus JKD6159 with both R9 and R10 chemistry long read sets from [Wick et al](https://doi.org/10.1128/mra.01129-22).
4.	Mycobacterium tuberculosis HR37v from [Chitale et al](https://www.nature.com/articles/s41467-022-34853-x).

## Wallclock Time Results

To compare the performance of Hybracter, we compared wall-clock runtime consumption on a machine with an Intel® Core™ i7-10700K CPU @ 3.80 GHz on a machine running Ubuntu 20.04.6 LTS with a total of 16 available threads (8 cores). We ran all tools with 8 and 16 threads to reflect performance and runtimes that are comparable to commonly available consumer hardware, with a cap at 32GB of RAM requested for each job. Hybracter hybrid and long were run with ‘hybracter hybrid-single’ and ‘hybracter long-single’ for each isolate in order to generate a comparable per sample runtime for comparison with the other tools. 

| **Tool**              | **Type** | **8 Threads**                                                           | **16 Threads**                                                          |
| --------------------- | -------- | ----------------------------------------------------------------------- | ----------------------------------------------------------------------- |
| **Hybracter hybrid**  | Hybrid   | Median = 00h:54m:23s Min = 00h:14m:22s Max = 02h:01m:37s | Median = 00h:40m:19s Min = 00h:12m:43s Max = 01h:21m:05s      |
| **Dragonflye hybrid** | Hybrid   | Median = 00h:10m:55s Min = 00h:04m:02s Max = 00h:13m:34s      | Median = 00h:07m:21s Min = 00h:03m:33s Max = 00h:09m:28s      |
| **Unicycler**         | Hybrid   | Median = 02h:03m:02s Min =  00h:39m:09s Max = 02h:48m:16s     | Median = 01h:06m:8s Min = 00h:21m:38s Max = 01h:30m:38s       |
| **Hybracter long**    | Long     | Median = 00h:45m:29s Min = 00h:10m:52s Max = 01h:23m:49s | Median = 00h:34m:56s Min = 00h:09m:49s Max = 00h:59m:21s |
| **Dragonflye long**   | Long     | Median = 00h:09m:24s Min = 00h:03m:52s Max = 00h:13m:32s      | Median = 00h:07m:00s Min = 00h:03m:22s Max = 00h:08m:56s |

## Chromosome Accuracy Results

The chromosome assemblies for all samples were complete and circularised, and were compared against the reference with [Dnadiff](https://github.com/marbl/MUMmer3/blob/master/docs/dnadiff.README) v1.3. Small (<60bp) insertions and deletions (indels), small nucleotide variants (SNVs) and large (>60bps) were compared. SNVs and small indels were normalised per 100kbp of reference sequence.

The summary results are presented below:

| **Tool**              | **Type** | **Small Indels per 100kbp**                         | **SNVs per 100kbp**                                 | **Small Indels + SNVs per 100kbp**               | **Large Indels**                                           |
| --------------------- | -------- | --------------------------------------------------- | --------------------------------------------------- | ------------------------------------------------ | ---------------------------------------------------------- |
| **Hybracter hybrid**  | Hybrid   | Median = 0.05 Min = 0 Max = 12.34    | Median = 0.16 Min = 0 Max = 3.52        | Median = 0.24 Min = 0 Max = 15.79 | Total = 59 Median = 0 Min = 0 Max = 26 |
| **Dragonflye hybrid** | Hybrid   | Median = 0.49 Min = 0 Max = 11.56   | Median = 0.03  Min = 0  Max = 2.21    | Median = 0.74  Min = 0  Max = 13.41     | Total = 91 Median = 1.5 Min = 1  Max = 29      |
| **Unicycler**         | Hybrid   | Median = 0.28 Min = 0  Max = 9.5      | Median = 1.25 Min = 0.25 Max = 4.13       | Median = 1.49 Min = 0.43 Max = 13.62   | Total = 134 Median = 2.5  Min = 0  Max = 41     |
| **Hybracter long**    | Long     | Median = 0.49 Min = 0.06 Max = 24.82 | Median = 1.07 Min = 0.07 Max = 10.46 | Median = 2.08 Min = 0.37 Max = 35.29   | Total = 66 Median = 1  Min = 0  Max = 27        |
| **Dragonflye long**   | Long     | Median = 3.01 Min = 1.61 Max = 43.8       | Median = 0.99 Min = 0.33 Max = 10.86 | Median = 3.81 Min = 2.01 Max = 53.1    | Total = 92 Median = 2 Min = 1 Max = 30    |


## Plasmid Recovery Results

All samples were analysed using the 4-step approach outlined below using summary length and GC% statistics for all contigs and the output of Dnadiff v1.3 comparisons generated for each sample and tool combination against the reference genome plasmids. 

1.	The number of circularised plasmid contigs recovered for each isolate was compared to the reference genome. If the tool recovered a circularised contig homologous to that in the reference (based on the genome lengths, GC% and Dnadiff output), it was denoted as completely recovered. For Dragonflye assemblies, some plasmids were multiplicated due to known issues with the long-read first assembly approach for small plasmids. Any circularised contigs that were multiplicated compared to the reference plasmid were therefore denoted as misassembled. 
2.	For additional circularised contigs not found in the reference recovered, these were tested for homology with NCBI nt database using the web version of blastn. If there was a hit to a plasmid, the Plassembler output within Hybracter was checked for Mash hits to the PLSDB. If there was a hit there too, the contig was denoted as an additional recovered plasmid. There were 2 in total, both from Lerminiaux Isolate G.
3.	Plasmids with contigs that were either not circularised but homologous to a reference plasmid, or circularised but incomplete (based on genome length and Dnadiff output) were denoted as partially recovered or misassembled. 
4.	Reference plasmids without any homologous contigs in the assembly were denoted as missed.

The summary results are presented below:

| **Tool**              | **Complete Plasmids Recovered** | **Total Plasmids Partially Recovered or Misassembled** | **Total Plasmids Missed** | **Additional Plasmids Recovered not in Reference** | **Additional Non-Plasmid Contigs Recovered** |
| --------------------- | ------------------------------- | ------------------------------------------------------ | ------------------------- | -------------------------------------------------- | -------------------------------------------- |
| **Hybracter hybrid**  | 57                              | 2                                                      | 0                         | 2                                                  | 6                                            |
| **Unicycler**         | 54                              | 2                                                      | 3                         | 1                                                  | 1                                            |
| **Dragonflye hybrid** | 34                              | 16                                                     | 9                         | 1                                                  | 7                                            |
| **Hybracter long**    | 54                              | 3                                                      | 2                         | 2                                                  | 3                                            |
| **Dragonflye long**   | 35                              | 16                                                     | 8                         | 1                                                  | 5                                            |

## Plasmid Accuracy Results

The plasmids assemblies for all samples were complete and circularised, and were compared against the reference with [Dnadiff](https://github.com/marbl/MUMmer3/blob/master/docs/dnadiff.README) v1.3. Small (<60bp) insertions and deletions (indels), small nucleotide variants (SNVs) and large (>60bps) were compared. SNVs and small indels were normalised per 100kbp of reference sequence.

| Tool              | Type   | Min SNVs per 100kbp | Median SNVs per 100kbp | Max SNVs per 100kbp | Min Small Indels per 100kbp | Median Small Indels per 100kbp | Max Small Indels per 100kbp | Min SNVs + Small Indels per 100kbp | Median SNVs + Small Indels per 100kbp | Max SNVs + Small Indels per 100kbp | Min  Large Indels | Median  Large Indels | Max  Large Indels | Total Large Indels |
| ----------------- | ------ | ----------------------- | ---------------------- | ----------------------- | ------------------------------- | ------------------------------ | ------------------------------- | -------------------------------------- | ------------------------------------- | -------------------------------------- | --------------------- | -------------------- | --------------------- | ------------------ |
| dragonflye_hybrid | hybrid | 0                       | 0                      | 9                       | 0                               | 0                              | 21.82                           | 0                                      | 0                                     | 21.82                                  | 1                     | 7.5                  | 16                    | 129                |
| dragonflye_long   | long   | 0                       | 1.38                   | 12.86                   | 1.3                             | 7.65                           | 33.77                           | 2.6                                    | 9.215                                 | 33.77                                  | 2                     | 6.5                  | 16                    | 123                |
| hybracter_hybrid  | hybrid | 0                       | 3.1                    | 46.19                   | 0                               | 0                              | 30.25                           | 0                                      | 4.15                                  | 47.34                                  | 0                     | 2                    | 9                     | 39                 |
| hybracter_long    | long   | 0                       | 4.635                  | 58.08                   | 0                               | 2.22                           | 125.42                          | 0                                      | 10.635                                | 130.25                                 | 0                     | 2                    | 7                     | 32                 |
| unicycler         | hybrid | 0                       | 2.405                  | 10.99                   | 0                               | 0                              | 3.45                            | 0                                      | 3.83                                  | 11.28                                  | 0                     | 2                    | 16                    | 51                 |
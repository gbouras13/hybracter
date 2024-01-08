# Benchmarking `hybracter` on Dorado v0.5.0 Data

## Background and Method

In December 2023, Ryan Wick updated his ongoing ONT-only accuracy analysis on his (blog](https://rrwick.github.io/2023/12/18/ont-only-accuracy-update.html) (you should read his blog if you are reading this and haven't yet!)

With the latest updates to Oxford Nanopore's basecaller Dorado, he found that the reads were were consistently Q20+ (median Q=20.5) and that when using Trycycler, there were significantly updates over earlier versions of Dorado thanks to new basecaller models - specifically 'dna_r10.4.1_e8.2_400bps_sup@v4.3.0'.

He found that - at least with Trycycler assemblies — we're pretty close to being able to generate perfect ONT-only assemblies - he found a total of 37 differences between the 'ground truth' and Trycycler (as an aside here, the closer we get to perfection, the more the ground-truths themselves come into question in my opinion - unless you use exhaustively generate the same assemblies using multiple sequencing modalities like this [paper](https://doi.org/10.1128/mra.01129-22) ).

I asked Ryan to run `hybracter` on these strains with `--no_medaka` (as he had previously found that Medaka polishing makes the assemblies worse with the newest ONT reads [here](https://rrwick.github.io/2023/10/24/ont-only-accuracy-update.html)). As the legend he is, he agreed kindly. 

He first tried to run `hybracter` on the entire raw read sets, but he found they failed the Flye step of Hybracter as there were too many long reads (~500x depth).

Accordingly he decided to run `hybracter` twice, with two different subsampling methods (both to 100x of the estimated genome size):

1. Subsampling with [Filtlong](https://github.com/rrwick/Filtlong) using '`--target_bases`. According to how Filtlong works, this will prefer both longer and higher-quality reads.
2. Subsampling randomly using [Rasusa](https://github.com/mbhall88/rasusa).

He then compared the results to the reference using the same methodology as he did for the Trycycler assemblies outlined in his [blog](https://rrwick.github.io/2023/12/18/ont-only-accuracy-update.html).

The full `hybracter` config parameters were:

```
  contaminants: none
  databases: null
  dnaapler_custom_db: none
  flyeModel: --nano-hq
  input: hybracter.csv
  log: hybracter/hybracter.log
  logic: best
  medakaModel: r1041_e82_400bps_sup_v4.2.0
  min_length: 1000
  min_quality: 9
  no_medaka: true
  output: hybracter
  single: false
  skip_qc: false
qc:
```

## Results

I wanted to look at the results from 2 angles:

1. The number of errors as calculated by Ryan.
2. The number of extra contigs - from some benchmarking of [Plassember](https://plassembler.readthedocs.io/en/latest/quality_control/), these usually indicate some impurities in the read sets of hybrid sequencing data. For long only, the same logic should apply - these would be areas of the chromosome that have low quality or other non-reference reads in the assemblies.

### Errors

The results of the errors are as follows:

| Genome                  | Trycycler errors | Hybracter (Filtlong) errors | Hybracter (Rasusa) errors |
| ----------------------- | ---------------- | --------------------------- | ------------------------- |
| _Campylobacter jejuni_    | 5                | 15                          | 16                        |
| _Campylobacter lari_      | 18               | 27                          | 40                        |
| _Escherichia coli_        | 1                | 2                           | _1328_ (see note)                           |
| _Listeria ivanovii_       | 5                | 7                           | 475                       |
| _Listeria monocytogenes_  | 0                | 2                           | 3                         |
| _Listeria welshimeri_     | 1                | 2                           | 1                         |
| _Salmonella enterica_     | 3                | 15                          | 18                        |
| _Vibrio cholerae_         | 2                | _31535_ (see note)                       | _28205_ (see note)                          |
| _Vibrio parahaemolyticus_ | 2                | 12                          | 7                         |

For some context of the _V. cholerae_  and _E. coli_ results in Ryan's words

> Importantly, the larger chromosome in _V. cholerae_ and the larger plasmid in _E. coli_ had some structural heterogeneity, so Hybracter's high error counts for those genomes aren't really a problem. For the Trycycler assemblies, I tried to use the majority variant, and Hybracter clearly settled on the other variant. I don't have an explanation for the glitch in the _L. ivanovii_ assembly (with Rasusa) - I think that may just be a Flye mistake.

Overall, excluding _Vibrio cholerae_, the 8 Trycyler assemblies had total errors of 35, while Hybracter with Filtlong had 82.

### Extra Contigs

Note: Ryan had pre-filtered the read sets to throw away all reads <10kp. Therefore, unsurprisingly, the 3 plasmids under 10kbp in E _E. coli_ and 1 in _S. enterica_ were not assembled so they are ignored. 

| Sample                    | Subsampling Method | Plasmids (Over 10kbp) in Reference | Number of Contigs | Plasmid Contigs | **Extra Non-Plasmid Contigs** |
| ------------------------- | ------------------ | ---------------------------------- | ----------------- | --------------- | ------------------------- |
| _Campylobacter jejuni_    | Filtlong           | 0                                  | 1                 | 0               | 0                         |
| _Campylobacter jejuni_    | Rasusa             | 0                                  | 3                 | 0               | 2                         |
| _Campylobacter lari_      | Filtlong           | 0                                  | 1                 | 0               | 0                         |
| _Campylobacter lari_      | Rasusa             | 0                                  | 1                 | 0               | 0                         |
| _Escherichia coli_        | Filtlong           | 2                                  | 3                 | 2               | 0                         |
| _Escherichia coli_        | Rasusa             | 2                                  | 3                 | 2               | 0                         |
| _Listeria ivanovii_       | Filtlong           | 0                                  | 1                 | 0               | 0                         |
| _Listeria ivanovii_       | Rasusa             | 0                                  | 6                 | 0               | 5                         |
| _Listeria monocytogenes_  | Filtlong           | 0                                  | 1                 | 0               | 0                         |
| _Listeria monocytogenes_  | Rasusa             | 0                                  | 4                 | 0               | 3                         |
| _Listeria welshimeri_     | Filtlong           | 0                                  | 1                 | 0               | 0                         |
| _Listeria welshimeri_     | Rasusa             | 0                                  | 5                 | 0               | 4                         |
| _Salmonella enterica_     | Filtlong           | 1                                  | 2                 | 1               | 0                         |
| _Salmonella enterica_     | Rasusa             | 1                                  | 2                 | 1               | 0                         |
| _Vibrio cholerae_         | Filtlong           | 0                                  | 2                 | 0               | 0                         |
| _Vibrio cholerae_         | Rasusa             | 0                                  | 13                | 0               | 11                        |
| _Vibrio parahaemolyticus_ | Filtlong           | 0                                  | 16                | 0               | 14                        |
| _Vibrio parahaemolyticus_ | Rasusa             | 0                                  | 2                 | 0               | 0                         |

Overall, 5 Rasusa assemblies had extra non-chromosome contigs, while only 1 Filtlong one did (_Vibrio parahaemolyticus_).

## Conclusions

1. Bacterial genome assembly is hard! While I think (along with the other benchmarking) that Hybracter is the best way of currently doing _automated_ long- and hybrid-read assemblies, it still doesn't always get it right. Structural variation and low-qualtiy reads that sneak into to the final read sets can cause some grief.
2. Nonetheless, most of the time, Hybracter is pretty good! 82 errors in the 8 genomes without structural variation is getting pretty close to perfect even with long reads only.
3. Subsampling deeply sequenced read sets with Filtlong outperforms Rasusa - both in terms of errors and also it was far less likely to assemble extra non-plasmid contigs.

Based on these results (and after getting requests for subsampling functionality), the `--subsample_depth` parameter has been added to Hybracter in v0.5.0. It defaults to 100. That value multiplied by the estimated chromosome size (`-c`) will be passed as `--target_bases` to Filtlong during Hybracter's QC phase. 

As a final word of warning: if you are really interested in small plasmid recovery and/or small plasmid copy number accuracy, you may want to disable subsampling, because `--target_bases` is not random - it prioritises the longest, highest quality reads (which may result in your small plasmid reads being dropped). 

You can turn off subsampling by specifying a very large number for `--subsample_depth` e.g. `--subsample_depth 100000`.

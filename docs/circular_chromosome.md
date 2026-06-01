# `--circular_chromosome`

## Overview

By default, `hybracter` classifies a sample as **complete** *only* if any Flye assembled contig exceeds the minimum chromosome length (`-c`) (it would probably have been more accurate if I had described this as "complete chromosome" in hindsight but alas). It *does not* enforce any circularity constraint. However, only circular contigs are reoriented with Dnaapler - linear contigs are still marked as a likely chromosome, but skipped in Dnaapler. 

This is the same default behaviour as all versions and was chosen to ensure users assembling organisms with **linear chromosomes** could still use Hybracter.

From v0.13.0 the `--circular_chromosome` flag can be specified to change this criterion such that a sample is only classified as complete if at least one contig is **both** above the minimum chromosome length **and** marked as circular by Flye. When `--circular_chromosome` is set:

- A contig that is long enough but **not** circularised by Flye → sample classified as **incomplete**.
- A contig that is long enough **and** circularised by Flye → sample classified as **complete** (and reoriented by dnaapler as normal).

## When should I use `--circular_chromosome`?

Use `--circular_chromosome` **if you know your organism has a circular chromosome and you want a stricter completeness check**.

This is a conservative flag. Its main benefit is as a safeguard against misassemblies caused by **structural heterogeneity** in your sample — they can be caused by a variety of causes, but most notably prophage induction.

### Prophage induction and assembly heterogeneity

Ryan Wick's [*magnus opus* on bacterial genome assembly](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010905)) describes in detail how structural heterogeneity at prophage integration sites can cause long-read assemblers like Flye to produce an incomplete and non-circular assembly.

The mechanism is as follows:

1. Many bacterial chromosomes carry one or more **prophages** — i.e. bacteriophage genomes integrated into the host chromosome.
2. In a bacterial culture, a subpopulation of cells may undergo **spontaneous prophage induction**, causing the phage to excise from the chromosome.
3. This creates **two versions** of the chromosome in your sequencing library:
   - The **intact** chromosome (prophage integrated).
   - The **excised** chromosome (prophage removed).
4. Because both versions are present in the long-read data, Flye sees structural variation at the prophage attachment site and can fail to circularise the chromosome, producing a linear instead, which can often still be more than the chromosome length value provided to Hybracter.

The `--circular_chromosome` flag makes `hybracter` treat such samples as **incomplete** rather than complete. This does not fix the underlying heterogeneity, but it flags the sample so you know to investigate further.

## When should I **not** use `--circular_chromosome`?

### Organisms with linear chromosomes

Do **not** use `--circular_chromosome` if you are assembling an organism with a **genuinely linear chromosome**. These include:

- **_Borrelia_ spp.** (e.g. _B. burgdorferi_, the Lyme disease agent) — have a linear chromosome and linear plasmids. Flye will never circularise these, so `--circular_chromosome` would classify every sample as incomplete.
- **_Streptomyces_ spp.** — same story

For these organisms the **default behaviour** (no flag) is desired: `hybracter` will classify the sample as complete based on length alone, and non-circular long contigs will be polished and included in the final assembly. Dnaapler will not reorient them.

### When the default is already sufficient

For the vast majority of users assembling common circular-chromosome bacteria (e.g. _E. coli_, _S. aureus_, _K. pneumoniae_, _P. aeruginosa_), the default behaviour (no flag) is perfectly reasonable and consistent with how `hybracter` has always worked - however, I now recommend using `--circular_chromosome` a strict flag if you are such a user




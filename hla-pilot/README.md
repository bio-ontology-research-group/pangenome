# First DR–DQ pilot — 16 September 2026

The pilot has run. Classical-type recovery and additional structural information were both examined. These results are exploratory; no novel HLA allele or disease association has been established.

## Classical-type recovery

The reference labels exclude the target donor and all linked pedigree members available in the metadata. Five duplicate JaSaPaGe assemblies were excluded in favour of HPRC. The Korean 076 trio and APR f/m/s trio are conservatively grouped by naming; their pedigrees have not been independently verified. Pedigree records cover 228 of 371 donors. Other undiscovered relatives remain a limitation.

The two methods are untuned nearest-neighbour assignment by Jaccard similarity. All tied nearest neighbours contribute to the candidate label set; disagreement produces an ambiguous call. No-calls count against the all-evaluated agreement denominator.

- **Bundles:** binary presence from the existing DRA–DMA decomposition. Labels are held out, but the decomposition used the whole panel. This is a transductive diagnostic, not an unseen-sample benchmark.
- **Sequence baseline:** canonical 31-mers sampled by a fixed content-hash rule, retaining exact strings. These features do not fit a vocabulary or weights on test labels. Inputs are existing annotated gene sequences plus their flanks, so this is conditional on the provided assemblies and annotations, not an end-to-end read-based typer.

Numeric two-field identity is the endpoint; expression suffixes and full-gene identity are not evaluated. Unresolved two-field labels are excluded per locus, while noncoding novelty can retain a usable two-field label. This makes the denominator smaller than the full panel. Training labels come from the old IPD-based Immuannot annotation.

| Gene | Bundle agreement, family excluded | Sequence agreement, family excluded | Bundle agreement, cohort excluded | Sequence agreement, cohort excluded |
|---|---:|---:|---:|---:|
| HLA-DRB1 | 595/724 (82.2%) | 689/724 (95.2%) | 528/724 (72.9%) | 644/724 (89.0%) |
| HLA-DQA1 | 699/739 (94.6%) | 727/740 (98.2%) | 665/739 (90.0%) | 728/740 (98.4%) |
| HLA-DQB1 | 687/736 (93.3%) | 729/737 (98.9%) | 643/736 (87.4%) | 715/737 (97.0%) |

The cohort exclusion uses the recorded analysis cohort (including HPRC ancestry subgroups); it is not a strict leave-study-out experiment. Family groups are also excluded in this comparison. Per-cohort results, call rates, rare alleles absent from training, and a training-majority baseline are in `results/metrics.tsv`.

## Experimental comparison

The [Gourraud et al. 2014 Sanger typing](https://doi.org/10.1371/journal.pone.0097282) covers antigen-recognition exons, not complete genes. We retain the published ambiguity lists and compare unordered diploid pairs. These are concordance results against an older nomenclature/database, not a perfect modern full-gene truth set. DQA1 is unavailable in this truth table.

| Gene | Assembly annotation | Bundles, family excluded | Sequence baseline, family excluded |
|---|---:|---:|---:|
| HLA-DRB1 | 59/64 (92.2%) | 42/64 (65.6%) | 55/64 (85.9%) |
| HLA-DQB1 | 58/62 (93.5%) | 52/62 (83.9%) | 58/62 (93.5%) |

The assembly baseline itself has discordances. These require inspection for assembly errors, resolution/nomenclature differences, and limitations of the experimental assay before being used to judge individual failures. The benchmark excludes unresolved assembly labels; it does not estimate performance across all novel or untyped haplotypes.

## Structural result

The signed, ordered gene-annotation screen found one three-gene combination with differing DRB5 content: **DRB1*15:03–DQA1*01:02–DQB1*06:02**. HG03139 haplotype 2 lacks a DRB5 annotation, whereas 15 other donors with that combination have one. All these entries originate from assemblies rather than upstream clipped graph paths. The screen only covers genes annotated by the existing pipeline; it is not a complete DRB pseudogene annotation.

A local assembly comparison against HG03195 haplotype 2 finds a disrupted alignment over the expected DRB5 region, with no N gaps in either extracted DR block. The control DRB5 sequence instead has a divergent main hit at the candidate’s DRB1 locus, illustrating why paralogous matches must not be counted as another DRB5 copy. This is assembly evidence warranting read validation, not a confirmed biological deletion or precise breakpoint call.

**This structural class is already known.** [Nesci et al. 1997](https://pubmed.ncbi.nlm.nih.gov/9027966/) reported DRB1*15:03 without detectable DRB5; a [2025 cohort study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12336363/) revisited it. It is a positive-control candidate for recovering gene content beyond DRB1–DQ labels, not a new structural class. Conventional HLA typing that already includes DRB5 can also capture this distinction.

Targeted competitive alignment of matched public reads is complete (DDBJ job 20604529). HG03139 has zero high-quality read depth across the control DRB5 gene, versus 20.10× mean depth and 99.86% of gene bases at ≥5 reads in HG03195. This supports absence alongside the assembly evidence. It is not a precise breakpoint, expression or general CNV validation. Recruitment bias and the use of a single control allele remain limitations; see `results/read_support_status.md`.

## Corrected plot and reproducibility

The dendrogram label-order bug was corrected locally and in both DDBJ workflow copies, with backups retained (`results/remote_patch_log.txt`). The regenerated plot is `results/corrected_dendrogram/figures/fig9_classII_dendrogram.png`. All 753 row labels agree with the reordered tree. The old remote figure/slide files have not been silently replaced.

Run from this folder:

```sh
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python3 run_pilot.py
python3 check_results.py
python3 write_report.py
```

To regenerate the corrected tree, run `python3 ../../classII_bundles.py --dendrogram-only` from `results/corrected_dendrogram/`. Input snapshots and SHA-256 checksums are in `source/manifest.json` and the earlier panel audit. Software versions are in `results/environment.json`.

The checks cover allele parsing, reverse-complement invariance, donor grouping, duplicate exclusion, denominators, experimental pair matching and tree label ordering. No parameter was tuned to these results. UKB access and phenotype testing remain post-hackathon.

## What the first run supports

Whole-class-II bundle similarity carries substantial HLA information, but locus-level sequence similarity is a stronger baseline in this pilot. A convincing bundle-based typing method must beat or complement that baseline. The structural screen can recover variation within a three-gene type, but its first hit is known biology; the next study stage needs finer structural/sequence definitions, stronger novelty checks and fully held-out graph construction or fixed-reference assignment.

# Structural MHC typing pilot

This is a technical pilot using held-out families and MHC-recruited public short reads. It does not establish clinical accuracy, complete SV reconstruction, a new structural class or disease relevance.

**Main result:** coarse DRB gene-content typing matches all 106 assembly labels. On 88 additional RCCX donors, pilot-calibrated total copy number matches 88/88, the full annotation-signature pair 73/88, and ordered C4 forms 74/88. Calibrated ordinary reference depth also recovers copy number in 88/88: no graph-specific copy-number advantage is established. Arrangement remains the harder, incompletely resolved endpoint.

## Catalogue

All 1,182 intervals (591 haplotypes at each locus) match the corresponding full-graph walk and assembly sequence exactly and contain no ambiguous bases. The frozen panel yields four coarse DRB gene-content types and 38 RCCX annotation signatures. RCCX signatures include nearest-prototype CYP21/TNX assignments, not validated functional alleles. The DRB interval omits a full pseudogene classification; annotation absence remains provisional.

2 classical-HLA label groups contain differing structural signatures in the strict subset; see `results/same_hla_structural_contrasts.tsv`. These are panel observations, not population frequencies or proof of disease effects.

`results/extended_HLA_panel.tsv` adds the structural layer to classical HLA labels for each of the 591 eligible assembled haplotypes: DRB content/order, RCCX component count, ordered C4 forms, stable structure IDs and sequence hashes. These are experimental descriptors, not new official HLA allele names. Assembly phase links these descriptors within a source haplotype; the short-read assay returns unordered structural pairs and does not phase them to classical HLA alleles.

## Prospective benchmark

Eighteen donors were selected deterministically before read results, enriching uncommon RCCX copy-number pairs and including known controls. Before computing predictions, the unchanged method was locked for the remaining 88 eligible donors. All 106 donors contribute two loci; per-set results are below. Candidate paths, marker family support and off-locus specificity exclude the query family. Unknown relatives outside available pedigree metadata may remain. The original graph alignment used all donors; this is a fixed-graph experiment, not a rebuilt-graph benchmark. Inference uses only retained training sequence profiles, not target-only nodes or target structural labels.

Exact canonical 31-mers are sampled deterministically at 1/64. Markers must occur in at least three training families, have no extra hits outside the locus in any retained training MHC assembly, and have at most eight occurrences per training path. This is MHC specificity, not whole-genome uniqueness. Reads exclude duplicate-marked, secondary, supplementary, unmapped and QC-failed records; all marker bases require Q20, and matching mates count once per fragment.

The prototype exhaustively fits unordered pairs of observed local paths by squared error against normalized marker dosage. Depth is estimated from training-defined near-universal single-copy markers. Scores are not calibrated genotype probabilities. Sequence-wide matching can impute arrangement through linked variation without a read spanning its defining junction; it is not independent physical phasing.

| Locus | Method | Copy number correct | Structural pair correct | Unique best pair |
|---|---|---:|---:|---:|
| RCCX | path_multiplicity | 52/106 | 38/106 | 106/106 |
| RCCX | binary_paths | 76/106 | 64/106 | 101/106 |
| RCCX | assembly_profile_oracle | 106/106 | 90/106 | 105/106 |
| DRB | path_multiplicity | 106/106 | 106/106 | 106/106 |
| DRB | binary_paths | 106/106 | 106/106 | 106/106 |
| DRB | assembly_profile_oracle | 106/106 | 106/106 | 106/106 |

The binary ablation discards within-path marker multiplicity. The flat-sequence control is analytically identical to multiplicity-aware path inference because verified path sequences and markers are identical; its copied result rows record that equivalence, not an independently tested third algorithm. No graph-storage advantage is demonstrated. The assembly-profile oracle supplies exact held-out assembly marker counts to the same reference search: it diagnoses panel/representation limits under perfect measurements and is not a deployable read-typing method.

| Locus | Evaluation set | Method | Copy number correct | Structural pair correct |
|---|---|---|---:|---:|
| RCCX | pilot | path_multiplicity | 8/18 | 4/18 |
| RCCX | pilot | binary_paths | 10/18 | 8/18 |
| RCCX | pilot | assembly_profile_oracle | 18/18 | 13/18 |
| RCCX | additional_validation | path_multiplicity | 44/88 | 34/88 |
| RCCX | additional_validation | binary_paths | 66/88 | 56/88 |
| RCCX | additional_validation | assembly_profile_oracle | 88/88 | 77/88 |
| DRB | pilot | path_multiplicity | 18/18 | 18/18 |
| DRB | pilot | binary_paths | 18/18 | 18/18 |
| DRB | pilot | assembly_profile_oracle | 18/18 | 18/18 |
| DRB | additional_validation | path_multiplicity | 88/88 | 88/88 |
| DRB | additional_validation | binary_paths | 88/88 | 88/88 |
| DRB | additional_validation | assembly_profile_oracle | 88/88 | 88/88 |

## Pilot-trained correction, independent additional evaluation

The pilot showed systematic overestimation of C4 dosage. A single factor (1.1526) was fitted from the 18 pilot donors as the median marker-dosage estimate divided by assembly copy number. The secondary method divides dosage by this factor and restricts path pairs to the rounded corrected copy number. It was frozen before any predictions on the remaining 88 donors were inspected; any overlapping query family is excluded from the calibration records. The original benchmark above is unchanged.

On the 88 additional donors, the corrected method recovers total RCCX copy number in **88/88** and the complete annotation-signature pair in **73/88**. These are assembly-label agreements, not independent breakpoint or functional validation. The copy-number result is supplied by calibrated marker dosage; graph path inference cannot claim credit for that improvement.

C4-only ordered forms are separately tabulated in `results/c4_order_pairs.tsv`, avoiding reliance on uncertain CYP21/TNX prototype distinctions. The structural pair still derives from sequence-panel inference and is not direct long-range phasing.

Validation copy-number distribution: {'3': 18, '4': 64, '5': 6}. The two- and six-copy donors occurred in the pilot, so this validation does not establish performance at those extremes. Each query family is excluded separately; other evaluation donors' assembled reference paths can remain in its reference panel.

As a fairness diagnostic, the reference-depth baselines receive the same pilot-only median-ratio calibration recipe. This comparison was added after observing the secondary method result; no validation outcomes were used to fit or choose its factors.

- Calibrated reference depth / WHR1_unique_region: 88/88 total-copy calls correct.
- Calibrated reference depth / TNXB_distal_region: 88/88 total-copy calls correct.
- The independent RCCX dosage sketch has no qualifying markers for the combined A/B-and-long/short forms. Total copy number alone does not uniquely resolve the signature pair. This is an identifiability diagnostic, not a competitive dedicated C4-typing baseline.

For the corrected method, Wilson 95% intervals are 95.8%–100.0% for total copy number and 73.8%–89.4% for the full signature pair. These intervals do not include annotation, recruitment or future-cohort uncertainty.


## Simpler baselines

| Locus | Baseline | Endpoint | Correct / attempted |
|---|---|---|---:|
| RCCX | Independent total-copy markers | Total copies | 36/106 |
| DRB | Sum of independent paralog dosages | Total copies | 105/106 |
| DRB | Independent component dosage | Structural pair (no forced ambiguity resolution) | 105/106 |
| RCCX | Reference depth / WHR1_unique_region | Total copies | 58/106 |
| RCCX | Reference depth / TNXB_distal_region | Total copies | 82/106 |

The depth baselines use original GRCh38 mappings, C4 exon means, two separately reported control regions and nearest-integer rounding, without tuning against these samples. Reference and paralog mapping bias remain. Independent dosage estimates require at least ten eligible component markers; missing/ambiguous calls count as failures in all-attempted accuracy.

## Limits and next decision

See `results/predictions.tsv` for represented/unrepresented types, margins and residuals, `results/uncertainty.tsv` for family-bootstrap intervals, and `results/dosage_baseline.tsv` for per-component estimates. A unique least-squares optimum is not a confident clinical call. Bootstrap intervals describe this enriched pilot only.

Neither an existing general-purpose SV caller nor graph read alignment is benchmarked here. Novel structures absent from training can be forced onto a known pair; residuals are recorded, but there is no validated novelty detector. Full chromosome phase, exact breakpoints, functional CYP21/TNX identity, gene conversion and UKB phenotype association are not established.

Before deployment, validate discordant structures with long reads or independent assemblies, extend beyond recruited reads, calibrate no-call thresholds on separate development donors, and test the locked assay on additional donors and ancestries. Any claim of graph-specific improvement requires a graph-alignment or topology-aware method against this flat-sequence baseline.

The [failure analysis](FAILURE_ANALYSIS.md) separates C4 marginal dosage errors, joint-form phase errors, chromosome-distribution errors and uncertain prototype distinctions. The 88 evaluation donors have now been inspected; further assay tuning requires a fresh validation set.

Methodological precedent: [PanGenie, Ebler et al. 2022](https://www.nature.com/articles/s41588-022-01043-w); direct MHC/RCCX precedent: [Logsdon et al. 2025](https://www.nature.com/articles/s41586-025-09140-6). This custom pilot does not implement PanGenie.

# Additional candidate validation

Flank read array 20605749 and specificity job 20605750 completed. Probes use exact 31-base genomic windows, Phred ≥20 across the window and deduplicated read names after removal of duplicate-marked/QC-failed/secondary/supplementary reads. Both donor MHC assemblies were searched for probe matches. Whole-allele phase and regulatory function are not tested.

## Flanking candidates

HLA-C contrast HLA-WITHIN-456a27fd141204be: the downstream position 1564 allele has 13 supporting fragments in HG03195 and 17 in HG03130. HG03195 has 11 alternative-probe fragments, consistent with its two different assembled MHC haplotypes. HG03130 has no alternative-probe fragments and its observed probe occurs on both assembled haplotypes. This supports the local sequence difference; it does not independently phase it to the classical HLA allele.

DRB1 contrast HLA-WITHIN-f7a4b2ab22ceb184: HG03130 has 12 fragments supporting its downstream position 1812 allele and zero alternative-probe fragments. HG03195 has only two assembly-probe fragments versus six alternative-probe fragments; the alternative probe is absent from both assembled MHC haplotypes. Retain this as unresolved read/assembly discordance, with a possible assembly or phasing artefact. Do not promote it as a validated regulatory variant.

## RCCX depth

Job 20605817 measured original GRCh38-aligned primary read depth across C4A/B coding exons and two nearby control regions. Base quality ≥20; mapping quality ≥0 intentionally retains repetitive-region alignments. The rough diploid copy estimate is four times the mean of C4A/B exon depths divided by the control depth, because GRCh38 represents two C4 genes per haploid reference. This estimator is exploratory and uncalibrated.

| Donor | Assembly diploid C4 copies | WHR1-normalised estimate | TNXB-normalised estimate |
|---|---:|---:|---:|
| HG03195 | 4 | 4.35 | 4.51 |
| NA21144 | 4 | 4.47 | 4.33 |

The control-dependent estimates are supporting context, not validated CN calls. No new structural classification is based on them. HG02735, the other donor in the one-versus-two matched-haplotype comparison, has no local CRAM; the Saudi pair also lacks a validated read-based structural assay. Module phase, duplication junctions and CYP21/TNX functional identity remain unresolved.

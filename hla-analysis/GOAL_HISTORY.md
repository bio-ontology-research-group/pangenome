# Full analysis acceptance ledger

Authoritative objective: the five numbered requirements in Robert's attachment
`/home/leechuck/.codex/attachments/54d0a140-db8d-41a7-bbb8-642c73e5a571/pasted-text-1.txt`.
This goal remains active until all five are evidenced. The previous pilot is evidence,
not a redefinition of completion. UKB execution is post-hackathon; the required UKB
deliverable here is candidate definitions and a realistic validation/association plan.

| Requirement | Existing evidence | Work still required |
|---|---|---|
| 1. Frozen QC panel | 754 path/input IDs reconciled; 371 name-reconciled donors; clipping, contig and gene flags; public pedigree; source hashes | Freeze final source/database/software versions and complete endpoint-specific eligibility; audit further donor/relative metadata and record unresolved identities explicitly |
| 2. Leakage-controlled benchmark | Fixed-feature sequence and transductive class-II bundle pilot; experimental diploid checks | Add gene-level bundles and gene-content controls, extend classical loci, report comparable denominators and uncertainty/rarity/cohort transfer; investigate experimental discordances; specify fixed-panel scope without claiming unseen-graph validation |
| 3. Beyond-label catalogue | DR–DQ ordered gene screen; known DRB5-absence positive control | Sequence and flanking variation; signatures conditional on all available classical labels; gene CN/order; distinct C4/RCCX analysis with limitations of annotation coverage |
| 4. Evidence-ranked candidates | HG03139 DRB5 absence supported by assembly/read contrast; literature identifies known class | Current frozen IPD comparison; catalogue-wide classification; inspect coding candidates and recurrent structural/flanking candidates; matched data checks where available |
| 5. Final report/workflow | Reproducible pilot and checks | Integrate all comparisons, stable candidate IDs, source hashes, evidence and unresolved cases; conclusions on incremental value; UKB genotyping and association plan; requirement-by-requirement final audit |

## Expanded analysis checkpoint

- Pinned IPD 3.65.0 and all required genomic/CDS/protein/G-group files to commit
  `5b915f27f7f620361cf83cb626eeac8a03c0247c`; URLs and checksums recorded.
- Extracted and audited 6,670 annotated gene entries across 11 loci. Current exact
  matches, annotation warnings and sequence hashes are in `sequence_catalogue.tsv`.
- Completed common-denominator fixed-panel benchmarks for all 11 loci, comparing
  gene sequence, sampled and unsampled coding sequence, flanks, gene bundles,
  regional bundles where relevant and gene content. Family/cohort exclusions,
  rarity, ambiguity and conditional family-bootstrap intervals are recorded.
- Expanded experimental comparison to five loci. Of 14 current-CDS versus old
  experimental disagreements, 10 have possible G-group compatibility; four remain
  unresolved. This is not final discordance validation.
- Defined 47 unresolved database-absent protein candidates with stable IDs and
  nearest-reference distances. No novel allele has been established.
- Generated conditional sequence/content/order/flank contrasts and a provisional
  C4 annotation screen. Full RCCX analysis and biological assessment remain open.
- Added a concrete post-hackathon `UKB_VALIDATION_PLAN.md` and generated intermediate
  `README.md`. The report explicitly preserves outstanding requirements.

Next: inspect/read-validate the four experimental discordances and leading coding
candidates; rank and check the within-type structural/flanking contrasts; establish
the added value, if any, of combining bundles with the unsampled sequence baseline;
finish panel/metadata audit and RCCX coverage before the final completion audit.

## Evaluation constraints

- No haplotype from the target donor, duplicate assembly, or linked known/provisional
  family may provide reference labels for that target.
- Existing bundles use a frozen complete panel, including the query assemblies.
  They are fixed-panel/transductive analyses, explicitly permitted by the objective;
  they do not measure generalisation to an unseen graph-building population.
- Retain ambiguity and no-calls; separate full database novelty from protein novelty.
- Missing annotations, upstream clipped paths and assembly gaps are not biological deletions.
- Evidence beyond DRB1–DQ alone is not automatically beyond full classical HLA typing.
- Flanking-sequence variation is not proof of regulatory function.
- Record known findings, artefacts and unresolved cases; discovering a new allele or
  phenotype association is not a completion requirement.

## Read-validation and literature checkpoint

- DDBJ array 20604953: six donors completed, three lack local CRAM files (HG01252, HG01358, HG01943). No live tasks remain. All logs and results archived.
- Specificity job 20605048 completed, checking probes against both donor MHC assemblies. Three database-absent protein candidates have unique local MHC probes with support for all distinguishing SNPs; whole-allele novelty remains unconfirmed.
- NA18943 A/DRB1 and NA19007 A have local read support for assembly bases. NA18943 also agrees across HPRC/JaSaPaGe coding sequences. NA18608 DRB1 remains unresolved because some probes are nonunique.
- Tested an exploratory bundle tie-breaker on the all-coding-kmer baseline: no correct calls rescued, with some false calls introduced. Increment tables preserve both exclusion schemes.
- Archived Logsdon et al. 2025 full text and supplementary tables: direct MHC PGR-TK/RCCX precedent. Published C4 signatures match four overlapping evaluable diploid donors. Full RCCX reconstruction in this panel is still outstanding.
- Cross-source assembly comparison and targeted read-count controls are reproducible; no UKB execution has occurred.

## Structural/QC checkpoint

- `freeze_panel.py` writes complete 754-entry donor/relationship/source QC and an 8,294-row gene-endpoint eligibility table. Unknown relationships outside the public pedigree remain explicit.
- DDBJ RCCX homology job 20605384 and assembled-span QC job 20605526 completed. There are 591 retained, non-clipped, single-contig, gap-free entries with balanced components and an order signature. CYP21/TNX functional paralog and fusion classification is not inferred.
- Two all-classical-two-field-matched backgrounds differ by RCCX component counts: ksa006#1/ksa008#1 (3/2) and HG02735#2/NA21144#2 (1/2). Coordinates, span hashes, component alignments and PNG/PDF figures are saved. These known module classes are candidate biological distinctions, not new structural discoveries.
- Strict source QC leaves 34 flanking, 24 unresolved genomic, 14 known-database genomic and six correlated structural contrasts within full classical label matches. Exclusions are explicit in `within_type_evidence.tsv`.
- Direct sequence comparison to the published HGSVC3 archive finds no exact ordered-CDS-block match for the 47 coding candidates; 1,005 known genomic matches serve as positive source controls. Negative matches do not establish unpublished novelty.
- Remaining: selected structural/flanking candidate validation where data permit; consolidate all evidence into the final report and runnable workflow; audit each acceptance requirement without treating unresolved biological identity as a positive finding.

# Completion audit against the requested five deliverables

This audit concerns the analysis of the agreed 754-haplotype panel. It does not
claim confirmed new alleles, validated functional RCCX configurations, disease
associations or executed UKB validation. Those outcomes were not prerequisites
for an evidence-ranked analysis that can return known, artefactual or unresolved
findings.

| Requirement | Evidence inspected | Scope and remaining scientific uncertainty |
|---|---|---|
| 1. Frozen, QC panel | `frozen_panel.tsv`, all 8,294 rows of `endpoint_eligibility.tsv`, original path/annotation reconciliation, `validation_groups.tsv`, IPD and source manifests, `remote_source_manifest.json` | All 754 sources recorded; 742 retained haplotypes/371 name-reconciled donors. Five duplicate donor sources and known/provisional family links handled. Fragmentation/clipping and unresolved identities/relationships remain explicit; no genome-wide kinship claim. |
| 2. Leakage-controlled benchmark | `benchmark.py` exclusion assertions, per-query predictions, common-method denominator checks, cohort/rarity/call-rate metrics, family-bootstrap intervals, five-locus experimental table, all 14 discordance reviews, bundle-increment results | Fixed-panel construction is explicitly allowed by the objective and disclosed throughout. Numeric two-field endpoint, eligible complete genes, ambiguity retained. Unknown relatives, experimental resolution and lack of unseen-graph testing limit generalisation. |
| 3. Beyond-label catalogue | Per-haplotype gene/CDS/protein/flank hashes and ordered annotation signatures; separate DR–DQ and full-classical label matching; `within_type_evidence.tsv`; all-panel RCCX PAF loci, span QC, eligibility and figure | Registered and unresolved sequence contrasts separated. C4/RCCX component counts/order assessed beyond inherited C4 labels. Functional paralog identity/fusions and regulatory effects are unidentified, not assumed. |
| 4. Evidence-ranked candidates | Frozen IPD exact matches; 47 candidate IDs and actual sequences; published HGSVC3 exact-block checks/positive controls; 2025 RCCX supplementary comparison; donor-matched probe counts and MHC specificity; cross-source assemblies; DRB5 positive control; flank and C4-depth checks | Three protein candidates receive local SNP support, none is an established new allele. A tested DRB1 flank is read-discordant and remains a possible artefact. Missing CRAMs and nonunique probes are retained as unresolved. Read depth does not validate module phase or whole structural alleles. |
| 5. Final report and workflow | `FINAL_REPORT.md`, `run_analysis.py`, full stage logs/return codes, `check_results.py`, stable candidate/signature IDs, source/software hashes, `REPRODUCIBILITY.md`, `UKB_VALIDATION_PLAN.md` | Report answers what tested bundles add, what remains unidentified and which WGS/WES validations are realistic. UKB genotyping, conditional association and replication are specified for post-hackathon execution. |

## Review decisions

- The sequence baseline outperforms the tested binary bundle features. The
  exploratory combination rescues no correct calls. A negative result is retained;
  performance claims are not expanded to all graph representations.
- Exact known protein/CDS sequence is not equivalent to exact genomic identity.
  Noncoding differences and database incompleteness are not labelled new proteins.
- A match at three DR–DQ loci is not described as a match at all classical loci.
  The two highlighted RCCX backgrounds meet the latter numeric two-field definition,
  including recorded DRB3/4/5 states; unknown absence and expression status remain
  limitations.
- Current prototype labels in RCCX mean sequence similarity, not functional
  CYP21A2 or pathogenic fusion identification. Clipped/fragmented examples are not
  promoted as deletions or novel module counts.
- Local read support does not establish complete allele phase, expression,
  unpublished novelty or disease relevance. A candidate with contradictory reads
  remains unresolved rather than being omitted.
- The additional prospective experiments listed in the final report are routes to
  resolve these scientific unknowns. They are not presented as having been run.

`results/completion_evidence.json` records the final cross-artifact checks. Its
presence is not by itself proof of scientific correctness: the table above and
the source-specific limitations define what the evidence actually supports.

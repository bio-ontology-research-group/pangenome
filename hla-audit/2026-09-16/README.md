# Panel audit and agreed scope

Subsequent execution: the [first pilot](../../hla-pilot/README.md) has now run, including the corrected dendrogram, family-excluded classical-type comparisons, a structural screen and targeted DDBJ read support. The inventory audit below remains the record of the preceding phase.

Robert confirmed on 16 September 2026:

- Use `/home/asianhla/data/upload/HLA/mhc_graph/`, the 754-entry MHC panel.
- Give recovery of classical HLA types and biological information beyond classical types equal priority.
- Individual-level UK Biobank WGS/WES and phenotypes are accessible, but technical access work puts UKB validation after the hackathon.

The hackathon therefore targets public-data validation, graph/bundle recovery of known labels, and supported candidate sequence/structural differences. Phenotype association is a later endpoint. No biobank access setup is required now.

## Reconciled inventory

`panel_manifest.tsv` is one row per graph input haplotype entry. `summary.json` provides counts. Inputs in `source/` were copied from DDBJ without modifying remote files; `source/manifest.json` records provenance and SHA-256 checksums. The two graph path lists were obtained with `vg paths -t 1 -E -x` from the existing GBZ indexes.

- All **754** graph input IDs match the annotation table and the actual full and clipped graph path inventories.
- **816** input contigs correspond to **816** full graph paths. The clipped graph has **819** paths and **144,220 fewer path bases across 198 entries**. These are bases summed over haplotype paths, not unique graph sequence or established gene losses.
- The cohort inventory is APR 106; HPRC 464; JaSaPaGe 38; CPC 116; K-PanRef 28; references 2.
- **376 biological sample labels map to 371 provisional donor IDs** after merging five documented HPRC/JaSaPaGe aliases: NA18940, NA18943, NA18945, NA18952 and NA18970. This is not a pedigree or genome-based identity check; relatedness and additional aliases remain to be assessed.
- **144 entries** originate from upstream clipped CPC/K-PanRef graphs. Using the downstream full graph avoids further clipping but cannot recover those upstream losses.
- **753** entries have existing DRA–DMA spans. HIFI032450D#2 lacks that span because its anchors lie on different contigs. Its DRB1, DQA1 and DQB1 annotations are on one contig, so it need not be excluded from a narrower DR–DQ analysis.
- **750/754** entries have exactly one annotated copy of each DRB1, DQA1 and DQB1 on the same contig. This is an annotation screen, not proof of biological single copy or correct phase.
- **727/754** entries have one copy each and numeric two-field labels for all three genes. This syntax check is not independent typing truth. Noncoding novelty can coexist with a known two-field type.
- **504 entries** have at least one DR–DQ annotation flagged novel by the old annotation pipeline. This mixes coding and noncoding novelty; it is not a count of newly discovered protein alleles. For example, 466 DRB1 annotations are classified as novel noncoding. The database release must be refreshed before a novelty claim.

## Confirmed visualization defect

The existing `classII_bundles.py` records the leaf list before calling `tree.ladderize()` and annotates the drawn tree with that old order. On the saved 753-leaf tree, ladderizing changes **752 tip positions**. Thus the positional mapping for the colour strips is wrong; some colours can nevertheless coincide because labels repeat. This defect concerns the displayed annotations, not the underlying bundle distance calculation or the saved purity calculation. Neither the original remote script nor its figures has been changed in this audit.

## Pilot design resulting from the audit

1. Keep the entire agreed panel in the inventory, and use full graph paths for sequence and structural extraction. Flag clipped-source and fragmented entries; do not interpret missing annotation as deletion.
2. Reconcile allele resolution, expression/novel suffixes, and independent truth provenance. Build a donor/relative grouping manifest before any claimed held-out evaluation.
3. Use a DR–DQ block with explicit ordered and oriented paths, alongside per-gene sequence and gene-copy annotations. The existing DRA–DMA bundle analysis remains exploratory and was fitted to the whole panel.
4. Evaluate known-type recovery with donor-held-out fitting and assignment, reporting accuracy, ambiguity, no-calls and ancestry/cohort performance. Distinguish agreement with assembly annotation from independent typing accuracy.
5. In parallel, rank sequence and structural differences within the same classical allele combination, assess robustness to representation and cohort, and seek donor-matched read/assembly evidence. Carry uncertain and novel candidates separately from known-label benchmarks.
6. Export stable candidate sequence/structural definitions and evidence for post-hackathon UKB genotyping, recurrence and conditional phenotype analyses.

Reproduce the metadata audit with:

```sh
python3 hla-audit/2026-09-16/audit_panel.py
```

The inventory uses Python's standard library. The tree-order check additionally uses Biopython. The script checks graph/input/annotation ID agreement. No sequence inference, graph rebuilding, model evaluation or cluster jobs were run in this audit.

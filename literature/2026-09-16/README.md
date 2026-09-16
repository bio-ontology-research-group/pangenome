# HLA graph literature and takeover assessment — 16 September 2026

Scope: literature search, downloaded references, and read-only inspection of the local repository and DDBJ. No cluster jobs submitted, graph rebuilt, or analysis code changed. The proposed analyses below await Robert's direction. This is a targeted literature search, not a systematic review or proof that an idea is unprecedented.

## Saved articles

XML files are full-text JATS from Europe PMC; matching TXT files are readable extractions. PDFs were checked for a PDF signature and converted with pdftotext. Publisher PDF downloads that returned HTML were rejected. `download-manifest.json` records the initial download attempts; `file-manifest.json` records the saved source files and hashes, including subsequent additions. The pangene PDF is a preprint; its XML is the published article. The pan-MHC PDF is explicitly a preprint.

| Article | Saved stem | Relevance |
| --- | --- | --- |
| [Chin et al., 2023, PGR-TK](https://doi.org/10.1038/s41592-023-01914-y) | `chin2023_pgrtk` (XML/TXT) | Direct precedent: principal-bundle decomposition, clustering and HLA class II gene combinations across 105 sequences. Merely making a larger bundle plot is not a new method. |
| [Li, Marin & Farhat, 2024, pangene](https://doi.org/10.1093/bioinformatics/btae456) | `li2024_pangene` (XML/PDF/TXT) | Gene order, orientation and copy-number graphs; HLA-DRB1 is an explicit example. |
| [Zhou, Song & Li, 2024, Immuannot](https://doi.org/10.1101/gr.278985.124) | `zhou2024_immuannot` (XML/TXT) | Assembly-based annotation, allele matching and candidate novel HLA/KIR sequences. Already used in this project. |
| [Dilthey et al., 2015, MHC population reference graph](https://doi.org/10.1038/ng.3257) | `dilthey2015_prg` (XML/TXT/metadata) | Early MHC graph-based genome inference. |
| [Dilthey et al., 2019, HLA*LA](https://pmc.ncbi.nlm.nih.gov/articles/PMC6821427/) | `dilthey2019_hlala` (XML/TXT) | Graph typing and an assembly-typing route. |
| [Lee & Kingsford, 2018, Kourami](https://doi.org/10.1186/s13059-018-1388-2) | `lee2018_kourami` (PDF/TXT) | Graph-guided novel HLA allele assembly, focused on typing exons. Novel discovery with graphs is established. |
| [Huijse et al., 2023, pan-MHC preprint](https://doi.org/10.1101/2023.09.01.555813) | `huijse2023_panmhc` (PDF/TXT) | Especially close: 246 phased MHC sequences, 1,246 putative novel allele sequences across the studied loci, C4 structural haplotypes and a reference graph. That count must not be described as 1,246 new classical protein types. |
| [Wang et al., 2026, 1000 Chinese Pangenome](https://doi.org/10.1038/s41586-026-10315-y) | `wang2026_1kcp` (XML/TXT/metadata) | Recent competing context: assembly HLA typing with HiFiHLA, HLA haplotype/LD analysis and multivariant imputation. Larger Asian sample size alone is a weak novelty argument. |
| [UKB HLA WES study, 2023](https://doi.org/10.1038/s42003-023-05496-5) | `ukb2023_hla` (XML/TXT) | HLA-HD calls in 454,824 participants and disease associations; distinguishes sequence-based calling from older imputed labels. |
| [Sekar et al., 2016, C4 and schizophrenia](https://doi.org/10.1038/nature16549) | `sekar2016_c4` (XML/TXT) | Biological precedent for MHC structural haplotypes affecting expression and disease. C4 is an MHC complement locus, not a classical HLA gene. |
| [Gourraud et al., 2014, experimental HLA typing](https://doi.org/10.1371/journal.pone.0097282) | `gourraud2014_experimental_hla` (XML/TXT) | Sanger antigen-recognition-exon typing; independent diploid comparison for the pilot, with ambiguity lists retained. |
| [DRB1*15:03–DRB5 haplotypes, 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12336363/) | `drb5_absence_2025` (XML/TXT) | Direct precedent for the pilot's structural hit. DRB5 absence on this background is known biology, not a novel structural class. |

Additional official references: [UKB HLA field 22182](https://biobank.ctsu.ox.ac.uk/ukb/field.cgi?id=22182), [IPD-IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/), and [novel-allele submission guidelines](https://www.ebi.ac.uk/ipd/imgt/hla/submission/guidelines/). Official allele naming is separate from assigning research cluster identifiers.

The structural-hit follow-up also identified [Nesci et al., 1997](https://pubmed.ncbi.nlm.nih.gov/9027966/), reporting DRB1*15:03 without detectable DRB5. Its abstract was inspected; full text has not been saved. The Europe PMC full-text endpoint for the related 2008 MS paper (PMC4346327) returned 404, so it is not listed as a saved article.

## What is already present

- Local tracked repository: older Saudi pangenome workflows and `hla-study-plan/`. Local `hla-viz/` is untracked and represents 610 haplotypes. Its figure captions contain prior biological interpretations, not independently confirmed findings from this inspection.
- DDBJ workflow: `/home/leechuck/hla/cwl/`, including `scripts/`, `analysis/`, and `nig/`. This directory did not expose Git metadata; its README refers to the Claude-era workflow. Do not infer a reproducible code revision from it yet.
- Shared graph: `/home/asianhla/data/upload/HLA/mhc_graph/`. Both `MHC.gbz` and `MHC.full.gbz` exist. The graph README reports Minigraph-Cactus 3.3.0 and a 15 September build. The input seqfile has 754 lines: README composition APR 106, HPRC r2 464, JaSaPaGe 38, K-PanRef 28, CPC 116, references 2. These are input haplotype entries, not 754 independent donors or verified intact MHC chromosomes.
- Newer results: `/home/leechuck/hla/results-stage5/`, `viz5/`, `classII5/`, and `slides5/`. These contain Immuannot calls, gene FASTAs, bundles, copy-number tables, candidate novel alleles and typing comparisons. Local 610-haplotype figures are stale relative to these files.
- K-PanRef and CPC input sequences were recovered from clipped source graphs. Even the downstream **full** MHC graph cannot restore sequence already removed upstream. Original assemblies or explicit missingness masks are needed before interpreting absences.
- Donor duplication across HPRC and JaSaPaGe is documented, including an apparent duplicated assembly haplotype for JaSaPaGe NA18952. Audit donor aliases and relatives before making train/test splits. Existing homozygosity and read-support claims have not been revalidated here.

## Specific analysis concerns found by reading existing code

`analysis/classII_bundles.py` and the saved purity table were inspected, not rerun.

- The saved DRB1 two-field cluster purities are 0.387, 0.451, 0.449 and 0.524 at requested cuts of 10, 20, 40 and 80. Singleton clusters are excluded and the remainder weighted by cluster size. These are descriptive in-panel purities, not held-out typing accuracy or evidence that finer HLA types have been recovered.
- The tree-cut routine greedily splits clades using child branch lengths; it is not the distance-threshold cut described by the nearby comment. Its choice needs justification and comparison with direct distance-based clustering.
- The script captures tip order before `tree.ladderize()` and uses that earlier order for annotations on the drawn ladderized tree. If ladderizing changes tip order, the colour strips are misassigned. Verify before interpreting the displayed dendrogram.
- The displayed bundle adjacency graph discards bundle orientation and suppresses consecutive repeated bundle IDs. Bundle presence alone does not encode gene copy number or ordered structural haplotypes. Revisit the actual ordered, oriented paths and coordinate-level gene annotations.
- `typing_concordance.py` compares Immuannot with FuFiHLA applied to pseudo-reads tiled from the same assemblies. This measures method agreement, not independent validation of assembly sequence. Its normalization strips novel/expression suffixes and it chooses the first annotated copy per haplotype; the resulting concordance cannot validate novel sequence, expression status or CNV. Published/experimental truth comparisons require a separate provenance and resolution audit.
- A candidate substitution occurring elsewhere in the panel is compatible with recombination, but also paralog confusion or misassembly. It does not prove the candidate allele is real. Short-read support at a novel base does not automatically establish phase across the whole allele.

## Proposed scientific direction (not implemented)

Maintain three explicit outputs: (1) conventional HLA allele calls with resolution and ambiguity, (2) candidate novel nucleotide/protein sequences, and (3) structural MHC haplotypes defined by ordered genes, orientation, copy number and surrounding sequence. A new bundle cluster is not itself a new HLA allele or a new functional type.

For known allele recovery, extract **observed haplotype paths**, identify gene/CDS sequences and compare with a frozen current IPD release. Arbitrary walks through the graph can create combinations never observed in any donor. Bundle compression may erase single-base distinctions required by HLA nomenclature. Existing sequence annotation is the baseline against which bundle-based inference should be tested.

Start with the DRB–DQA1–DQB1 block. Compare coding-only sequence, complete-gene sequence, flanking sequence, ordered bundles, and gene content/order/CN. Learn representations and tune parameters on training donors only, then assign held-out donors with uncertainty and an unknown state. Keep both haplotypes, duplicate assemblies and close relatives together; add leave-cohort-out tests. Report coverage, ambiguity, per-gene/per-resolution accuracy and ancestry-stratified performance. If only clustering is done on the full panel, label it transductive exploration rather than unseen-donor prediction.

The strongest biological question is whether reliable structural distinctions exist **within the same classical allele combination**, and whether they add expression, functional or disease information after accounting for classical HLA alleles, local sequence variation, ancestry, relatedness and technical batch. A bundle-versus-label association alone may simply recover linkage disequilibrium. Study local blocks because recombination makes one whole-MHC cluster label too coarse. C4/RCCX is a useful separate extension and positive-control locus if scope allows.

For novelty, distinguish protein-changing alleles, synonymous CDS alleles, noncoding variants, extensions of incompletely sequenced known alleles, and new combinations of known alleles. Refresh the old IPD 3.55 comparison before calling anything novel. Require intact locus assignment and phase, inspect homopolymers, paralogs and assembly switches, and seek donor-matched reads, another technology/assembly or inheritance evidence. Use independent donors for replication. Maintain candidate status until these checks pass.

UKB imputed HLA labels can compare common known types but cannot establish a new allele absent from their reference panel. WES can support covered coding sites; WGS provides a better opportunity for structural genotyping, still with short-read identifiability limits. Validate junctions/CN and uncertainty before disease association. Prefer an initial public donor-matched validation panel; use UKB or another accessible cohort for recurrence and phenotype testing. Failure to find an Asian rare allele in UKB would not by itself refute it.

## Decisions resolved after the initial assessment

1. Robert confirmed the 754-entry shared MHC graph as the starting panel.
2. Conventional-type recovery and additional biology have equal priority.
3. Individual-level UKB WGS/WES and phenotypes are available, but validation there is deferred until after the hackathon because technical access requires additional work.

The subsequent [panel audit](../../hla-audit/2026-09-16/README.md) reconciles the actual graph paths and annotations, records donor aliases and clipping flags, and confirms the tree-annotation defect. The initial assessment above is retained as the record of what was known before that audit.

First proposed implementation milestone after confirmation: reconcile the sample/path/assembly manifest; audit existing annotation and truth labels; repair the verified visualization issues; then produce a donor-held-out DR–DQ pilot and an evidence-ranked candidate list. Do not rebuild the entire graph or launch a biobank scan before this pilot establishes what is identifiable.

## Additional direct precedent located during validation

Logsdon et al. (2025), *Complex genetic variation in nearly complete human genomes*, Nature 644:430–441, [doi:10.1038/s41586-025-09140-6](https://doi.org/10.1038/s41586-025-09140-6). Full text: `logsdon2025_complex_genetic_variation.xml` and `.txt`; supplementary tables: `logsdon2025_supplementary_tables.xlsx`. This is direct prior art for MHC PGR-TK visualisation and phased gene/pseudogene analysis of RCCX. Its Supplementary Table 53 supplies 130 RCCX haplotype annotations. Our unordered diploid C4 signatures agree in all four overlapping evaluable donors, without establishing whole-RCCX equivalence or independent assembly provenance. The generic approach and known module classes must not be claimed as novel.

Current gene naming: [NCBI Gene 8859](https://www.ncbi.nlm.nih.gov/gene/8859) names the former STK19/RP1 gene WHR1. Preserve STK19 as a historical alias when comparing RCCX papers. Current UCSC RefSeq annotation snapshot is in `hla-analysis/source/rccx_refseq_ucsc.json`; it is preparatory reference material, not a completed panel annotation.

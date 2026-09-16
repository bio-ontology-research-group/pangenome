# C4/RCCX structural screen

Reference prototypes: GRCh38 UCSC RefSeq gene spans plus the Ensembl WHR1B pseudogene span. WHR1/WHR1B are the current names corresponding to historical STK19/STK19B. Prototype coordinates, accessions, sequence and source snapshots are frozen in `source/rccx_*`. DDBJ minimap2 job 20605384 and span-QC job 20605526 completed.

This screen locates homologous components, retaining ≥85% prototype-span coverage and ≥95% alignment identity and collapsing overlapping alternative-prototype matches. C4 A/B and long/short labels are inherited from the original annotations. CYP21A1P-like versus CYP21A2-like is nearest-prototype similarity, not a functional or gene-conversion call. Shared TNX/WHR1 alignments contained in a longer homolog are not counted twice.

Across 754 entries, 750 have equal observed component counts on one contig. The conservative subset has 591 retained, non-clipped, single-contig, gap-free entries with matching counts, a usable order signature and ≥1 kb distance from contig edges. Component-count distribution: {1: 65, 2: 466, 3: 56, 4: 4}. These are panel counts, not population frequencies.

![RCCX within matching HLA types](rccx_matched_hla.png)

The ksa006#1 versus ksa008#1 contrast has three versus two C4/CYP21/TNX/WHR1 components. HG02735#2 versus NA21144#2 has one versus two. Both pairs match all available classical numeric two-field calls, including annotated DRB3/4/5 states. These distinctions therefore extend beyond the three-gene DR–DQ label comparison. The spans have no ambiguous bases. Their module-count classes are known biology; neither a new structural class nor a phenotype effect is established.

Four entries need caution: HIFI032007D#2 and HIFI032462D#2 have unbalanced components from upstream-clipped graph inputs; HIFI032164D#1 has five of each component but lacks a consistent C4 orientation signature and is also from clipped inputs; ksa004#2 has no RCCX calls in a fragmented ten-contig MHC extraction. These are not validated deletions, inversions or five-module discoveries.

Four four-component entries (HG02392#1, NA18948#2, NA18952#1, apr048#1) merit targeted structural validation. The reference-based screen does not resolve chimeric CYP21/TNX genes, functional status, pathogenicity, or exact duplication breakpoints. Short-read depth/junction or independent long-read evidence is still needed before genotype deployment.

Logsdon et al. 2025 provides direct prior art and a separate published RCCX catalogue. Its C4 patterns agree with all four overlapping evaluable diploid donors; our broader analysis should not claim the method or module classes as novel. Stable screen-signature IDs and per-haplotype coordinates are in `rccx_eligibility.tsv` and `rccx_homology_loci.tsv`.

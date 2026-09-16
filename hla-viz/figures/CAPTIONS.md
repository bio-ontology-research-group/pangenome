# HLA pangenome figures (610 haplotypes: APR 106, HPRC r2 464, JaSaPaGe Saudi 18 / Japanese 20, GRCh38, CHM13)

All calls are Immuannot (IPD-IMGT/HLA 3.55) on MHC regions cut from the assemblies with pgr-tk; per-gene
sequences are cut at the Immuannot coordinates; bundle analyses use pgr-tk 0.5.1. Scripts: `hla/analysis/` in
`leechuck/pangenome-bh26`; data tables in `data/`.

**fig3_mhc_homozygosity.png - Two HPRC assemblies are genuinely MHC-homozygous; one JaSaPaGe assembly is a duplicated haplotype.**
(a) Heterozygous SNPs per 100 kb across the MHC from 1000G high-coverage Illumina genotypes (assembly-independent).
NA18976 (JPT) and NA19909 (ASW) carry <100 heterozygous sites outside the paralog-rich HLA-A/-H and DRB windows,
versus >10,000 for typical individuals. (b) The same individuals' 1000G reads realigned to their own assembled
haplotype 1: sites with 25-75% alternative base per 50 kb. Both haplotypes of JaSaPaGe NA18952 (dashed) show the
full heterozygosity of the individual, i.e. the second haplotype is missing from that assembly, whereas the HPRC
assembly of the same individual separates the two. (c) Substitutions between the two assembled haplotypes of each
individual (minimap2 asm20): NA18976 202 and NA19909 38 differences over 5.2 Mb; JaSaPaGe NA18952 hap1 = hap2 (10 differences).

**fig4_novel_coding_alleles.png - Candidate novel coding alleles and their support.**
(a) Immuannot ":new" alleles with CDS differences in classical genes (10 haplotypes). "private" = number of
substitutions whose 31-mer occurs in no other of the 610 haplotypes; most substitutions of the multi-change alleles
(apr003 A*24, apr011 C*03, apr001 DRB1*16) are seen in hundreds of other haplotypes, i.e. these look like recombinant
alleles of known segments rather than errors. Rows in red differ from the closest allele only by 1-bp indels (HiFi
homopolymer artefacts). (b, c) Base-level pileups of the individual's own 1000G Illumina reads realigned to the
assembled contig around the novel codon: HG02717 DQB1 Ala>Asp (8/14 reads carry the assembly base, the rest the other
haplotype) and NA20346 DPA1 Ala>Met (24/51). APR/JaSaPaGe reads are not public, so those rows are assembly-only.

**fig5_population_hla.png - HLA allele landscape by cohort.**
Top: two-field allele frequencies (top 6 per cohort) for HLA-A, -B, -DRB1; e.g. A*24:02 35% and DRB1*09:01 25% in
Japanese, B*51:01 28% in Saudi, DRB1*03:01 18% and B*40:06 11% in APR. Bottom: secondary DRB gene (DR haplogroup)
frequencies, C4A/C4B long/short forms (Saudi C4B mostly long, other cohorts mostly short), and the fraction of gene
copies whose full-length sequence is absent from IPD-IMGT/HLA (60-78% for DRB1, <5% for A/C/DQB1).

**fig6_classII_haplotype_flow.png / fig7_classII_flow_by_cohort.png - Gene-level class II haplotype flow (after Chin, ASHI 2023).**
Alluvial plot of two-field alleles along DRB3/4/5 - DRB1 - DQA1 - DQB1 - DQA2 - DQB2 - TAP2 - TAP1 for all 608
haplotypes (ribbons coloured by secondary DRB gene), and per cohort for the DR-DQ block. 425 distinct gene-level
strings among 608 haplotypes; the DRB1*13-DQA1*01-DQB1*05 combination Chin flagged in HG03516/NA18906 recurs in 9
HPRC r2 haplotypes.

**fig8_classII_diplotype_pca.png - Diplotype PCA on class II two-field alleles.** Each point is a haplotype, grey lines
join the two haplotypes of an individual; colour = secondary DRB gene, labels = DRB1 allele group.

**fig9_classII_dendrogram.png - pgr-tk bundle-distance dendrogram of the DRA..DMA region (610 haplotypes).**
Colour strips: DRB1 allele group and secondary DRB gene. Sequence-level clusters follow DR haplogroups; the
purity table `data/classII_cluster_purity.tsv` gives the fraction of a cluster sharing the DRB1 allele / the
DRB345-DRB1-DQA1-DQB1 string at 10-80 clusters.

**fig10_classII_bundle_graph.png - Principal-bundle graph of the class II region.** 863 bundles from
pgr-pbundle-decomp (w=48, k=56, r=2, min_span=8); nodes = bundles (size = length, colour = number of haplotypes),
edges = consecutive bundles on a haplotype. The DRB block forms the tangled part; DQ/DO/TAP/DM is nearly linear.

**fig11_classII_bundle_pca.png - Diplotype PCA on principal-bundle presence (sequence level).** The 608 x 863 bundle
presence matrix separates DRB4 (DR7/DR9/DR4), DRB5 (DR15/DR16), DRB3 (DR3/11/12/13/14) and DRB1-only (DR1/DR8/DR10)
haplogroups on PC1/PC2 (43% of variance), reproducing Chin's slide 19 with 6x more haplotypes; grey lines join the
two haplotypes of each individual.

**Interactive / per-gene:** `gene_bundles/<GENE>.html` (pgr-tk bundle plots, 610 rows), `classII_bundle/classII.html`
(DRA..DMA), `mhc_bundle/MHC.html` (whole MHC), `gene_graphs/<GENE>.viz.png` (pggb + odgi), `plots/` (Immuannot summaries).

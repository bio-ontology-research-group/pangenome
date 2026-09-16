# HLA region extraction, typing and pangenome visualisation (CWL)

`hla_pangenome.cwl` takes haplotype-resolved human assemblies and produces, for
each haplotype, the MHC region, Immuannot HLA/KIR/C4 gene annotation with
full-resolution allele calls, per-gene sequence sets (every annotated gene
copy cut out at its Immuannot coordinates, or fetched by homology with pgr-tk),
and two visualisations per gene: pgr-tk principal-bundle plots and pggb graphs
rendered with odgi. pgr-tk (PanGenome Research Tool Kit, Chin et al. 2023) is
used wherever it fits: region extraction (`pgr-query`), optional gene fetching
(`pgr-query`) and structural visualisation (`pgr-pbundle-decomp`,
`pgr-pbundle-bed2dist`, `pgr-pbundle-bed2svg`).

## Steps

| step | tool | what it does |
|------|------|--------------|
| `prepare_reference` | samtools faidx, minimap2 -d | pull GRCh38 chr6 out of the reference, write the extended MHC query (GRCh38 chr6:28,510,120-33,480,577 +/-100 kb), build an asm20 index for the minimap2 alternative |
| `extract_pgrtk` (scatter, default) | pgr-tk `pgr-query --fastx-file` + `scripts/pgr_query_rename.py` | SHIMMER-index the whole assembly in memory (~30 s, ~9 GB per 3 Gb haplotype), fetch segments homologous to the MHC query, drop spurious hits (<20 kb or <20 anchors), name them `sample#hap#contig:start-end` (minus-strand hits reverse-complemented, `_rc`) |
| `extract_minimap2` (scatter, `extractor: minimap2`) | minimap2 asm20 + `scripts/extract_mhc.py` | alternative: map contigs to chr6 and project the MHC interval onto the contigs (about 20x slower) |
| `immuannot` (scatter; skipped when `precomputed_gtfs` is given) | Immuannot v3 (minimap2, IPD-IMGT/HLA 3.55, IPD-KIR 2.13, RefSeq C4) | gene structure and allele calls (GTF) |
| `aggregate` | pandas / matplotlib (`scripts/aggregate_immuannot.py`) | call table, haplotype x gene matrix, copy numbers, allele-frequency / diversity / novel-allele plots, one query fasta per gene (gene +/- 2 kb from GRCh38, else CHM13, else the first carrier) |
| `concat_mhc` | cat | all MHC haplotypes in one PanSN fasta (the pgr-tk database) |
| `gene_extract` (default, `gene_source: immuannot`) | samtools faidx (`scripts/extract_genes.py`) | cut every annotated copy of each gene (+/- 2 kb, gene orientation) out of `MHC.fa` at the Immuannot coordinates; records named `sample#hap#GENE[_n]` |
| `gene_fetch` (scatter over genes, `gene_source: pgr-query`) | pgr-tk `pgr-query` | alternative: fetch each gene from every MHC haplotype by homology to the reference gene. Defaults w=24/k=32 recover all 610 HLA-A/-C haplotypes; HLA-DRB1 haplogroups diverge so far in the introns that only w=16/k=24 finds them all (w=48/k=56 found 223 HLA-A and 113 DRB1) |
| `gene_bundle` (scatter) | pgr-tk `pgr-pbundle-decomp`, `bed2dist`, `bed2sorted`, `bed2svg` | MAP-graph GFA, principal-bundle bed, bundle distance + dendrogram, SVG/HTML bundle plot per gene |
| `gene_graph` (scatter) | pggb 0.7.4 | one graph per HLA gene (`-p 90 -s 2000 -n <sequences>`) |
| `gene_viz`, `gene_draw` | odgi viz / odgi layout + draw | 1D (rows = haplotypes) and 2D renderings |
| `mhc_bundle` | pgr-tk bundle tools (w=80, k=56, r=6, min_span=28) | principal-bundle decomposition and plot of the whole MHC across all haplotypes |
| `mhc_graph`, `mhc_viz` | pggb, odgi | optional whole-MHC graph (`build_mhc_graph: true`, sparsified with `mhc_graph_n`) |

Use the Immuannot `consensus` field (column `consensus` in `hla_calls.tsv`) as
the allele call; `:new` marks alleles absent from IPD-IMGT/HLA.

## Inputs

`scripts/make_inputs.py` scans assembly folders and writes the job YAML
(parallel arrays `assemblies`, `samples`, `haplotypes`, `cohorts`). Naming
rules are in the script docstring; HPRC `pat`/`mat` map to haplotype 1/2.
The GRCh38 and CHM13 PanSN fastas are picked up as cohort `REF`.
JaSaPaGe is split into `JaSaPaGe-Saudi` (ksa*) and `JaSaPaGe-Japanese` (1000G JPT NA*). With `--populations
data/hprc_r2_populations.tsv` HPRC r2 is split by population (1000G sequence indexes + HPRC r2 sample metadata):
`HPRC-Japanese` (JPT), `HPRC-Jewish` (HG002, Ashkenazi; the only Jewish individual in HPRC r2), `HPRC-EastAsian`
(CHB, CHS, CDX, KHV, HG005) and `HPRC-Rest`.

K-PanRef (Korean pangenome, 16 samples, Minigraph-Cactus graph; assembly FASTAs not released) is read from the graph:
`nig/kpanref_mhc.sbatch` pulls every haplotype path of each sample from `KPanRef.gbz` (`vg paths -F`), keeps the path
segments with an alignment of >= 50 kb matching bases (MAPQ >= 20) to the GRCh38 MHC, and writes
`<sample>.hap<N>.mhcctg.fa` (cohort `KPanRef-Korean`). `nig/add_kpanref.sbatch` then runs extraction and Immuannot on
them, adds them to the precomputed MHC/GTF folders and reruns aggregation, gene graphs and the class II span.

```bash
python3 scripts/make_inputs.py \
  --dir APR=/home/asianhla/data/upload/APR/assemblies \
  --dir HPRC_r2=/home/asianhla/data/HPRC_r2/fasta \
  --dir JaSaPaGe=/home/asianhla/data/JaSaPaGe/assembly_clean/fasta \
  --dir KPanRef=~/hla/kpanref/fasta --populations data/hprc_r2_populations.tsv \
  --reference /home/asianhla/data/HPRC_r2/fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set.PanSN.fa \
  --immuannot-dir ~/hla/Immuannot --immuannot-ref ~/hla/Data-2024Feb02 \
  --out nig/inputs-nig.yml            # add --per-cohort-limit 2 for a smoke test
```

## Running

Anywhere with Docker/Singularity (tools carry `DockerRequirement` hints):

```bash
cwltool --parallel hla_pangenome.cwl inputs.yml
```

On the NIG BioHackathon node there is no container runtime available to users
(no docker daemon; unprivileged user namespaces are blocked), so tools come
from a micromamba environment and cwltool runs with `--no-container`:

```bash
bash nig/setup_env.sh                                  # once, on a login node
sbatch nig/run_nig.sbatch nig/inputs-nig.yml ~/hla/results
```

`nig/run_nig.sbatch` requests 28 cores / 200 GB on partition `asianhla-c32`
and runs `toil-cwl-runner --batchSystem single_machine --no-container`.
cwltool's `--parallel` executor deadlocked after ~390 of the 610 extraction
jobs (process alive, no children), so Toil is used; pass `restart` as third
argument to resume from the job store after a failure. The Toil venv is first
on PATH, so its `python3` runs the CWL python scripts: `setup_env.sh` installs
pandas and matplotlib into it.

The 610-haplotype run on the NIG node was done in stages, each a plain run of
the same workflow with precomputed intermediates passed back in:

1. extraction only (`~/hla/mhc_all`: `<sample>_<hap>.mhc.fa/.tsv`, 8 cores per
   pgr-query, ~30 s each);
2. `nig/stage2.sh` - `mhc_fastas`/`mhc_tsvs` given, Immuannot on 610 haplotypes
   (4 threads each, 7 in parallel, ~2 h) plus everything downstream;
3. `nig/stage3.sh <mhc_dir> <gtf_dir>` - `precomputed_gtfs` given too, so only
   aggregation, gene extraction, bundles and graphs rerun (~1.5 h, dominated by
   the whole-MHC bundle decomposition of 3 Gb).

Samples present in more than one cohort (the 1000G JPT trios NA18940/43/45/52/70
are in both HPRC r2 and JaSaPaGe) keep the first label and get `<sample>.<cohort>`
for the later cohort. Assemblies whose contigs are already PanSN-named (HPRC r2,
the references, JaSaPaGe NA*) are not prefixed a second time.

## Outputs

- `*.mhc.fa`: extracted MHC segments per haplotype; `MHC.fa`: all of them
- `*.gtf.gz`: Immuannot annotation per haplotype
- `hla_calls.tsv`, `hla_calls_matrix.tsv`, `gene_copy_number.tsv`, `mhc_extraction_summary.tsv`
- `plots/`: `allele_freq_<GENE>.png`, `allele_diversity.png`, `novel_allele_rate.png`, `gene_copy_number.png`, `mhc_extraction_coverage.png`
- per gene: `<GENE>.fa` (gene sequences), `<GENE>.regions.tsv` (Immuannot cut) or `<GENE>.hits.tsv` (pgr-query fetch), `<GENE>.svg` / `<GENE>.html` (pgr-tk bundle plot), `<GENE>.bed`, `<GENE>.pmapg.gfa`, `<GENE>.nwk`, `<GENE>.ctg.summary.tsv`, `<GENE>.gfa` / `<GENE>.og` (pggb), `<GENE>.stats.tsv`, `<GENE>.viz.png`, `<GENE>.viz_depth.png`, `<GENE>.draw.png`
- whole MHC: `MHC.svg` / `MHC.html`, `MHC.bed`, `MHC.pmapg.gfa`, `MHC.nwk`, `MHC.ctg.summary.tsv`; plus `MHC.gfa` / `MHC.og` / `MHC.viz.png` when `build_mhc_graph` is set

## Analyses (`analysis/`, `results/`)

`results/figures/CAPTIONS.md` describes the figures; `results/tables/` holds the call table (`hla_calls.tsv.gz`,
cohorts APR / JaSaPaGe-Saudi / HPRC-Jewish / JaSaPaGe-Japanese / HPRC-Japanese / KPanRef-Korean / HPRC-EastAsian /
HPRC-Rest / REF; order, labels and colours in `analysis/cohorts.py`), copy numbers, class II gene-level haplotype strings,
novel-allele support and 1000G heterozygosity tracks. Scripts (run inside a directory holding the workflow outputs):

- `analysis/fig5_population.py` - allele-frequency landscape by cohort, DR haplogroups, C4 forms, novel-allele fraction
- `scripts/read_support.sh` + `analysis/fig3_mhc_homozygosity.py` - MHC-region reads from the public 1000G high-coverage
  CRAMs realigned to the individual's own assembly (pileup heterozygosity), hap1-vs-hap2 divergence, 1000G SNP density
- `scripts/novel_alleles.py` + `analysis/fig4_novel_alleles.py` - Immuannot `cds_mut` to contig coordinates, 31-mer
  search of each novel variant across all haplotypes, base-level read pileups at the novel codons
- `analysis/classII_flow_pca.py` - gene-level class II haplotype strings, alluvial haplotype flow (all / per cohort),
  diplotype PCA on two-field alleles (after Chin, ASHI 2023, speakerdeck.com/jchin/human-pangenome-graph-analysis-for-the-mhc)
- `scripts/extract_span.py` + `tools/pgr_pbundle.cwl` + `analysis/classII_bundles.py` - DRA..DMA span from every
  haplotype, pgr-tk bundle decomposition, dendrogram coloured by DRB1 group, bundle graph, bundle-presence diplotype PCA

## References

- pgr-tk: Chin et al., "Multiscale analysis of pangenomes enables improved representation of genomic diversity for repetitive and clinically relevant genes", Nature Methods 2023; https://github.com/GeneDx/pgr-tk (v0.5.1 binaries)
- Immuannot: Zhou et al., Genome Research 2024, doi:10.1101/gr.278985.124; data bundle Zenodo 10.5281/zenodo.10948964
- pggb: Garrison et al., Nature Methods 2024; odgi: Guarracino et al., Bioinformatics 2022
- minimap2: Li, Bioinformatics 2018

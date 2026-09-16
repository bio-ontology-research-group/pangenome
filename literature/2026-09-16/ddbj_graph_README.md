# Whole-MHC pangenome graph of 754 haplotypes (`/home/asianhla/data/upload/HLA/mhc_graph/`)

Built 15 September 2026 with Minigraph-Cactus (Cactus 3.3.0 static release, `cactus-pangenome`) on the NIG supercomputer
(32 cores, 240 GB, 3 h 43 min). Job script: `workflow/nig/mhc_mc_graph.sbatch`.

## Input

- The 754 MHC haplotype fastas of the workflow (`../mhc_fastas/`): APR 106, HPRC release 2 464, JaSaPaGe 38, K-PanRef 28,
  CPC 116, GRCh38 and CHM13 (MHC segments from `extract_mhc`, GRCh38 region chr6:28,510,120-33,480,577 plus flanks)
- Sample names: `<sample>.<haplotype>` (e.g. `HG00096.1`), references `GRCh38` and `CHM13`
- Contigs were renamed `ctg<N>` because the PanSN names (`sample#hap#contig:start-end`) clash with Cactus/vg path naming;
  `contig_names.tsv` maps sample, haplotype, new name, original record name

## Command

```bash
cactus-pangenome js seqfile.txt --outDir out --outName MHC --reference GRCh38 CHM13 \
  --gfa clip full --gbz clip full --giraffe clip --vcf --odgi full --viz
```

## Files

- `MHC.gbz` - clipped graph with haplotypes (GRCh38 reference; sequence not aligned in minigraph is clipped): use for mapping
- `MHC.dist`, `MHC.shortread.withzip.min`, `MHC.shortread.zipcodes` - vg giraffe indexes for `MHC.gbz`
- `MHC.gfa.gz` - clipped graph as GFA; `MHC.full.gbz`, `MHC.full.gfa.gz`, `MHC.full.og` - unclipped graph (GBZ, GFA, odgi)
- `MHC.vcf.gz` (+ `.tbi`) - variants against GRCh38 (decomposed); `MHC.raw.vcf.gz` - raw `vg deconstruct` output
- `MHC.full.hal` - Cactus alignment; `MHC.snarls`, `MHC.full.snarls`; `MHC.stats/`, `MHC.viz/` - statistics and odgi views
- `MHC.sv.gfa.gz` - minigraph SV graph; `MHC.paf`, `MHC.gaf.gz` - haplotype-to-graph mappings
- Graph size: 417,896 nodes, 579,137 edges, 5.89 Mb of sequence

## Mapping short reads (example: one 1000G sample from `../../1000G_MHC/cram/`)

```bash
samtools collate -u -O --reference GRCh38_full_analysis_set_plus_decoy_hla.fa HG00096.mhc.cram tmp \
  | samtools fastq -F 0x900 -1 r1.fq.gz -2 r2.fq.gz -0 /dev/null -s /dev/null -n -
vg giraffe -t 16 -Z MHC.gbz -d MHC.dist -m MHC.shortread.withzip.min -z MHC.shortread.zipcodes \
  -f r1.fq.gz -f r2.fq.gz -o gaf > HG00096.gaf
```

Use the vg bundled with Cactus 3.3.0 (`~leechuck/hla/cactus/cactus-bin-v3.3.0/bin/vg`) or a vg version that reads these
indexes. Reference paths are `GRCh38#0#ctg<N>` / `CHM13#0#ctg<N>`; `contig_names.tsv` gives the GRCh38 MHC coordinates.

## Caveats

- The graph covers only the MHC: reads from paralogous or HLA-like sequence elsewhere in the genome can be forced onto it.
  The 1000G CRAM extracts contain only reads that bwa placed in the MHC, the chr6 alternative haplotypes or HLA contigs.
- K-PanRef and CPC haplotypes come from Minigraph-Cactus graphs themselves (clipped paths), not from the original assemblies.
- JaSaPaGe NA18952 carries one haplotype twice (see the main README).

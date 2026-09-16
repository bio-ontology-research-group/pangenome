# Reproducing the HLA analysis

The entry point is `FINAL_REPORT.md`. The original pilot and intermediate README
are historical checkpoints; neither overrides the final evidence classifications.

## Local analytic workflow

From the repository root:

```sh
python3 hla-analysis/run_analysis.py
```

This executes the ordered stages listed in `run_analysis.py`, stops on the first
failure and writes per-stage logs and return codes to `results/reproduction/`.
It recomputes allele extraction, candidate comparisons, the 11-locus benchmark,
uncertainty, within-type contrasts, eligibility, evidence summaries and final report.
It consumes archived remote measurements; it does not pretend to regenerate reads,
assemblies or the supplied pgr-tk decomposition.

Dependencies: Python 3.13, Biopython, NumPy, SciPy, Matplotlib and openpyxl. Exact
observed versions and platform are recorded in `results/artifact_manifest.json`
and the earlier pilot environment record. No scikit-learn or GPU is required.
The driver fixes BLAS/OpenMP threads at two. Family-bootstrap seed is 20260916.

The following sibling folders must be retained:

- `hla-audit/2026-09-16/`: panel, original calls, source metadata and provenance.
- `hla-pilot/source/`: regional bundles, pedigree and experimental typing table.
- `hla-pilot/results/validation_groups.tsv`: reconciled donor/family exclusions,
  reproducible with the pilot's `run_pilot.py` and audited in the expanded workflow.
- `literature/2026-09-16/`: full texts and the published HGSVC3 supplementary workbook.

All 33 gene FASTA/region/bundle files and 754 GTFs are archived under `source/`.
IPD files are pinned to commit `5b915f27f7f620361cf83cb626eeac8a03c0247c` (3.65.0).
Source manifests record hashes, URLs and remote paths. The downloaded HGSVC3 allele
archive is a published comparison set, not a replacement truth set.

## Archived remote stages

Remote commands use SSH alias `ddbj`, with Slurm submission through its `a001`
interactive gateway. Compute partition/account are `asianhla-c32`/`asianhla-group`.
The remote environment is `/home/leechuck/hla/mm/envs/hla/bin`; minimap2 2.28 and
samtools 1.21 were used. The working CRAM decode reference is the actual file
`/home/leechuck/hla/ref1kg/local/GRCh38_full_analysis_set_plus_decoy_hla.fa`.
The similarly named symlink outside `local/` was unusable on the compute node.

| Stage | Script and remote work directory | Recorded job(s) | Local evidence |
|---|---|---|---|
| DRB5 positive control | `../hla-pilot/drb5_read_support.sbatch`; `codex-pilot/drb5/` | 20604529 (20604506 failed reference lookup) | Pilot read-support archive/report |
| Coding/experimental SNP probes | `read_markers.sbatch`, `count_read_markers.py`; `codex-analysis/read-markers/` | 20604953 array | `source/read_support_20604953.tar.gz` |
| Probe MHC specificity | `marker_specificity.py`, same directory | 20605048 | `results/read_support/marker_specificity.tsv` |
| RCCX homology | `rccx_map.sbatch`; `codex-analysis/rccx/` | 20605384 | `source/rccx_map_20605384.tar.gz` |
| RCCX span gaps/edges | `rccx_span_qc.py`, same directory | 20605526 | `results/rccx_span_qc.tsv` and log |
| Flanking probes | `flank_markers.sbatch`; `codex-analysis/flank-validation/` | 20605749 array; specificity 20605750 | `source/flank_validation_20605749.tar.gz` |
| Exploratory C4 depth | `rccx_depth.sbatch`; `codex-analysis/rccx/` | 20605817 | `source/rccx_depth_20605817.tar.gz` |

Remote directories above are beneath `/home/leechuck/hla/`. The original MHC FASTAs
are in `mhc_all/`; public recruited CRAMs are in
`/home/asianhla/data/upload/1000G_MHC/cram/`. Their paths, availability and SHA-256
hashes are recorded in `source/remote_source_manifest.json`. Three requested coding
candidate CRAMs are absent; their failed array tasks are retained, not hidden or
reported as successful validation.

To recompute a remote stage, copy its script and explicitly referenced inputs into
the recorded work directory, verify source hashes, and submit through Slurm. For
example, from the DDBJ gateway:

```sh
ssh -o BatchMode=yes a001 sbatch /home/leechuck/hla/codex-analysis/rccx/rccx_map.sbatch
```

The coding array uses indices 0–8 against `read_marker_donors.txt`; the flank array
uses 0–1 in its separate directory. `prepare_read_markers.py` and
`prepare_flank_markers.py` regenerate the probes. Specificity jobs run
`python3 marker_specificity.py` from the corresponding directory. The RCCX span
job runs `python3 rccx_span_qc.py` after copying the generated `rccx_screen.tsv`.
Only transfer completed measurement tables/logs into the corresponding local
directories; do not overwrite one donor's `support.tsv` with another's.

The local symlink `source/rccx_paf -> rccx_run/paf` exposes the archived PAF files.
Preserve it or recreate it after extracting `rccx_map_20605384.tar.gz` into
`source/rccx_run/`. The decompositions themselves are fixed inputs: rebuilding the
graph or changing pgr-tk parameters defines a new experiment.

## Verification contracts

`check_results.py` exercises strand coordinate extraction, exon/stop interval union,
IPD hashes, a known partial-CDS exclusion, a newly registered allele control,
candidate sequence IDs, paired-read counting/quality filtering, common benchmark
denominators, reference RCCX component counts and conservative structural QC.
Training-family exclusion is asserted while generating each benchmark prediction.
The full workflow log complements these focused tests; neither proves clinical
accuracy, novelty or biological function.

Final file hashes are recorded after reproduction and review. Re-running figure
generation may change PDF metadata without changing the underlying data. Treat TSV
tables, sequence definitions and source hashes as the quantitative record.
